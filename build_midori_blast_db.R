## build_midori_blast_db.R
## Rebuilds a MIDORI2 COI BLAST database with NCBI taxonomy embedded, so that
## blastn output fields sscinames, scomnames, and staxids are populated.
##
## MIDORI2 headers already contain full NCBI taxids in the lineage string, e.g.:
##   >AB000001.1.<1.>648###root_1;Eukaryota_2759;...;Gadus_morhua_8049
## The last _<number> in the lineage is the species-level NCBI taxid.
##
## Steps:
##   1. Extract FASTA from the existing (v4, no-taxonomy) BLAST database
##   2. Parse seqid and taxid from every header
##   3. Write a taxid_map (seqid<TAB>taxid)
##   4. Rebuild with makeblastdb -blastdb_version 5 -taxid_map
##
## External dependencies:
##   blastdbcmd and makeblastdb must be on PATH (NCBI BLAST+)

## ============================================================================
## Section 1: Configuration
## ============================================================================

SOURCE_DB  <- "/Volumes/Clupea/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST"
OUTPUT_DIR <- "/Volumes/Clupea/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST"   # rebuild in-place
DB_NAME    <- "MIDORI2_UNIQ_NUC_GB270_CO1_BLAST"

## ============================================================================
## Section 2: Validation
## ============================================================================

stopifnot(
  "blastdbcmd not on PATH"  = nchar(Sys.which("blastdbcmd"))  > 0,
  "makeblastdb not on PATH" = nchar(Sys.which("makeblastdb")) > 0
)

DB_FASTA   <- file.path(OUTPUT_DIR, paste0(DB_NAME, "_clean.fasta"))
TAXID_MAP  <- file.path(OUTPUT_DIR, paste0(DB_NAME, "_taxid_map.txt"))
DB_OUT     <- file.path(OUTPUT_DIR, DB_NAME)

message("Source DB  : ", SOURCE_DB)
message("Output FASTA: ", DB_FASTA)
message("Output DB  : ", DB_OUT)

## ============================================================================
## Section 3: Extract FASTA from existing database
## ============================================================================

message("\nExtracting FASTA from existing database...")
t0 <- proc.time()

extract_exit <- system2("blastdbcmd",
  args   = c("-db", shQuote(SOURCE_DB), "-entry", "all",
             "-outfmt", "%f", "-out", shQuote(DB_FASTA)),
  stdout = TRUE, stderr = TRUE
)
exit_code <- attr(extract_exit, "status")
if (!is.null(exit_code) && exit_code != 0) {
  stop("blastdbcmd failed:\n", paste(extract_exit, collapse = "\n"))
}

elapsed <- proc.time() - t0
message(sprintf("  Done in %.0f s. Size: %.0f MB",
  elapsed["elapsed"],
  file.info(DB_FASTA)$size / 1e6))

## ============================================================================
## Section 4: Sanitize FASTA and build taxid map
## ============================================================================

# MIDORI2 seqids embed the full taxonomy lineage, making them >50 chars and
# containing <, >, # which makeblastdb -parse_seqids rejects.
# Strategy: strip everything from ### onward (that's the taxonomy), remove
# < and > from the coordinate fields, collapse any resulting double-dots.
# Taxid is extracted from the last _<number> in the taxonomy string.

message("\nSanitizing FASTA headers and building taxid map...")

RAW_FASTA <- file.path(OUTPUT_DIR, paste0(DB_NAME, ".fasta"))  # blastdbcmd output

py_script <- tempfile(fileext = ".py")
writeLines(c(
  "import re, sys",
  paste0("fasta_in  = '", RAW_FASTA,  "'"),
  paste0("fasta_out = '", DB_FASTA,   "'"),
  paste0("tmap_out  = '", TAXID_MAP,  "'"),
  "MAX_LEN = 50",
  "n_seqs = n_taxid = n_trunc = n_dup = 0",
  "seen = {}",
  "def sanitize(raw):",
  "    s = raw.replace('<','').replace('>','').replace(',','_')",
  "    while '..' in s: s = s.replace('..','.')",
  "    return s.rstrip('.')",
  "with open(fasta_in) as fin, open(fasta_out, 'w') as fout, open(tmap_out, 'w') as tmap:",
  "    for line in fin:",
  "        if line.startswith('>'):",
  "            header = line[1:].rstrip()",
  "            if '###' in header:",
  "                seqid_raw, tax = header.split('###', 1)",
  "                m = re.search(r'_(\\d+)$', tax)",
  "                taxid = m.group(1) if m else None",
  "            else:",
  "                seqid_raw, taxid = header, None",
  "            seqid = sanitize(seqid_raw)",
  "            if len(seqid) > MAX_LEN:",
  "                n_trunc += 1",
  "                seqid = seqid[:MAX_LEN - 4]",
  "            base = seqid",
  "            if base in seen:",
  "                seen[base] += 1",
  "                seqid = f'{base}_{seen[base]:03d}'",
  "                n_dup += 1",
  "            else:",
  "                seen[base] = 0",
  "            n_seqs += 1",
  "            if taxid:",
  "                n_taxid += 1",
  "                tmap.write(f'{seqid}\\t{taxid}\\n')",
  "            fout.write(f'>{seqid}\\n')",
  "        else:",
  "            fout.write(line)",
  "print(f'Sequences : {n_seqs:,}', flush=True)",
  "print(f'With taxid: {n_taxid:,}', flush=True)",
  "print(f'Truncated : {n_trunc:,}', flush=True)",
  "print(f'Duplicates: {n_dup:,}',  flush=True)"
), py_script)

py_result <- system2("python3", args = py_script, stdout = TRUE, stderr = TRUE)
message(paste(py_result, collapse = "\n"))
unlink(py_script)

if (!file.exists(DB_FASTA) || file.info(DB_FASTA)$size == 0) {
  stop("Sanitized FASTA not created: ", DB_FASTA)
}
message("  Taxid map written: ", TAXID_MAP)

## ============================================================================
## Section 5: Rebuild with makeblastdb
## ============================================================================

message("\nRunning makeblastdb...")

db_result <- system2("makeblastdb",
  args = c(
    "-in",           shQuote(DB_FASTA),
    "-input_type",   "fasta",
    "-dbtype",       "nucl",
    "-taxid_map",    shQuote(TAXID_MAP),
    "-parse_seqids",
    "-blastdb_version", "5",
    "-title",        DB_NAME,
    "-out",          shQuote(DB_OUT)
  ),
  stdout = TRUE, stderr = TRUE
)
db_exit <- attr(db_result, "status")
message(paste(db_result, collapse = "\n"))
if (!is.null(db_exit) && db_exit != 0) {
  stop("makeblastdb failed (exit ", db_exit, ")")
}

message("\nDone. New taxonomy-aware database: ", DB_OUT)
message("The existing taxdb.btd / taxdb.bti files in the same directory")
message("will be used by BLAST for sscinames/scomnames lookups.")
