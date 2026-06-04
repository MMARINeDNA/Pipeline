## build_bold_blast_db.R
## Converts the BOLD COI reference database into a taxonomy-aware local BLAST database.
##
## Source data: /Volumes/Clupea/20514192/ (coidb v0.6.0, BOLD Data Package 2026-01-30)
##   - coidb.clustered.fasta.gz : one sequence per BOLD BIN, headers: >{processid} bin_uri:{BIN}
##   - coidb.BOLD_BIN.consensus_taxonomy.inclNA.tsv.gz : BIN-level taxonomy
##   - coidb.BOLD_BIN.consensus_taxonomy.exclNA.tsv.gz : BIN-level taxonomy (excludes NAs)
##
## External dependencies (must be on PATH or set below):
##   - taxonkit  : https://bioinf.shenwei.me/taxonkit/
##   - makeblastdb : part of NCBI BLAST+
##   - python3   : for streaming FASTA filter
##
## NCBI taxonomy data for taxonkit must be in ~/.taxonkit/ (or TAXONKIT_DB below).
## Download: https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Configuration — edit these before each run
# ==============================================================================

SOURCE_DIR    <- "/Volumes/Clupea/20514192"       # folder with the downloaded coidb files
OUTPUT_DIR    <- "/Volumes/Clupea/bold_coi"       # where the BLAST DB files will be written
DB_NAME       <- "bold_coi"                       # prefix for all output files

TAXONKIT_PATH <- "/Users/rpk/taxonkit"            # path to taxonkit binary
TAXONKIT_DB   <- "~/.taxonkit"                    # NCBI taxonomy data dir for taxonkit

USE_INCL_NA   <- TRUE    # TRUE = inclNA taxonomy (more conservative); FALSE = exclNA

# ==============================================================================
# Section 2: Derived paths + validation
# ==============================================================================

tax_suffix  <- if (USE_INCL_NA) "inclNA" else "exclNA"
FASTA_GZ    <- file.path(SOURCE_DIR, "coidb.clustered.fasta.gz")
TAX_TSV_GZ  <- file.path(SOURCE_DIR, paste0("coidb.BOLD_BIN.consensus_taxonomy.", tax_suffix, ".tsv.gz"))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

DB_FASTA    <- file.path(OUTPUT_DIR, paste0(DB_NAME, ".fasta"))
TAXID_MAP   <- file.path(OUTPUT_DIR, paste0(DB_NAME, "_taxid_map.txt"))
NOHIT_LOG   <- file.path(OUTPUT_DIR, paste0(DB_NAME, "_no_taxid.tsv"))
DB_OUT      <- file.path(OUTPUT_DIR, DB_NAME)

stopifnot(
  "SOURCE FASTA not found"     = file.exists(FASTA_GZ),
  "TAXONOMY TSV not found"     = file.exists(TAX_TSV_GZ),
  "taxonkit binary not found"  = file.exists(TAXONKIT_PATH),
  "makeblastdb not on PATH"    = nchar(Sys.which("makeblastdb")) > 0,
  "python3 not on PATH"        = nchar(Sys.which("python3")) > 0
)

message("Source FASTA : ", FASTA_GZ)
message("Source taxonomy : ", TAX_TSV_GZ)
message("Output directory: ", OUTPUT_DIR)

# ==============================================================================
# Section 3: Parse FASTA headers  →  processid / bin_uri table
# ==============================================================================

message("\nParsing FASTA headers...")

# Stream headers only; grep returns lines like ">AANIC001-09 bin_uri:BOLD:AAA0008"
raw_headers <- system2("bash",
  args = c("-c", shQuote(paste0("gunzip -c '", FASTA_GZ, "' | grep '^>' || true"))),
  stdout = TRUE, stderr = FALSE
)

header_df <- tibble(raw = raw_headers) %>%
  mutate(
    processid = str_extract(raw, "(?<=^>)\\S+"),
    bin_uri   = str_extract(raw, "(?<=bin_uri:)\\S+")
  ) %>%
  select(processid, bin_uri) %>%
  filter(!is.na(bin_uri))

message("  Sequences with bin_uri: ", nrow(header_df))

# ==============================================================================
# Section 4: Load BIN taxonomy; clean "unresolved.*" placeholders
# ==============================================================================

message("\nLoading BIN taxonomy (", tax_suffix, ")...")

bin_tax <- read_tsv(TAX_TSV_GZ, show_col_types = FALSE) %>%
  rename_with(tolower)

# coidb uses "unresolved.Rank" strings where rank is unknown — treat as NA
clean_unresolved <- function(x) {
  ifelse(grepl("^unresolved", x, ignore.case = TRUE), NA_character_, x)
}
bin_tax <- bin_tax %>%
  mutate(across(c(kingdom, phylum, class, order, family, genus, species), clean_unresolved))

message("  BINs loaded: ", nrow(bin_tax))

# ==============================================================================
# Section 5: Map taxonomy names to NCBI taxids via taxonkit name2taxid
# ==============================================================================

lookup_taxids <- function(names_vec, taxonkit_path, taxonkit_db) {
  unique_names <- unique(na.omit(names_vec))
  if (length(unique_names) == 0) return(setNames(character(0), character(0)))

  tmp_in  <- tempfile(fileext = ".txt")
  tmp_out <- tempfile(fileext = ".txt")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  writeLines(unique_names, tmp_in)

  system2(taxonkit_path,
    args   = c("name2taxid", tmp_in, "-o", tmp_out, "--data-dir", taxonkit_db),
    stdout = FALSE, stderr = FALSE
  )

  if (!file.exists(tmp_out) || file.info(tmp_out)$size == 0) {
    return(setNames(rep(NA_character_, length(unique_names)), unique_names))
  }

  res <- read_tsv(tmp_out, col_names = c("name", "taxid"), col_types = "cc", progress = FALSE) %>%
    filter(!is.na(taxid) & taxid != "") %>%
    group_by(name) %>% slice(1) %>% ungroup()

  out <- setNames(rep(NA_character_, length(unique_names)), unique_names)
  out[res$name] <- res$taxid
  out
}

# Hierarchical fallback: species → genus → family → order → class → phylum
# We batch-lookup each rank level, then fill in gaps with coarser ranks.

ranks <- c("species", "genus", "family", "order", "class", "phylum")

message("\nLooking up NCBI taxids (6 rank levels, batched)...")
taxid_maps <- list()
for (rk in ranks) {
  message("  rank: ", rk, " (", sum(!is.na(bin_tax[[rk]])), " non-NA names)")
  taxid_maps[[rk]] <- lookup_taxids(bin_tax[[rk]], TAXONKIT_PATH, TAXONKIT_DB)
}

# Build per-BIN taxid with fallback
bin_tax <- bin_tax %>%
  mutate(taxid = NA_character_)

for (rk in ranks) {
  tmap <- taxid_maps[[rk]]
  still_missing <- is.na(bin_tax$taxid)
  bin_tax$taxid[still_missing] <- tmap[bin_tax[[rk]][still_missing]]
}

n_resolved   <- sum(!is.na(bin_tax$taxid))
n_unresolved <- sum( is.na(bin_tax$taxid))
message(sprintf("  Resolved: %d / %d BINs (%.1f%%); unresolved: %d",
  n_resolved, nrow(bin_tax), 100 * n_resolved / nrow(bin_tax), n_unresolved))

# Log unresolved BINs for inspection
if (n_unresolved > 0) {
  bin_tax %>%
    filter(is.na(taxid)) %>%
    write_tsv(NOHIT_LOG)
  message("  Unresolved BINs written to: ", NOHIT_LOG)
}

# ==============================================================================
# Section 6: Join processid → taxid; write taxid_map
# ==============================================================================

message("\nBuilding taxid map...")

seq_taxids <- header_df %>%
  left_join(bin_tax %>% select(bin_uri, taxid), by = "bin_uri") %>%
  filter(!is.na(taxid))

message("  Sequences with taxid: ", nrow(seq_taxids), " / ", nrow(header_df))

# makeblastdb -taxid_map expects: seqid<tab>taxid, no header
write_tsv(seq_taxids %>% select(processid, taxid), TAXID_MAP, col_names = FALSE)
message("  Taxid map written: ", TAXID_MAP)

# ==============================================================================
# Section 7: Filter FASTA to sequences with taxids (streaming via python3)
# ==============================================================================

message("\nFiltering FASTA to ", nrow(seq_taxids), " sequences with taxids...")

keep_ids_file <- tempfile(fileext = ".txt")
writeLines(seq_taxids$processid, keep_ids_file)

py_script <- tempfile(fileext = ".py")
writeLines(c(
  "import gzip, sys",
  paste0("keep_file = '", keep_ids_file, "'"),
  paste0("fasta_in  = '", FASTA_GZ, "'"),
  paste0("fasta_out = '", DB_FASTA, "'"),
  "keep = set(open(keep_file).read().strip().split('\\n'))",
  "written = 0",
  "with gzip.open(fasta_in, 'rt') as fasta, open(fasta_out, 'w') as out:",
  "    keep_seq = False",
  "    for line in fasta:",
  "        if line.startswith('>'):",
  "            seqid = line[1:].split()[0]",
  "            keep_seq = seqid in keep",
  "            if keep_seq: written += 1",
  "        if keep_seq:",
  "            out.write(line)",
  "print(f'Sequences written: {written}', flush=True)"
), py_script)

py_result <- system2("python3", args = py_script, stdout = TRUE, stderr = TRUE)
message(paste(py_result, collapse = "\n"))
unlink(c(keep_ids_file, py_script))

if (!file.exists(DB_FASTA) || file.info(DB_FASTA)$size == 0) {
  stop("Filtered FASTA was not created or is empty: ", DB_FASTA)
}

# ==============================================================================
# Section 8: Build BLAST database with makeblastdb
# ==============================================================================

message("\nRunning makeblastdb...")

makeblastdb_args <- c(
  "-in",           DB_FASTA,
  "-input_type",   "fasta",
  "-dbtype",       "nucl",
  "-taxid_map",    TAXID_MAP,
  "-parse_seqids",
  "-blastdb_version", "5",
  "-title",        paste0(DB_NAME, ".", tax_suffix),
  "-out",          DB_OUT
)

db_result   <- system2("makeblastdb", args = makeblastdb_args, stdout = TRUE, stderr = TRUE)
db_exitcode <- attr(db_result, "status")

message(paste(db_result, collapse = "\n"))

if (!is.null(db_exitcode) && db_exitcode != 0) {
  stop("makeblastdb failed (exit ", db_exitcode, ")")
}

message("\nBOLD COI BLAST database created successfully.")
message("Database prefix: ", DB_OUT)
message("Use this path as BLAST_DB in pipeline_code_5.0.R")
