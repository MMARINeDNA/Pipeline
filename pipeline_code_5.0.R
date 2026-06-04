##############################################################################
##                         Welcome to The Kelly Lab's                       ##
##              Metabarcoding Taxonomic Assignment Pipeline                 ##
## This pipeline will take you from de-multiplexed fastq files to a matrix ##
## with taxon names and counts. Works on Linux or Mac (PC: use WSL2).      ##
##############################################################################

here::i_am("pipeline_code_5.0.R")

suppressMessages({
  library(tidyverse)
  library(dada2)
  library(digest)
  library(seqinr)
  library(sys)
  library(ShortRead)
  library(here)
})

source(here("local_blast.R"))
source(here("LCA_function_4.0.R"))

## ============================================================================
## SECTION 1: CONFIGURATION — edit only this section for each run
## ============================================================================

# Run identity
RUN_NAME   <- "20260530_WOAC_COI"
PRIMERNAME <- "LRY"   # must match a row name in primer.data.csv
SP_THOLD   <- 98      # % identity floor for species-level assignment

# Directories
PARENT_DIR    <- "/Users/rpk/Downloads/20260530_COI"  # contains Fastq/ subfolder with raw reads
PROCESSED_DIR <- "/Users/rpk/Desktop"                   # final run folder is created here
DB_DIR        <- "/Users/rpk/Desktop/lca_dbs"           # per-primer _database.csv files live here

# Local BLAST database
BLAST_DB <- "/Volumes/Clupea/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST"

# Local tool paths (use just the command name if it is on PATH)
CUTADAPT      <- "/Library/Frameworks/Python.framework/Versions/3.10/bin/cutadapt"
TAXONKIT_PATH <- "/Users/rpk/taxonkit"

# Optional: path to taxids-to-exclude file; set to NULL to disable
NEGATIVE_TAXIDS_FILE <- here("MURI_taxids_to_exclude.txt")

# Set TRUE to skip DADA2 and reload a previous run's ASV table for re-annotation
SKIP_DADA2 <- FALSE

# DADA2 size filter overrides — set to NA to use values from primer.data.csv
MAX_AMPLICON_OVERRIDE <- 350
MIN_AMPLICON_OVERRIDE <- 200
OVERLAP_OVERRIDE      <- 20

# Set TRUE to stop after DADA2 (skips BLAST/LCA — useful for testing)
STOP_AFTER_DADA2 <- FALSE

## ============================================================================
## SECTION 2: SETUP — all paths derived here; no path construction downstream
## ============================================================================

# Temporary working directories under PARENT_DIR (large, removed on cleanup)
FASTQ_DIR    <- file.path(PARENT_DIR, "Fastq")
TRIMMED_DIR  <- file.path(PARENT_DIR, "for_dada2")
FILTERED_DIR <- file.path(PARENT_DIR, "filtered")

# Permanent output directories written directly to the final run folder
RUN_DIR    <- file.path(PROCESSED_DIR, paste(RUN_NAME, PRIMERNAME, sep = "_"))
CODE_DIR   <- file.path(RUN_DIR, "code_etc")
OUTPUT_DIR <- file.path(RUN_DIR, "outputs")
LOG_DIR    <- file.path(RUN_DIR, "logs")

# Key file paths referenced throughout the script
RUN_PREFIX       <- paste(RUN_NAME, PRIMERNAME, sep = "_")
HASH_KEY_CSV     <- file.path(OUTPUT_DIR, paste0(RUN_PREFIX, "_hash_key.csv"))
HASH_KEY_FASTA   <- file.path(OUTPUT_DIR, paste0(RUN_PREFIX, "_hash_key.fasta"))
ASV_CSV          <- file.path(OUTPUT_DIR, paste0(RUN_PREFIX, "_ASV_table.csv"))
SEQS_TO_ANNOTATE <- file.path(OUTPUT_DIR, "seqs_to_annotate.fasta")
BLAST_RESULTS    <- file.path(OUTPUT_DIR, "new_annotations.txt")
DATABASE_FILE    <- file.path(DB_DIR, paste0(PRIMERNAME, "_database.csv"))

# Create all directories (silently if they already exist)
for (d in c(TRIMMED_DIR, FILTERED_DIR, CODE_DIR, OUTPUT_DIR, LOG_DIR, DB_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Verify external tool availability
if (system2(CUTADAPT, "--version", stdout = FALSE, stderr = FALSE) != 0) {
  stop("cutadapt not found at: ", CUTADAPT)
}
if (!file.exists(TAXONKIT_PATH)) {
  stop("taxonkit not found at: ", TAXONKIT_PATH)
}

# Load and validate primer metadata
primer.data <- read.csv(here("primer.data.csv"))
primer.row  <- filter(primer.data, name == PRIMERNAME)
if (nrow(primer.row) == 0) stop("PRIMERNAME '", PRIMERNAME, "' not found in primer.data.csv")
if (nrow(primer.row) > 1) stop("Multiple rows match PRIMERNAME '", PRIMERNAME, "' in primer.data.csv")

PRIMERSEQ_F         <- primer.row$seq_f
PRIMERSEQ_R         <- primer.row$seq_r
MAX_AMPLICON_LENGTH <- if (!is.na(MAX_AMPLICON_OVERRIDE)) MAX_AMPLICON_OVERRIDE else primer.row$max_amplicon_length
MIN_AMPLICON_LENGTH <- if (!is.na(MIN_AMPLICON_OVERRIDE)) MIN_AMPLICON_OVERRIDE else primer.row$min_amplicon_length
OVERLAP             <- if (!is.na(OVERLAP_OVERRIDE))      OVERLAP_OVERRIDE      else primer.row$overlap

# Archive pipeline scripts for reproducibility
file.copy(here("primer.data.csv"),       CODE_DIR, overwrite = TRUE)
file.copy(here("pipeline_code_5.0.R"),   CODE_DIR, overwrite = TRUE)
file.copy(here("local_blast.R"),         CODE_DIR, overwrite = TRUE)
file.copy(here("LCA_function_4.0.R"),    CODE_DIR, overwrite = TRUE)
if (!is.null(NEGATIVE_TAXIDS_FILE) && file.exists(NEGATIVE_TAXIDS_FILE)) {
  file.copy(NEGATIVE_TAXIDS_FILE, CODE_DIR, overwrite = TRUE)
}

## ============================================================================
## SECTION 3: PRIMER TRIMMING
## ============================================================================

# Helper: extract sample name from a file path.
# Strips the primer prefix (if present, with - or _ separator), then the
# Illumina standard suffix (_S{n}_L{lane}_R{1/2}...), falling back to
# stripping from _R1/_R2 for non-Illumina filenames.
sample_name_from_file <- function(paths, primername) {
  x <- basename(paths)
  x <- str_remove(x, paste0("^", primername, "[-_]?"))  # primer prefix (optional)
  x <- str_remove(x, "_S\\d+_L\\d+_R[12].*$")           # Illumina standard suffix
  x <- str_remove(x, "_R[12].*$")                        # fallback for non-Illumina naming
  x
}

# Write the primer-trimming shell script inline (no dependency on make_primer_shell_script.R).
# Uses find instead of ls so the primer-separator character (- or _) doesn't matter.
# R2 is derived by replacing only the last _R1 occurrence to avoid corrupting sample names.
trim_script <- file.path(CODE_DIR, "trim_primers.sh")
writeLines(paste0(
'#!/bin/bash
set -euo pipefail

PRIMERNAME=', PRIMERNAME, '
PRIMERSEQ_F=', PRIMERSEQ_F, '
PRIMERSEQ_R=', PRIMERSEQ_R, '
FASTQ_DIR=', shQuote(FASTQ_DIR), '
TRIMMED_DIR=', shQuote(TRIMMED_DIR), '
LOG_DIR=', shQuote(LOG_DIR), '
CUTADAPT=', shQuote(CUTADAPT), '

FILES=$(find "${FASTQ_DIR}" -maxdepth 1 -name "*_R1*.fastq.gz" | sort)
if [ -z "$FILES" ]; then
    echo "ERROR: No *_R1*.fastq.gz files found in ${FASTQ_DIR}" >&2
    exit 1
fi

for R1 in $FILES; do
    # Derive R2 by replacing only the final _R1 in the filename, not the path
    DIRNAME=$(dirname "$R1")
    BASENAME=$(basename "$R1")
    R2_BASENAME="${BASENAME%_R1*}_R2${BASENAME#*_R1}"
    R2="${DIRNAME}/${R2_BASENAME}"
    if [ ! -f "$R2" ]; then
        echo "WARNING: No R2 file found for $BASENAME — skipping" >&2
        continue
    fi
    echo "Trimming: $BASENAME"
    ${CUTADAPT} \\
        -g ${PRIMERSEQ_F} \\
        -G ${PRIMERSEQ_R} \\
        -o "${TRIMMED_DIR}/${BASENAME}" \\
        -p "${TRIMMED_DIR}/${R2_BASENAME}" \\
        --discard-untrimmed \\
        -j 0 \\
        "$R1" "$R2" >> "${LOG_DIR}/cutadapt_trim_report.txt"
done
'), trim_script)
system2("chmod", c("+x", shQuote(trim_script)))
system2("sh", shQuote(trim_script))

## ============================================================================
## SECTION 4: DADA2
## ============================================================================

if (SKIP_DADA2) {

  message("SKIP_DADA2 = TRUE: loading existing ASV table and hash key from ", OUTPUT_DIR)
  current_asv  <- read_csv(ASV_CSV,      show_col_types = FALSE)
  conv_table   <- read_csv(HASH_KEY_CSV, show_col_types = FALSE)
  sample.names <- sort(unique(current_asv$Sample_name))

} else {

  ### Discover trimmed files --------------------------------------------------
  fnFs <- sort(list.files(TRIMMED_DIR, pattern = "_R1.*\\.fastq\\.gz$", full.names = TRUE))
  fnRs <- sort(list.files(TRIMMED_DIR, pattern = "_R2.*\\.fastq\\.gz$", full.names = TRUE))

  if (length(fnFs) == 0) stop("No trimmed R1 files found in ", TRIMMED_DIR)
  if (length(fnFs) != length(fnRs)) {
    stop("Unequal number of R1 (", length(fnFs), ") and R2 (", length(fnRs), ") files in ", TRIMMED_DIR)
  }

  sample.names <- sample_name_from_file(fnFs, PRIMERNAME)

  if (anyDuplicated(sample.names)) {
    stop("Duplicate sample names after parsing filenames: ",
         paste(unique(sample.names[duplicated(sample.names)]), collapse = ", "), "\n",
         "Check that PRIMERNAME is correct and filenames follow <PRIMER>_<sample>_R1_*.fastq.gz")
  }

  ### Set up filtered file paths (full paths throughout — no setwd needed) -----
  filtFs <- file.path(FILTERED_DIR, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(FILTERED_DIR, paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  ### Optional: inspect quality before committing to truncLen
  # plotQualityProfile(fnFs[1:4])
  # plotQualityProfile(fnRs[1:4])

  ### Quality filtering -------------------------------------------------------
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       truncLen   = MIN_AMPLICON_LENGTH,
                       maxN       = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                       compress   = TRUE, multithread = FALSE, matchIDs = TRUE)

  # Drop samples that produced no output (keeps sample.names complete for tracking)
  exists <- file.exists(filtFs) & file.exists(filtRs)
  filtFs <- filtFs[exists]
  filtRs <- filtRs[exists]

  ### Error model + denoising -------------------------------------------------
  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)

  dadaFs <- dada(filtFs, err = errF, selfConsist = TRUE, multithread = TRUE, MAX_CONSIST = 20)
  dadaRs <- dada(filtRs, err = errR, selfConsist = TRUE, multithread = TRUE, MAX_CONSIST = 20)

  ### Merge + chimera removal -------------------------------------------------
  mergers       <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                              minOverlap = OVERLAP, verbose = TRUE, trimOverhang = TRUE)
  seqtab        <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = TRUE, verbose = TRUE)

  ### Size filter -------------------------------------------------------------
  keep   <- nchar(colnames(seqtab.nochim)) >= MIN_AMPLICON_LENGTH &
            nchar(colnames(seqtab.nochim)) <= MAX_AMPLICON_LENGTH
  cleaned.seqtab <- seqtab.nochim[, keep,  drop = FALSE]
  write.csv(seqtab.nochim[, !keep, drop = FALSE], file.path(LOG_DIR, "filtered_out_asv.csv"))

  ### Read tracking -----------------------------------------------------------
  getN <- function(x) sum(getUniques(x))
  track <- cbind(
    out,
    denoisedF = sapply(sample.names, \(s) if (s %in% names(dadaFs))        getN(dadaFs[[s]])         else 0L),
    denoisedR = sapply(sample.names, \(s) if (s %in% names(dadaRs))        getN(dadaRs[[s]])         else 0L),
    merged    = sapply(sample.names, \(s) if (s %in% names(mergers))       getN(mergers[[s]])        else 0L),
    nonchim   = sapply(sample.names, \(s) if (s %in% rownames(seqtab.nochim)) rowSums(seqtab.nochim)[s] else 0L)
  )
  rownames(track) <- sample.names
  write.csv(track, file.path(LOG_DIR, "tracking_reads.csv"))

  ### Hashing -----------------------------------------------------------------
  message("Creating ASV table and hash key... ", Sys.time())
  seqtab.df  <- as.data.frame(cleaned.seqtab)
  conv_table <- tibble(
    Hash     = map_chr(colnames(seqtab.df),
                       ~ digest(.x, algo = "sha1", serialize = FALSE, skip = "auto")),
    Sequence = colnames(seqtab.df)
  )

  write_csv(conv_table, HASH_KEY_CSV)
  write.fasta(as.list(conv_table$Sequence), as.list(conv_table$Hash), file.out = HASH_KEY_FASTA)

  current_asv <- rownames_to_column(seqtab.df, "Sample_name") |>
    mutate(Label = PRIMERNAME, .after = "Sample_name") |>
    pivot_longer(cols = -c(Sample_name, Label), names_to = "Sequence", values_to = "nReads") |>
    filter(nReads > 0) |>
    left_join(conv_table, by = "Sequence") |>
    select(-Sequence) |>
    relocate(Hash, .after = Label)

  write_csv(current_asv, ASV_CSV)

} # end SKIP_DADA2 block

## ============================================================================
## SECTION 5: ANNOTATION (BLAST + LCA)
## ============================================================================

if (isTRUE(STOP_AFTER_DADA2)) {
  message("STOP_AFTER_DADA2 = TRUE — pipeline stopped after hashing. Outputs: ", OUTPUT_DIR)
  quit("no")
}

# Determine which hashes need BLASTing
if (file.exists(DATABASE_FILE)) {
  db            <- read.csv(DATABASE_FILE, row.names = 1)
  seqs_to_blast <- filter(conv_table, !Hash %in% db$Hash)
} else {
  db            <- NULL
  seqs_to_blast <- conv_table
}

if (nrow(seqs_to_blast) == 0) {
  message("All hashes already in database — skipping BLAST")
} else {
  write.fasta(as.list(seqs_to_blast$Sequence), as.list(seqs_to_blast$Hash),
              file.out = SEQS_TO_ANNOTATE)

  local_blast(
    PATH_TO_FASTA         = SEQS_TO_ANNOTATE,
    PATH_FOR_RESULTS      = BLAST_RESULTS,
    BLAST_DB              = BLAST_DB,
    PATH_FOR_SHELL_SCRIPT = file.path(CODE_DIR, "blast_new_seqs.sh"),
    NEGATIVE_TAXIDS_FILE  = if (!is.null(NEGATIVE_TAXIDS_FILE))
                              file.path(CODE_DIR, basename(NEGATIVE_TAXIDS_FILE))
                            else NULL,
    CUL_N = "50"
  )

  if (!file.exists(BLAST_RESULTS) || file.size(BLAST_RESULTS) == 0L) {
    message("BLAST returned no results — no new annotations added")
  } else if (!is.null(db)) {
    # Append new annotations to existing database
    LCA(BLASTOUTPUT = BLAST_RESULTS,
        FASTA       = SEQS_TO_ANNOTATE,
        DB_PATH_IN  = DATABASE_FILE,
        DB_PATH_OUT = DATABASE_FILE,
        SP_THOLD    = SP_THOLD)
  } else {
    # Create the database from scratch
    db_new <- LCA(BLASTOUTPUT = BLAST_RESULTS,
                  FASTA       = SEQS_TO_ANNOTATE,
                  SP_THOLD    = SP_THOLD)
    write.csv(distinct(db_new), DATABASE_FILE)
  }

  db <- read.csv(DATABASE_FILE, row.names = 1)
}

if (is.null(db) || nrow(db) == 0) {
  stop("No annotation database available. BLAST may have returned no results and no prior ",
       "database exists at:\n  ", DATABASE_FILE)
}

## ============================================================================
## SECTION 6: OUTPUT TABLES
## ============================================================================

# All samples, including those that produced zero reads, appear as zero columns
full_samples <- sort(sample.names)

pad_missing_samples <- function(wide_tbl, all_samples) {
  missing <- setdiff(all_samples, colnames(wide_tbl))
  for (s in missing) wide_tbl[[s]] <- 0L
  wide_tbl
}

### Taxon table (long) --------------------------------------------------------
tax_table <- current_asv |>
  left_join(select(db, Hash, BestTaxon, Class), by = "Hash") |>
  group_by(Sample_name, BestTaxon, Class) |>
  summarise(nReads = sum(nReads), .groups = "drop")
write_csv(tax_table, file.path(OUTPUT_DIR, paste0(RUN_PREFIX, "_taxon_table.csv")))

### Taxon table (wide) --------------------------------------------------------
tax_wide <- tax_table |>
  pivot_wider(id_cols = c(BestTaxon, Class),
              names_from = Sample_name, values_from = nReads, values_fill = 0) |>
  pad_missing_samples(full_samples) |>
  select(BestTaxon, Class, all_of(full_samples))
write_csv(tax_wide, file.path(OUTPUT_DIR, paste0(RUN_PREFIX, "_taxon_table_wide.csv")))

### Haplotype table (wide) ----------------------------------------------------
haplotype_table <- current_asv |>
  left_join(db, by = "Hash") |>
  mutate(BestTaxon_Haplotype = paste(BestTaxon, HaplotypeNumber, sep = "_")) |>
  group_by(BestTaxon_Haplotype, Sample_name, Class) |>
  summarise(nReads = sum(nReads, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(id_cols = c(BestTaxon_Haplotype, Class),
              names_from = Sample_name, values_from = nReads, values_fill = 0) |>
  pad_missing_samples(full_samples) |>
  select(BestTaxon_Haplotype, Class, all_of(full_samples))
write.csv(haplotype_table, file.path(OUTPUT_DIR, paste0(RUN_PREFIX, "_haplotype_table.csv")),
          row.names = FALSE)

message("Pipeline complete. Outputs: ", OUTPUT_DIR)

## ============================================================================
## SECTION 7: CLEANUP — remove large intermediate directories from PARENT_DIR
## ============================================================================

CLEANUP <- TRUE
if (CLEANUP) {
  for (d in c(TRIMMED_DIR, FILTERED_DIR)) {
    unlink(d, recursive = TRUE)
  }
}
