local_blast <- function(
    PATH_TO_FASTA,
    PATH_FOR_RESULTS,
    BLAST_DB              = "/Volumes/Clupea/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST/MIDORI2_UNIQ_NUC_GB270_CO1_BLAST",
    PATH_FOR_SHELL_SCRIPT = NULL,      # if NULL, sits next to PATH_FOR_RESULTS
    NEGATIVE_TAXIDS_FILE  = NULL,
    CUL_N                 = "50",
    PERCENT_IDENTITY      = "90",
    WORD_SIZE             = "30",
    EVALUE                = "1e-30",
    MAX_TARGET_SEQS       = "50",
    NUM_THREADS           = 8
) {
  if (!file.exists(PATH_TO_FASTA)) stop("Query FASTA not found: ", PATH_TO_FASTA)
  if (!dir.exists(dirname(BLAST_DB)))  stop("BLAST_DB directory not found: ", dirname(BLAST_DB))

  # Optional negative taxid list: strip comments and blank lines, write clean temp file
  neg_arg <- ""
  if (!is.null(NEGATIVE_TAXIDS_FILE) && file.exists(NEGATIVE_TAXIDS_FILE)) {
    raw       <- readLines(NEGATIVE_TAXIDS_FILE, warn = FALSE)
    clean_ids <- sub("\\s.*$", "", trimws(raw))           # drop everything after first whitespace
    clean_ids <- clean_ids[grepl("^[0-9]+$", clean_ids)] # keep only bare integers
    if (length(clean_ids) > 0) {
      clean_taxid_file <- tempfile(fileext = ".txt")
      writeLines(clean_ids, clean_taxid_file)
      neg_arg <- sprintf("  -negative_taxidlist '%s' \\\n", clean_taxid_file)
    }
  }

  # Shell script location
  if (is.null(PATH_FOR_SHELL_SCRIPT)) {
    PATH_FOR_SHELL_SCRIPT <- sub("\\.txt$", "_blast.sh", PATH_FOR_RESULTS)
  }

  sh_content <- sprintf(paste0(
    "#!/bin/bash\n",
    "cd '%s'\n\n",
    "blastn \\\n",
    "  -query '%s' \\\n",
    "  -db '%s' \\\n",
    "  -perc_identity %s \\\n",
    "  -word_size %s \\\n",
    "  -evalue %s \\\n",
    "  -max_target_seqs %s \\\n",
    "  -culling_limit %s \\\n",
    "%s",
    "  -num_threads %d \\\n",
    "  -outfmt '6 sscinames scomnames qseqid sseqid pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs' \\\n",
    "  -out '%s'\n"
  ),
    dirname(BLAST_DB),
    PATH_TO_FASTA, BLAST_DB,
    PERCENT_IDENTITY, WORD_SIZE, EVALUE,
    MAX_TARGET_SEQS, CUL_N,
    neg_arg,
    NUM_THREADS,
    PATH_FOR_RESULTS
  )

  writeLines(sh_content, PATH_FOR_SHELL_SCRIPT)
  Sys.chmod(PATH_FOR_SHELL_SCRIPT, mode = "0755")

  n_seqs <- length(grep("^>", readLines(PATH_TO_FASTA, warn = FALSE)))
  message("Running local BLAST (", n_seqs, " sequences)... ", Sys.time())

  result    <- system2("bash", args = shQuote(PATH_FOR_SHELL_SCRIPT), stdout = TRUE, stderr = TRUE)
  exit_code <- attr(result, "status")

  if (!is.null(exit_code) && exit_code != 0) {
    stop("BLAST failed (exit ", exit_code, "):\n", paste(result, collapse = "\n"))
  }
  if (length(result) > 0) {
    if (any(grepl("requires additional data files", result, fixed = TRUE)) && neg_arg != "") {
      warning(
        "BLAST database lacks taxonomy index — -negative_taxidlist was ignored.\n",
        "  Taxid filtering requires a database built with makeblastdb -taxid_map.\n",
        "  Set NEGATIVE_TAXIDS_FILE = NULL in the pipeline config to suppress this warning."
      )
    }
    message(paste(result, collapse = "\n"))
  }

  if (!file.exists(PATH_FOR_RESULTS) || file.info(PATH_FOR_RESULTS)$size == 0) {
    warning("BLAST output missing or empty: ", PATH_FOR_RESULTS)
  } else {
    message("BLAST complete: ", PATH_FOR_RESULTS)
  }
}
