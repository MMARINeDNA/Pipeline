ceg_blast_submit <- function(
    PATH_TO_FASTA = NULL,
    PATH_TO_BLAST_TEMPLATE = "code/MURIblast_template.sh",     # local
    PATH_FOR_CREATED_SHELL_SCRIPT = "code_etc/MURIblast.sh",    # local (edited template)
    PATH_FOR_RESULTS = "blast_results/MURI_blast_results.txt",  # local desired path
    CUL_N = "100",
    SSH_SERVER = "KellyCEG@frustule.ocean.washington.edu:3004",
    NEGATIVE_TAXIDS_FILE = NULL,  # local path (optional)
    
    # ---------- Slurm controls ----------
    SLURM_TIME = "02:00:00",
    SLURM_CPUS = 16,
    SLURM_MEM  = "16G",
    SLURM_JOB_NAME = "blastn",
    SLURM_PARTITION = NULL,       # e.g., "main" / "short" / NULL
    SLURM_LOG_DIR = "slurm-logs", # remote (under $HOME)
    
    # ---------- Metadata file (local) ----------
    JOB_META_FILE = file.path(dirname(PATH_FOR_CREATED_SHELL_SCRIPT), "slurm_job_meta.txt")
) {
  suppressMessages(require(ssh))
  
  add_ssh_key <- function() system("eval $(ssh-agent -s) && ssh-add ~/.ssh/ceg_rsa")
  sh_quote <- function(x) paste0("'", gsub("'", "'\"'\"'", x, fixed = TRUE), "'")
  rbasename <- function(x) basename(normalizePath(x, mustWork = FALSE))
  
  add_ssh_key()
  session <- ssh_connect(SSH_SERVER)
  on.exit({ try(ssh_disconnect(session), silent = TRUE) }, add = TRUE)
  
  # Ensure remote dirs
  ssh_exec_wait(session, command = c(
    "mkdir -p tmp/raw tmp/processed tmp/scripts",
    paste("mkdir -p", sh_quote(SLURM_LOG_DIR))
  ))
  
  # Optional negative taxids
  negative_taxidlist_path_remote <- ""
  if (!is.null(NEGATIVE_TAXIDS_FILE) && file.exists(NEGATIVE_TAXIDS_FILE)) {
    scp_upload(session, files = NEGATIVE_TAXIDS_FILE, to = "tmp/raw")
    negative_taxidlist_path_remote <- paste0("/mnt/nfs/home/KellyCEG/tmp/raw/", rbasename(NEGATIVE_TAXIDS_FILE))
  }
  
  # Edit your template (keep structure)
  f <- readLines(PATH_TO_BLAST_TEMPLATE, warn = FALSE)
  f <- sub("^QUERY_FASTA=.*",
           paste0("QUERY_FASTA='/mnt/nfs/home/KellyCEG/tmp/raw/", rbasename(PATH_TO_FASTA), "'"), f)
  f <- sub("^CULLING=.*", paste0("CULLING=\"", CUL_N, "\""), f)
  f <- sub("^BLAST_OUTPUT=.*",
           paste0("BLAST_OUTPUT=\"/mnt/nfs/home/KellyCEG/tmp/processed/", rbasename(PATH_FOR_RESULTS), "\""), f)
  f <- sub("^NEGATIVE_TAXIDLIST=.*",
           paste0("NEGATIVE_TAXIDLIST=\"", negative_taxidlist_path_remote, "\""), f)
  # Make num_threads follow Slurm allocation even if hard-coded in template
  f <- gsub("-num_threads\\s+\\d+",
            "-num_threads $SLURM_CPUS_PER_TASK", f, perl = TRUE)
  
  writeLines(f, con = PATH_FOR_CREATED_SHELL_SCRIPT)
  
  message("Uploading FASTA and script to CEG server…")
  scp_upload(session, files = PATH_TO_FASTA, to = "tmp/raw")
  scp_upload(session, files = PATH_FOR_CREATED_SHELL_SCRIPT, to = "tmp/scripts")
  
  # Slurm submit wrapper
  submit_name   <- paste0("submit_", rbasename(PATH_FOR_CREATED_SHELL_SCRIPT))
  submit_remote <- file.path("tmp/scripts", submit_name)
  submit_local  <- tempfile(pattern = "sbatch_", fileext = ".sh")
  
  sb_lines <- c(
    "#!/bin/bash",
    paste0("#SBATCH --job-name=", SLURM_JOB_NAME),
    paste0("#SBATCH -o ", SLURM_LOG_DIR, "/%x-%j.out"),
    paste0("#SBATCH --cpus-per-task=", SLURM_CPUS),
    paste0("#SBATCH --time=", SLURM_TIME),
    paste0("#SBATCH --mem=", SLURM_MEM)
  )
  if (!is.null(SLURM_PARTITION) && nzchar(SLURM_PARTITION)) {
    sb_lines <- c(sb_lines, paste0("#SBATCH --partition=", SLURM_PARTITION))
  }
  sb_lines <- c(
    sb_lines,
    "",
    "set -euo pipefail",
    "echo \"$(date): starting $SLURM_JOB_NAME on $HOSTNAME with $SLURM_CPUS_PER_TASK cores; job $SLURM_JOB_ID\"",
    paste0("bash ", sh_quote(file.path("tmp/scripts", rbasename(PATH_FOR_CREATED_SHELL_SCRIPT)))),
    "rc=$?",
    "echo \"$(date): finished $SLURM_JOB_NAME; exit code $rc\"",
    "exit $rc"
  )
  writeLines(sb_lines, submit_local)
  scp_upload(session, files = submit_local, to = submit_remote)
  ssh_exec_wait(session, command = c(paste0("chmod +x ", sh_quote(submit_remote))))
  
  message("Submitting Slurm job…")
  sb_out <- ssh_exec_internal(session, command = c(paste("sbatch", sh_quote(submit_remote))))
  if (sb_out$status != 0) stop("Failed to submit sbatch job: ", rawToChar(sb_out$stderr))
  sbatch_msg <- rawToChar(sb_out$stdout)
  job_id <- sub(".*Submitted batch job ([0-9]+).*", "\\1", sbatch_msg)
  if (!grepl("^[0-9]+$", job_id)) stop("Could not parse Slurm job ID from sbatch output: ", sbatch_msg)
  
  # Write local meta to find the job later without objects
  remote_result <- paste0("tmp/processed/", rbasename(PATH_FOR_RESULTS))
  remote_log    <- file.path(SLURM_LOG_DIR, paste0(SLURM_JOB_NAME, "-", job_id, ".out"))
  
  meta <- c(
    paste0("JOB_ID=", job_id),
    paste0("REMOTE_RESULT=", remote_result),
    paste0("REMOTE_LOG=", remote_log),
    paste0("LOCAL_RESULT_DIR=", dirname(PATH_FOR_RESULTS)),
    paste0("LOCAL_RESULT_BASENAME=", rbasename(PATH_FOR_RESULTS)),
    paste0("SLURM_JOB_NAME=", SLURM_JOB_NAME),
    paste0("SLURM_LOG_DIR=", SLURM_LOG_DIR),
    paste0("SSH_SERVER=", SSH_SERVER)
  )
  dir.create(dirname(JOB_META_FILE), showWarnings = FALSE, recursive = TRUE)
  writeLines(meta, JOB_META_FILE)
  
  message("Slurm job ID: ", job_id)
  message("Job metadata written to: ", JOB_META_FILE)
  message("To check later, just call check_ceg_blast_quick(...) with the same run inputs.")
  invisible(job_id)
}

check_ceg_blast <- function(
    PROCESSED_LOCATION,
    RUN_NAME,
    PRIMERNAME,
    META_BASENAME = "slurm_job_meta.txt",
    print_command = TRUE
) {
  suppressMessages(require(ssh))
  
  add_ssh_key <- function() system("eval $(ssh-agent -s) && ssh-add ~/.ssh/ceg_rsa")
  sh_quote <- function(x) paste0("'", gsub("'", "'\"'\"'", x, fixed = TRUE), "'")
  add_ssh_key()
  
  # ---- read meta ----
  meta_file <- file.path(PROCESSED_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME), "code_etc", META_BASENAME)
  if (!file.exists(meta_file)) {
    message("No job metadata found at: ", meta_file)
    return(invisible(FALSE))
  }
  kv <- readLines(meta_file, warn = FALSE)
  kv <- kv[nzchar(kv)]
  parse_kv <- function(x) { p <- strsplit(x, "=", fixed = TRUE)[[1]]; setNames(list(paste(p[-1], collapse="=")), p[1]) }
  meta <- do.call(c, lapply(kv, parse_kv))
  
  need <- c("JOB_ID","REMOTE_RESULT","REMOTE_LOG","LOCAL_RESULT_DIR","LOCAL_RESULT_BASENAME","SSH_SERVER")
  miss <- setdiff(need, names(meta))
  if (length(miss)) {
    message("Meta file missing entries: ", paste(miss, collapse = ", "))
    return(invisible(FALSE))
  }
  
  job_id     <- meta[["JOB_ID"]]
  remote_res <- meta[["REMOTE_RESULT"]]    # may be relative
  remote_log <- meta[["REMOTE_LOG"]]       # may be relative
  local_dir  <- meta[["LOCAL_RESULT_DIR"]]
  local_name <- meta[["LOCAL_RESULT_BASENAME"]]
  ssh_server <- meta[["SSH_SERVER"]]
  
  # ---- connect (agent, then keyfile fallback) ----
  session <- tryCatch(ssh_connect(ssh_server),
                      error = function(e) try(ssh_connect(ssh_server, keyfile = "~/.ssh/ceg_rsa"), silent = TRUE))
  if (inherits(session, "try-error")) stop("Authentication with ssh server failed (agent + keyfile).")
  on.exit({ try(ssh_disconnect(session), silent = TRUE) }, add = TRUE)
  
  # Discover remote $HOME (for absolute’ing relative paths)
  rem_home <- rawToChar(ssh_exec_internal(session, command = "printf $HOME")$stdout)
  rem_home <- trimws(rem_home)
  
  # ---- resolve authoritative stdout path via scontrol ----
  sc <- ssh_exec_internal(session, command = paste("scontrol show job", job_id, "-o"))
  if (sc$status == 0) {
    line <- rawToChar(sc$stdout)
    # StdOut=..., WorkDir=...
    m_out <- regexpr("StdOut=([^ ]+)", line)
    if (m_out > 0) {
      path <- sub("^StdOut=", "", regmatches(line, m_out))
      path <- trimws(path)
      if (nchar(path)) remote_log <- path
    }
    m_wd <- regexpr("WorkDir=([^ ]+)", line)
    workdir <- if (m_wd > 0) trimws(sub("^WorkDir=", "", regmatches(line, m_wd))) else ""
    # If remote_log is relative, anchor it to WorkDir or $HOME
    if (!grepl("^/", remote_log)) {
      remote_log <- if (nzchar(workdir)) file.path(workdir, remote_log) else file.path(rem_home, remote_log)
    }
  } else {
    # fallback: if relative, make absolute under $HOME
    if (!grepl("^/", remote_log)) remote_log <- file.path(rem_home, remote_log)
  }
  
  # ---- status via sacct (preferred), else squeue ----
  state <- NA_character_; exitc <- NA_character_; maxrs <- NA_character_
  acct <- ssh_exec_internal(session, command = paste("sacct -j", job_id, "--format=JobID,State,ExitCode,MaxRSS -P -n"))
  if (acct$status == 0) {
    lines <- strsplit(rawToChar(acct$stdout), "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(lines)]
    main  <- grep(paste0("^", job_id, "\\|"), lines, value = TRUE)
    if (length(main)) {
      parts <- strsplit(main[1], "\\|", perl = TRUE)[[1]]
      state <- parts[2]; exitc <- parts[3]; if (length(parts) >= 4) maxrs <- parts[4]
    }
  }
  if (is.na(state)) {
    q <- ssh_exec_internal(session, command = paste("squeue -j", job_id, "-h -o %T"))
    qstate <- trimws(rawToChar(q$stdout))
    state <- if (nchar(qstate) > 0) qstate else "COMPLETED?"
    exitc <- if (nchar(qstate) > 0) "0:0" else exitc
  }
  
  message("Job ", job_id, " state: ", state, if (!is.na(exitc)) paste0(" (ExitCode: ", exitc, ")") else "")
  message("Log path: ", remote_log)
  
  # ---- print the exact blastn command from the slurm log (RUNNING or COMPLETED) ----
  if (print_command) {
    # Pull up to 5 matching lines to be safe; no file download
    cmd <- paste("grep -E", sh_quote("\\bblastn[[:space:]]"), "-m 5", sh_quote(remote_log))
    g <- ssh_exec_internal(session, command = cmd)
    if (g$status == 0 && nchar(rawToChar(g$stdout)) > 0) {
      cat("\n=== BLAST COMMAND (from slurm log) ===\n",
          rawToChar(g$stdout),
          "======================================\n", sep = "")
    } else {
      # fallback: try to print the definition line if present
      alt <- ssh_exec_internal(session, command = paste("grep -m 1 -n", sh_quote("BLAST_CMD="), sh_quote(remote_log)))
      if (alt$status == 0 && nchar(rawToChar(alt$stdout)) > 0) {
        cat("\n=== BLAST_CMD line (from slurm log) ===\n",
            rawToChar(alt$stdout),
            "=======================================\n", sep = "")
      } else {
        message("Blast command not found in log yet (log not flushed or echo not reached).")
      }
    }
  }
  
  # ---- on completion, download result (handle relative/absolute) ----
  if (grepl("^COMPLETED", state)) {
    if (!dir.exists(local_dir)) dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
    if (!grepl("^/", remote_res)) remote_res <- file.path(rem_home, remote_res)
    scp_download(session, files = remote_res, to = local_dir)  # 'to' must be a directory
    message("✅ Results downloaded to: ", file.path(local_dir, local_name))
    return(invisible(TRUE))
  } else {
    message("Still running or pending; re-run this check later.")
    return(invisible(FALSE))
  }
}
