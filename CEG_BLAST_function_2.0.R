ceg_blast <- function(
    PATH_TO_FASTA = NULL,
    PATH_TO_BLAST_TEMPLATE = "code/MURIblast_template.sh", # on local drive
    PATH_FOR_CREATED_SHELL_SCRIPT = "code_etc/MURIblast.sh", # on local drive
    PATH_FOR_RESULTS = "blast_results/MURI_blast_results.txt", # on local drive
    CUL_N = "50",
    SSH_SERVER = "KellyCEG@frustule.ocean.washington.edu:3004",
    NEGATIVE_TAXIDS_FILE = NULL  # Path to the negative taxid file
) {
  
  suppressMessages(require(ssh))
  suppressMessages(require(here))
  
  # Function to add SSH key
  add_ssh_key <- function() {
    system("eval $(ssh-agent -s) && ssh-add ~/.ssh/ceg_rsa")
  }
  
  # Add SSH key
  add_ssh_key()
  
  # Initialize variables
  negative_taxidlist_path_remote <- ""
  
  # Connect to SSH server
  session <- ssh_connect(SSH_SERVER)
  
  # Read the template BLAST command to assign files and params
  f <- readLines(PATH_TO_BLAST_TEMPLATE, warn = FALSE)
  
  # Modify the shell script to assign variables
  # Replace QUERY_FASTA
  f <- sub("^QUERY_FASTA=.*", paste0("QUERY_FASTA='/mnt/nfs/home/KellyCEG/tmp/raw/", basename(PATH_TO_FASTA), "'"), f)
  # Replace CULLING
  f <- sub("^CULLING=.*", paste0("CULLING=\"", CUL_N, "\""), f)
  # Replace BLAST_OUTPUT
  f <- sub("^BLAST_OUTPUT=.*", paste0("BLAST_OUTPUT=\"/mnt/nfs/home/KellyCEG/tmp/processed/", basename(PATH_FOR_RESULTS), "\""), f)
  
  # Set NEGATIVE_TAXIDLIST variable
  f <- sub("^NEGATIVE_TAXIDLIST=.*", paste0("NEGATIVE_TAXIDLIST=\"", negative_taxidlist_path_remote, "\""), f)
  
  # Write the modified shell script
  writeLines(f, con = PATH_FOR_CREATED_SHELL_SCRIPT)
  
  print("Uploading Files to CEG Server")
  
  # Upload the FASTA file
  scp_upload(session,
             files = PATH_TO_FASTA,
             to = "tmp/raw")
  
  # Process NEGATIVE_TAXIDS_FILE
  if (!is.null(NEGATIVE_TAXIDS_FILE) && file.exists(NEGATIVE_TAXIDS_FILE)) {
    # Upload the negative taxid file to the server
    scp_upload(session, 
               files = NEGATIVE_TAXIDS_FILE,
               to = "tmp/raw")
    
    # Set the path to the negative taxid file on the remote server
    negative_taxidlist_path_remote <- paste0("/mnt/nfs/home/KellyCEG/tmp/raw/", basename(NEGATIVE_TAXIDS_FILE))
  }
  
  # Upload the shell script
  scp_upload(session,
             files = PATH_FOR_CREATED_SHELL_SCRIPT,
             to = "tmp/scripts")
  
  print("Blasting on CEG Server")
  
  # Execute the shell script on the server
  ssh_exec_wait(session, command = c(paste0("sh tmp/scripts/", basename(PATH_FOR_CREATED_SHELL_SCRIPT))))
  
  # Download the results
  scp_download(session,
               files = paste0("tmp/processed/", basename(PATH_FOR_RESULTS)),
               to = dirname(PATH_FOR_RESULTS))
  
  # Disconnect from SSH session
  ssh_disconnect(session)
  
  print("Completed")
}
