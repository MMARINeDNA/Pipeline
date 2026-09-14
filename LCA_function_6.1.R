LCA <- function(
    BLASTOUTPUT = here("blast_results/testout.txt"), #path to tabular blast output in format "6 sscinames scomnames qseqid sseqid pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs"  
    FASTA = here("test_fastas/temp2.fasta"), #from which to get the actual sequences, because for some reason they aren't in the blast results table
    DB_PATH_IN = NULL, # here("results/MURI_MFU_combinedFish_Annotation.csv") #existing database to which to add these sequence annotations, if desired
    DB_PATH_OUT = NULL, #if desired, name of database file to write out; may be the same as `DB_PATH_IN`, if appending to existing db
    SP_MAX_THOLD = 96, #Threshold below which species assignments are ignored, and only higher taxonomic assignments are considered
    SP_MIN_THOLD = 98.7 #pident floor for the second tier of hits when there are no 100% matches
    ){
  
  suppressMessages(require(here))
  suppressMessages(require(tidyverse))
  suppressMessages(require(sys)) 
  suppressMessages(require(ShortRead)) 
  
  taxonkit_path <- TAXONKIT_PATH
  
  blastoutput <- read_delim(BLASTOUTPUT, delim = "\t", col_names = F, show_col_types = FALSE)
  fasta <- readFasta(FASTA)
  
  # Initialize `db` as an empty dataframe with correct columns if no database exists
  if (!is.null(DB_PATH_IN) && file.exists(DB_PATH_IN)){
    db <- read.csv(DB_PATH_IN, row.names = 1) %>% 
      mutate(
        BlastTaxIDs = as.character(BlastTaxIDs),  # Ensure BlastTaxIDs is character
        HaplotypeNumber = as.numeric(HaplotypeNumber)  # Ensure HaplotypeNumber is numeric
      )
  } else {
    db <- data.frame(
      Hash = character(), 
      BestTaxon = character(), 
      HaplotypeNumber = numeric(), 
      Tax_list = character(),
      Max_pident = character(), 
      Mismatches = character(), 
      BlastTaxIDs = character(),
      AccessionNumbers= character(), 
      LCA_taxid = character(), 
      Kingdom = character(), 
      Phylum = character(), 
      Class = character(), 
      Order = character(), 
      Family = character(), 
      Genus = character(), 
      Species = character(), 
      Sequence = character(), 
      DateAdded = character(), 
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(DB_PATH_OUT)){
    db_outfile <- DB_PATH_OUT
  }
      
  # Process blastoutput
  names(blastoutput) <- c("Taxon","CommonName",
                          "Hash","Accession",
                          "pident", "length", "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "qlen", "qcovs")
  
  ##Remove unconfirmed sequence hits
  blastoutput <- blastoutput %>%
    filter(
      !str_detect(Accession, "XM_") &
        !str_detect(Accession, "XP_") &
        !str_detect(Accession, "XR_")
    )
  
  result_list <- list()
  
  # Apply decision tree for each unique Hash
    for (hash_id in unique(blastoutput$Hash)) {
    hash_hits <- blastoutput %>%
      filter(Hash == hash_id)
    
    # Step 1: Check for pident == 100
    selected_hits <- hash_hits %>%
      filter(pident == 100)
    
    # Step 3: If no hits with pident == 100, check for pident > 98 (3 SNPs)
    if (nrow(selected_hits) == 0) {
      selected_hits <- hash_hits %>%
        filter(pident > SP_MIN_THOLD)
    }
    
    # Step 4: If no hits with pident > 98.5, check for pident > 96
    if (nrow(selected_hits) == 0) {
      selected_hits <- hash_hits %>%
        filter(pident > 96)
    }
    
    # Step 5: If no hits with pident > 97, check for pident > 94
    if (nrow(selected_hits) == 0) {
      selected_hits <- hash_hits %>%
        filter(pident > 94)
    }
    
    # Step 6: If no hits meet the thresholds, select all available hits
    if (nrow(selected_hits) == 0) {
      selected_hits <- hash_hits
    }
    
    # Store the selected hits for this hash_id
    result_list[[hash_id]] <- selected_hits
  }
  
  
  selected_data <- do.call(rbind, result_list)
  
  # *Compute max pident per Hash*
  Max_pident_df <- selected_data %>%
    group_by(Hash) %>%
    summarise(Max_pident = max(pident, na.rm = TRUE))

  
  #Build data `a` to run LCA
  a <- selected_data %>% 
    select(Hash, staxids) %>% 
    drop_na() %>% 
    unique() %>% 
    group_by(Hash) %>% 
    summarise(tax_list = paste(staxids, collapse = " ")) %>% 
    mutate(tax_list = gsub(";", " ", tax_list))
  
  #summarize accession numbers for each Hash ----
  acc_df <- selected_data %>%       
    select(Hash, Accession) %>%
    group_by(Hash) %>%
    summarise(AccessionNumbers = paste(Accession, collapse = ", "), .groups = "drop")
  
  # Filter out hashes already in the db
  if (!is.null(DB_PATH_IN) && nrow(db) > 0){
    a <- a %>% 
      filter(!Hash %in% db$Hash)
  }
  
  if (dir.exists(here("tmp")) == FALSE){
    system2("mkdir", shQuote(here("tmp")))  
  }
  
  #write out to run taxonkit LCA
  filepath <- here("tmp/taxids.txt")
  fileout <- here("tmp/taxids_lca.txt")
  lineagesout <- here("tmp/taxids_lineages.txt")
  
  write.table(a, file = filepath, row.names = F, col.names = F, quote = F, sep = "\t")
  exec_internal(taxonkit_path, args = c("lca", filepath, "-i", "2", "-U", "-D", "-o", fileout))
  
  exec_internal(
    taxonkit_path,
    args = c("reformat", fileout, "-I", "3",
             "-f", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}",
             "-o", lineagesout)
  )
  
  # Read columns deterministically
  b <- read.table(lineagesout, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, header = FALSE)
  names(b)[1:10] <- c("Hash","BlastTaxIDs","LCA_taxid",
                      "Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
  d <- blastoutput %>% 
    select(Hash, Taxon) %>% 
    drop_na() %>% 
    unique() %>% 
    group_by(Hash) %>% 
    summarise(Tax_list = paste(Taxon, collapse = ", ")) 
  
  e <- blastoutput %>% 
    select(Hash, mismatch) %>% 
    drop_na() %>% 
    unique() %>% 
    group_by(Hash) %>% 
    summarise(Mismatches = paste(mismatch, collapse = ", "))
  
  f <- d %>%
    left_join(e, by = "Hash") %>%
    left_join(b, by = "Hash") %>%
    left_join(acc_df,  by = "Hash") %>%
    # Keep every hash that actually got a lineage back from taxonkit.
    # (Hashes already present in `db` were filtered out of `a` upstream, so they
    #  come back as all-NA here and are the ones we want to drop.)
    # NB: do NOT filter on Order != "" -- NCBI leaves the order rank blank for
    # many fish families (Embiotocidae, Moronidae, Sciaenidae, Pomacentridae...),
    # which silently discarded species-level annotations.
    filter(!is.na(LCA_taxid))
  
  # ** Join the max pident to f **
  f <- f %>%
    left_join(Max_pident_df, by = "Hash")
  
  # ** If Max_pident < Species threshold, blank out Species **
  f <- f %>%
    mutate(Max_pident = suppressWarnings(as.numeric(Max_pident))) %>%  # ensure numeric
    mutate(Species = if_else(!is.na(Max_pident) & Max_pident < SP_MAX_THOLD, "", Species))

  # ** Now run rank-choosing logic **
  # Descend the lineage and take the finest rank that is actually populated.
  # Blank ("") and NA are treated identically, and gaps at any rank are skipped
  # rather than aborting the descent.
  blank_to_na <- function(x) if_else(is.na(x) | x == "", NA_character_, x)

  f <- f %>%
    mutate(BestTaxon = coalesce(
      blank_to_na(Species),
      blank_to_na(Genus),
      blank_to_na(Family),
      blank_to_na(Order),
      blank_to_na(Class),
      blank_to_na(Phylum),
      blank_to_na(Kingdom)
    ))
  
  # Haplotype assignment logic
  f <- f %>%
    group_by(BestTaxon) %>%
    mutate(
      MaxHaplotypeNumber = if_else(
        BestTaxon %in% db$BestTaxon,
        coalesce(
          as.numeric(max(db$HaplotypeNumber[db$BestTaxon == BestTaxon], na.rm = TRUE)),
          0
        ),
        0
      ),
      MaxHaplotypeNumber = replace_na(MaxHaplotypeNumber, 0),
      HaplotypeNumber = row_number() + MaxHaplotypeNumber
    ) %>%
    mutate(HaplotypeNumber = if_else(MaxHaplotypeNumber == 0, row_number(), HaplotypeNumber)) %>%
    select(-MaxHaplotypeNumber) %>%
    relocate(Hash, BestTaxon, HaplotypeNumber, Tax_list, Mismatches)
  
  #ensure format is right
  f <- f %>%
    mutate(
      HaplotypeNumber = as.numeric(HaplotypeNumber),
      AccessionNumbers = as.character(AccessionNumbers) , 
      BlastTaxIDs = as.character(BlastTaxIDs),
      HaplotypeNumber = if_else(is.na(HaplotypeNumber) | HaplotypeNumber == "#NAME?", NA_real_, HaplotypeNumber)
    ) %>%
    filter(
      !is.na(HaplotypeNumber) &
        !is.na(BestTaxon) &
        !is.na(Hash) &
        !is.na(BlastTaxIDs)
    )
  
  ids_in_fasta <- as.character(ShortRead::id(fasta))
  f            <- f %>% filter(Hash %in% ids_in_fasta)          # drop rows with no sequence
  
  f$Sequence   <- sread(fasta[ match(f$Hash, ids_in_fasta) ]) %>% as.character()
  
  f$DateAdded <- format(Sys.Date(), "%Y-%m-%d")
  
  if (!is.null(DB_PATH_IN) && nrow(db) > 0){
    g <- db %>%
      mutate(Mismatches = as.character(Mismatches)) %>%  
      bind_rows(f %>% mutate(Mismatches = as.character(Mismatches))) 
  } else {
    g <- f
  }
  
  if (!is.null(DB_PATH_OUT)){
    write.csv(g, file = DB_PATH_OUT)
  } else {
    return(f)
  }
}
