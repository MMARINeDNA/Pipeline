##############################################################################
##                         Welcome to The Kelly Lab's                       ##
##              Metabarcoding Taxonomic Assignment Pipeline                 ##
## This pipeline will take you from de-multiplexed fastq files to a matrix  ##
## with taxon names and counts. It should work for linux or mac. PC users,  ##
## look into WSL2.                                                          ##
##############################################################################
##    Before you start, ensure you have access to the following files.      ##
##                                                                          ##
## Required Files:                                                          ##
## 1. Your Fastqs, of course. Their names must start with the primer name   ##
## 2. CEG_BLAST_function.R - Function to run BLAST on the CEG cluster       ##
## 3. LCA_function.R - Function to run LCA                                  ##
## 4. Pipeline_code.R - This very code                                      ##
## 5. MURI_taxids_to_exclude.txt #optional                                  ##
## 6. MURIblast_*_template.sh - Code for taxonomic assignment               ##
##    depending on delimitation                                             ##
## 7. make_primer_shell_script.R - A script that writes another script      ##
## 8. primer_data.csv - Sheet with primer information like sequence, etc    ##
##                                                                          ##
## Required Programs:                                                       ##
## 1. Taxonkit: https://bioinf.shenwei.me/taxonkit/                         ##
## 2. Cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html ##
##############################################################################

# Restart your R session and clear the environment.
rm(list = ls())

## put the path to this very code here.
here::i_am("pipeline_code_4.0.R")

# Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(dada2))
suppressMessages(library(digest))
suppressMessages(library(seqinr))
suppressMessages(library(ssh))
suppressMessages(library(sys))
suppressMessages(library(ShortRead))
suppressMessages(library(here))

# Set working directory (this will change again downstream)
setwd(here())

# Source external functions for BLAST and LCA
source("CEG_BLAST_function_2.0.R")
source("LCA_function_4.0.R")
# Define paths and parameters
PARENT_LOCATION <- "/mnt/c/Users/pedro/OneDrive/Documents/UW/Pipeline/test"  # Pedro's PC
if (strsplit(PARENT_LOCATION, "")[[1]][nchar(PARENT_LOCATION)] != "/") {
  PARENT_LOCATION <- paste0(PARENT_LOCATION, "/")
}

RUN_NAME <- "Michelle"
PRIMERNAME <- "MV1"  # must match one of the known primers
BLASTSCOPE <- "eukaryote"  # either "vertebrate" or "eukaryote"
SP_THOLD <- 96  # sequence similarity threshold

NEGATIVE_TAXIDS_FILE <- "MURI_taxids_to_exclude.txt"  # optional

# Where the hash database is
DATABASE_LOCATION <- "/mnt/c/Users/pedro/OneDrive/Documents/UW/Pipeline/LCA3_databases"
if (strsplit(DATABASE_LOCATION, "")[[1]][nchar(DATABASE_LOCATION)] != "/") {
  DATABASE_LOCATION <- paste0(DATABASE_LOCATION, "/")
}

# Define the output location for processed files
PROCESSED_LOCATION <- "/mnt/c/Users/pedro/OneDrive/Documents/UW/Pipeline"  # Pedro's PC

system2("mkdir", shQuote(paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME)))  # create folder for code, etc
FASTQ_LOCATION <- paste0(PARENT_LOCATION, "Fastq")          # folder within Parent Location
CODE_LOCATION <- paste0(PARENT_LOCATION, "code_etc")         # folder within Parent Location
system2("mkdir", args = shQuote(CODE_LOCATION))              # create folder for code, etc
TRIMMED_LOCATION <- paste0(PARENT_LOCATION, "for_dada2")       # folder within Parent Location
FILTERED_LOCATION <- paste0(PARENT_LOCATION, "filtered")       # folder within Parent Location
OUTPUT_LOCATION <- paste0(PARENT_LOCATION, "outputs")          # folder within Parent Location
system2("mkdir", args = shQuote(FILTERED_LOCATION))           # create folder for processed reads
system2("mkdir", args = shQuote(OUTPUT_LOCATION))             # create folder for pipeline outputs

# Check dependencies and paths for cutadapt and taxonkit
CUTADAPT <- "/home/pedrobdfp/.local/bin/cutadapt"  # Pedro's laptop
system2(CUTADAPT, args = shQuote("--version"))

TAXONKIT_PATH <- "/usr/local/bin/taxonkit"  # Pedro's laptop
system2(TAXONKIT_PATH, args = shQuote("version"))

# Read in primer data
primer.data <- read.csv("/mnt/c/Users/pedro/OneDrive/Documents/UW/Pipeline/primer.data.csv")
PRIMERSEQ_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_f)
PRIMERSEQ_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_r)
PRIMERLENGTH_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_f)
PRIMERLENGTH_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_r)
MAX_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(max_amplicon_length)
MIN_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(min_amplicon_length)
OVERLAP <- primer.data %>% filter(name == PRIMERNAME) %>% pull(overlap)

### COPY FILES TO SCRIPTS FOLDER and write primer-trimming script via make_primer_shell_script.R
write.csv(primer.data, paste0(CODE_LOCATION, "/primer.data.csv"))
system2("cp", args = c(shQuote(here("make_primer_shell_script.R")), shQuote(CODE_LOCATION)))
system2("cp", args = c(shQuote(here(NEGATIVE_TAXIDS_FILE)), shQuote(CODE_LOCATION)))  # taxa masked from blast search
source(paste0(CODE_LOCATION, "/make_primer_shell_script.R"))  # create shell script for primer trimming
system2("cp", args = c(shQuote(here("pipeline_code_3.1.R")), shQuote(CODE_LOCATION)))  # copy this pipeline code

if (BLASTSCOPE == "vertebrate") {
  system2("cp", args = c(shQuote(here("MURIblast_vertebrate_template.sh")), shQuote(CODE_LOCATION)))  # Copy vertebrate template  
}
if (BLASTSCOPE == "eukaryote") {
  system2("cp", args = c(shQuote(here("MURIblast_eukaryote_template.sh")), shQuote(CODE_LOCATION)))  # Copy eukaryote template  
}

## Set this to true if you are rerunning the same data and already have the hash_key and ASV_table
SKIP_DADA2 <- FALSE

#-----------------------------------------------------------------------------------------------
### RUN PRIMER-TRIMMING -----------------------------------------------------------
system2("sh", args = shQuote(paste0(CODE_LOCATION, "/trim_primers.sh")))

### RUN DADA2 ---------------------------------------------------------------------
filelist <- system2("ls", args = shQuote(TRIMMED_LOCATION), stdout = TRUE)

# Identify forward and reverse files
fnFs <- filelist[str_detect(filelist, "_R1")]
fnRs <- filelist[str_detect(filelist, "_R2")]

# Extract sample names (assuming they are the first two fields separated by "_")
sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")

### Name filtered files in filtered/subdirectory ----------------------------------
filtFs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### Filter out Empty Samples (Modified to retain empty samples as zeros) ---------------
setwd(TRIMMED_LOCATION)
file.empty <- function(filenames) file.info(filenames)$size == 20
empty_files <- file.empty(fnFs) | file.empty(fnRs)
# Instead of dropping empty samples, we record the empty indices but keep all samples.
# (Do NOT subset away empty samples from fnFs, fnRs, filtFs, filtRs, or sample.names)

### (Optional) Manually check quality
plotQualityProfile(fnFs[6:10])
plotQualityProfile(fnRs[2:5])

##########
# Run filtering and trimming. For empty files, filterAndTrim() should return zeros.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     #truncLen = 200,
                     #trimRight = c(PRIMERLENGTH_R,PRIMERLENGTH_F),
                     #truncLen = round(MAX_AMPLICON_LENGTH/2)+OVERLAP,
                     truncLen = MIN_AMPLICON_LENGTH,
                     #maxLen = MAX_AMPLICON_LENGTH+OVERLAP,
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = FALSE, matchIDs = TRUE)

### Filter ---------------------------------------------------------------
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

## Learn error rates ---------------------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

### Dereplicate and Learn Error Rates ---------------------------------------------------------------
dadaFs <- dada(filtFs, err = errF, selfConsist = TRUE, multithread = TRUE, MAX_CONSIST = 20)
dadaRs <- dada(filtRs, err = errR, selfConsist = TRUE, multithread = TRUE, MAX_CONSIST = 20)

### Merge Paired Reads ---------------------------------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = OVERLAP, verbose = TRUE, trimOverhang = TRUE)

### Construct sequence table ---------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)

### Remove chimeras ------------------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
freq.nochim <- sum(seqtab.nochim) / sum(seqtab)

### Filter by Size -------------------------------------------------------------------------
indexes.to.keep <- which((nchar(colnames(seqtab.nochim)) <= MAX_AMPLICON_LENGTH) & 
                           (nchar(colnames(seqtab.nochim)) >= MIN_AMPLICON_LENGTH))
cleaned.seqtab.nochim <- seqtab.nochim[, indexes.to.keep]
filteredout.seqtab.nochim <- seqtab.nochim[, !indexes.to.keep]
write.csv(filteredout.seqtab.nochim, paste0(FASTQ_LOCATION, "/../logs/", "filtered_out_asv.csv"))

### Track reads through pipeline (Modified to include empty samples as zeros) -----------
getN <- function(x) sum(getUniques(x))
denoisedF <- sapply(sample.names, function(s) if(s %in% names(dadaFs)) getN(dadaFs[[s]]) else 0)
denoisedR <- sapply(sample.names, function(s) if(s %in% names(dadaRs)) getN(dadaRs[[s]]) else 0)
merged_vals <- sapply(sample.names, function(s) if(s %in% names(mergers)) getN(mergers[[s]]) else 0)
nonchim   <- sapply(sample.names, function(s) if(s %in% rownames(seqtab.nochim)) rowSums(seqtab.nochim)[s] else 0)

track <- cbind(out, denoisedF, denoisedR, merged_vals, nonchim)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, paste0(FASTQ_LOCATION, "/../logs/", "tracking_reads.csv"))

### Create Hashing  ------------------------------------------------------------------------
conv_file <- paste0(OUTPUT_LOCATION, "/", paste0(RUN_NAME, "_", PRIMERNAME, "_hash_key.csv"))
conv_file.fasta <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_hash_key.fasta"))
ASV_file <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_ASV_table.csv"))
taxonomy_file <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_taxonomy_output.csv"))
bootstrap_file <- file.path(OUTPUT_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME, "_tax_bootstrap.csv"))

print(paste0("creating ASV table and hash key...", Sys.time()))
seqtab.nochim.df <- as.data.frame(cleaned.seqtab.nochim)
conv_table <- tibble(Hash = "", Sequence = "")
Hashes <- map_chr(colnames(seqtab.nochim.df), ~ digest(.x, algo = "sha1", serialize = F, skip = "auto"))
conv_table <- tibble(Hash = Hashes, Sequence = colnames(seqtab.nochim.df))

write_csv(conv_table, conv_file)
write.fasta(sequences = as.list(conv_table$Sequence),
            names = as.list(conv_table$Hash),
            file.out = conv_file.fasta)
sample.df <- tibble::rownames_to_column(seqtab.nochim.df, "Sample_name")
sample.df <- data.frame(append(sample.df, c(Label = PRIMERNAME), after = 1))
current_asv <- bind_cols(sample.df %>% dplyr::select(Sample_name, Label), seqtab.nochim.df)
current_asv <- current_asv %>%
  pivot_longer(cols = c(-Sample_name, -Label), names_to = "Sequence", values_to = "nReads") %>%
  filter(nReads > 0)
current_asv <- merge(current_asv, conv_table, by = "Sequence") %>%
  select(-Sequence) %>%
  relocate(Hash, .after = Label)

write_csv(current_asv, ASV_file)

### Move files to processed folder
system2("mv", args = c(shQuote(paste0(PARENT_LOCATION, "code_etc")),
                       shQuote(paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME))))
system2("mv", args = c(shQuote(paste0(PARENT_LOCATION, "outputs")),
                       shQuote(paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME))))
system2("mv", args = c(shQuote(paste0(PARENT_LOCATION, "logs")),
                       shQuote(paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME))))

## Optional: If SKIP_DADA2 is TRUE, load previous ASV and hash key
if (SKIP_DADA2) {
  existing_asv_table <- file.path(PROCESSED_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME), "outputs", 
                                  paste0(RUN_NAME, "_", PRIMERNAME, "_ASV_table.csv"))
  existing_hash_key <- file.path(PROCESSED_LOCATION, paste0(RUN_NAME, "_", PRIMERNAME), "outputs", 
                                 paste0(RUN_NAME, "_", PRIMERNAME, "_hash_key.csv"))
  library(readr)
  current_asv <- read_csv(existing_asv_table)
  conv_table <- read_csv(existing_hash_key)
  sample.names <- sort(unique(current_asv$Sample_name))
}

### ANNOTATION ------------------------------------------------------------------------------------------
if (file.exists(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))) {
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
  seen <- which(conv_table$Hash %in% db$Hash)
  notseen <- which(!conv_table$Hash %in% db$Hash)
  
  write.fasta(sequences = as.list(conv_table$Sequence[notseen]),
              names = as.list(conv_table$Hash[notseen]),
              file.out = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/seqs_to_annotate.fasta"))
} else {
  write.fasta(sequences = as.list(conv_table$Sequence),
              names = as.list(conv_table$Hash),
              file.out = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/seqs_to_annotate.fasta"))
}

# Blast new sequences on CEG server
ceg_blast(
  PATH_TO_FASTA = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/seqs_to_annotate.fasta"),
  PATH_TO_BLAST_TEMPLATE = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/code_etc/MURIblast_", BLASTSCOPE, "_template.sh"),
  PATH_FOR_CREATED_SHELL_SCRIPT = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/code_etc/MURIblast_new_seqs.sh"),
  PATH_FOR_RESULTS = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/new_annotations.txt"),
  NEGATIVE_TAXIDS_FILE = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/code_etc/", NEGATIVE_TAXIDS_FILE),
  CUL_N = "50"
)

if (file.size(paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/new_annotations.txt")) == 0L) {
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
} else if (file.exists(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))) {
  LCA(BLASTOUTPUT = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/new_annotations.txt"),
      FASTA = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/seqs_to_annotate.fasta"),
      DB_PATH_IN = paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"),
      DB_PATH_OUT = paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"),
      SP_THOLD = SP_THOLD)
  
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
} else {
  db <- LCA(BLASTOUTPUT = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/new_annotations.txt"),
            FASTA = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/seqs_to_annotate.fasta"),
            SP_THOLD = SP_THOLD)
  write.csv(db %>% distinct(), paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))
}

## Writing the output taxon tables
tax_table <- current_asv %>%
  left_join(db %>% dplyr::select(Hash, BestTaxon, Class)) %>%
  group_by(Sample_name, BestTaxon, Class) %>%
  summarise(nReads = sum(nReads))
write_csv(tax_table, file = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/", RUN_NAME, "_", PRIMERNAME, "_", "taxon_table.csv"))

# Create taxon_table_wide with all samples (pad missing ones with zeros)
full_samples <- sample.names
tax_wide <- tax_table %>%
  pivot_wider(id_cols = c(BestTaxon, Class), names_from = Sample_name, values_from = nReads, values_fill = 0)
missing_samples <- setdiff(full_samples, colnames(tax_wide))
if(length(missing_samples) > 0) {
  for(s in missing_samples) {
    tax_wide[[s]] <- 0
  }
}
tax_wide <- tax_wide %>% dplyr::select(BestTaxon, Class, all_of(sort(full_samples)))
write_csv(tax_wide, file = paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/", RUN_NAME, "_", PRIMERNAME, "_", "taxon_table_wide.csv"))

# A taxon table but that lists variants/haplotypes separately
merged_data <- current_asv %>%
  left_join(db, by = "Hash")
merged_data <- merged_data %>%
  mutate(BestTaxon_Haplotype = paste(BestTaxon, HaplotypeNumber, sep = "_"))
summarized_data <- merged_data %>%
  group_by(BestTaxon_Haplotype, Sample_name, Class) %>%
  summarise(nReads = sum(nReads, na.rm = TRUE), .groups = 'drop')
reshaped_data <- summarized_data %>%
  pivot_wider(id_cols = c(BestTaxon_Haplotype, Class), names_from = Sample_name, values_from = nReads, values_fill = 0)

# Pad missing samples with zeros in the haplotype table as well
missing_samples_haplo <- setdiff(full_samples, colnames(reshaped_data))
if(length(missing_samples_haplo) > 0) {
  for(s in missing_samples_haplo) {
    reshaped_data[[s]] <- 0
  }
}
reshaped_data <- reshaped_data %>% dplyr::select(BestTaxon_Haplotype, Class, all_of(sort(full_samples)))
write.csv(reshaped_data, paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME, "/outputs/", RUN_NAME, "_", PRIMERNAME, "_", "haplotype_table.csv"), row.names = FALSE)

message(paste0("Pipeline complete. Outputs are now available in ", paste0(PROCESSED_LOCATION, "/", RUN_NAME, "_", PRIMERNAME)))

### CLEANUP
CLEANUP = TRUE
if (CLEANUP) {
  for (j in c("code_etc", "outputs", "filtered", "for_dada2", "logs")) {
    system2("rm", args = c("-r", shQuote(paste0(PARENT_LOCATION, j))))
  }
}
