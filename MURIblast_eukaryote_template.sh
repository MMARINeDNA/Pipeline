#!/bin/bash

# Change to the directory containing the BLAST database
cd /mnt/nfs/home/KellyCEG/blastdb_euk/

# Update the PATH environment variable to include the BLAST+ binaries
PATH=$PATH:/mnt/nfs/home/KellyCEG/ncbi-blast-2.15.0+/bin
export PATH=${PATH}:${HOME}/edirect
export BLASTDB=/mnt/nfs/home/KellyCEG/blastdb_euk/
# Define BLAST database and query FASTA file
BLAST_DB='/mnt/nfs/home/KellyCEG/blastdb_euk/nt_euk'  # The BLAST database
QUERY_FASTA='/mnt/nfs/home/KellyCEG/tmp/raw/seqs_to_annotate.fasta'  # The FASTA file to BLAST

# BLAST PARAMETERS
PERCENT_IDENTITY="90"
WORD_SIZE="30"  
EVALUE="1e-30"  
MAXIMUM_MATCHES="100"
CULLING="50"
BLAST_OUTPUT="/mnt/nfs/home/KellyCEG/tmp/processed/new_annotations.txt"

# Initialize NEGATIVE_TAXIDLIST
NEGATIVE_TAXIDLIST='/mnt/nfs/home/KellyCEG/tmp/raw/MURI_taxids_to_exclude.txt'  # The exclude list

# Handle NEGATIVE_TAXIDLIST_OPTION
# Read taxids to exclude into a variable
if [ -n "$NEGATIVE_TAXIDLIST" ]; then
    EXCLUDE_TAXIDS=$(grep -o '^[^#]*' "$NEGATIVE_TAXIDLIST" \
        | sed '/^\s*$/d;s/[[:blank:]]//g' \
        | paste -sd ',' -)
    NEGATIVE_TAXIDS_OPTION="-negative_taxids \"$EXCLUDE_TAXIDS\""
else
    NEGATIVE_TAXIDS_OPTION=""
fi

# BLAST command
BLAST_CMD="/mnt/nfs/home/KellyCEG/ncbi-blast-2.15.0+/bin/blastn \
    -query \"${QUERY_FASTA}\" \
    -db \"${BLAST_DB}\" \
    -num_threads 16 \
    -perc_identity \"${PERCENT_IDENTITY}\" \
    -word_size \"${WORD_SIZE}\" \
    -evalue \"${EVALUE}\" \
    -max_target_seqs \"${MAXIMUM_MATCHES}\" \
    -culling_limit \"${CULLING}\" \
    ${NEGATIVE_TAXIDS_OPTION} \
    -outfmt \"6 sscinames scomnames qseqid sseqid pident length mismatch gapopen \
qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs\" \
    -out \"${BLAST_OUTPUT}\""
echo $BLAST_CMD
# Execute the BLAST command
eval $BLAST_CMD
