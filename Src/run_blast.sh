#!/bin/bash

# Set paths
MOUSE_DB="blast_db/mouse_splice_db"
QUERY="human_urine_peptides.fasta"
DB_INPUT="mouse_splice_peptides.fasta"
RESULTS_DIR="blast_results"
OUT_FILE="$RESULTS_DIR/blastp_results.tsv"

mkdir -p blast_db $RESULTS_DIR

# Step 1: Make BLAST database
makeblastdb -in $DB_INPUT -dbtype prot -out $MOUSE_DB

# Step 2: Run BLASTP (optimized for short peptides)
blastp -query $QUERY -db $MOUSE_DB -out $OUT_FILE -evalue 1000 -word_size 2 -outfmt 6

echo "BLASTP finished. Results saved to $OUT_FILE"
