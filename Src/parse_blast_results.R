# Load Required Libraries
library(dplyr)
library(readr)

# Load BLASTP results (outfmt 6 columns)
blast_cols <- c("query_id", "subject_id", "perc_identity", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", 
                "evalue", "bit_score")

blast_results <- read_tsv("blast_results/blastp_results.tsv", col_names = blast_cols)

# Filter criteria
min_identity <- 90  # Percentage Identity Cutoff
min_alignment_length <- 6  # Minimum peptide match length
max_evalue <- 1e-3  # Stringent e-value cutoff

filtered_hits <- blast_results %>%
  filter(perc_identity >= min_identity,
         alignment_length >= min_alignment_length,
         evalue <= max_evalue) %>%
  arrange(desc(perc_identity), evalue)

# Output filtered results
write_tsv(filtered_hits, "blast_results/filtered_peptide_hits.tsv")

cat("Filtered BLAST hits saved to blast_results/filtered_peptide_hits.tsv\n")
