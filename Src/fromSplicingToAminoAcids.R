# Load necessary libraries
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(GenomicRanges)

# --- Step 1: Define Exon Coordinates ---
# Example exons (chromosome, start, end, strand) for mouse
exons <- GRanges(seqnames = c("chr1", "chr1"),
                 ranges = IRanges(start = c(1000, 2000),
                                  end = c(1100, 2100)),
                 strand = c("+", "+"))

# --- Step 2: Extract Exon Sequences from Mouse Genome ---
genome <- BSgenome.Mmusculus.UCSC.mm10
exon_seqs <- getSeq(genome, exons)

# --- Step 3: Concatenate Exon Sequences in Transcript Order ---
full_transcript_seq <- paste0(as.character(exon_seqs), collapse = "")
full_transcript_seq <- DNAString(full_transcript_seq)

# --- Step 4: Adjust Reading Frame (Assuming Starts at First Base) ---
# If you have CDS start position, trim sequence accordingly
# Here, we assume translation starts at position 1 (adjust if needed)

# --- Step 5: Translate to Amino Acid Sequence ---
protein_seq <- translate(full_transcript_seq)

# --- Output ---
cat("Amino Acid Sequence:\n")
cat(as.character(protein_seq), "\n")
