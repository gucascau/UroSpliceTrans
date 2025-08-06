# UroSpliceTrans
UroSpliceTrans is a bioinformatics pipeline designed to translate alternative splicing events into predicted urinary peptide sequences. By integrating exon-level splicing information with genomic sequences, UroSpliceTrans reconstructs transcript isoforms, identifies open reading frames (ORFs), and generates corresponding peptide sequences that reflect splicing-specific alterations. This tool enables researchers to connect transcriptomic splicing changes to the urinary peptidome, facilitating the discovery of novel biomarkers and providing mechanistic insights into splicing-driven disease processes.

# Availability 
1. Splice-Aware Isoform Reconstruction based on user-defined exon combinations.
2. Genomic Sequence Extraction for accurate exon retrieval.
3. Translation of Splice Isoforms into Peptides with reading frame correction.
4. Urinary Peptidome Mapping for biomarker discovery.
5. Define Urinary peptide biomarkers of disease specific splicing alternations


# Dependencies

. blast (ncbi-blast-2.8.1+)(https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

. bedtools (bedtools-2.25.0) (https://bedtools.readthedocs.io/en/latest/)

. R packages: BSgenome.Hsapiens.UCSC.hg38, Biostrings, GenomicRanges

# Install

```
    cd ~
    git clone https://github.com/gucascau/UroSpliceTrans.git
```   

# Usage

### 1. Splice-Aware Isoform Reconstruction by rMATS-turbo. 
This software will identify differential alternative splicing events and quantifies percent spliced in (PSI) values. It analyzes five basic AS event types: skipped exon (SE), alternative 5′ splice sites (A5SS), alternative 3′ splice sites (A3SS), mutually exclusive exons (MXE) and retained intron (RI).

### 2. Transfer the splicing events into peptides 
```
Rscript Src/fromSplicingToAminoAcids.R
```
### 3. Blast human urinary peptides against these mouse splicing-derived peptides.
```
bash Src/run_blast.sh # Run Shell BLAST script

Rscript Src/parse_blast_results.R # Run R script to filter results
```

# Workflow Overview
**1. Generation of Mouse Splicing-Derived Peptides:**
Alternative splicing events were translated into peptide sequences using the UroSpliceTrans pipeline.
The resulting peptide sequences were compiled into a FASTA file representing the mouse splicing-derived peptidome.  
**2. Human Urinary Peptide Dataset Preparation:**
Human urinary peptides identified from mass spectrometry were curated and formatted as a FASTA query file.  
**3. BLASTP Search (Peptide-to-Peptide Alignment):**
A BLAST protein database was created from the mouse splice-derived peptides.
Human urinary peptides were aligned against this database using BLASTP, optimized for short peptide alignments.
The search parameters included a high e-value threshold and small word size to capture short peptide matches.  
**4. Filtering and Hit Selection:**
BLAST alignments were filtered based on percentage identity (≥90%), minimum alignment length (≥6 amino acids), and stringent e-value thresholds (≤1e-3).
The filtered high-confidence peptide matches were compiled for downstream analysis.  

# Contact

For each insertion event, we selected the most reliable representative read that showed the highest read quality, highest donor identity and the highest read count support. 

For more detail information, please feel free to contact: xin.wang@nationwidechildrens.org

This project is licensed under the terms of the MIT license.
