# -----------------------------------------------------------
# Script: ensg_to_bed.R
# Purpose: Convert a list of Ensembl Gene IDs (ENSG) into a
#          6-column directional BED file using Ensembl biomaRt.
# David Shechter - 6-11-25
# -----------------------------------------------------------

# Load required libraries
library(biomaRt)   # for querying Ensembl annotations
library(readr)     # for reading/writing files
library(dplyr)     # for data manipulation

# --------------------------
# Step 1: Load Input IDs
# --------------------------

filename <- "histoneH1_genes_HDB20"
workingdir <- "C:/Users/david/OneDrive/Bioinformatics/PRMTi/HistoneGeneAnalysis/"
setwd(workingdir)

# Input: CSV file with a single column of ENSG IDs (no header)
input_file <- paste0(workingdir, filename, ".csv")
ensg_ids <- read_csv(input_file, col_names = FALSE)$X1

# --------------------------
# Step 2: Connect to Ensembl
# --------------------------

# Use the Ensembl Genes dataset for human
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# --------------------------
# Step 3: Query Gene Coordinates
# --------------------------

# Define attributes to fetch from BioMart
attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand"
)

# Query Ensembl using our ENSG ID list
gene_coords <- getBM(
  attributes = attributes,
  filters = "ensembl_gene_id",
  values = ensg_ids,
  mart = mart
)

# --------------------------
# Step 4: Clean and Format BED
# --------------------------

# Restrict to standard chromosomes (1–22, X, Y, MT)
standard_chr <- c(as.character(1:22), "X", "Y", "MT")
gene_coords <- gene_coords %>%
  filter(chromosome_name %in% standard_chr)

# Construct BED fields:
# - BED is 0-based, so subtract 1 from start
# - Use gene symbol or ENSG as name
# - Use strand to get '+' or '-'
bed_df <- gene_coords %>%
  mutate(
    chrom = paste0("chr", chromosome_name),
    start = start_position - 1,
    end = end_position,
    name = external_gene_name,
    name = ifelse(is.na(name) | name == "", ensembl_gene_id, name),
    score = 0,
    strand = ifelse(strand == 1, "+", "-")
  ) %>%
  select(chrom, start, end, name, score, strand) %>%
  arrange(chrom, start)

# --------------------------
# Step 5: Write BED File
# --------------------------

# Output file
output_file <- paste0(workingdir, filename, ".bed")

write_tsv(bed_df, output_file, col_names = FALSE)

# Done!
message("BED file written to: ", output_file)