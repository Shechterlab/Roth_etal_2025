#!/usr/bin/env Rscript
# Maxim Maron & David Shechter
# Updated 2024-11-20

setwd(


# Load necessary libraries
library(DESeq2)
library(tibble)
library(dplyr)

# Function to run DESeq2 analysis
run_deseq2 <- function(countData, colData, condition1, condition2) {
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = condition2)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", condition1, condition2))
  res <- as.data.frame(res)
  return(res)
}

# Main script starts here
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No file provided", call. = FALSE)
}

# Read the matrix file
featureCountsMatrix <- read.csv(args[1], sep = "\t", row.names = 1, header = TRUE)

# Extract condition names and prepare metadata
conditions <- colnames(featureCountsMatrix)
conditions <- gsub("_REP[12]$", "", conditions)
colData <- data.frame(sampleName = colnames(featureCountsMatrix), condition = conditions)

# Iterate through each condition and run DESeq2 analysis against the control
for (cond in unique(conditions)) {
  if (cond != "CONTROL") {
    print(paste("Analyzing", cond, "vs CONTROL"))
    res <- run_deseq2(countData = featureCountsMatrix, colData = colData, condition1 = cond, condition2 = "CONTROL")
    
    # Add GeneID and GeneName
    res$GeneID <- rownames(res)
    res$GeneName <- featureCountsMatrix$GeneName[match(rownames(res), featureCountsMatrix$GeneID)]
    
    # Write results to file
    write.csv(res, file = paste0(cond, "_vs_CONTROL_proSeq-genes_DESeq2.csv"), row.names = FALSE)
  }
}

print("DESeq2 analysis completed.")
