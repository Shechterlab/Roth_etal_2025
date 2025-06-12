library(DESeq2)
library(data.table)

# Load the count matrix with gene names (assuming Geneid is the first column)
countData <- fread("updated_featureCounts_matrix_with_gene_names.txt")
rownames(countData) <- countData$Geneid
countData <- countData[, -c(1, 2), with = FALSE]  # Remove Geneid and GeneName columns for DESeq2

# Define the metadata for your samples (this needs to be customized)
# For example:
# colData <- data.frame(
#   condition = factor(c("control", "treatment", "control", "treatment")),
#   row.names = colnames(countData)
# )

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Run the differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Ordering results by p-value
resOrdered <- res[order(res$pvalue),]

# Write results to a file
write.csv(as.data.frame(resOrdered), "DESeq2_results.csv")

