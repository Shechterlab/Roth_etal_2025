library(data.table)
library(org.Hs.eg.db)

# Read your existing matrix
matrix_file <- "merged_featureCounts_matrix.txt"
dt <- fread(matrix_file)

# Strip version numbers from Ensembl IDs
gene_ids <- gsub("\\..*", "", dt$Geneid)

# Map Ensembl IDs (without version numbers) to gene symbols
gene_names <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Map Ensembl IDs to gene symbols, handling 1:many mappings by listing all symbols
gene_names <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "list")

# Convert list to a comma-separated string for each Ensembl ID
gene_names_combined <- sapply(gene_names, function(x) paste(x, collapse = ", "))

# Add gene names to the matrix
dt[, GeneName := gene_names_combined[gene_ids]]





# Add gene names to the matrix
dt[, GeneName := gene_names[gene_ids]]

# Optionally, reorder columns to have GeneName after Geneid
setcolorder(dt, c("Geneid", "GeneName", setdiff(names(dt), c("Geneid", "GeneName"))))

# Save the updated matrix
fwrite(dt, "merged_featureCounts-geneNames_matrix.txt", sep = "\t")

