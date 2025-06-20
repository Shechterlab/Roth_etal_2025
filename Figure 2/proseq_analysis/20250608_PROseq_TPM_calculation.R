setwd('E:\\PRMT\\proseq/peppro/hg38/')
library(Rsubread)
library(tidyverse)
sampleTable <- read.table("2020_arpeggio_metadata.txt", sep = "\t", header = T)
sampleTable <- sampleTable[c(1:2,11:12,15:16,19:20,31:34,41:42),]
sampleTable$name <- sub("_R1_001", "_sort.bam", sampleTable$fastq_file)
filenames <- sampleTable$name
gtffile <- 'hg38_ncbiRefSeq_noTSS.gtf'
fc <- featureCounts(files = filenames, annot.ext=gtffile, isGTFAnnotationFile=TRUE, isPairedEnd=FALSE, GTF.attrType="gene_name", GTF.featureType = 'transcript', strandSpecific = 1, nthreads = 2, ignoreDup = TRUE, countMultiMappingReads = FALSE)
colnames(fc$counts) <- sampleTable$description
countdata <- fc$counts
count_df <- as.data.frame(countdata)
count_df$gene_name <- rownames(count_df)
gtf <- read.table(gtffile, header = F, sep = '\t')
gtf <- gtf %>%
  mutate(V9 = str_split(V9, ";\\s*")) %>%
  unnest(V9) %>%
  separate(V9, into = c("key", "value"), sep = " ") %>%
  pivot_wider(names_from = key, values_from = value)
gene_start_end <- gtf %>%
  group_by(gene_name) %>%
  summarise(start = min(V4), end = max(V5), .groups = "drop")
gene_start_end$width <- (gene_start_end$end - gene_start_end$start)
count_df <- merge(count_df, gene_start_end, by='gene_name')

count_df <- count_df[c(grep('gene_name', colnames(count_df)),grep('width', colnames(count_df)),1:nrow(sampleTable)+1)]

count_df <- count_df[count_df$width > 50,]

countToTpm <- function(counts, effLen) {
  counts <- as.numeric(counts)
  effLen <- as.numeric(effLen)
  
  #Normalize counts by length (counts per base)
  rpk <- counts / effLen
  
  #Calculate scaling factor (sum of all RPKs)
  scaling_factor <- sum(rpk)
  
  #Scale to 1 million
  tpm <- (rpk / scaling_factor) * 1e6
  
  return(tpm)
}

#This is a loop to apply the above function to  all columns containing counts in a dataframe
tpm <- data.frame(matrix(ncol = (ncol(count_df)-2), nrow = nrow(count_df)))
for (i in 3:(ncol(count_df))) {
  tpm[,i-2] <- countToTpm(counts = count_df[,i], effLen = count_df$width)
}
colnames(tpm) <- colnames(count_df[1:nrow(sampleTable)+2])
tpm$gene_name <- count_df$gene_name

min_count <- as.data.frame((tpm[, -ncol(tpm)] > 10) * 1)


min_count$gene_name <- tpm$gene_name
min_count <- min_count[rowSums(min_count[,1:(ncol(min_count)-1)]) > 0, ]

tpm <- tpm[(tpm$gene_name %in% min_count$gene_name),]

write.table(tpm, '20250608_GSK591_PROseq_TPM.txt', col.names = TRUE, row.names = FALSE)
#Breaking the code apart and performing it with 3 separate loops
#First step is to adjust the reads for the gene length
#df <- data.frame(matrix(ncol = (ncol(count_df)-2), nrow = nrow(count_df)))
#for (i in 2:(ncol(count_df)-1)) {
#  for (j in 1:nrow(count_df)) {
#    df[j,i-1] <- count_df[j,i] / count_df[j,6]
#  }
#}
#second step is to sum all of the adjusted reads per condition
#denom <- data.frame(matrix(ncol = (ncol(df)), nrow = 1))
#for (i in 1:(ncol(df))) {
#    denom[,i] <- sum(df[,i])
#}
#next step is to divide each adjust read by the sum of all reads and multiply by 1e6
#tpm <- data.frame(matrix(ncol = ncol(df), nrow = nrow(df)))
#for (i in 1:(ncol(df))) {
#  for (j in 1:nrow(df)) {
#    tpm[j,i] <- df[j,i] / denom[,i] * 1e6
#  }
#}