# load the DESeq2 package
library(DESeq2)
# load the Rsubread package
library(Rsubread)
# load the apeglm package
library(apeglm)
# load the EnhancedVolcano package
library(EnhancedVolcano)
# load the RColorBrewer package
library(RColorBrewer)
# load the pheatmap package
library(pheatmap)

##3 args in total
#1 -> path to samplesheet e.g. "/gs/gsfs0/users/shechter-lab/data/NGS/pro-seq/2020-07_proseq/Processed/results_pipeline/DESEQ2/samplesheet.csv"
#2 -> path to "results_pipeline" e.g. "/gs/gsfs0/users/shechter-lab/data/NGS/pro-seq/2020-07_proseq/Processed/results_pipeline/"
#3 -> path to gtf file e.g. "/gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38_ncbiRefSeq_noTSS.gtf"

#make script accept positiional command line arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#get the name of the analysis from the directory name
analysis_name <- basename(getwd())

#read in the sample table
sampleTable <- read.csv(args[1], row.names=1)
#create a vector of file names
filenames <- paste0(args[2],sampleTable$run,'/aligned_hg38/',sampleTable$run,'_sort.bam')
#read in the gtf file
gtffile <- args[3]
#counts the number of reads that map to each gene
fc <- featureCounts(files = filenames, annot.ext=gtffile, isGTFAnnotationFile=TRUE, isPairedEnd=FALSE, GTF.attrType="gene_name", GTF.featureType = 'transcript', strandSpecific = 1, nthreads = 40, ignoreDup = TRUE, countMultiMappingReads = FALSE)
#assign the column names of the count data to the sample names
colnames(fc$counts) <- sampleTable$run
#assign the count data to a variable
countdata <- fc$counts
#assign the condition data to a variable
coldata <- as.data.frame(sampleTable$condition)
#assign the column name of the condition data to 'condition'
colnames(coldata) <- 'condition'
#create a DESeqDataSet object from the count data, condition data, and design formula
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
#keep only genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts(dds)) >= 10
#filter the DESeqDataSet object to only include genes with at least 10 counts in at least 3 samples
dds <- dds[keep,]

#relevel the condition column so that the control is the reference
dds$condition <- relevel(dds$condition, ref = "Untreated")

dds <- DESeq(dds)

# create a vst object from the dds object
vsd <- vst(dds, blind=FALSE)

# load the column data
cdata <- colData(dds)


# create a heatmap of the count matrix
pdf(paste0(analysis_name, "_Heatmap_CountMatrix.pdf"), height = 8, width = 8)
pheatmap(assay(vsd),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = as.data.frame(cdata[,"condition"], row.names=rownames(cdata)))
dev.off()

# calculate the distance between each sample
sampleDists <- dist(t(assay(vsd)))

# create a matrix of the sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste0(analysis_name, "_Heatmap_SampleDistances.pdf"), height = 8, width = 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#Create a PCA of the sample data
pdf(paste0(analysis_name, "_PCA.pdf"), height = 8, width = 8)
plotPCA(vsd, intgroup=c("condition"))
dev.off()


# for loop to iterate through each comparison
for (i in 2:length(unique(resultsNames(dds)))){
  # name the comparison
  name <- unique(resultsNames(dds))[i]
  # run DESeq2 analysis
  res <- lfcShrink(dds, coef=name, type="apeglm")
  # add gene_id column to results
  res$gene_id <- rownames(res)
  # convert results to dataframe
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj),]
  write.csv(res_df, file=paste0(name,"_DESEQ2.csv"), row.names=F)
  pdf(paste0(name,"_MAplot.pdf"), height = 8, width = 8)
  plotMA(res)
  dev.off()
  res_df <- res_df[!is.na(res_df$padj),]
  pdf(paste0(name,'_Volcano_plot.pdf'), height = 8, width = 8)
  print(EnhancedVolcano(res_df,
                        lab = res_df$gene_id,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = name,
                        pCutoff = 0.05,
                        FCcutoff = 0.2,
                        pointSize = 2.0,
                        labSize = 5.0,
                        col=c('#d4d4d4', '#d4d4d4', '#fdbb84', '#e34a33'),
                        colAlpha = 0.5))
  dev.off()
}
