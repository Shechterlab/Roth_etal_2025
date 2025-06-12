#set working directory
setwd('C:\\Users\\maxim\\OneDrive - Montefiore Medicine\\PROseq\\gtf\\hg38')

#read in bed file
hg38<-rtracklayer::import('hg38.ncbiRefSeq.gtf')

df <- as.data.frame(hg38)
chromosomes <- c('1', '2', '3',  '4',  '5',  '6',  '7',  'X',  '8',  '9',  '11', '10', '12', '13', '14', '15', '16', '17', '18', '20', '19', 'Y', '22', '21')
df <-df[df$seqnames %in% chromosomes,]
df <- df[df$type == 'transcript',]



#remove duplicates
df <- df[!duplicated(cbind(df$seqnames,df$start, df$end)),]

#create a new bed file with the TSS of each gene
for (i in 1:nrow(df)){
  if (df[i,5] == '-') {
    df[i,2] <- df[i,3] - 500
    df[i,3] <- df[i,3] + 500
  } else {
    df[i,3] <- df[i,2] + 500
    df[i,2] <- df[i,2] - 500
  }
}

#keep only the first 4 columns of the bed file
hg38_tss_bed <- df[c(1:3)]


# write the TSS region of each gene to a file
write.table(hg38_tss_bed, 'hg38_refseq_tss.bed', sep = '\t', row.names = F, col.names = F, quote = F)


##Below commands to also annotate TES, however, will not use them to maximize number of short genes that are included
#df <- read.table('hg38_refseq_all_genes.bed', sep = '\t', header = F)
#df <- df[!duplicated(cbind(df$V1,df$V2, df$V3)),]
#for (i in 1:nrow(df)){
#  if (df[i,6] == '-') {
#    df[i,3] <- df[i,2] + 500
#    df[i,2] <- df[i,2] - 500
#  } else {
#    df[i,2] <- df[i,3] - 500
#    df[i,3] <- df[i,3] + 500
#  }
#}

#hg38_tes_bed <- df[c(1:4)]

#hg38_tes_bed <- hg38_tes_bed[!duplicated(cbind(hg38_tes_bed$V1,hg38_tes_bed$V2, hg38_tes_bed$V3)),]

#write.table(hg38_tes_bed, 'hg38_tes.bed', sep = '\t', row.names = F, col.names = F, quote = F)


#hg38_tss_tes_bed <- rbind(hg38_tss_bed, hg38_tes_bed)

#hg38_tss_tes_bed$V1<-gsub("chr","",as.character(hg38_tss_tes_bed$V1))

chromosomes <- c('1', '2', '3',  '4',  '5',  '6',  '7',  'X',  '8',  '9',  '11', '10', '12', '13', '14', '15', '16', '17', '18', '20', '19', 'Y', '22', '21')
hg38_tss_tes_bed <-hg38_tss_tes_bed[hg38_tss_tes_bed$V1 %in% chromosomes,]

#write.table(hg38_tss_tes_bed, 'hg38_tss_tes.bed', sep = '\t', row.names = F, col.names = F, quote = F)

#Take this bed file and use bedtools subtract to remove these intervals from the desired GTF file, then print gene level features to use for donwstream differential expression analysis

# bedtools subtract -a gtf/hg38/hg38.ncbiRefSeq.gtf -b bed/hg38/hg38_tss.bed  > gtf/hg38/hg38_ncbiRefSeq_noTSS.gtf

#Then go back in R and filter the GTF file for gene only types as below


library(rtracklayer)
library(devtools)
library(GenomicRanges)
setwd('C:\\Users\\maxim\\OneDrive - Montefiore Medicine\\PROseq\\gtf\\hg38')

#hg38 <-rtracklayer::import('hg38.ncbiRefSeq.gtf')
#df <- as.data.frame(hg38)
#df$seqnames<-gsub("chr","",as.character(df$seqnames))
#rtracklayer::export(df, 'hg38.ncbiRefSeq.gtf')

hg38<-rtracklayer::import('hg38_ncbiRefSeq_noTSS.gtf')
df <- as.data.frame(hg38)
df_gene <- df[df$type == 'transcript',]
chromosomes <- c('1', '2', '3',  '4',  '5',  '6',  '7',  'X',  '8',  '9',  '11', '10', '12', '13', '14', '15', '16', '17', '18', '20', '19', 'Y', '22', '21')
df_gene <-df_gene[df_gene$seqnames %in% chromosomes,]
df_gene <- df_gene[!duplicated(cbind(df_gene$seqnames,df_gene$start,df_gene$end)),]


#rtracklayer::export(df_gene, 'hg38_ncbiRefSeq_noTSS.gtf')
#test <- rtracklayer::import('2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_noTSS_gene.gtf')
#df<-as.data.frame(test)

#df <- df[c(1:4)]
#write.table(df, '2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_noTSS_gene.bed', sep = '\t', row.names = F, col.names = F, quote = F)


