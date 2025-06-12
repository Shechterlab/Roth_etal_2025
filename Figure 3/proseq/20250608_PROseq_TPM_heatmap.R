library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
df <- read.table('E:/PRMT\\proseq\\peppro/hg38/20250608_GSK591_PROseq_TPM.txt', header = T, row.names = NULL)

df.m <- df[,1:ncol(df)-1]

rownames(df.m) <- df$gene_name
df.m <- as.matrix(df.m)

transposed_matrix <- t(df.m)

z_tr_mt <- scale(transposed_matrix)

my_matrix <- t(z_tr_mt)

#colors <- brewer.pal(3, 'RdYlBu')

ht_opt("heatmap_row_names_gp" = gpar(fontsize = 4), 'heatmap_column_names_gp' = gpar(fontsize = 8))
col_fun = colorRamp2(c(-2,0,2), c('#55439B', 'grey90', '#E76F00'))

hmap <- Heatmap(my_matrix, name = 'z-score', col = col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, row_names_side = 'left', show_row_names = FALSE, row_km = 4, column_km = NULL)

pdf('proseq_heatmap_zscore_GSK591.pdf')
hmap <- draw(hmap)
dev.off()
