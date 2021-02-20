library(ggplot2)
library(viridis)


setwd("")

cell_by_gene_df = read.table("gene_neighborhood_mean_exp_fold_change_3.0_um_5_max_neighbors_percentile_99.0.txt", header = TRUE, sep = "\t",row.names = NULL)

gene_names = cell_by_gene_df$Row
cell_by_gene_df$Row <- NULL 

row.names(cell_by_gene_df) <- gene_names

enrichment_matrix <- as.matrix(cell_by_gene_df)
correlation_mat <- cor(enrichment_matrix,method = 'pearson')

pdf("neighborhood_correlation_matrix.pdf",width=24,height=24)
out <- heatmap.2(correlation_mat,
                 Rowv = T,Colv=T,
                 distfun = dist,hclustfun = hclust,
                 symm=FALSE, symbreaks = T, symkey=T,scale="none",dendrogram="both",
                 col= viridis(200),
                 density.info = "none",
                 trace = "none",
                 sepcolor="black",
                 sepwidth=c(5,5),
                 margins = c(5, 5),
                 key =T
)
dev.off()