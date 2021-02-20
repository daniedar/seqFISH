cell_by_gene_umap_annotated.txt 
 - mRNA counts per gene in each cell analyzed in the umap
 - conditon in LB growth curve
 - umap cluster # per cell (starting at 0 instead of 1 as in table S4).

spatial_params_and_umap_table.txt
 - spatial params per cell (cell length, DAPI and 16S rRNA signal).

Files for running scanpy:
 - matrix.mtx
 - barcodes
 - genes
 - LB_scanpy_analysis.ipynb

umap_cluster_expression_wilcoxon.txt
 - differential expression output for all clusters and genes with log fold-change and p-values
