# Load data

log_counts <- read.csv("/Users/andresnoe/obds_sep20/obds_training/r_pca_and_clustering/data/logcounts.csv", row.names = 1)
# View(log_counts)

cell_md <- read.csv("/Users/andresnoe/obds_sep20/obds_training/r_pca_and_clustering/data/cell_metadata.csv", row.names = 1)
# View(cell_md)

gene_md <- read.csv("/Users/andresnoe/obds_sep20/obds_training/r_pca_and_clustering/data/gene_metadata.csv", row.names = 1)
# View(gene_md)
