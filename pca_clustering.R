library(tidyverse)

# Load data

log_counts <- read.csv("/Users/andresnoe/obds_sep20/obds_training/r_pca_and_clustering/data/logcounts.csv", row.names = 1)
# View(log_counts)

cell_md <- read.csv("/Users/andresnoe/obds_sep20/obds_training/r_pca_and_clustering/data/cell_metadata.csv", row.names = 1)
# View(cell_md)

gene_md <- read.csv("/Users/andresnoe/obds_sep20/obds_training/r_pca_and_clustering/data/gene_metadata.csv", row.names = 1)
# View(gene_md)


# Convert data frame to matrix
class(log_counts)
log_count_matrix <- as.matrix(log_counts)

# Good idea to check dimensions of matrix before PCA to estimate computation time
str(log_count_matrix)

# ========================================================================
# Perform PCA. How many principal components do you think you should keep for
# follow up analysis?
count_pca_scaled <- prcomp(t(log_count_matrix), center = TRUE, scale. = TRUE)
# EACH SAMPLE IN THE PCA MUST BE A ROW (therefore we have transposed the matrix)
summary(count_pca_scaled)
str(count_pca_scaled)
View(count_pca_scaled$x)

count_pca_scaled$x

count_pca_not_scaled <- prcomp(t(log_count_matrix), center = TRUE, scale. = FALSE)
summary(count_pca_not_scaled)
str(count_pca_not_scaled)

# As general rule, scaling is recommended for PCA
#Plot with ggplot
pca_df <- as.data.frame(count_pca_scaled$x)
pca_df %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point()

# Screeplot to choose how many PCs - take about 10-12 PCs
screeplot(count_pca_scaled, npcs = 25, type="lines") # base R version

screeplot_df <- data.frame(
    'var' = (count_pca_scaled$sdev)^2,
    'PC' = paste("PC",seq_along(count_pca_scaled$sdev),
                 sep=""))
screeplot_df$PC <- factor(screeplot_df$PC,screeplot_df$PC)

screeplot_df <- screeplot_df %>%
    mutate(percent_var = 100*(screeplot_df$var/sum(screeplot_df$var)),
           cumulative_var = cumsum(percent_var))
head(screeplot_df)

screeplot_df[1:20,] %>%
    ggplot(aes(x = PC,  y = percent_var)) +
    geom_point()
