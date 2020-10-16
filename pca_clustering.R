library(tidyverse)
install.packages("ggthemes")
library(ggthemes)
library(cowplot)

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


# Use the experimental metadata to visualise which cell types visually cluster
# together in principal component space.
pca_df <- as.data.frame(count_pca_scaled$x)
pca_df <- rownames_to_column(pca_df)
cell_md <- rownames_to_column(cell_md)
pca_df <- pca_df %>%
    left_join(cell_md)

str(pca_df)


pca_df %>%
    ggplot(aes(x = PC1, y = PC2, color = rowname)) +
    geom_point()


# https://dplyr.tidyverse.org/reference/tidyeval.html
for (metadata in factor(colnames(cell_md)[-1])){
    p <- pca_df %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = metadata)) +
        geom_point()
    print(p)
}

ggplot(pca_df) +
    geom_density(aes(PC1, fill = Status), color = "black", alpha = 0.5) +
    facet_grid(Time~Infection) +
    theme_cowplot()


# Find the top genes associated with the top principal components.
# Visualise gene expression values against PC coordinates.
View(count_pca_scaled$rotation)
top_genes_PC1 <- count_pca_scaled$rotation[,"PC1"]

top_genes_PC1 <- sort(top_genes_PC1, decreasing = TRUE)

sort(top_genes_PC1, decreasing = TRUE) %>%
    head

sort(top_genes_PC1, decreasing = TRUE) %>%
    tail



top_genes_PC2 <- count_pca_scaled$rotation[,"PC2"]

top_genes_PC2 <- sort(top_genes_PC2, decreasing = TRUE)

sort(top_genes_PC2, decreasing = TRUE) %>%
    head

sort(top_genes_PC2, decreasing = TRUE) %>%
    tail


#  ============================================
#  CLUSTERING #
#  ============================================
kmean <- kmeans(t(log_count_matrix), centers = 4)
kmean
str(kmean)

head(kmeans$cluster, 3)

# View(pca_df)
# View(kmean)

pca_df$cluster <- as.factor(kmean$cluster[pca_df$rowname])
pca_df$cluster

pca_df %>%
    ggplot(aes(x = PC1, y = PC2, color = cluster)) +
    geom_point()

kmean$withinss
kmean$betweenss

candidate_k <- 2:20
km <- sapply(candidate_k, function(i){
    km <- kmeans(t(log_count_matrix), centers = i)
    sum(km$withinss)
})
str(km)

kmeans_plot <- data.frame(sum_of_withinss = km, k=candidate_k)
kmeans_plot
ggplot(kmeans_plot, aes(x=k, y=sum_of_withinss))+
    geom_point()
