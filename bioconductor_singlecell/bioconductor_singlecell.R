# renv::restore()

library(tidyverse)
library(scater)
library(scran)
library(DropletUtils)
library(DelayedMatrixStats)
library(DelayedArray)
library(Matrix)
# BiocManager::install('sparseMatrixStats')
# BiocManager::install('iSEE')
library(sparseMatrixStats)
library(iSEE)


#Loading data and creating single cell experiment object
sce <- read10xCounts(
    c(filtered="/Users/andresnoe/obds_sep20/working_directory/bioc/filtered_feature_bc_matrix"),
    col.names = TRUE)

# Exploring single cell experiment object (sce)
sce
dim(sce)
str(sce)
slotNames(sce)
metadata(sce) # Miscellaneous list of extra metadata

colnames(colData(sce)) # Cell metadata
colData(sce)$Barcode

colnames(rowData(sce)) # Gene metadata
rowData(sce)

assayNames(sce)
assay(sce)


# QC -------
# Per cell QC ----
# THIS DATA SET IS ALREADY FILTERED

barcode <- barcodeRanks(assay(sce, "counts"))
# Be explicit about which assay you're referring to with "counts" in the assay function
barcode
ggplot(data.frame(barcode), aes(x=rank, y=total, col=fitted))+
    geom_point()+
    scale_x_log10()+
    scale_y_log10()
# THIS DATA SET IS ALREADY FILTERED, so there is a steep drop off


# THIS DATA SET IS ALREADY FILTERED
per_cellQC <- perCellQCMetrics(sce, exprs_values="counts")
per_cellQC
# Values:
# sum: numeric, the sum of counts for each cell.
# detected: numeric, the number of GENES (observations) above detection_limit.
# # should be less than sum because some genes will have multiple counts to them
# # detection_limit is defined in the function (a limit of 1 is usually the lower
# # limit of detection for single cell rna-seq)
# percent_top: numeric matrix, the percentage of counts assigned to the percent_topage of most highly expressed genes. Each column of the matrix corresponds to an entry of the sorted percent_top, in increasing order.
# total: numeric, the total sum of counts for each cell across main and alternative Experiments.

# THIS DATA SET IS ALREADY FILTERED
# BUT: There are a few cells that aren't very clean (look at plot 2)

ggplot(data.frame(per_cellQC), aes(x=percent_top_50))+
    geom_histogram(bins=100, colour="black")

ggplot(data.frame(per_cellQC), aes(x=percent_top_50, y=total))+
    geom_point()

ggplot(data.frame(per_cellQC))+
    geom_col(aes(y=sum,x=rownames(data.frame(per_cellQC))), bins=100, colour="black") +
    geom_col(aes(y=detected, x=rownames(data.frame(per_cellQC))), bins=100, colour="red")
rownames(data.frame(per_cellQC))

per_cellQC %>% as.data.frame %>%
    pivot_longer(-c(sum, detected, total)) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    ggplot(aes(x=value, y = total)) +
    geom_point() +
    facet_wrap(~name)


#Add the perCellQC metric to the sce
sce <- addPerCellQC(sce, exprs_values="counts")


# Per feature QC ----
per_featureQC <- perFeatureQCMetrics(sce)
per_featureQC
ggplot(data.frame(per_featureQC), aes(detected))+
    geom_histogram(bins=1000)

# mean: numeric, the mean counts for each feature.
# detected: numeric, the percentage of observations above detection_limit.

sce <- addPerFeatureQC(sce)
rowData(sce)


# Convert counts to normalised counts
sce <- logNormCounts(sce)
range(assay(sce,"logcounts")) # 0.00000 12.68275

mean_var_log_counts <- data.frame(
    mean=rowMeans(assay(sce,"logcounts")),
    variance=sparseMatrixStats::rowVars(assay(sce,"logcounts")))

ggplot(data.frame(mean_var_log_counts),
       aes(x=mean,y=variance))+
    geom_point()

mean_var_counts <- data.frame(
    mean=rowMeans(log(assay(sce,"counts")+1)),
    variance=sparseMatrixStats::rowVars(log(as.matrix(assay(sce,"counts")+1))))

ggplot(data.frame(mean_var_counts),
       aes(x=mean,y=variance))+
    geom_point()

# Model the variance of the log-expression profiles for each gene,
# decomposing it into technical and biological components based on
# a fitted mean-variance trend.
decomposed_log <- modelGeneVar(sce)
decomposed_log

plot(decomposed_log$mean, decomposed_log$total)
curve(metadata(decomposed_log)$trend(x), add=TRUE, col="dodgerblue")
points(metadata(decomposed_log)$mean, metadata(decomposed_log)$var, col="red", pch=16)
# Each dot in the plot is a gene
# Line = trend of total variation
# HVGs the top most dots are driving heterogeneity but want to grab enough genes
# However, we want to grab ENOUGH of these genes,
# that are above technical variation (represented here by the blue line)

# Other plots
# 1. Library size on x vs. raw counts on for given gene where you should see before normalisation that the counts of that gene are affected by library size. Then the same plot after normalisation where the counts of the gene should not be correlated with library size anymore.
# 2. The other plot is the mean on x and variance on y for normalised counts (not logged) and normalised counts (logged). Log should help to deal with the mean-variance relationship


# Select features for downstream analyses, e.g. highly variable genes; use scran.
HVG <- getTopHVGs(decomposed_log,
           var.field = "bio",
           1000)
rowData(sce)$HVG <- rownames(sce) %in% HVG
rowData(sce)

# Apply dimensionality reduction to compact the data and further reduce noise; use scater.----
# PCA
set.seed(123)
sce <- runPCA(sce,
              name='RNA_PCA',
              ntop=Inf, # Use all the HVGs I give you
              subset_row=HVG)
str(reducedDims(sce)$RNA_PCA)
colData(sce)
plotReducedDim(sce,
               dimred = "RNA_PCA",
               ncomponents = 2,
               colour_by="percent_top_50")

# Scree plot to choose PCs for UMAP

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")

sce <- runUMAP(sce,
        name="RNA_UMAP",
        dimred="RNA_PCA",
        n_dimred=30)
# 30 Looks good on the scree plot

plotReducedDim(sce,
               dimred = "RNA_UMAP",
               ncomponents = 2,
               colour_by="percent_top_50")

# Clustering -----
# Cluster cells; use scran.
colData(sce)$quickCluster <- quickCluster(sce, assay.type="logcounts")
plotReducedDim(sce,
               dimred = "RNA_UMAP",
               ncomponents = 2,
               colour_by="quickCluster")


# Using SNN graph - suggested method
SNNGraph <- buildSNNGraph(sce, use.dimred = 'RNA_PCA')
colData(sce)[["cluster_louvain"]] <- factor(igraph::cluster_louvain(SNNGraph)$membership)
colData(sce)[["cluster_louvain"]]

plotReducedDim(sce,
               dimred = "RNA_UMAP",
               ncomponents = 2,
               colour_by="cluster_louvain",
               shape_by = "quickCluster",
               size_by = "quickCluster")

# Identify cluster markers ------
# Identify markers for each cluster; use scran.
markers <- findMarkers(sce,
            assay.type="logcounts",
            groups=sce$cluster_louvain)
markers[[1]] # The result of comparing cluster 1 to each other cluster


o
#USING iSEE to explore data
# # Really only used for exploration and visualisation, not analysis
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)
rownames(sce) # All rownames should be gene symbols unless NA,
# in which case it is a gene ID
iSEE(sce)
# Can export code required to make plots. WICKED!


# https://www.nature.com/articles/s41587-019-0071-9?platform=hootsuite
# Pseudotime analysis review

