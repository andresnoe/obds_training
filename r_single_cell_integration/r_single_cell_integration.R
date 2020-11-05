renv::restore()
BiocManager::install("batchelor")
BiocManager::install("Seurat")
library(tidyverse)
library(batchelor)
library(scran)
library(Seurat)
library(scater)
library(DropletUtils)

# https://osca.bioconductor.org/integrating-datasets.html#linear-regression
# Load in data/create sce -----

#Loading data and creating single cell experiment object
sce <- read10xCounts(
    c(v2="/Users/andresnoe/obds_sep20/working_directory/pbmc_1k_v2",
      v3="/Users/andresnoe/obds_sep20/working_directory/pbmc_1k_v3"),
    col.names = TRUE)


# Exploring single cell experiment object (sce)
sce
dim(sce)
str(sce)
slotNames(sce)
metadata(sce) # Miscellaneous list of extra metadata

colnames(colData(sce)) # Cell metadata
colData(sce)

colnames(rowData(sce)) # Gene metadata
rowData(sce)

assayNames(sce)
assay(sce)

# QC -----
is.mito <- grepl("^MT-", rowData(sce)$Symbol)
table(is.mito)
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
# quickPerCellQC uses the isOutlier function
# # Convenience function to determine which values in a numeric vector are
# # outliers based on the median absolute deviation (MAD).

filtered
table(filtered)
sce <- sce[, !filtered$discard]


# Normalisation post-filtering ----

sce <- logNormCounts(sce)
range(assay(sce,"logcounts")) # 0.00000 11.28685

decomposed_var <- modelGeneVar(sce)
HVG <- getTopHVGs(decomposed_var, prop=0.1)

# Quick exploration -----
# These datasets are pre-filtered

set.seed(123)
sce <- runPCA(sce,
              name='PCA',
              ntop=Inf, # Use all the HVGs I give you
              subset_row=HVG)
str(reducedDims(sce)$RNA_PCA)
colData(sce)
plotReducedDim(sce,
               dimred = "PCA",
               ncomponents = 2,
               colour_by = "Sample")
# Can we decide whether a linear batch correction is sufficient just off PCA?
# Probably not

# Scree plot to choose PCs for UMAP

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")

sce <- runUMAP(sce,
               name="UMAP",
               dimred="PCA",
               n_dimred=20)
# 20 Looks good on the scree plot

plotReducedDim(sce,
               dimred = "UMAP",
               ncomponents = 2)

# Integrate v2 and v3 datasets -----
# We use the rescaleBatches() function from the batchelor package to remove
# the batch effect. This is roughly equivalent to applying a linear regression
# to the log-expression values per gene, with some adjustments to improve
# performance and efficiency. For each gene, the mean expression in each batch
# is scaled down until it is equal to the lowest mean across all batches.

# Split v2 and v3 datasets

sce_v2 <- sce[,sce$Sample=="v2"]
dim(sce_v2)
sce_v3 <- sce[,sce$Sample=="v3"]
dim(sce_v3)

# LINEAR CORRECTION

rescaled <- rescaleBatches(sce_v2, sce_v3)
rescaled
dim(rescaled)
colData(rescaled)

set.seed(123) # To ensure reproducibility of IRLBA. NOTE: USING 123
rescaled <- runPCA(rescaled,
                   subset_row=HVG,
                   exprs_values="corrected")

str(reducedDim(rescaled))
plotReducedDim(rescaled,
               dimred = "PCA",
               ncomponents = 2,
               colour_by = "batch")
# In this case, the linear batch correction method employed in rescaleBatches
# was insufficient to fully correct for the differences across all 'clusters'
# (eg the right-most cluster)
# However, there is no linear method to FULLY correct for all differences among
# batches, so this isn't a terrible outcome.

# Correction diagnostics

# One way of checking batch correction is to use SNN
# If integration works well, then cells from all batches should be present in
# all clusters
snn.gr <- buildSNNGraph(rescaled, use.dimred="PCA")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
colData(rescaled)$cluster_walktrap_SNN <- as.factor(clusters.resc)
colData(rescaled)$batch <- as.factor(colData(rescaled)$batch)
tab.resc <- table(Cluster=clusters.resc, Batch=rescaled$batch)
tab.resc

plotReducedDim(rescaled,
               dimred = "PCA",
               ncomponents = 2,
               colour_by = "cluster_walktrap_SNN",
               shape_by = "batch")

rescaled <- runUMAP(rescaled,
               name="UMAP",
               dimred="PCA",
               n_dimred=20)
# Haven't checked scree plot, assume 20 is fine from previous check
library(patchwork)
p1 <- plotReducedDim(rescaled,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by = "cluster_walktrap_SNN")
p2 <- plotReducedDim(rescaled,
                     dimred = "UMAP",
                     ncomponents = 2,
                     colour_by = "batch")

p1+p2


# MNN CORRECTION

set.seed(123)
mnn.out <- fastMNN(sce_v2, sce_v3,
                   d=50, k=20,
                   subset.row=HVG,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
str(BiocSingular::RandomParam(deferred=TRUE))
# BiocSingular::RandomParam(deferred=TRUE) is a complex object that contains
# a bunch of random parameters. This object specifies the algorithm to use for PCA.
# This uses a fast approximate algorithm from irlba by default
# k = 	An integer scalar specifying the number of nearest
# neighbors to consider when identifying MNNs. The k parameter also depends
# on the size of the clusters you are expecting. If desired cell cluster is very
# small, then k parameter will need to be small (as number of neighbours will be smaller)
# d = Numeric scalar specifying the number of dimensions to use for
# dimensionality reduction in multiBatchPCA. If NA, no dimensionality
# reduction is performed and any matrices in ... are used as-is.

mnn.out
# Contains a reconstructed matrix in the assays slot, containing the low-rank reconstruction of the expression matrix. This can be interpreted as per-gene corrected log-expression values (after cosine normalization, if cos.norm=TRUE) but should not be used for quantitative analyses. This has number of rows equal to the number of input genes if subset.row=NULL or correct.all=TRUE, otherwise each row corresponds to a gene in subset.row.
# Can get back corrected matrix for all genes if specify the argument
# 'correct.all' = TRUE


# The corrected matrix in the reducedDims() contains the low-dimensional corrected coordinates for all cells, which we will use in place of the PCs in our downstream analyses.
dim(reducedDim(mnn.out, "corrected"))

# A reconstructed matrix in the assays() contains the corrected
# expression values for each gene in each cell, obtained by projecting
# the low-dimensional coordinates in corrected back into gene expression space.
# IT IS NOT RECOMMENDED to use this for anything other than visualization.
assay(mnn.out, "reconstructed")

mnn.out <- runUMAP(mnn.out,
                    name="UMAP",
                    dimred="corrected",
                    n_dimred=20)
mnn.out$batch <- as.factor(mnn.out$batch)
p3 <- plotReducedDim(mnn.out,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by = "batch")
p2+p3

# Correction diagnostics

snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.wt.mnn <- igraph::cluster_walktrap(snn.gr)$membership
colData(mnn.out)$cluster_walktrap_SNN <- as.factor(clusters.wt.mnn)
tab.wt.mnn <- table(Cluster=clusters.wt.mnn, Batch=mnn.out$batch)
tab.wt.mnn


clusters.louvain.mnn <- igraph::cluster_louvain(snn.gr)$membership
colData(mnn.out)$cluster_louvain_SNN <- as.factor(clusters.louvain.mnn)
tab.louvain.mnn <- table(Cluster=clusters.louvain.mnn, Batch=mnn.out$batch)
tab.louvain.mnn

p4 <- plotReducedDim(mnn.out,
                     dimred = "UMAP",
                     ncomponents = 2,
                     colour_by = "cluster_walktrap_SNN")

(p1+p2)/(p4+p3)

p5 <- plotReducedDim(mnn.out,
                     dimred = "UMAP",
                     ncomponents = 2,
                     colour_by = "cluster_louvain_SNN")

# Compare both clustering methods (just for fun)
tab.cluster <- table(Walktrap=mnn.out$cluster_walktrap_SNN,
                     Louvain=mnn.out$cluster_louvain_SNN)
tab.cluster

heatmap(log1p(tab.cluster))



# For fastMNN(), one useful diagnostic is the proportion of variance within
# each batch that is lost during MNN correction.
# Specifically, this refers to the within-batch variance that is
# removed during orthogonalization with respect to the average
# correction vector at each merge step.
metadata(mnn.out)$merge.info$lost.var
# Large proportions of lost variance (>10%) suggest that correction
# is removing genuine biological heterogeneity.
# This would occur due to violations of the assumption of orthogonality
# between the batch effect and the biological subspace (Haghverdi et al. 2018).
# In this case, the proportion of lost variance is small,
# indicating that non-orthogonality is not a major concern.



# ALTERNATIVE TO USING HVG for CORRECTION ------------------------------
# 13.6.2 Encouraging consistency with marker genes
# https://osca.bioconductor.org/integrating-datasets.html

# In some situations, we will already have performed within-batch
# analyses to characterize salient aspects of population heterogeneity.
# This is not uncommon when merging datasets from different sources where
# each dataset has already been analyzed, annotated and interpreted separately.
# It is subsequently desirable for the integration procedure to retain these
# “known interesting” aspects of each dataset in the merged dataset.
# We can encourage this outcome by using the marker genes within each
# dataset as our selected feature set for fastMNN() and related methods.
# This focuses on the relevant heterogeneity and represents a semi-supervised
# approach that is a natural extension of the strategy described in Section 8.4.


# It is preferable to perform DE analyses using the uncorrected
# expression values with blocking on the batch (holding batch variability constant)
# Looks for cluster variability, while holding batch variability constant
m.out <- findMarkers(sce,
                     groups=clusters.louvain.mnn,
                     block=sce$Sample,
                     direction="up",
                     test.type="t", # t-test comparing cluster 1 to the average of all other clusters
                     lfc=1,
                     row.data=rowData(sce)[,3,drop=FALSE])
# direction: A string specifying the direction of log-fold changes
# to be considered in the alternative hypothesis.

# This strategy is based on the expectation that any genuine DE between
# clusters should still be present in a within-batch comparison where batch effects are absent.
# It penalizes genes that exhibit inconsistent DE across batches, thus protecting
# against misleading conclusions when a population in one batch is aligned
# to a similar-but-not-identical population in another batch
str(clusters.mn)

m.out[[1]]
# For each gene and cluster, the summary effect size is defined as the effect size from the pairwise comparison with the largest p-value. This reflects the fact that, with this approach, a gene is only as significant as its weakest DE. Again, this value is not directly used for ranking and are only reported for the sake of the user.

# Subset rownames sce to get the symbols of the markers
m.out[[2]]$Symbol <- rowData(sce)[rownames(m.out[[2]]), "Symbol"]
View(as.data.frame(m.out[[2]]))
