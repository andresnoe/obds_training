renv::restore()

library(tidyverse)
library(Seurat)
library(patchwork)
library(clustree)

# features <- read_tsv("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/features.tsv.gz")
# matrix <- read_tsv("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/matrix.mtx.gz", comment = "%")
# barcodes <- read_tsv("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/barcodes.tsv.gz")
# tail(features)
# head(matrix)
# head(barcodes)


# scRNA-SEQ ANALYSES
# https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html


# Give folder with the:
# 1) barcodes.tsv.gz
# 2) matrix.mtx.gz
# 3) features.tsv.gz
pbmc.data <- Read10X("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/")

pbmc.data
View(pbmc.data)
class(pbmc.data) #This is a list 
summary(pbmc.data$`Antibody Capture`)
rownames(pbmc.data$`Antibody Capture`)
colnames(pbmc.data$`Antibody Capture`)

# Let's now create the Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`,
                               project = "pbmc3k")
pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc.data$`Antibody Capture`)

Assays(pbmc)
DefaultAssay(pbmc) #RNA is the default already, don't need to change it yet

# Initial quality control (just looking)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Store experimental variables on a per cell basis
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# grep for any genes that start with "MT-".
# Can also do for any other genes that contain a string in each cell
pbmc[["percent.mt"]]

pbmc[["percent.rib"]] <- PercentageFeatureSet(pbmc, pattern = "^RPS|^RPL")
pbmc[["percent.rib"]]

View(pbmc[[]])

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA")) +
    geom_hline(yintercept=1250)
# (total number of RNAs)

VlnPlot(pbmc, features = c("nCount_RNA"))
# (total number of counts)
VlnPlot(pbmc, features = c("percent.mt")) +
    geom_hline(yintercept=14)

VlnPlot(pbmc, features = c("percent.rib"))

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC plot of nCount_RNA by nFeature_RNA coloured by percent.mt (use ggplot2)
ggplot(pbmc[[]], aes(x=nCount_RNA, y=nFeature_RNA, col=percent.mt))+
    geom_point()
# Suggestion: if just looking at this plot, can cut off at nFeature_RNA = 6000
# But this is an iterative process

# Cutoffs:
# # 1000 for nFeature_RNA 
# # 12.5 for percent.mt
# # Leave nCount_RNA for now
# # Leave percent.rib for now (but some people say that if it's above 50% in PBMCs then they should be filtered out, but be careful)

# Not yet filtered --------- 
# WILL NOT FILTER YET - get up to clustering step and look for cluster of dead cells

# Remember that SCTransform performs these three functions:
# NormalizeData
# FindVariableFeatures
# ScaleData

pbmc <- SCTransform(pbmc,
                    assay = 'RNA',
                    seed.use = 1448145,
                    verbose = TRUE)
# This produces a new assay which has all the sc transform stuff
# Look at documentation of SCTransform for defaults, particularly:
# new.assay.name = "SCT"
# variable.features.n = 3000
# return.only.var.genes = TRUE - only returns genes that are HVGs
# # Might be important if gene OUTSIDE of this HVG list is diff. exp.
# vars.to.regress = NULL (can regress out percent.mt if wanted)
# # BUT do not try to regress out BIG batch effects as this is only meant for small-y effecting variables

# WARNING: Don't worry about this one "1: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached" - seems to work fine

VariableFeatures(pbmc) # find the HVGs and extract them as characters.
# identifies which genes have been pulled by SCTransform


# Perform linear dimensional reduction
DefaultAssay(pbmc) # SCTransform automatically changes default assay to SCT

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc, reduction = "pca", group.by = "percent.mt") +
    NoLegend()

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# Print the PC loadings: genes that most contribute to PC

VizDimLoadings(pbmc, dims = 1:3, reduction = "pca")
# Plots above print statement essentially
# balanced = FALSE by default
# "balanced"  = FALSE: Return an equal number of genes with + and - scores. If FALSE (default), returns the top genes ranked by the scores absolute values


DimHeatmap(pbmc, dims = 1, cells = 2000, balanced = TRUE)
# DimHeatmap can give an indication of which PCs we should
# include in the dim red (similar to Scree plot)

# JackStraw cannot be run on SCTransform-normalized data.
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

ElbowPlot(pbmc, ndims = 50)
# Using SCTransform can include a few more PCs than otherwise, so err on side of increasing
# In this case we will choose between 20-30 (25)



# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20)
# dims = 25 was decided in the elbow plot above
# k.param play around with this (default is 20), depending on size of clusters
# annoy.metric: Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
# # annoy = approximate  nearest neighbours (oh yeah) is another R package
# # Just a distance metric for measuring distance between cells

pbmc <- FindClusters(pbmc, resolution = c(0.1,0.2,0.3,0.5, 0.8, 1, 1.3, 1.5, 2))
# From  vignette: We find that setting this parameter [resolution]
# between 0.4-1.2 typically returns good results for single-cell
# datasets of around 3K cells
# # Affected by the NUMBER OF CELLS you have
# # Usually use cluster tree to determine more definitively, but this is just a start
View(pbmc[[]])
Idents(pbmc) <- "SCT_snn_res.0.8"
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc,cells.highlight = WhichCells(pbmc, expression = nCount_RNA > 11000))
DimPlot(pbmc,cells.highlight = WhichCells(pbmc, expression = nFeature_RNA < 1000))
DimPlot(pbmc,cells.highlight = WhichCells(pbmc, expression = percent.mt > 12.5))
# The thresholds look quite reasonable at this stage

ggplot(pbmc[[]], aes(x=nFeature_RNA))+
    geom_histogram(bins = 100)

pbmc <- subset(pbmc,
               subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & percent.mt < 12.5)

# Filtered data --------- 
# Below I have copied the code from above (in section titled "not yet filtering")
pbmc <- SCTransform(pbmc,
                    assay = 'RNA',
                    seed.use = 1448145,
                    verbose = TRUE)

VariableFeatures(pbmc) 


# Perform linear dimensional reduction
DefaultAssay(pbmc)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc, reduction = "pca", group.by = "percent.mt") +
    NoLegend()

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(pbmc, dims = 1:3, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 2000, balanced = TRUE)


ElbowPlot(pbmc, ndims = 50)

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20)


pbmc <- FindClusters(pbmc, resolution = c(0.1,0.2,0.3,0.5, 0.8, 1, 1.3, 1.5, 2))

View(pbmc[[]])
Idents(pbmc) <- "SCT_snn_res.0.8"
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")
# B cells bottom left?
# Top monocytes? Lucy's career is on the line

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))


#  Post-filtering steps ---------

# Let's find out how many clusters we should have
clustree(pbmc)
# Let's choose a resolution of 0.8
# But would usually go through and look at markers of the clusters

Idents(pbmc) <- "SCT_snn_res.0.8"


# MULTI-MODAL ANALYSES ------
# https://satijalab.org/seurat/v3.2/multimodal_vignette.html

# Now we can repeat the preprocessing (normalization and scaling) steps that we typically run
# with RNA, but modifying the 'assay' argument.  For CITE-seq data, we do not recommend typical
# LogNormalization. Instead, we use a centered log-ratio (CLR) normalization, computed
# independently for each feature.  This is a slightly improved procedure from the original look up DBS
pbmc <- NormalizeData(pbmc,
                      assay = "ADT",
                      normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay = "ADT")

DefaultAssay(pbmc) <- "RNA"
rownames(pbmc)

# https://www.biorxiv.org/content/10.1101/2020.02.24.963603v1
# DSB normalisation method

FeaturePlot(pbmc,
            features = c("CD3-TotalSeqB", "CD4-TotalSeqB", "CD8a-TotalSeqB", "CD11b-TotalSeqB", "CD14-TotalSeqB", "CD15-TotalSeqB", "CD16-TotalSeqB", "CD19-TotalSeqB", "CD4", "CD3"),
            min.cutoff = "q05",
            max.cutoff = "q95")

RidgePlot(pbmc, features = c("adt_CD3-TotalSeqB", "adt_CD8a-TotalSeqB", "adt_CD45RA-TotalSeqB", "adt_CD45RO-TotalSeqB"), ncol = 4)

isotype_controls <- c("IgG1-control-TotalSeqB", "IgG2a-control-TotalSeqB", "IgG2b-control-TotalSeqB")
DefaultAssay(pbmc) <- "ADT"
pbmc <- RunPCA(pbmc,
               features = rownames(pbmc)[!rownames(pbmc) %in% isotype_controls],
               reduction.name = "pca_adt",
               reduction.key = "pca_adt_",
               npcs = 15,
               verbose = FALSE)
Reductions(pbmc)
# Shows the dimensional reductions that you have
DimPlot(pbmc, reduction = "pca_adt")

ElbowPlot(pbmc, ndims=15)
cbmc <- RunUMAP(pbmc,
                dims = 1:12,
                reduction = "pca_adt",
                reduction.key = "adtTSNE_",
                reduction.name = "tsne_adt")
pbmc <- FindNeighbors(pbmc, features = rownames(cbmc), dims = NULL)
pbmc <- FindClusters(pbmc, resolution = c(0.2,0.4,0.6,0.8,1,1.2,1.4), graph.name = "ADT_snn")

# We can compare the RNA and protein clustering
# and use this to annotate the protein clustering
# (we could also of course use FindMarkers)
DefaultAssay(pbmc) <- "ADT"
Idents(pbmc) <- "ADT_snn_res.0.8"
clustering.table <- table(Idents(pbmc), pbmc@meta.data$SCT_snn_res.0.8)
clustering.table
