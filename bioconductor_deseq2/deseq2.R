#DESeq2 exercises

library(tidyverse)
library(biomaRt)
library(pheatmap)
# install.packages("ggthemes")
# library(ggthemes)
library(DESeq2)

# Read data and metadata in 
sample_table <- read_tsv("data/obds_sampletable.tsv")
counts_table <- read_tsv("data/obds_countstable.tsv.gz")
View(sample_table)
View(counts_table)

# Convert the counts table (obds_countstable.tsv.gz) and the sample information table (obds_sampletable.tsv)
# into a suitable format for generating a DESeqDataSet object
class(counts_table)
counts_table <- column_to_rownames(counts_table,"Geneid")
counts_table <- as.matrix(counts_table)

class(sample_table)
sample_table <- column_to_rownames(sample_table,"Sample_accession")
rownames_sample_table <- sample_table

# Do the colnames of counts_table == to rownames of sample_table ?
table(colnames(counts_table)==rownames(sample_table))

# Set Egr2/3 DKO CD8 cells as the reference level
# Separate up sample_title column
sample_table <- sample_table %>% 
    separate(sample_title, c("egr_locus", "genotype", "cell_type", 'replicate'), sep = "_") %>%
    unite(col = "condition", egr_locus, genotype, cell_type, sep = "_") %>%
    dplyr::select(-c(species,library_layout)) %>%
    mutate(condition=factor(condition, levels=c("Egr2/3_DKO_CD8", "Egr2/3_DKO_CD4", "Egr2_Kin_CD4", "Egr2_Kin_CD8")))
sample_table
levels(sample_table$condition)

# Make DESeq Data Set (dds) from matrix
dds <- DESeqDataSetFromMatrix(counts_table,
                       sample_table,
                       ~ condition)
colData(dds)
rowRanges(dds)
design(dds)
counts(dds) # Is the equivalent to assays(dds)$counts


# # Calculate the size factors for each sample – estimateSizeFactors()
# dds <- estimateSizeFactors(dds)
# sizeFactors(dds)
# 
# # sizefactorsdf <- data.frame(sizeFactors(dds))
# 
# sizefactorsdf <- data.frame(sample = names(sizeFactors(dds)),
#                            size_factor = sizeFactors(dds),
#                            sample_group = colData(dds)$condition)
# 
# 
# ggplot(sizefactorsdf, aes(y=size_factor, x=sample, fill=sample_group))+
#     geom_col()+
#     theme(axis.text.x = element_text(angle=45, hjust=1))
#     
# # Obtain dispersion estimates for each gene – estimateDispersions()
# dds <- estimateDispersions(dds)
# dispersions(dds)
# # Plot the per-gene dispersion estimates
# # (DESeq2 has a helper function for this)
# plotDispEsts(dds) # Looks similar to example data, good quality
# 
# # Perform the Wald test – nbinomWaldTest()
# dds <- nbinomWaldTest(dds)

# Use the DESeq() function to perform steps 5-7 in one go
# DESeq() Runs all three functions (estimateSizeFactors, estimateDispersions, nbinomWalTest) at the same time
dds <- DESeq(dds)

# Access the coefficients of the NB GLM
# NAs may represent independently filtered out genes, or several other reasons
# (outliers,)
head(coef(dds),3)
    

# Access the results table for the comparison
# between CD8+ and CD4+ T cells from Egr2/3 DKO mice
res <- results(dds)
View(res)
# This object contains the results columns:
# baseMean, log2FoldChange, lfcSE, stat, pvalue and padj
# and also includes metadata columns of variable information.

res_df <- as.data.frame(res, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'))
resultsNames(dds)
# res_df <- as.data.frame(res, name = "condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8")

View(res_df)


# Plot a histogram of the raw and BH-adjusted p-values
# do they look as expected?
install.packages("patchwork")
install.packages("cowplot")
library(patchwork)
library(cowplot)

pval_plot <- ggplot(res_df, aes(x=pvalue))+
    geom_histogram()

padj_plot <- ggplot(res_df, aes(x=padj))+
    geom_histogram()

plot_grid(pval_plot, padj_plot, align = "h")

pval_plot + padj_plot

# Generate an MA plot of the log2 FC values for all genes
# https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot
plotMA(res, ylim=c(-8,10))

resultsNames(dds)

# Shrink the log2 FC values using the normal, apeglm and ashr methods

resNorm <- lfcShrink(dds, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type = "normal")
resAsh <- lfcShrink(dds, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type = "ashr")
resapeglm <- lfcShrink(dds, coef='condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8', type ='apeglm')
# Compare methods with different MA plots
plotMA(res)
plotMA(resNorm)
norm_plot <- recordPlot()
plotMA(resAsh)
ashr_plot <- recordPlot()
plotMA(resapeglm)
apeglm_plot <- recordPlot()

plot_grid(norm_plot,ashr_plot,apeglm_plot, ncol = 1)
#  recordPlot allows us to save base R plots to an object
# This is useful to then facet/panel figures

# Generate a results table (one shrinkage method) containing mgi symbols
# Use the EnsDb.Mmusculus.v79 package
res_df <- as.data.frame(resapeglm)
head(res_df)
library(EnsDb.Mmusculus.v79)
# columns(EnsDb.Mmusculus.v79)
# keytypes(EnsDb.Mmusculus.v79)
keys <- rownames(res_df)
res_symbols <- select(EnsDb.Mmusculus.v79,
       keys=keys,
       columns=c("SYMBOL", "GENEID"),
       keytype="GENEID")
head(res_symbols)
class(res_symbols)

# Count number of duplicates
table(duplicated(res_symbols$SYMBOL))
# Say which are duplicates
res_symbols$SYMBOL[duplicated(res_symbols$SYMBOL)]

duplicated_genes <- dplyr::filter(res_symbols,SYMBOL %in%  res_symbols$SYMBOL[duplicated(res_symbols$SYMBOL)])

# Count how many NAs
dim(res_df[is.na(res_df$SYMBOL),])

#Left join w symbols
res_df <- rownames_to_column(res_df, "GENEID")

res_df <- left_join(res_df,res_symbols)
head(res_df)

# Remove all genes with a padj of NA
head(res_df)
res_df <- res_df[!is.na(res_df$padj),]


# Write the results table to a CSV file
head(res_df)
dim(res_df)

write.csv(res_df, file = "results/res_df.csv", quote=FALSE, row.names = FALSE)

# Filter the results table for padj < 0.05 & logFC >1, and write to a CSV file
res_df_p_05 <- dplyr::filter(res_df, padj<0.05)  %>% 
    dplyr::filter(abs(log2FoldChange)>1)
dim(res_df_p_05)
write.csv(res_df_p_05, file = "results/res_df_p_05.csv", quote=FALSE, row.names = FALSE)

# Generate VST and rlog transformed counts:

# The point of these two transformations, the VST and the rlog,
# is to remove the dependence of the variance on the mean, particularly
# the high variance of the logarithm of count data when the mean is low.
# Both VST and rlog use the experiment-wide trend of variance over mean,
# in order to transform the data to remove the experiment-wide trend.

# Generally recommend blind = FALSE
dds_vst <- vst(dds, blind = FALSE)
dds_rlog <- rlog(dds, blind = FALSE)

vst_plot <- vsn::meanSdPlot(assay(dds_vst))
rlog_plot <- vsn::meanSdPlot(assay(dds_rlog))

# Using both sets of transformed counts:
# Generate a PCA plot either using all genes, or top 500 most variable genes
dim(dds_vst)
DESeq2::plotPCA(dds_vst, ntop=nrow(dds_vst))
DESeq2::plotPCA(dds_vst, ntop=500)

# Generate a heatmap of the top 20 (by shrunken FC)
# differentially-expressed genes – label samples by
# condition and genes by mgi symbol

library("pheatmap")

top_20 <- res_df_p_05[order(-abs(res_df_p_05$log2FoldChange)),] %>% 
    slice_head(n=20) %>% 
    pull(GENEID,SYMBOL)
heatmap_genes <- as.data.frame(assay(dds_vst)) %>%
    dplyr::filter(rownames(.) %in% top_20)

anno <- data.frame("Condition" = colData(dds)[,c("condition")])
rownames(anno) <- colnames(assay(dds_vst)[top_20,])

pheatmap(assay(dds_vst)[top_20,], annotation_col=anno, scale = "row")


# Plot a volcano plot:
# Ø Highlight significantly differentially expressed genes (p adj < 0.05, log2FC > 1) in red
dim(res_df)
head(res_df)
colour <- c('black','red')
sig_deg <- with(res_df, factor((padj<0.05) & abs(log2FoldChange)>1))
colour <- colour[sig_deg]

vp_df <- res_df %>% 
    dplyr::mutate(Log10_padj=-log10(padj))

ggplot(vp_df, aes(x=log2FoldChange, y=Log10_padj)) +
    geom_point(colour=colour) +
    ylim(0, 310)


# Ø Add labels to highlight the location of some interesting genes
top_genes <- names(top_20)[!is.na(names(top_20))]
vp_df_subset <- vp_df %>% 
    dplyr::filter(SYMBOL %in% top_genes)

labels <- data.frame(c('Rgs20', 'Oprk1', 'Tram1', 'Kcnq5'))

ggplot(vp_df, aes(x=log2FoldChange, y=Log10_padj)) +
    geom_point(colour=colour) +
    ylim(0, 310) +
    geom_label(data=vp_df_subset, aes(label=SYMBOL)) # Can change aesthetics throughout the ggplot and feed new data (aes usually has to correspond to data for given geom)



