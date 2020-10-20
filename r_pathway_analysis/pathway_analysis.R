renv::restore()
install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
# BiocManager::install("GO.db")
options(timeout=300)
# BiocManager::install("EnsDb.Mmusculus.v79")
# BiocManager::install("tidyverse")
library(clusterProfiler)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(org.Mm.eg.db)

data <- read.csv("data/res_df.csv")
dim(data)
# View(data)


# columns(EnsDb.Mmusculus.v79)
# keytypes(EnsDb.Mmusculus.v79)
keys <- data$GENEID
gene_id <- select(EnsDb.Mmusculus.v79,
                      keys=keys,
                      columns=c("ENTREZID", "GENEID"),
                      keytype="GENEID")
head(gene_id)
class(gene_id)
data <- data %>% 
    left_join(gene_id)


# Upregulated_genes - where padj < 0.05 & log2FoldChange > 1
upregulated_genes <- data[data$padj < 0.05 & data$log2FoldChange > 1,]

# Downregulated genes - where padj < 0.05 & log2FoldChange < -1
downregulated_genes <- data[data$padj < 0.05 & data$log2FoldChange < -1,]

# Remove NAs

data_noNA <- data %>% 
    dplyr::filter(!is.na(ENTREZID))
upregulated_genes_noNA <- upregulated_genes %>% 
    dplyr::filter(!is.na(ENTREZID))
downregulated_genes_noNA <- downregulated_genes %>% 
    dplyr::filter(!is.na(ENTREZID))

dim(upregulated_genes_noNA)
dim(downregulated_genes_noNA)

# Make genelists
# The geneList contains three features
# 1) numeric vector: fold change or other type of numerical variable
# 2) named vector: every number was named by the corresponding gene ID
# 3) sorted vector: number should be sorted in decreasing order

genelist <- data_noNA[,3]
names(genelist) <- data_noNA$ENTREZID
genelist <- sort(genelist,decreasing=TRUE)
length(genelist) #13968

genelist_up <- upregulated_genes_noNA[,3]
names(genelist_up) <- upregulated_genes_noNA$ENTREZID
genelist_up <- sort(genelist_up,decreasing=TRUE)
length(genelist_up) #1244

genelist_down <- downregulated_genes_noNA[,3]
names(genelist_down) <- downregulated_genes_noNA$ENTREZID
genelist_down <- sort(genelist_down,decreasing=TRUE)
length(genelist_down) #879

# Perform ORA with GO gene sets
ego <- enrichGO(gene=names(genelist_up),
                universe=names(genelist),
                OrgDb = org.Mm.eg.db)

# Perform ORA with KEGG gene sets
ekegg <- enrichKEGG(gene=names(genelist_up),
                    universe=names(genelist),
                    organism='mmu',
                    pvalueCutoff=0.05,
                    pAdjustMethod='BH')
head(ekegg)
class(ekegg)

# Plot

barplot(ego)
dotplot(ego)
ego_readable <- setReadable(ego, "org.Mm.eg.db", "ENTREZID")
cnetplot(ego_readable, foldChange = genelist, circular=TRUE, colorEdge=TRUE)
heatplot(ego_readable, foldChange = genelist)

barplot(ekegg)
dotplot(ekegg)
ekegg_readable <- setReadable(ekegg, "org.Mm.eg.db", "ENTREZID")
cnetplot(ekegg_readable, foldChange = genelist, circular=TRUE, colorEdge=TRUE)
heatplot(ekegg_readable, foldChange = genelist)

plot <- dotplot(ego, showCategory=20)
class(plot)
plot + 
    ggtitle("dotplot") +
    theme(axis.text.y = element_text(size=2))


# GSEA

#before doing GSEA, remove duplicate gene names
genelist_no_dups <- genelist[!duplicated(names(genelist))]
length(genelist)
length(genelist_no_dups)

gsa <- gseGO(geneList     = genelist_no_dups,
             OrgDb        = org.Mm.eg.db,
             ont          = "ALL",
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = TRUE)
gsa$Description[1:3]

dotplot(gsa)
barplot(gsa)
gsa_readable <- setReadable(gsa, "org.Mm.eg.db", "ENTREZID")
cnetplot(gsa_readable, foldChange = genelist, circular=TRUE, colorEdge=TRUE)
heatplot(gsa_readable, foldChange = genelist)
enrichplot::gseaplot2(gsa, geneSetID = c("GO:0045087","GO:0002768"), pvalue_table = TRUE, ES_geom = "dot")

terms <- gsa$Description[1:3]
p <- enrichplot::pmcplot(terms, 2010:2020)
p2 <- enrichplot::pmcplot(terms, 2010:2017, proportion=FALSE)
plot_grid(p, p2, ncol=2)
