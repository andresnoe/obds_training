# Code clinic 12/11/20

library(tidyverse)
library(Seurat)
library(patchwork)
library(DropletUtils)
library(BiocGenerics)
library(scDblFinder)

pbmc1k <- DropletUtils::read10xCounts(
  c(pbmc1k= "/Users/andresnoe/obds_sep20/working_directory/pbmc_1k_v3_raw"),
  col.names = TRUE
)
pbmc1k # 33538 x 6794880

pbmc1k_filtered <- DropletUtils::read10xCounts(
  c(pbmc1k= "/Users/andresnoe/obds_sep20/working_directory/pbmc_1k_v3_filtered"),
  col.names = TRUE
)
pbmc1k_filtered

set.seed(100)
e.out <- emptyDrops(assay(pbmc1k))
# lower = 100 by default
# 	A numeric scalar specifying the lower bound on the total UMI count,
# at or below which all barcodes are assumed to correspond to empty droplets.
e.out$rank<- rank(-e.out$Total)
e.out
dim(e.out)

plotdata <- dplyr::filter(as.data.frame(e.out), Total>0)
dim(plotdata)

ggplot(plotdata, aes(Total, -LogProb))+
  geom_point()+
  labs(x="Total UMI count", y="-Log Probability")

ggplot(plotdata, aes(x=rank, y=Total, colour=FDR>0.01))+
  geom_point()+
  labs(x="Cell rank according to count", y="Total UMI count")+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")
# blue  = true cells; red are most likely empty droplets


is.cell <- e.out$FDR <= 0.01
# FDR that's small is a cell, FDR that's larger likely means an empty droplet
sum(is.cell, na.rm=TRUE)
table(Limited=e.out$Limited, Significant=is.cell)

pbmc1k_droplets_filtered <- pbmc1k[, which(e.out$FDR <= 0.01)]
pbmc1k_droplets_filtered

dim(pbmc1k_droplets_filtered) # 33538  1207
sum(is.cell, na.rm=TRUE) # 1207

dim(pbmc1k_filtered) # cellRanger filtered dataset

# table(emptyDroplets=colnames(pbmc1k_droplets_filtered),
      # cellRanger=colnames(pbmc1k_filtered))


# DOUBLET DETECTION

doublet.out <- scDblFinder(pbmc1k_droplets_filtered)

colData(doublet.out)
doublet.out
table(call=doublet.out$scDblFinder.class)

ggplot(as.data.frame(colData(doublet.out)), aes(scDblFinder.score))+
  geom_histogram()


pbmc1k_droplets_dbl_filtered <-
  doublet.out[, which(doublet.out$scDblFinder.class=="singlet")]

  