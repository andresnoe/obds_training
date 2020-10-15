#renv::install("tidyverse")
#renv::install("BiocManager")
#BiocManager::install("cowplot")
#BiocManager::install("gridExtra")
#BiocManager::install("patchwork")
#BiocManager::install("biomaRt")
#BiocManager::install("pheatmap")

library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)
#library(biomaRt)
#library(pheatmap)


# PLOTTING IN BASE R

data("mtcars")

par()$lwd
par(mar=c(6,6,5,5), lwd=5, mfrow=c(1,2)) # mfrow sets the plotting area into a 1*2 array


str(mtcars)
table(mtcars$gear)


barplot(table(mtcars$gear),
        xlab = "Number of gears",
        ylab = "Number of cars",
        main = "A main title\nspread across two lines",
        col = "red",
        border = "cadetblue")
abline(h=6, lwd=70)

# Making vector of colours, named by gears
colours <- c("red","green","blue")
names(colours) <- c(3,4,5)
colours
gears_column <- as.character(mtcars$gear)
gears_column
colours[gears_column]

par(mar=c(6,6,5,5), lwd=1)
plot(mtcars$mpg, mtcars$hp,
     col = colours[gears_column],
     pch = 16,
     cex = 6,
     xlab = "Miles per gallon",
     ylab = "Horse power",
     cex.lab = 2)
legend(legend = sort(unique(mtcars$gear)),
       x = "right",
       fill = c("red","green","blue"))

# PLOTTING IN ggplot2

#import file
getwd()
cgr <- read.table("data/coding_gene_region.bed")

str(cgr)
head(cgr)
tail(cgr)
colnames(cgr) <- c("chr","start","end","ensemblID","score","strand")
head(cgr)
class(cgr)


cgr$length <- cgr$end-cgr$start
cgr

ggplot(cgr, aes(x = length)) +
  geom_histogram(bins = 100, fill = "red", colour = "green") +
  labs(title = "Histogram of interval length",
       x = "Interval length",
       y = "Count") +
  theme(axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 10, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 135, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust=0.5)) +
  xlim(0,1500000)
  
