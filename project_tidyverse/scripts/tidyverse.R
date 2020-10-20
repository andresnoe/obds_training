#TIDYVERSE exercises

library(tidyverse)
library(biomaRt)
library(pheatmap)
# install.packages("ggthemes")
library(ggthemes)

# Read data and metadata in 
sample_table <- read_tsv("data/obds_sampletable.tsv")
counts_table <- read_tsv("data/obds_countstable.tsv.gz")
sample_table
counts_table


# Tidy and annotate to 3 columns: Geneid, sample, count
test_with_na <- counts_table %>% 
  pivot_longer(-Geneid, names_to = "sample_name", values_to = "count")

prepared_counts <- counts_table %>% 
  pivot_longer(-Geneid, names_to = "sample_name", values_to = "count", values_drop_na = TRUE)

# Join with gene info to get mgi_symbol: use biomaRt package
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
ensembl <- useMart("ensembl")
searchDatasets(mart = ensembl, pattern = "mus")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
# Using mmusculus_gene_ensembl version GRCm38.p6
View(ensembl)
listAttributes(ensembl)
symbols <- getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'),
                                mart = ensembl,
                                values = unique(prepared_counts$Geneid))
symbol_counts <- left_join(prepared_counts,symbols, by = c("Geneid" = "ensembl_gene_id"))

# Tidy metadata file: one variable per column and don't need species and library_layout columns
sample_table
# Separate up sample_title column
sample_table <- sample_table %>% 
  separate(sample_title, c("egr_locus", "genotype", "cell_type", 'replicate'), sep = "_") %>%
  unite(col = "genotype", egr_locus, genotype, sep = "_") %>%
  dplyr::select(-c(species,library_layout))
sample_table

# Add metadata to table with counts and gene info
symbol_counts
sample_table
processed_joined <- left_join(symbol_counts,  sample_table, by=c('sample_name'='Sample_accession'))
processed_joined

# Calculate counts per million (CPM) - use group_by() and mutate()
grouped <- processed_joined %>%
  group_by(sample_name) %>%
  mutate(total_count_per_sample = sum(count)) %>%
  mutate(total_count_in_million = total_count_per_sample/1e6) %>%
  mutate(cpm = count/total_count_in_million) %>%
  mutate(log_cpm = log2(cpm+1))
View(grouped)

# Plot read depth per sample
# Use group_by() and summarise()

grouped %>%
  ggplot(aes(x=sample_name, y=total_count_per_sample, colour = replicate)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # theme is the thing we use to change the non-data components of the plot

# How many genes have no counts for any sample?

zero_count_genes <- group_by(grouped,mgi_symbol) %>%
  summarise(total_gene_count = sum(count)) %>% # ungroups data and creates new column with name "total_gene..."
  filter(total_gene_count == 0) %>%
  pull(mgi_symbol)
length(zero_count_genes)
length(unique(zero_count_genes))
  
# Draw a density plot of log2(CPM + 1) for all genes

grouped %>%
  ggplot(aes(x= log_cpm,  colour = sample_name)) +
           geom_density()

# Filter out genes that have low expression in 3 or fewer samples
# For low expression use CPM < 0.5

# filter (i.e. keep) genes that have detectable expression (> 0.5 CPM) in at least 3 samples
# Keep genes that are expressed in number of samples that form smallest experimental unit
grouped_filtered <- grouped %>% 
  group_by(Geneid) %>%
  mutate(high_expression = sum(cpm>0.5)) %>%
  filter(high_expression>=3)

# What proportion of genes are lowly expressed?
gene_list <- unique(grouped$Geneid)
gene_list
length(gene_list)
length(unique(grouped_filtered$Geneid))

proportion_lowly_expressed <- (100-(100*(length(unique(grouped_filtered$Geneid))) / length(gene_list)))

proportion_lowly_expressed

# Make a density plot of log2(CPM + 1) with the filtered data

grouped_filtered %>%
  ggplot(aes(x= log_cpm,  colour = sample_name)) +
  geom_density()

grouped_filtered %>%
  ggplot(aes(x= log_cpm,  colour = sample_name)) +
  geom_density() +
  facet_wrap(~sample)

# Plot CD4 and CD8 expression for all samples - does it make sense? # Colour by replicate and facet by genotype against cell type
filtered_tcell_genes <- grouped_filtered %>%
  filter(mgi_symbol %in% c("Cd4","Cd8a")) %>%
  filter(mgi_symbol == "Cd4" | mgi_symbol == "Cd8a")
filtered_tcell_genes

filtered_tcell_genes %>%
  ggplot(aes(x= mgi_symbol, y = log_cpm, fill=replicate)) +
  geom_col(position = position_dodge()) +
  facet_grid(genotype ~ cell_type)
# If using geom_bar only provide x aesthetic
# If using geom_col can provide both x and y

# Choose 8 biologically relevant genes and plot a heatmap using the pheatmap package
hvg <- grouped_filtered %>% 
  group_by(Geneid) %>% 
  summarise(variance=var(log_cpm)) %>% 
  arrange(desc(variance)) %>% 
  slice_head(n=20)
hvg

hvg_counts <- grouped_filtered %>%
  filter(Geneid %in% hvg$Geneid)

#Filter is for rows, select is for columns
count_table_heatmap <- hvg_counts %>%
  ungroup() %>% 
  dplyr::select(mgi_symbol, sample_name, count) %>% 
  pivot_wider(names_from = sample_name, values_from = count) %>% 
  column_to_rownames(var="mgi_symbol")
count_table_heatmap
pheatmap(count_table_heatmap, scale = 'row')
# Scale = remove mean and divide by standard deviation
