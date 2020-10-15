#TIDYVERSE exercises

library(tidyverse)
library(biomaRt)
library(pheatmap)

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

# mutate(count_per_total = count/total_count_per_sample) %>%
#   mutate(cpm = count_per_total*1e6)
