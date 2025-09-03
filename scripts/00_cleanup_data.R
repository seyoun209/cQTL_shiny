setwd("/work/users/s/e/seyoun/cQTL_shiny")

# Load libraries
library(tidyverse)
library(plotly)
library(jsonlite)

# Read data
counts <- read_csv("data/raw_counts_with_coords.csv")
meta <- read_tsv("data/donor_samples.txt")
rna <- read_tsv("data/rna_extraction.txt")
de <- read_csv("data/de_genes_fromAll_fnf_results.csv")

# Clean and merge metadata
meta_data <- rna %>%
  select(Sample, Donor, Condition, RIN, Tech_Rep) %>%
  left_join(meta %>% select(Donor, Sex, Age, FragmentBatch), by = "Donor") |> 
  mutate(
    donor = Donor,
    condition = Condition,
    sample = paste(Donor, Condition, Tech_Rep , Sex, sep = "_")
  ) %>%
  rename(age = Age, sex = Sex) %>%
  select(-Sample)

meta_data <- meta_data %>%
  distinct(sample, donor, condition, sex, .keep_all = TRUE)

# Reshape counts to long format
counts_long <- counts %>%
  pivot_longer(cols = starts_with("AM"), names_to = "sample", values_to = "count") %>%
  separate(sample, into = c("donor", "condition", "rep", "sex"), sep = "_", remove = FALSE)

# Add the base Ensembl ID (without version) to DE results for joining
de <- de %>%
  mutate(base_gene_id = str_replace(gene_id, "\\.[0-9]+$", ""))

# Merge with metadata 
counts_annot <- counts_long %>%
  left_join(meta_data, by = c("sample", "donor", "condition", "sex")) %>%
  left_join(de, by = c(
    "ensgID" = "base_gene_id",
    "seqnames" = "seqnames",
    "start" = "start",
    "end" = "end"
  ))



save(counts_annot, file="data/counts_annot.Rdata")

load("data/counts_annot.Rdata")

# Filter significant genes
sig_genes <- counts_annot %>% 
  filter(!is.na(padj) & padj < 0.05) %>%
  pull(ensgID) %>%
  unique()

print(paste("Number of significant genes:", length(sig_genes)))

# Prepare gene list for JSON
genes_for_json <- counts_annot %>%
  #filter(ensgID %in% sig_genes) %>%
  distinct(ensgID, SYMBOL, log2FoldChange, padj, baseMean) %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), ensgID, SYMBOL)) %>%
  rename(gene_id = ensgID)

# Export gene list
write_json(genes_for_json, "data/genes.json")

# Create expression data subset
counts_subset <- counts_annot %>% 
  #filter(ensgID %in% sig_genes) %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), ensgID, SYMBOL)) %>%
  select(ensgID, SYMBOL, sample, donor, condition, sex, count, age) %>%
  rename(Age = age)
  
print(paste("Number of rows in counts_subset:", nrow(counts_subset)))

# Export expression data
write_json(counts_subset, "data/expression_data.json")
