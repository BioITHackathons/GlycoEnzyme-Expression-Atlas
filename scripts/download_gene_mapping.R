#!/usr/bin/env Rscript

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)

# Get all Ensembl IDs and their corresponding gene symbols
message("Getting gene mappings from org.Hs.eg.db...")
ensembl_to_symbol <- select(org.Hs.eg.db, 
                          keys=keys(org.Hs.eg.db, keytype="ENSEMBL"),
                          columns=c("SYMBOL", "ENSEMBL"),
                          keytype="ENSEMBL")

# Save the mapping
message("Saving gene ID mapping...")
write.csv(ensembl_to_symbol, "data/processed/gene_id_mapping.csv", row.names = FALSE)
message("Gene ID mapping saved to data/processed/gene_id_mapping.csv") 