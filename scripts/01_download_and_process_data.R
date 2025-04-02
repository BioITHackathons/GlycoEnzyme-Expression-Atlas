#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("tibble", quietly = TRUE))
  install.packages("tibble")

library(GEOquery)
library(limma)
library(dplyr)
library(tibble)

# Create necessary directories
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Download FSHD dataset (GSE140261)
message("Downloading FSHD dataset (GSE140261)...")
gset <- getGEO("GSE140261", GSEMatrix = TRUE, getGPL = TRUE)

# Extract expression data and phenotype information
expr_matrix <- exprs(gset[[1]])
pheno_data <- pData(gset[[1]])

# Save raw data
message("Saving raw data...")
saveRDS(expr_matrix, "data/raw/fshd_expr_matrix.rds")
saveRDS(pheno_data, "data/raw/fshd_pheno_data.rds")

# Process phenotype data
message("Processing phenotype data...")
# Extract relevant columns and rename them
pheno_processed <- pheno_data %>%
  select(
    title,
    geo_accession,
    characteristics_ch1,
    characteristics_ch1.1,
    characteristics_ch1.2,
    characteristics_ch1.3
  ) %>%
  rename(
    sample_name = title,
    sample_id = geo_accession,
    tissue = characteristics_ch1,
    condition = characteristics_ch1.1,
    age = characteristics_ch1.2,
    gender = characteristics_ch1.3
  )

# Clean up the phenotype data
pheno_processed <- pheno_processed %>%
  mutate(
    tissue = gsub("tissue: ", "", tissue),
    condition = gsub("condition: ", "", condition),
    age = gsub("age: ", "", age),
    gender = gsub("gender: ", "", gender)
  )

# Save processed phenotype data
message("Saving processed phenotype data...")
saveRDS(pheno_processed, "data/processed/fshd_pheno_processed.rds")

# Create a summary of the dataset
message("Creating dataset summary...")
dataset_summary <- list(
  total_samples = nrow(pheno_processed),
  total_genes = nrow(expr_matrix),
  conditions = table(pheno_processed$condition),
  tissues = table(pheno_processed$tissue),
  gender_distribution = table(pheno_processed$gender)
)

# Save dataset summary
saveRDS(dataset_summary, "data/processed/fshd_dataset_summary.rds")

message("Data processing completed successfully!") 