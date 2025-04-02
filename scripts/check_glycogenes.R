#!/usr/bin/env Rscript

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(biomaRt)

# Read the significant genes and all genes
sig_genes <- read.csv("data/processed/fshd_significant_genes.csv", row.names = 1)
all_genes <- read.csv("data/processed/fshd_counts_matrix.csv", row.names = 1)

# Read the glycogenes list
glycogenes <- readLines("data/isha_data/all_glycogenes.txt")

# Function to clean Ensembl IDs (remove version numbers)
clean_ensembl_ids <- function(ids) {
  gsub("\\.[0-9]+$", "", ids)
}

# Clean the rownames of our data
rownames(all_genes) <- clean_ensembl_ids(rownames(all_genes))
rownames(sig_genes) <- clean_ensembl_ids(rownames(sig_genes))

# Read the gene ID mapping
message("Reading gene ID mapping...")
mapping <- read.csv("data/processed/gene_id_mapping.csv")

# Create mapping dictionaries
symbol_to_ensembl <- setNames(mapping$ENSEMBL, mapping$SYMBOL)
ensembl_to_symbol <- setNames(mapping$SYMBOL, mapping$ENSEMBL)

# Find which glycogenes are present in the dataset
glycogene_ensembl_ids <- unique(symbol_to_ensembl[glycogenes])
glycogene_ensembl_ids <- glycogene_ensembl_ids[!is.na(glycogene_ensembl_ids)]
present_glycogenes <- intersect(glycogene_ensembl_ids, rownames(all_genes))

# Print summary
cat("\nAnalysis Results:\n")
cat("----------------\n")
cat("Total glycogenes in list:", length(glycogenes), "\n")
cat("Glycogenes with Ensembl IDs:", length(glycogene_ensembl_ids), "\n")
cat("Glycogenes found in dataset:", length(present_glycogenes), "\n\n")

if(length(present_glycogenes) > 0) {
  cat("Expression analysis of glycogenes:\n")
  cat("--------------------------------\n")
  
  # Create a data frame for results
  results <- data.frame(
    ensembl_id = present_glycogenes,
    gene_symbol = ensembl_to_symbol[present_glycogenes],
    log2FoldChange = NA,
    padj = NA,
    baseMean = NA,
    stringsAsFactors = FALSE
  )
  
  # Get expression data for each glycogene
  for(i in 1:nrow(results)) {
    ensembl_id <- results$ensembl_id[i]
    if(ensembl_id %in% rownames(sig_genes)) {
      # If gene is in significant genes, get its stats
      results$log2FoldChange[i] <- sig_genes[ensembl_id, "log2FoldChange"]
      results$padj[i] <- sig_genes[ensembl_id, "padj"]
    }
    # Get base mean expression
    results$baseMean[i] <- mean(as.numeric(all_genes[ensembl_id,]))
  }
  
  # Sort by absolute log2FoldChange
  results <- results[order(abs(results$log2FoldChange), decreasing = TRUE, na.last = TRUE),]
  
  # Print results
  for(i in 1:nrow(results)) {
    gene <- results$gene_symbol[i]
    if(!is.na(results$log2FoldChange[i])) {
      cat(sprintf("%s (%s): log2FoldChange = %.2f, padj = %.2e, baseMean = %.2f\n", 
                  gene, results$ensembl_id[i],
                  results$log2FoldChange[i], results$padj[i], results$baseMean[i]))
    } else {
      cat(sprintf("%s (%s): Not significantly different, baseMean = %.2f\n", 
                  gene, results$ensembl_id[i], results$baseMean[i]))
    }
  }
  
  # Save results to file
  write.csv(results, "data/processed/glycogene_expression_analysis.csv", row.names = FALSE)
  cat("\nResults saved to data/processed/glycogene_expression_analysis.csv\n")
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  cat("------------------\n")
  cat("Total glycogenes analyzed:", nrow(results), "\n")
  cat("Significantly different glycogenes:", sum(!is.na(results$padj)), "\n")
  if(sum(!is.na(results$padj)) > 0) {
    cat("Up-regulated glycogenes:", sum(results$log2FoldChange > 0, na.rm = TRUE), "\n")
    cat("Down-regulated glycogenes:", sum(results$log2FoldChange < 0, na.rm = TRUE), "\n")
  }
} 