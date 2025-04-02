#!/usr/bin/env Rscript

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
if (!requireNamespace("enrichplot", quietly = TRUE))
  BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(enrichplot)

# Read DESeq2 results
results <- read.csv("data/Deseq2_results/fshd_significant_genes.csv", row.names = 1)

# Clean Ensembl IDs (remove version numbers)
clean_ensembl_ids <- function(ids) {
  return(gsub("\\.[0-9]+$", "", ids))
}

# Clean the Ensembl IDs
clean_ids <- clean_ensembl_ids(rownames(results))
rownames(results) <- clean_ids

# Convert Ensembl IDs to Entrez IDs
entrez_ids <- bitr(clean_ids, 
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Create a mapping between Ensembl and Entrez IDs
id_mapping <- setNames(entrez_ids$ENTREZID, entrez_ids$ENSEMBL)

# Create gene list with Entrez IDs
gene_list <- results$log2FoldChange[clean_ids %in% entrez_ids$ENSEMBL]
names(gene_list) <- id_mapping[clean_ids[clean_ids %in% entrez_ids$ENSEMBL]]
gene_list <- sort(gene_list, decreasing = TRUE)

# KEGG Pathway Analysis using Entrez IDs
kegg_result <- gseKEGG(geneList = gene_list,
                       organism = 'hsa',
                       keyType = 'kegg',
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       nPermSimple = 10000)

# Reactome Pathway Analysis
reactome_result <- gsePathway(geneList = gene_list,
                            organism = 'human',
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 0.05,
                            nPermSimple = 10000)

# Calculate term similarity
kegg_similarity <- pairwise_termsim(kegg_result)
reactome_similarity <- pairwise_termsim(reactome_result)

# Save pathway analysis results
write.csv(as.data.frame(kegg_result), "results/pathway/kegg_pathways.csv", row.names = FALSE)
write.csv(as.data.frame(reactome_result), "results/pathway/reactome_pathways.csv", row.names = FALSE)

# Create pathway visualization plots
pdf("results/figures/kegg_pathway_dotplot.pdf", width = 12, height = 8)
dotplot(kegg_result, showCategory = 20) + 
  theme_bw() +
  ggtitle("Top 20 KEGG Pathways")
dev.off()

pdf("results/figures/reactome_pathway_dotplot.pdf", width = 12, height = 8)
dotplot(reactome_result, showCategory = 20) + 
  theme_bw() +
  ggtitle("Top 20 Reactome Pathways")
dev.off()

# Create enrichment maps
pdf("results/figures/kegg_enrichment_map.pdf", width = 12, height = 10)
emapplot(kegg_similarity, showCategory = 30)
dev.off()

pdf("results/figures/reactome_enrichment_map.pdf", width = 12, height = 10)
emapplot(reactome_similarity, showCategory = 30)
dev.off()

# Create category netplots
pdf("results/figures/kegg_cnetplot.pdf", width = 14, height = 10)
cnetplot(kegg_result, 
         categorySize = "pvalue", 
         showCategory = 10,
         foldChange = gene_list)
dev.off()

pdf("results/figures/reactome_cnetplot.pdf", width = 14, height = 10)
cnetplot(reactome_result, 
         categorySize = "pvalue", 
         showCategory = 10,
         foldChange = gene_list)
dev.off()

# Prepare data for STRING/Cytoscape analysis
# Create a file with significant genes and their fold changes
string_input <- data.frame(
  ensembl_id = clean_ids,
  entrez_id = id_mapping[clean_ids],
  log2FC = results$log2FoldChange,
  padj = results$padj
)

# Save STRING input file
write.csv(string_input, "results/network/string_input.csv", row.names = FALSE)

# Create a summary of the analysis
cat("\nPathway and Network Analysis Summary:\n------------------------------------\n")
cat("Number of significant genes analyzed:", nrow(results), "\n")
cat("Number of genes with Entrez IDs:", length(gene_list), "\n")
cat("Number of KEGG pathways found:", nrow(kegg_result), "\n")
cat("Number of Reactome pathways found:", nrow(reactome_result), "\n\n")

cat("Output files created:\n")
cat("1. results/pathway/kegg_pathways.csv\n")
cat("2. results/pathway/reactome_pathways.csv\n")
cat("3. results/figures/kegg_pathway_dotplot.pdf\n")
cat("4. results/figures/reactome_pathway_dotplot.pdf\n")
cat("5. results/figures/kegg_enrichment_map.pdf\n")
cat("6. results/figures/reactome_enrichment_map.pdf\n")
cat("7. results/figures/kegg_cnetplot.pdf\n")
cat("8. results/figures/reactome_cnetplot.pdf\n")
cat("9. results/network/string_input.csv\n") 