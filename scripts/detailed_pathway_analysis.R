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
if (!requireNamespace("pheatmap", quietly = TRUE))
  BiocManager::install("pheatmap")
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(pheatmap)
library(igraph)

# Create necessary directories
dir.create("results/figures/detailed_pathways", recursive = TRUE, showWarnings = FALSE)
dir.create("results/network", recursive = TRUE, showWarnings = FALSE)

# Read pathway analysis results
kegg_results <- read.csv("results/pathway/kegg_pathways.csv")
reactome_results <- read.csv("results/pathway/reactome_pathways.csv")

# Read gene expression data
results <- read.csv("data/Deseq2_results/fshd_significant_genes.csv", row.names = 1)

# Clean Ensembl IDs
clean_ensembl_ids <- function(ids) {
  return(gsub("\\.[0-9]+$", "", ids))
}
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

# Function to create detailed pathway visualizations
create_pathway_visualizations <- function(pathway_id, pathway_name, gene_list, type = "kegg") {
  tryCatch({
    # Get pathway genes
    if (type == "kegg") {
      pathway_genes <- kegg_results$core_enrichment[kegg_results$ID == pathway_id]
      pathway_genes <- unlist(strsplit(pathway_genes, "/"))
    } else {
      pathway_genes <- reactome_results$core_enrichment[reactome_results$ID == pathway_id]
      pathway_genes <- unlist(strsplit(pathway_genes, "/"))
    }
    
    # Get matching genes
    matching_indices <- which(clean_ids %in% entrez_ids$ENSEMBL[entrez_ids$ENTREZID %in% pathway_genes])
    
    if (length(matching_indices) > 0) {
      # Create data frame for visualization
      pathway_expr <- data.frame(
        log2FoldChange = results$log2FoldChange[matching_indices],
        padj = results$padj[matching_indices],
        row.names = rownames(results)[matching_indices]
      )
      
      # Create expression matrix for heatmap
      expr_matrix <- as.matrix(pathway_expr$log2FoldChange)
      rownames(expr_matrix) <- rownames(pathway_expr)
      colnames(expr_matrix) <- "log2FoldChange"
      
      # Create heatmap
      pdf(paste0("results/figures/detailed_pathways/", pathway_id, "_heatmap.pdf"), width = 10, height = 8)
      pheatmap(expr_matrix,
               scale = "none",
               show_rownames = TRUE,
               main = paste(pathway_name, "Gene Expression"),
               cluster_cols = FALSE)
      dev.off()
      
      # Create volcano plot
      pathway_volcano <- ggplot(pathway_expr, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = padj < 0.05)) +
        theme_bw() +
        ggtitle(paste(pathway_name, "Volcano Plot")) +
        scale_color_manual(values = c("grey", "red"), name = "Significant") +
        xlab("log2 Fold Change") +
        ylab("-log10(adjusted p-value)")
      
      ggsave(paste0("results/figures/detailed_pathways/", pathway_id, "_volcano.pdf"), 
             pathway_volcano, width = 8, height = 6)
      
      cat("Created visualizations for pathway:", pathway_name, "\n")
    } else {
      cat("No matching genes found for pathway:", pathway_name, "\n")
    }
  }, error = function(e) {
    cat("Error processing pathway:", pathway_name, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
}

# Create detailed visualizations for top pathways
cat("\nProcessing KEGG pathways...\n")
top_kegg_pathways <- head(kegg_results, 5)
for (i in 1:nrow(top_kegg_pathways)) {
  create_pathway_visualizations(
    pathway_id = top_kegg_pathways$ID[i],
    pathway_name = top_kegg_pathways$Description[i],
    gene_list = gene_list,
    type = "kegg"
  )
}

cat("\nProcessing Reactome pathways...\n")
top_reactome_pathways <- head(reactome_results, 5)
for (i in 1:nrow(top_reactome_pathways)) {
  create_pathway_visualizations(
    pathway_id = top_reactome_pathways$ID[i],
    pathway_name = top_reactome_pathways$Description[i],
    gene_list = gene_list,
    type = "reactome"
  )
}

# Prepare data for Cytoscape visualization
cat("\nPreparing Cytoscape files...\n")

# Create a file with pathway-gene relationships and gene attributes
pathway_gene_relationships <- data.frame()

# Add KEGG pathway relationships
for (i in 1:nrow(kegg_results)) {
  genes <- unlist(strsplit(kegg_results$core_enrichment[i], "/"))
  pathway_gene_relationships <- rbind(pathway_gene_relationships,
                                    data.frame(
                                      pathway_id = kegg_results$ID[i],
                                      pathway_name = kegg_results$Description[i],
                                      gene_id = genes,
                                      pathway_type = "KEGG",
                                      NES = kegg_results$NES[i],
                                      pvalue = kegg_results$pvalue[i]
                                    ))
}

# Add Reactome pathway relationships
for (i in 1:nrow(reactome_results)) {
  genes <- unlist(strsplit(reactome_results$core_enrichment[i], "/"))
  pathway_gene_relationships <- rbind(pathway_gene_relationships,
                                    data.frame(
                                      pathway_id = reactome_results$ID[i],
                                      pathway_name = reactome_results$Description[i],
                                      gene_id = genes,
                                      pathway_type = "Reactome",
                                      NES = reactome_results$NES[i],
                                      pvalue = reactome_results$pvalue[i]
                                    ))
}

# Save pathway-gene relationships for Cytoscape
write.csv(pathway_gene_relationships, 
          "results/network/pathway_gene_relationships.csv", 
          row.names = FALSE)

# Create gene attribute file for Cytoscape
gene_attributes <- data.frame(
  gene_id = names(gene_list),
  log2FoldChange = gene_list,
  padj = results$padj[match(names(gene_list), id_mapping[clean_ids])]
)

write.csv(gene_attributes,
          "results/network/gene_attributes.csv",
          row.names = FALSE)

# Create a summary of the analysis
cat("\nDetailed Pathway Analysis Summary:\n------------------------------------\n")
cat("Number of detailed pathway visualizations created:", 
    nrow(top_kegg_pathways) + nrow(top_reactome_pathways), "\n")
cat("\nOutput files created:\n")
cat("1. Pathway visualizations in results/figures/detailed_pathways/\n")
cat("2. Cytoscape network file: results/network/pathway_gene_relationships.csv\n")
cat("3. Cytoscape gene attributes: results/network/gene_attributes.csv\n")

cat("\nNext steps for Cytoscape visualization:\n")
cat("1. Open Cytoscape\n")
cat("2. Import network from results/network/pathway_gene_relationships.csv\n")
cat("3. Import node attributes from results/network/gene_attributes.csv\n")
cat("4. Use 'pathway_type' and 'NES' for pathway node styling\n")
cat("5. Use 'log2FoldChange' and 'padj' for gene node styling\n") 