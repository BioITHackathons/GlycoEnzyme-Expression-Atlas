#!/usr/bin/env Rscript

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("pheatmap", quietly = TRUE))
  BiocManager::install("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
if (!requireNamespace("plotly", quietly = TRUE))
  BiocManager::install("plotly")
if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  BiocManager::install("RColorBrewer")

library(pheatmap)
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer)

# Read the data
results <- read.csv("data/processed/glycogene_expression_analysis.csv", row.names = 1)
all_glycogenes <- readLines("data/isha_data/all_glycogenes.txt")
counts_matrix <- read.csv("data/processed/fshd_counts_matrix.csv", row.names = 1)

# Check for matching genes
matching_genes <- intersect(results$gene_symbol, all_glycogenes)
cat("\nMatching genes found in both datasets:\n-------------------------------------\n")
cat(paste(matching_genes, collapse = ", "))
cat("\n\nTotal matching genes:", length(matching_genes), "\n\n")

# Filter significant results
sig_results <- results[results$padj < 0.05, ]

# Create publication-ready volcano plot
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, 
                               ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"),
                               "Not significant"),
             alpha = ifelse(padj < 0.05, 0.8, 0.3)),
             size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "#E41A1C", 
                               "Down-regulated" = "#377EB8",
                               "Not significant" = "#999999")) +
  geom_text(data = subset(results, abs(log2FoldChange) > 2 & padj < 0.05),
            aes(label = gene_symbol),
            size = 3,
            hjust = 0.5,
            vjust = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "top") +
  labs(x = "log2 Fold Change",
       y = "-log10(adjusted p-value)")

# Save volcano plot
ggsave("data/processed/glycogene_volcano_plot.pdf", volcano_plot, width = 8, height = 6)

# Create interactive volcano plot
interactive_plot <- ggplotly(volcano_plot)
htmlwidgets::saveWidget(interactive_plot, "data/processed/glycogene_volcano_plot_interactive.html")

# Prepare data for heatmap
sig_genes <- rownames(sig_results)
sig_counts <- counts_matrix[sig_genes, ]
sig_counts[is.na(sig_counts)] <- 0
scaled_counts <- t(scale(t(sig_counts)))

# Remove rows with NA or Inf values
scaled_counts <- scaled_counts[complete.cases(scaled_counts), ]
scaled_counts[is.infinite(scaled_counts)] <- 0

# Create publication-ready heatmap
sample_groups <- rep(c("Control", "FSHD"), c(8, 27))  # Adjusted to match actual sample numbers
annotation_col <- data.frame(
  Condition = factor(sample_groups),
  row.names = colnames(scaled_counts)
)

annotation_colors <- list(
  Condition = c(Control = "#377EB8", FSHD = "#E41A1C")
)

heatmap <- pheatmap(scaled_counts,
                    scale = "none",
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    clustering_method = "ward.D2",
                    annotation_col = annotation_col,
                    annotation_colors = annotation_colors,
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    fontsize = 8,
                    border_color = NA,
                    main = "Differentially Expressed Glycogenes",
                    legend = TRUE,
                    legend_breaks = c(-2, -1, 0, 1, 2),
                    legend_labels = c("-2", "-1", "0", "1", "2"),
                    color = colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"))(100))

# Save heatmap
pdf("data/processed/glycogene_heatmap.pdf", width = 10, height = 12)
heatmap
dev.off()

# Print summary
cat("\nVisualization Summary:\n---------------------\n")
cat("Total glycogenes analyzed:", length(matching_genes), "\n")
cat("Significant glycogenes (padj < 0.05):", nrow(sig_results), "\n")
cat("Up-regulated significant glycogenes:", sum(sig_results$log2FoldChange > 0), "\n")
cat("Down-regulated significant glycogenes:", sum(sig_results$log2FoldChange < 0), "\n\n")

# Save differentially expressed glycogenes to a CSV file with additional information
sig_results <- sig_results[order(abs(sig_results$log2FoldChange), decreasing = TRUE), ]
sig_results$regulation <- ifelse(sig_results$log2FoldChange > 0, "Up-regulated", "Down-regulated")
sig_results$ensembl_id <- rownames(sig_results)
write.csv(sig_results, "data/processed/differentially_expressed_glycogenes.csv", row.names = FALSE)

cat("Output files created:\n")
cat("1. data/processed/glycogene_volcano_plot.pdf\n")
cat("2. data/processed/glycogene_volcano_plot_interactive.html\n")
cat("3. data/processed/glycogene_heatmap.pdf\n")
cat("4. data/processed/differentially_expressed_glycogenes.csv\n") 