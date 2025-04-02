#!/usr/bin/env Rscript

# FSHD RNAseq DESeq2 Analysis for GlycoEnzyme-Expression-Atlas
# 
# This script:
# 1. Loads FSHD count matrix and sample metadata from processed data directory
# 2. Performs differential expression analysis using DESeq2
# 3. Filters results to focus on glycogenes
# 4. Creates visualizations for the differential expression results
# 5. Analyzes correlation between glycogenes and DUX4 expression
#
# Requirements:
# - R 4.0+
# - Bioconductor packages: DESeq2, ggplot2, pheatmap, biomaRt

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(biomaRt)
  library(RColorBrewer)
  library(ggrepel)
})

# Set up directories
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
processed_data_dir <- file.path(data_dir, "processed")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(results_dir, "figures")

# Create directories if they don't exist
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Function to load FSHD dataset
load_fshd_dataset <- function() {
  # Check if the required files exist
  counts_file <- file.path(processed_data_dir, "fshd_counts_matrix.csv")
  metadata_file <- file.path(processed_data_dir, "fshd_sample_metadata.csv")
  
  if (!file.exists(counts_file) || !file.exists(metadata_file)) {
    stop("FSHD dataset files not found. Please run download_rnaseq_data.R first.")
  }
  
  # Load data
  message("Loading FSHD dataset...")
  counts <- read.csv(counts_file, row.names = 1)
  metadata <- read.csv(metadata_file, row.names = 1)
  
  # Check for DUX4 scores
  dux4_scores_file <- file.path(processed_data_dir, "fshd_dux4_scores.csv")
  if (file.exists(dux4_scores_file)) {
    message("Loading DUX4 scores...")
    dux4_scores <- read.csv(dux4_scores_file)
  } else {
    dux4_scores <- NULL
  }
  
  return(list(counts = counts, metadata = metadata, dux4_scores = dux4_scores))
}

# Function to load glycogenes list
load_glycogenes <- function() {
  # Try to load from the all_glycogenes.txt file if it exists
  glycogenes_file <- file.path(data_dir, "isha_data", "all_glycogenes.txt")
  if (file.exists(glycogenes_file)) {
    glycogenes <- readLines(glycogenes_file)
    return(glycogenes)
  }
  
  # Alternative: try to load from the CSV file we might have created
  glycogenes_file_csv <- file.path(data_dir, "isha_data", "glycogenes_list.csv")
  if (file.exists(glycogenes_file_csv)) {
    glycogenes_df <- read.csv(glycogenes_file_csv)
    return(glycogenes_df$gene_symbol)
  }
  
  message("No glycogenes list found. Will analyze all genes.")
  return(NULL)
}

# Function to perform DESeq2 analysis
run_deseq2_analysis <- function(counts, metadata) {
  message("Running DESeq2 analysis for FSHD vs control...")
  
  # Ensure counts are integers (required by DESeq2)
  counts <- round(as.matrix(counts))
  
  # Make sure sample names match between counts and metadata
  shared_samples <- intersect(colnames(counts), rownames(metadata))
  if (length(shared_samples) == 0) {
    stop("No matching sample names found between count matrix and metadata")
  }
  
  # Filter to keep only samples present in both datasets
  counts <- counts[, shared_samples]
  metadata <- metadata[shared_samples, ]
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
  )
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Get results (FSHD vs control)
  res <- results(dds, contrast = c("condition", "FSHD", "control"))
  
  # Add gene names as a column
  res$gene <- rownames(res)
  
  # Convert to dataframe
  res_df <- as.data.frame(res)
  
  # Save results
  write.csv(res_df, file.path(processed_data_dir, "deseq2_results_fshd_vs_control.csv"))
  
  # Create MA plot
  pdf(file.path(figures_dir, "ma_plot_fshd_vs_control.pdf"))
  plotMA(res)
  dev.off()
  
  # Generate PCA plot
  vsd <- vst(dds, blind = FALSE)
  pdf(file.path(figures_dir, "pca_plot_fshd_vs_control.pdf"))
  plotPCA(vsd, intgroup = "condition")
  dev.off()
  
  # Return results and DESeq dataset
  return(list(results = res_df, dds = dds))
}

# Function to filter results for glycogenes
filter_for_glycogenes <- function(results_df, dds) {
  # Load glycogenes list
  glycogenes <- load_glycogenes()
  
  if (is.null(glycogenes) || length(glycogenes) == 0) {
    message("No glycogenes list available. Analyzing all significant genes.")
    # Filter for significant genes (adjusted p-value < 0.05)
    sig_genes <- results_df[results_df$padj < 0.05 & !is.na(results_df$padj), ]
    sig_genes <- sig_genes[order(sig_genes$padj), ]
    message(paste0("Found ", nrow(sig_genes), " significant genes (padj < 0.05)"))
    return(list(all_sig_genes = sig_genes, glyco_sig_genes = NULL))
  }
  
  # Filter for glycogenes
  glycogene_results <- results_df[rownames(results_df) %in% glycogenes, ]
  message(paste0("Found ", nrow(glycogene_results), " glycogenes in the dataset"))
  
  # Filter for significant glycogenes
  sig_glycogenes <- glycogene_results[glycogene_results$padj < 0.05 & !is.na(glycogene_results$padj), ]
  sig_glycogenes <- sig_glycogenes[order(sig_glycogenes$padj), ]
  message(paste0("Found ", nrow(sig_glycogenes), " significantly differentially expressed glycogenes (padj < 0.05)"))
  
  # Filter for all significant genes for comparison
  sig_genes <- results_df[results_df$padj < 0.05 & !is.na(results_df$padj), ]
  sig_genes <- sig_genes[order(sig_genes$padj), ]
  
  # Save filtered results
  write.csv(glycogene_results, file.path(processed_data_dir, "all_glycogenes_results.csv"))
  write.csv(sig_glycogenes, file.path(processed_data_dir, "significant_glycogenes.csv"))
  
  return(list(all_sig_genes = sig_genes, glyco_sig_genes = sig_glycogenes, all_glycogenes = glycogene_results))
}

# Function to analyze correlation with DUX4 scores
analyze_dux4_correlation <- function(dds, dux4_scores, filtered_results) {
  if (is.null(dux4_scores)) {
    message("DUX4 scores not available. Skipping correlation analysis.")
    return(NULL)
  }
  
  message("Analyzing correlation between gene expression and DUX4 scores...")
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Extract DUX4 scores and ensure they match the samples
  scores <- dux4_scores$dux4_score
  names(scores) <- dux4_scores$sample
  
  # Only use samples with DUX4 scores that are also in our count data
  shared_samples <- intersect(colnames(norm_counts), names(scores))
  valid_samples <- !is.na(scores[shared_samples])
  valid_samples <- shared_samples[valid_samples]
  
  if (length(valid_samples) < 3) {
    message("Too few samples with valid DUX4 scores. Skipping correlation analysis.")
    return(NULL)
  }
  
  # Function to compute correlation for a set of genes
  compute_correlations <- function(gene_set, label) {
    if (is.null(gene_set) || nrow(gene_set) == 0) {
      message(paste0("No ", label, " to analyze for DUX4 correlation."))
      return(NULL)
    }
    
    # Compute correlation for each gene
    cor_results <- data.frame(
      gene = rownames(gene_set),
      log2FoldChange = gene_set$log2FoldChange,
      padj = gene_set$padj,
      correlation = NA
    )
    
    for (i in 1:nrow(cor_results)) {
      gene <- cor_results$gene[i]
      if (gene %in% rownames(norm_counts)) {
        expr <- norm_counts[gene, valid_samples]
        cor_results$correlation[i] <- cor(expr, scores[valid_samples], method = "spearman")
      }
    }
    
    # Remove NA correlations
    cor_results <- cor_results[!is.na(cor_results$correlation), ]
    
    # Sort by correlation (absolute value)
    cor_results <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]
    
    # Save correlation results
    write.csv(cor_results, file.path(processed_data_dir, paste0("dux4_correlated_", label, ".csv")))
    
    message(paste0("Correlation analysis with DUX4 scores completed for ", label))
    message(paste0("Results saved to ", file.path(processed_data_dir, paste0("dux4_correlated_", label, ".csv"))))
    
    return(cor_results)
  }
  
  # Compute correlations for different gene sets
  all_sig_corr <- compute_correlations(filtered_results$all_sig_genes, "significant_genes")
  glyco_sig_corr <- compute_correlations(filtered_results$glyco_sig_genes, "significant_glycogenes")
  all_glyco_corr <- compute_correlations(filtered_results$all_glycogenes, "all_glycogenes")
  
  return(list(
    all_sig_corr = all_sig_corr,
    glyco_sig_corr = glyco_sig_corr,
    all_glyco_corr = all_glyco_corr
  ))
}

# Function to create enhanced visualizations
create_enhanced_visualizations <- function(deseq_results, filtered_results, correlation_results) {
  message("Creating enhanced visualizations...")
  
  # 1. Create volcano plot highlighting glycogenes
  create_volcano_plot(deseq_results$results, filtered_results)
  
  # 2. Create heatmap of significant glycogenes
  if (!is.null(filtered_results$glyco_sig_genes) && nrow(filtered_results$glyco_sig_genes) > 0) {
    create_glycogene_heatmap(deseq_results$dds, filtered_results$glyco_sig_genes)
  }
  
  # 3. Create DUX4 correlation plots if available
  if (!is.null(correlation_results)) {
    create_correlation_plots(correlation_results)
  }
  
  # 4. Create bar plots for top genes
  create_bar_plots(filtered_results)
}

# Function to create volcano plot highlighting glycogenes
create_volcano_plot <- function(results_df, filtered_results) {
  # Filter NA values
  results_filtered <- results_df[!is.na(results_df$padj), ]
  
  # Create classification for plotting
  results_filtered$category <- "Not Significant"
  results_filtered$category[results_filtered$padj < 0.05 & abs(results_filtered$log2FoldChange) >= 1] <- "Significant"
  
  # Mark glycogenes if available
  if (!is.null(filtered_results$all_glycogenes)) {
    results_filtered$category[rownames(results_filtered) %in% rownames(filtered_results$all_glycogenes)] <- "Glycogene"
    results_filtered$category[rownames(results_filtered) %in% rownames(filtered_results$glyco_sig_genes)] <- "Significant Glycogene"
  }
  
  # Create color palette
  color_palette <- c(
    "Not Significant" = "gray",
    "Significant" = "blue",
    "Glycogene" = "orange",
    "Significant Glycogene" = "red"
  )
  
  # Create volcano plot
  p <- ggplot(results_filtered, aes(x = log2FoldChange, y = -log10(pvalue), color = category)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = color_palette) +
    theme_minimal() +
    labs(
      title = "Volcano Plot: FSHD vs Control",
      subtitle = "Highlighting Glycogenes",
      x = "log2 Fold Change",
      y = "-log10 p-value",
      color = "Category"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.8) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # Add labels for top glycogenes if available
  if (!is.null(filtered_results$glyco_sig_genes) && nrow(filtered_results$glyco_sig_genes) > 0) {
    # Label top 10 significant glycogenes
    top_n <- min(10, nrow(filtered_results$glyco_sig_genes))
    top_glycogenes <- filtered_results$glyco_sig_genes[1:top_n, ]
    
    # Add labels
    p <- p + ggrepel::geom_text_repel(
      data = results_filtered[rownames(results_filtered) %in% rownames(top_glycogenes), ],
      aes(label = gene),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 20,
      size = 3
    )
  }
  
  # Save the plot
  ggsave(file.path(figures_dir, "volcano_plot_fshd_glycogenes.png"), p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(figures_dir, "volcano_plot_fshd_glycogenes.pdf"), p, width = 10, height = 8)
  
  message("Volcano plot created and saved")
}

# Function to create heatmap of significant glycogenes
create_glycogene_heatmap <- function(dds, sig_glycogenes) {
  # Extract normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Extract counts for significant glycogenes
  glyco_counts <- norm_counts[rownames(norm_counts) %in% rownames(sig_glycogenes), ]
  
  # Z-score transformation for better visualization
  z_scores <- t(scale(t(glyco_counts)))
  
  # Get sample metadata
  sample_data <- as.data.frame(colData(dds))
  
  # Create annotation for samples
  annotation_col <- data.frame(
    Condition = sample_data$condition,
    row.names = rownames(sample_data)
  )
  
  # Color palette for annotation
  ann_colors <- list(
    Condition = c(control = "#1B9E77", FSHD = "#D95F02")
  )
  
  # Create heatmap
  pdf(file.path(figures_dir, "heatmap_significant_glycogenes.pdf"), width = 10, height = 12)
  pheatmap(
    z_scores,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = "Expression of Significant Glycogenes in FSHD",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    fontsize_row = 8
  )
  dev.off()
  
  # Also save as PNG for easier viewing
  png(file.path(figures_dir, "heatmap_significant_glycogenes.png"), width = 1000, height = 1200, res = 120)
  pheatmap(
    z_scores,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = "Expression of Significant Glycogenes in FSHD",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    fontsize_row = 8
  )
  dev.off()
  
  message("Heatmap of significant glycogenes created and saved")
}

# Function to create correlation plots
create_correlation_plots <- function(correlation_results) {
  # Create scatter plot for glycogene correlation with DUX4 scores
  if (!is.null(correlation_results$glyco_sig_corr) && nrow(correlation_results$glyco_sig_corr) > 0) {
    # Get top correlated glycogenes
    top_n <- min(20, nrow(correlation_results$glyco_sig_corr))
    top_corr <- correlation_results$glyco_sig_corr[1:top_n, ]
    
    # Create color based on correlation direction
    top_corr$color <- ifelse(top_corr$correlation > 0, "Positive", "Negative")
    
    # Create bar plot of correlations
    p1 <- ggplot(top_corr, aes(x = reorder(gene, correlation), y = correlation, fill = color)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Positive" = "#D55E00", "Negative" = "#0072B2")) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "Glycogenes Correlated with DUX4 Expression",
        x = "Gene",
        y = "Spearman Correlation",
        fill = "Correlation"
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    # Save the plot
    ggsave(file.path(figures_dir, "glycogene_dux4_correlation.png"), p1, width = 10, height = 8, dpi = 300)
    ggsave(file.path(figures_dir, "glycogene_dux4_correlation.pdf"), p1, width = 10, height = 8)
    
    message("Glycogene-DUX4 correlation plot created and saved")
  }
  
  # Create scatter plot comparing fold change vs correlation
  if (!is.null(correlation_results$all_glyco_corr) && nrow(correlation_results$all_glyco_corr) > 0) {
    # Create scatter plot
    p2 <- ggplot(correlation_results$all_glyco_corr, aes(x = log2FoldChange, y = correlation)) +
      geom_point(aes(color = padj < 0.05), alpha = 0.7, size = 2) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
      theme_minimal() +
      labs(
        title = "Glycogene Expression: Fold Change vs DUX4 Correlation",
        x = "log2 Fold Change (FSHD vs Control)",
        y = "Correlation with DUX4 Score",
        color = "Significant"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    # Label top correlated genes
    top_n_pos <- correlation_results$all_glyco_corr %>% 
      filter(correlation > 0) %>% 
      arrange(desc(correlation)) %>% 
      head(5)
    
    top_n_neg <- correlation_results$all_glyco_corr %>% 
      filter(correlation < 0) %>% 
      arrange(correlation) %>% 
      head(5)
    
    top_to_label <- rbind(top_n_pos, top_n_neg)
    
    # Add labels
    p2 <- p2 + ggrepel::geom_text_repel(
      data = top_to_label,
      aes(label = gene),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 20,
      size = 3
    )
    
    # Save the plot
    ggsave(file.path(figures_dir, "glycogene_fold_change_vs_correlation.png"), p2, width = 10, height = 8, dpi = 300)
    ggsave(file.path(figures_dir, "glycogene_fold_change_vs_correlation.pdf"), p2, width = 10, height = 8)
    
    message("Fold change vs correlation plot created and saved")
  }
}

# Function to create bar plots for top genes
create_bar_plots <- function(filtered_results) {
  # Create bar plot for top significant glycogenes
  if (!is.null(filtered_results$glyco_sig_genes) && nrow(filtered_results$glyco_sig_genes) > 0) {
    # Get top glycogenes (limit to 20)
    top_n <- min(20, nrow(filtered_results$glyco_sig_genes))
    top_genes <- filtered_results$glyco_sig_genes[1:top_n, ]
    
    # Create direction variable for coloring
    top_genes$direction <- ifelse(top_genes$log2FoldChange > 0, "Up in FSHD", "Down in FSHD")
    
    # Create the plot
    p <- ggplot(top_genes, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = direction)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Up in FSHD" = "#D55E00", "Down in FSHD" = "#0072B2")) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste0("Top ", top_n, " Differentially Expressed Glycogenes in FSHD"),
        x = "Gene",
        y = "log2 Fold Change",
        fill = "Direction"
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      geom_text(
        aes(label = sprintf("p=%.2e", padj), hjust = ifelse(log2FoldChange > 0, -0.1, 1.1)),
        size = 3
      )
    
    # Save the plot
    ggsave(file.path(figures_dir, "top_glycogenes_fold_change.png"), p, width = 10, height = 8, dpi = 300)
    ggsave(file.path(figures_dir, "top_glycogenes_fold_change.pdf"), p, width = 10, height = 8)
    
    message("Bar plot of top glycogenes created and saved")
  }
}

# Main function
main <- function() {
  message("Starting FSHD RNAseq analysis workflow...")
  
  # Load FSHD dataset
  dataset <- load_fshd_dataset()
  
  # Run DESeq2 analysis
  deseq_results <- run_deseq2_analysis(dataset$counts, dataset$metadata)
  
  # Filter for glycogenes
  filtered_results <- filter_for_glycogenes(deseq_results$results, deseq_results$dds)
  
  # Analyze correlation with DUX4 scores if available
  correlation_results <- analyze_dux4_correlation(deseq_results$dds, dataset$dux4_scores, filtered_results)
  
  # Create enhanced visualizations
  create_enhanced_visualizations(deseq_results, filtered_results, correlation_results)
  
  message("FSHD RNAseq analysis workflow completed!")
  message("Results have been saved to the processed_data and figures directories.")
}

# Run the main function
main() 