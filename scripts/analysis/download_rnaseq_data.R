#!/usr/bin/env Rscript

# Download FSHD RNAseq dataset for GlycoEnzyme-Expression-Atlas
# 
# This script:
# 1. Downloads the FSHD (facioscapulohumeral dystrophy) RNAseq dataset GSE140261 from GEO
# 2. Processes it into a format suitable for DESeq2 analysis
# 3. Saves count matrix and sample metadata to the raw data directory
#
# Requirements:
# - R 4.0+
# - Bioconductor packages: GEOquery, Biobase

# Set up directories
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
raw_data_dir <- file.path(data_dir, "raw")
processed_data_dir <- file.path(data_dir, "processed")

# Create directories if they don't exist
dir.create(raw_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages if not already installed
required_packages <- c("GEOquery", "Biobase", "DESeq2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing package: ", pkg))
    BiocManager::install(pkg)
  }
}

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(DESeq2)
})

# Function to download and process FSHD GEO dataset
download_fshd_dataset <- function(geo_accession = "GSE140261") {
  message(paste0("Processing FSHD dataset ", geo_accession, "..."))
  
  # Set direct file paths based on the actual location
  geo_dir <- file.path(raw_data_dir, geo_accession)
  counts_file <- file.path(geo_dir, paste0(geo_accession, "_year2_genecounts.csv.gz"))
  metadata_file <- file.path(geo_dir, paste0(geo_accession, "_year2_column_data.csv.gz"))
  
  # Check if files already exist
  if (!dir.exists(geo_dir)) {
    dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Only download if needed
  if (!file.exists(counts_file) || !file.exists(metadata_file)) {
    message("Files not found. Downloading from GEO...")
    tryCatch({
      supp_files <- getGEOSuppFiles(geo_accession, fetch_files = TRUE, baseDir = raw_data_dir)
    }, error = function(e) {
      message("Error downloading files: ", e$message)
      message("Using existing files if available...")
    })
  } else {
    message("Files already downloaded. Using existing files.")
  }
  
  # Verify files now exist
  if (!file.exists(counts_file) || !file.exists(metadata_file)) {
    stop("Required files still not available. Please check your internet connection or manually download the files.")
  }
  
  # Read files directly with gzfile to examine the raw format
  message("Checking file structure...")
  counts_conn <- gzfile(counts_file, "rt")
  counts_header <- readLines(counts_conn, n = 2)
  close(counts_conn)
  
  metadata_conn <- gzfile(metadata_file, "rt")
  metadata_header <- readLines(metadata_conn, n = 2)
  close(metadata_conn)
  
  message("Count matrix header: ", counts_header[1])
  message("Metadata header: ", metadata_header[1])
  
  # Read the count matrix with check.names=FALSE to preserve original column names
  message("Reading count matrix...")
  counts <- read.csv(gzfile(counts_file), row.names = 1, check.names = FALSE)
  
  # Read the metadata
  message("Reading sample metadata...")
  metadata <- read.csv(gzfile(metadata_file), row.names = 1, check.names = FALSE)
  
  # Fix: Print some debug information
  message("Count matrix has ", nrow(counts), " genes and ", ncol(counts), " samples")
  message("Count matrix column names: ", paste(head(colnames(counts), 5), collapse=", "), "...")
  
  message("Metadata has ", nrow(metadata), " samples and ", ncol(metadata), " columns")
  message("Metadata columns: ", paste(head(colnames(metadata), 5), collapse=", "), "...")
  
  # The sample names in the count matrix columns should match the sample_name column in metadata
  # Let's create a mapping to check this
  message("Fixing sample names to ensure metadata and counts match...")
  if("sample_name" %in% colnames(metadata)) {
    # Check which sample_name values match the count matrix column names
    count_samples <- colnames(counts)
    metadata_samples <- metadata$sample_name
    
    message("First 5 count_samples: ", paste(head(count_samples, 5), collapse=", "))
    message("First 5 metadata_samples: ", paste(head(metadata_samples, 5), collapse=", "))
    
    shared_samples <- intersect(count_samples, metadata_samples)
    
    message("Count matrix has ", length(count_samples), " sample columns")
    message("Metadata has ", length(metadata_samples), " sample rows")
    message("Shared samples: ", length(shared_samples))
    
    if (length(shared_samples) == 0) {
      # Try to clean up sample names to match
      message("No direct match found. Trying to clean up sample names...")
      
      # Clean up count sample names (remove any 'X' prefix that R might add)
      clean_count_samples <- gsub("^X", "", count_samples)
      clean_count_samples <- gsub("\\.", "-", clean_count_samples)
      
      message("First 5 clean_count_samples: ", paste(head(clean_count_samples, 5), collapse=", "))
      
      # Check for matches again
      shared_samples <- intersect(clean_count_samples, metadata_samples)
      message("After cleanup, shared samples: ", length(shared_samples))
      
      if (length(shared_samples) > 0) {
        # Create a mapping between original column names and cleaned names
        name_mapping <- setNames(count_samples, clean_count_samples)
        
        # Filter metadata to only include samples that exist in the cleaned count matrix names
        metadata <- metadata[metadata$sample_name %in% shared_samples, ]
        
        # Set row names to sample_name for easier matching
        rownames(metadata) <- metadata$sample_name
        
        # Rename count matrix columns to match metadata
        matching_cols <- name_mapping[shared_samples]
        counts <- counts[, matching_cols]
        colnames(counts) <- names(matching_cols)
      }
    } else {
      # Filter metadata to only include samples that exist in the count matrix
      metadata <- metadata[metadata$sample_name %in% shared_samples, ]
      
      # Set row names to sample_name for easier matching
      rownames(metadata) <- metadata$sample_name
      
      # Filter counts to only include samples in metadata
      counts <- counts[, shared_samples]
    }
  }
  
  # Clean up metadata for DESeq2
  # Extract condition information (control vs FSHD)
  # In this dataset, we'll use pheno_type column to identify controls/FSHD
  if("pheno_type" %in% colnames(metadata)) {
    metadata$condition <- ifelse(metadata$pheno_type == "Control", "control", "FSHD")
    message("Using pheno_type to determine condition. Found ", sum(metadata$condition == "control"), 
            " control and ", sum(metadata$condition == "FSHD"), " FSHD samples")
  } else {
    # Fallback to using sample name
    metadata$condition <- ifelse(grepl("control", rownames(metadata), ignore.case = TRUE), 
                                "control", "FSHD")
  }
  
  # Ensure sample names match between counts and metadata
  shared_samples <- intersect(colnames(counts), rownames(metadata))
  if (length(shared_samples) == 0) {
    stop("No matching sample names found between count matrix and metadata")
  }
  
  message(paste0("Found ", length(shared_samples), " matching samples between counts and metadata"))
  
  # Filter to keep only samples present in both datasets
  counts <- counts[, shared_samples]
  metadata <- metadata[shared_samples, ]
  
  # Save processed files
  write.csv(counts, file.path(processed_data_dir, "fshd_counts_matrix.csv"))
  write.csv(metadata, file.path(processed_data_dir, "fshd_sample_metadata.csv"))
  
  message("FSHD dataset processed and saved:")
  message(paste0("  - Count matrix: ", file.path(processed_data_dir, "fshd_counts_matrix.csv")))
  message(paste0("  - Sample metadata: ", file.path(processed_data_dir, "fshd_sample_metadata.csv")))
  
  # Also download the series matrix file to get additional metadata
  message("Downloading series matrix for additional metadata...")
  gse <- getGEO(geo_accession, GSEMatrix = TRUE)
  if (length(gse) > 0) {
    eset <- gse[[1]]
    pheno_data <- pData(eset)
    write.csv(pheno_data, file.path(processed_data_dir, "fshd_complete_metadata.csv"))
    message(paste0("  - Complete metadata: ", file.path(processed_data_dir, "fshd_complete_metadata.csv")))
  }
  
  # Extract DUX4 gene expression data
  # The DUX4 score is included in the metadata as "dux4.rlogsum"
  if ("dux4.rlogsum" %in% colnames(metadata)) {
    message("DUX4 score data found in metadata")
    # Create a file with DUX4 scores for samples
    dux4_data <- data.frame(
      sample = rownames(metadata),
      condition = metadata$condition,
      dux4_score = metadata$dux4.rlogsum,
      dux4_group = metadata$dux4.group
    )
    write.csv(dux4_data, file.path(processed_data_dir, "fshd_dux4_scores.csv"), row.names = FALSE)
    message(paste0("  - DUX4 scores: ", file.path(processed_data_dir, "fshd_dux4_scores.csv")))
  }
  
  # Run initial exploratory analysis to identify potential glycogenes
  explore_dux4_related_genes(counts, metadata)
  
  return(list(counts = counts, metadata = metadata))
}

# Function to identify potential glycogenes related to DUX4 expression
explore_dux4_related_genes <- function(counts, metadata) {
  message("Performing exploratory analysis to identify DUX4-related genes...")
  
  # Check if DUX4 score is available
  if (!"dux4.rlogsum" %in% colnames(metadata)) {
    message("DUX4 score not found in metadata. Skipping exploratory analysis.")
    return(NULL)
  }
  
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
  
  # Convert to dataframe
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Save results
  write.csv(res_df, file.path(processed_data_dir, "fshd_initial_deseq_results.csv"))
  
  # Filter for significant genes
  sig_genes <- res_df[which(res_df$padj < 0.05 & !is.na(res_df$padj)), ]
  sig_genes <- sig_genes[order(sig_genes$padj), ]
  
  # Save significant genes
  write.csv(sig_genes, file.path(processed_data_dir, "fshd_significant_genes.csv"))
  
  message(paste0("Found ", nrow(sig_genes), " significantly differentially expressed genes"))
  message(paste0("Results saved to ", file.path(processed_data_dir, "fshd_significant_genes.csv")))
  
  # Look for correlation with DUX4 scores
  if (nrow(sig_genes) > 0) {
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Calculate correlation with DUX4 scores
    dux4_scores <- metadata$dux4.rlogsum
    
    # Only use samples with DUX4 scores
    valid_samples <- !is.na(dux4_scores)
    if (sum(valid_samples) < 3) {
      message("Too few samples with valid DUX4 scores. Skipping correlation analysis.")
      return(NULL)
    }
    
    # Compute correlation for significant genes
    cor_results <- data.frame(
      gene = rownames(sig_genes),
      log2FoldChange = sig_genes$log2FoldChange,
      padj = sig_genes$padj,
      correlation = NA
    )
    
    for (i in 1:nrow(cor_results)) {
      gene <- cor_results$gene[i]
      expr <- norm_counts[gene, valid_samples]
      cor_results$correlation[i] <- cor(expr, dux4_scores[valid_samples], method = "spearman")
    }
    
    # Sort by correlation (absolute value)
    cor_results <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]
    
    # Save correlation results
    write.csv(cor_results, file.path(processed_data_dir, "fshd_dux4_correlated_genes.csv"))
    
    message(paste0("Correlation analysis with DUX4 scores completed"))
    message(paste0("Results saved to ", file.path(processed_data_dir, "fshd_dux4_correlated_genes.csv")))
    
    # Look for known glycogenes among the top correlated genes
    # This requires a separate reference list of glycogenes
    glycogenes_file <- file.path(data_dir, "isha_data", "all_glycogenes.txt")
    if (file.exists(glycogenes_file)) {
      glycogenes <- readLines(glycogenes_file)
      
      # Find which significant genes are glycogenes
      glyco_sig <- cor_results[cor_results$gene %in% glycogenes, ]
      
      if (nrow(glyco_sig) > 0) {
        # Save glycogene results
        write.csv(glyco_sig, file.path(processed_data_dir, "fshd_significant_glycogenes.csv"))
        
        message(paste0("Found ", nrow(glyco_sig), " significantly differentially expressed glycogenes"))
        message(paste0("Results saved to ", file.path(processed_data_dir, "fshd_significant_glycogenes.csv")))
      } else {
        message("No known glycogenes found among significant genes")
      }
    } else {
      message("Glycogenes list not found. To identify glycogenes, please create a list at:")
      message(paste0(glycogenes_file))
    }
  }
}

# Main function
main <- function() {
  message("Starting download and processing of FSHD RNAseq dataset (GSE140261)...")
  
  # Download and process the dataset
  data <- download_fshd_dataset("GSE140261")
  
  message("FSHD dataset download and processing completed!")
  message("You can now run the rnaseq_deseq_analysis.R script to analyze this data.")
  message("The analysis will include:")
  message("  - Differential expression between FSHD and control samples")
  message("  - Correlation with DUX4 expression")
  message("  - Identification of glycogenes associated with FSHD")
}

# Run the main function
main() 