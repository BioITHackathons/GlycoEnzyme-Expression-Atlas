#!/usr/bin/env Rscript

# Install required R packages for FSHD RNAseq Analysis
# 
# This script installs all the R packages required for:
# 1. Downloading and processing FSHD RNAseq data (download_rnaseq_data.R)
# 2. Performing DESeq2 analysis (rnaseq_deseq_analysis.R)
# 3. Running pathway enrichment analysis (glycogene_pathway_analysis.R)
#
# Usage: Rscript scripts/setup/install_r_packages.R

# Function to install CRAN packages if not available
install_cran_packages <- function(packages) {
  new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  if (length(new_packages) > 0) {
    message(paste0("Installing the following CRAN packages: ", paste(new_packages, collapse = ", ")))
    install.packages(new_packages, repos = "https://cloud.r-project.org")
  } else {
    message("All required CRAN packages are already installed.")
  }
}

# Function to install Bioconductor packages if not available
install_bioc_packages <- function(packages) {
  # Check if BiocManager is installed, and install if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  
  # Check which packages need to be installed
  new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  if (length(new_packages) > 0) {
    message(paste0("Installing the following Bioconductor packages: ", paste(new_packages, collapse = ", ")))
    BiocManager::install(new_packages)
  } else {
    message("All required Bioconductor packages are already installed.")
  }
}

# Main function
main <- function() {
  message("Installing required R packages for FSHD RNAseq Analysis...")
  
  # Define required CRAN packages
  cran_packages <- c(
    "ggplot2",
    "dplyr", 
    "tidyr",
    "ggrepel",
    "RColorBrewer",
    "pheatmap"
  )
  
  # Define required Bioconductor packages
  bioc_packages <- c(
    "GEOquery",
    "Biobase",
    "DESeq2",
    "biomaRt",
    "clusterProfiler", 
    "org.Hs.eg.db", 
    "enrichplot", 
    "DOSE", 
    "ReactomePA",
    "pathview"
  )
  
  # Install packages
  install_cran_packages(cran_packages)
  install_bioc_packages(bioc_packages)
  
  # Check if any packages failed to install
  missing_packages <- c(
    cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)],
    bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
  )
  
  if (length(missing_packages) > 0) {
    warning(paste0("The following packages could not be installed: ", 
                   paste(missing_packages, collapse = ", ")))
  } else {
    message("All required packages have been successfully installed!")
  }
  
  message("Setup complete. You can now run the following scripts:")
  message("1. scripts/analysis/download_rnaseq_data.R - Download FSHD RNAseq data")
  message("2. scripts/analysis/rnaseq_deseq_analysis.R - Perform DESeq2 analysis")
  message("3. scripts/analysis/glycogene_pathway_analysis.R - Run pathway analysis")
}

# Run the main function
main() 