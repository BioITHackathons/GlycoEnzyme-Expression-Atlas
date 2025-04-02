#!/usr/bin/env Rscript

# Glycogene Pathway Analysis for FSHD RNAseq Data
# 
# This script:
# 1. Loads the significant glycogenes from the DESeq2 analysis
# 2. Performs pathway enrichment analysis using clusterProfiler
# 3. Creates visualizations of enriched pathways
# 4. Investigates glycosylation-specific pathways
#
# Requirements:
# - R 4.0+
# - Bioconductor packages: clusterProfiler, org.Hs.eg.db, enrichplot, DOSE

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ReactomePA)
})

# Set up directories
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
processed_data_dir <- file.path(data_dir, "processed")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(results_dir, "figures")
pathway_dir <- file.path(results_dir, "pathway_analysis")

# Create directories if they don't exist
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

# Function to load glycogene results
load_glycogene_results <- function() {
  # Check if the significant glycogenes file exists
  sig_glycogenes_file <- file.path(processed_data_dir, "significant_glycogenes.csv")
  if (!file.exists(sig_glycogenes_file)) {
    stop("Significant glycogenes file not found. Please run rnaseq_deseq_analysis.R first.")
  }
  
  # Load significant glycogenes
  sig_glycogenes <- read.csv(sig_glycogenes_file, row.names = 1)
  message(paste0("Loaded ", nrow(sig_glycogenes), " significant glycogenes"))
  
  # Also load all glycogenes for comparison
  all_glycogenes_file <- file.path(processed_data_dir, "all_glycogenes_results.csv")
  if (file.exists(all_glycogenes_file)) {
    all_glycogenes <- read.csv(all_glycogenes_file, row.names = 1)
    message(paste0("Loaded ", nrow(all_glycogenes), " total glycogenes"))
  } else {
    all_glycogenes <- NULL
  }
  
  # Also load DUX4 correlated glycogenes if available
  corr_glycogenes_file <- file.path(processed_data_dir, "dux4_correlated_all_glycogenes.csv")
  if (file.exists(corr_glycogenes_file)) {
    corr_glycogenes <- read.csv(corr_glycogenes_file, row.names = 1)
    message(paste0("Loaded ", nrow(corr_glycogenes), " DUX4-correlated glycogenes"))
  } else {
    corr_glycogenes <- NULL
  }
  
  return(list(
    sig_glycogenes = sig_glycogenes,
    all_glycogenes = all_glycogenes,
    corr_glycogenes = corr_glycogenes
  ))
}

# Function to convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_list) {
  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  message(paste0("Converted ", nrow(gene_ids), " genes to Entrez IDs out of ", length(gene_list), " total"))
  
  return(gene_ids)
}

# Function to perform pathway enrichment analysis
run_pathway_analysis <- function(gene_ids, label) {
  message(paste0("Running pathway analysis for ", label, "..."))
  
  # Set up results list
  enrichment_results <- list()
  
  # GO Biological Process enrichment
  go_bp <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  enrichment_results$GO_BP <- go_bp
  
  if (nrow(go_bp) > 0) {
    message(paste0("Found ", nrow(go_bp), " enriched GO biological processes"))
    # Save results
    write.csv(as.data.frame(go_bp), file.path(pathway_dir, paste0(label, "_GO_BP.csv")))
    
    # Create dotplot
    p1 <- dotplot(go_bp, showCategory = 20, title = paste0(label, ": GO Biological Process"))
    ggsave(file.path(figures_dir, paste0(label, "_GO_BP_dotplot.png")), p1, width = 10, height = 8, dpi = 300)
    
    # Create cnetplot
    if (nrow(go_bp) >= 5) {
      p2 <- cnetplot(go_bp, categorySize = "pvalue", foldChange = NULL, showCategory = 5)
      ggsave(file.path(figures_dir, paste0(label, "_GO_BP_cnetplot.png")), p2, width = 12, height = 10, dpi = 300)
    }
  } else {
    message("No enriched GO biological processes found")
  }
  
  # GO Molecular Function enrichment
  go_mf <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  enrichment_results$GO_MF <- go_mf
  
  if (nrow(go_mf) > 0) {
    message(paste0("Found ", nrow(go_mf), " enriched GO molecular functions"))
    # Save results
    write.csv(as.data.frame(go_mf), file.path(pathway_dir, paste0(label, "_GO_MF.csv")))
    
    # Create dotplot
    p3 <- dotplot(go_mf, showCategory = 20, title = paste0(label, ": GO Molecular Function"))
    ggsave(file.path(figures_dir, paste0(label, "_GO_MF_dotplot.png")), p3, width = 10, height = 8, dpi = 300)
  } else {
    message("No enriched GO molecular functions found")
  }
  
  # KEGG pathway enrichment
  kegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  enrichment_results$KEGG <- kegg
  
  if (nrow(kegg) > 0) {
    message(paste0("Found ", nrow(kegg), " enriched KEGG pathways"))
    # Convert to readable gene symbols
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    
    # Save results
    write.csv(as.data.frame(kegg), file.path(pathway_dir, paste0(label, "_KEGG.csv")))
    
    # Create dotplot
    p4 <- dotplot(kegg, showCategory = 20, title = paste0(label, ": KEGG Pathways"))
    ggsave(file.path(figures_dir, paste0(label, "_KEGG_dotplot.png")), p4, width = 10, height = 8, dpi = 300)
    
    # Create pathway plot for top pathway
    if (nrow(kegg) > 0) {
      tryCatch({
        top_pathway_id <- kegg$ID[1]
        p5 <- pathview::pathview(
          gene.data = gene_ids$ENTREZID,
          pathway.id = top_pathway_id,
          species = "hsa",
          out.suffix = paste0(label, "_top_pathway"),
          kegg.dir = pathway_dir
        )
        message(paste0("Created pathway visualization for top KEGG pathway: ", kegg$Description[1]))
      }, error = function(e) {
        message("Error creating pathway visualization: ", e$message)
      })
    }
  } else {
    message("No enriched KEGG pathways found")
  }
  
  # Reactome pathway enrichment
  reactome <- enrichPathway(
    gene = gene_ids$ENTREZID,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  enrichment_results$Reactome <- reactome
  
  if (nrow(reactome) > 0) {
    message(paste0("Found ", nrow(reactome), " enriched Reactome pathways"))
    # Convert to readable gene symbols
    reactome <- setReadable(reactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    
    # Save results
    write.csv(as.data.frame(reactome), file.path(pathway_dir, paste0(label, "_Reactome.csv")))
    
    # Create dotplot
    p6 <- dotplot(reactome, showCategory = 20, title = paste0(label, ": Reactome Pathways"))
    ggsave(file.path(figures_dir, paste0(label, "_Reactome_dotplot.png")), p6, width = 10, height = 8, dpi = 300)
  } else {
    message("No enriched Reactome pathways found")
  }
  
  # Return all enrichment results
  return(enrichment_results)
}

# Function to analyze glycosylation-specific pathways
analyze_glycosylation_pathways <- function(glycogene_results, enrichment_results) {
  message("Analyzing glycosylation-specific pathways...")
  
  # Define glycosylation-related terms to search for in pathway results
  glyco_terms <- c(
    "glyco", "glycan", "glycosylation", "carbohydrate", "sugar", "sialic", 
    "fucosylation", "mannose", "galactose", "GalNAc", "GlcNAc", "lectin"
  )
  
  # Function to filter enrichment results for glycosylation-related pathways
  filter_glyco_pathways <- function(enrichment_result, result_name) {
    if (is.null(enrichment_result) || nrow(enrichment_result) == 0) {
      return(NULL)
    }
    
    # Search for glycosylation terms in pathway descriptions
    glyco_pathways <- enrichment_result[grep(paste(glyco_terms, collapse = "|"), 
                                           enrichment_result$Description, 
                                           ignore.case = TRUE), ]
    
    if (nrow(glyco_pathways) > 0) {
      message(paste0("Found ", nrow(glyco_pathways), " glycosylation-related pathways in ", result_name))
      write.csv(as.data.frame(glyco_pathways), 
               file.path(pathway_dir, paste0("glycosylation_pathways_", result_name, ".csv")))
      
      # Create a plot if there are results
      if (nrow(glyco_pathways) > 0) {
        p <- dotplot(glyco_pathways, showCategory = nrow(glyco_pathways), 
                    title = paste0("Glycosylation-Related Pathways (", result_name, ")"))
        ggsave(file.path(figures_dir, paste0("glycosylation_pathways_", result_name, ".png")), 
              p, width = 10, height = 8, dpi = 300)
      }
      
      return(glyco_pathways)
    } else {
      message(paste0("No glycosylation-related pathways found in ", result_name))
      return(NULL)
    }
  }
  
  # Extract glycosylation-related pathways from each enrichment result
  glyco_pathways <- list()
  
  for (result_name in names(enrichment_results)) {
    for (db in c("GO_BP", "GO_MF", "KEGG", "Reactome")) {
      if (!is.null(enrichment_results[[result_name]][[db]])) {
        glyco_pathways[[paste0(result_name, "_", db)]] <- 
          filter_glyco_pathways(enrichment_results[[result_name]][[db]], paste0(result_name, "_", db))
      }
    }
  }
  
  # Create a combined visualization of all glycosylation pathways found
  combined_results <- data.frame()
  
  for (name in names(glyco_pathways)) {
    if (!is.null(glyco_pathways[[name]]) && nrow(glyco_pathways[[name]]) > 0) {
      # Extract the data we need
      temp <- as.data.frame(glyco_pathways[[name]])
      temp$Source <- name
      temp <- temp[, c("Description", "pvalue", "p.adjust", "Count", "Source")]
      combined_results <- rbind(combined_results, temp)
    }
  }
  
  if (nrow(combined_results) > 0) {
    # Save combined results
    write.csv(combined_results, file.path(pathway_dir, "all_glycosylation_pathways.csv"))
    
    # Create a combined plot
    combined_results$NegLog10PAdj <- -log10(combined_results$p.adjust)
    combined_results$Pathway <- paste0(combined_results$Description, " (", combined_results$Source, ")")
    
    # Keep only top 30 pathways if there are more
    if (nrow(combined_results) > 30) {
      combined_results <- combined_results[order(combined_results$p.adjust), ][1:30, ]
    }
    
    # Reorder for plotting
    combined_results$Pathway <- factor(combined_results$Pathway, 
                                     levels = combined_results$Pathway[order(combined_results$NegLog10PAdj)])
    
    p <- ggplot(combined_results, aes(x = NegLog10PAdj, y = Pathway, size = Count, color = Source)) +
      geom_point() +
      theme_minimal() +
      labs(
        title = "Glycosylation-Related Pathways",
        x = "-log10(adjusted p-value)",
        y = NULL,
        size = "Gene Count",
        color = "Source"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave(file.path(figures_dir, "all_glycosylation_pathways.png"), p, width = 12, height = 10, dpi = 300)
    ggsave(file.path(figures_dir, "all_glycosylation_pathways.pdf"), p, width = 12, height = 10)
    
    message("Created combined visualization of all glycosylation-related pathways")
  } else {
    message("No glycosylation-related pathways found in any enrichment results")
  }
  
  return(glyco_pathways)
}

# Main function
main <- function() {
  message("Starting glycogene pathway analysis for FSHD...")
  
  # Load glycogene results
  glycogene_results <- load_glycogene_results()
  
  # Initialize enrichment results list
  enrichment_results <- list()
  
  # Run pathway analysis for significant glycogenes
  if (!is.null(glycogene_results$sig_glycogenes) && nrow(glycogene_results$sig_glycogenes) > 0) {
    sig_entrez <- convert_to_entrez(rownames(glycogene_results$sig_glycogenes))
    enrichment_results$sig <- run_pathway_analysis(sig_entrez, "significant_glycogenes")
  }
  
  # Run pathway analysis for DUX4-correlated glycogenes if available
  if (!is.null(glycogene_results$corr_glycogenes) && nrow(glycogene_results$corr_glycogenes) > 0) {
    # Filter for significantly correlated genes
    sig_corr <- glycogene_results$corr_glycogenes[abs(glycogene_results$corr_glycogenes$correlation) > 0.5, ]
    
    if (nrow(sig_corr) > 0) {
      corr_entrez <- convert_to_entrez(sig_corr$gene)
      enrichment_results$corr <- run_pathway_analysis(corr_entrez, "dux4_correlated_glycogenes")
    }
  }
  
  # Analyze glycosylation-specific pathways
  glyco_pathways <- analyze_glycosylation_pathways(glycogene_results, enrichment_results)
  
  message("Glycogene pathway analysis completed!")
  message("Results have been saved to the pathway_analysis directory.")
}

# Run the main function
main() 