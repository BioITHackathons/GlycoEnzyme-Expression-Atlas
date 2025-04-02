# Load required packages
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(plotly)

# Function to load and prepare data
load_app_data <- function() {
  # Load DESeq2 results
  deseq_results <- read_csv("data/Deseq2_results/fshd_significant_genes.csv")
  
  # Load glycogene list
  glycogenes <- read_csv("data/isha_data/all_glycogenes.txt", col_names = "gene_symbol")
  
  # Load pathway analysis results
  kegg_pathways <- read_csv("results/pathway/kegg_pathways.csv")
  reactome_pathways <- read_csv("results/pathway/reactome_pathways.csv")
  
  # Load network data
  network_data <- read_csv("results/network/pathway_gene_relationships.csv")
  gene_attributes <- read_csv("results/network/gene_attributes.csv")
  
  # Prepare glycogene classification
  glycogene_types <- data.frame(
    type = c(
      "Glycosyltransferase",
      "Glycosidase",
      "Glycan-binding protein",
      "Glycan biosynthesis enzyme",
      "Glycan modification enzyme",
      "Glycan transport protein",
      "Glycan recognition protein",
      "Other"
    ),
    description = c(
      "Enzymes that transfer sugar moieties to acceptor molecules",
      "Enzymes that break down glycans",
      "Proteins that bind to specific glycan structures",
      "Enzymes involved in glycan biosynthesis",
      "Enzymes that modify existing glycan structures",
      "Proteins involved in glycan transport",
      "Proteins that recognize specific glycan structures",
      "Other glycan-related proteins"
    )
  )
  
  # Create gene type mapping (example - you'll need to update this with actual data)
  gene_type_mapping <- data.frame(
    gene_symbol = glycogenes$gene_symbol,
    type = sample(glycogene_types$type, nrow(glycogenes), replace = TRUE)
  )
  
  # Prepare expression data for heatmap
  expression_matrix <- deseq_results %>%
    select(gene_symbol, log2FoldChange, padj) %>%
    spread(gene_symbol, log2FoldChange)
  
  # Create volcano plot data
  volcano_data <- deseq_results %>%
    mutate(
      significant = padj < 0.05,
      regulation = case_when(
        log2FoldChange > 1 & padj < 0.05 ~ "Up-regulated",
        log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Prepare pathway data
  pathway_data <- bind_rows(
    kegg_pathways %>% mutate(source = "KEGG"),
    reactome_pathways %>% mutate(source = "Reactome")
  ) %>%
    arrange(p.adjust)
  
  # Return list of data frames
  list(
    deseq_results = deseq_results,
    glycogenes = glycogenes,
    glycogene_types = glycogene_types,
    gene_type_mapping = gene_type_mapping,
    expression_matrix = expression_matrix,
    volcano_data = volcano_data,
    pathway_data = pathway_data,
    network_data = network_data,
    gene_attributes = gene_attributes
  )
}

# Function to create heatmap
create_heatmap <- function(expression_matrix, selected_genes) {
  if (is.null(selected_genes) || length(selected_genes) == 0) {
    return(NULL)
  }
  
  # Filter matrix for selected genes
  matrix_subset <- expression_matrix %>%
    filter(gene_symbol %in% selected_genes)
  
  # Create heatmap
  pheatmap(matrix_subset,
           scale = "row",
           show_rownames = TRUE,
           show_colnames = TRUE,
           clustering_method = "ward.D2",
           main = "Gene Expression Heatmap")
}

# Function to create volcano plot
create_volcano_plot <- function(volcano_data) {
  ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not significant" = "grey")) +
    theme_minimal() +
    labs(title = "Volcano Plot of Differential Expression",
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)") +
    theme(legend.title = element_blank())
}

# Function to create pathway enrichment plot
create_pathway_plot <- function(pathway_data, source = "KEGG", top_n = 10) {
  pathway_data %>%
    filter(source == !!source) %>%
    head(top_n) %>%
    ggplot(aes(x = reorder(Description, NES), y = NES)) +
    geom_bar(stat = "identity", aes(fill = p.adjust)) +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue") +
    theme_minimal() +
    labs(title = paste("Top", top_n, source, "Pathways"),
         x = "Pathway",
         y = "Normalized Enrichment Score") +
    theme(axis.text.y = element_text(size = 8))
}

# Function to prepare network data for Cytoscape
prepare_network_data <- function(network_data, gene_attributes, selected_pathway = NULL) {
  if (!is.null(selected_pathway)) {
    network_data <- network_data %>%
      filter(pathway_id == selected_pathway)
  }
  
  # Prepare node attributes
  node_attributes <- gene_attributes %>%
    filter(gene_symbol %in% unique(c(network_data$source, network_data$target)))
  
  # Prepare edge attributes
  edge_attributes <- network_data %>%
    select(source, target, pathway_id, pathway_name)
  
  list(
    nodes = node_attributes,
    edges = edge_attributes
  )
} 