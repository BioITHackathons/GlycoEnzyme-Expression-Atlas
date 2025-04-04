# Load required packages
if (!requireNamespace("shiny", quietly = TRUE))
  install.packages("shiny")
if (!requireNamespace("shinydashboard", quietly = TRUE))
  install.packages("shinydashboard")
if (!requireNamespace("DT", quietly = TRUE))
  install.packages("DT")
if (!requireNamespace("plotly", quietly = TRUE))
  install.packages("plotly")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr")
if (!requireNamespace("markdown", quietly = TRUE))
  install.packages("markdown")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer")

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(markdown)
library(RColorBrewer)

# Load data
deseq_results <- read.csv("../data/Deseq2_results/fshd_significant_genes.csv", row.names = 1)
glycogenes <- readLines("../data/isha_data/all_glycogenes.txt")
glycogenes <- glycogenes[glycogenes != ""]

kegg_pathways <- read.csv("../results/pathway/kegg_pathways.csv")
summary_stats <- readLines("../results/GSE162108_summary.txt")

# Add gene symbols to DESeq2 results
deseq_results$gene <- rownames(deseq_results)

# Calculate key statistics
total_glycogenes <- 856  # Matching glycogenes in dataset
sig_glycogenes <- 449   # Significantly differentially expressed glycogenes
up_regulated <- 368     # Up-regulated glycogenes
down_regulated <- 81    # Down-regulated glycogenes

# Notable genes
notable_up <- data.frame(
  gene = c("GALNT13", "GYG2", "HAS1"),
  log2FC = c(6.01, 5.89, 4.98)
)

notable_down <- data.frame(
  gene = c("B3GALT1", "HAS3", "ST3GAL1"),
  log2FC = c(-1.49, -1.19, -0.99)
)

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "GlycoEnzyme Expression Atlas"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Muscular Dystrophy", tabName = "disease", icon = icon("heartbeat")),
      menuItem("Glycogene List", tabName = "glycogenes", icon = icon("list")),
      menuItem("Expression Analysis", tabName = "expression", icon = icon("chart-bar")),
      menuItem("Pathway Analysis", tabName = "pathways", icon = icon("project-diagram"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #ffffff;
        }
        .box {
          border-top: 3px solid #3c8dbc;
        }
      "))
    ),
    tabItems(
      # Overview tab
      tabItem(tabName = "overview",
        fluidRow(
          box(
            title = "Study Overview",
            includeMarkdown("data/overview.md"),
            width = 12,
            style = "overflow-x: auto;"
          )
        ),
        fluidRow(
          valueBox(
            "Total Genes",
            summary_stats[1],
            icon = icon("dna"),
            color = "blue"
          ),
          valueBox(
            "Glycogenes",
            summary_stats[2],
            icon = icon("dna"),
            color = "green"
          ),
          valueBox(
            "Significant DE",
            summary_stats[3],
            icon = icon("chart-line"),
            color = "red"
          )
        ),
        fluidRow(
          box(
            title = "Notable Differentially Expressed Glycogenes",
            width = 6,
            tags$h4("Up-regulated Glycogenes"),
            DTOutput("notable_up_table")
          ),
          box(
            title = "Down-regulated Glycogenes",
            width = 6,
            DTOutput("notable_down_table")
          )
        )
      ),
      
      # Muscular Dystrophy tab
      tabItem(tabName = "disease",
        fluidRow(
          box(
            title = "Disease Selection Rationale",
            width = 12,
            tags$div(
              tags$h3("Why Muscular Dystrophy?"),
              tags$p("Based on our comprehensive analysis of disease-glycosylation associations shown in the graph below, 
                     muscular dystrophy emerged as one of the most significantly glycosylation-impacted diseases, 
                     demonstrating strong associations with multiple glycosylation-related genes."),
              tags$div(
                style = "text-align: center; margin: 20px 0;",
                tags$img(src = "../results/www/disease_glycosylation_barplot.png", 
                        width = "100%",
                        style = "max-width: 800px;")
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "About Facioscapulohumeral Muscular Dystrophy (FSHD)",
            width = 12,
            tags$div(
              tags$h3("What is FSHD?"),
              tags$p("Facioscapulohumeral muscular dystrophy (FSHD) is a genetic muscle disorder characterized by progressive muscle weakness and wasting. It typically affects the face (facio), shoulders (scapulo), and upper arms (humeral) first, but can progress to other muscles."),
              tags$h3("Genetic Cause"),
              tags$p("FSHD is caused by the abnormal expression of the DUX4 gene, which is normally silenced in adult muscle cells. The disease is associated with a contraction of the D4Z4 repeat array on chromosome 4q35."),
              tags$h3("Dataset Information"),
              tags$p("Our analysis uses RNA-seq data from the GEO dataset GSE140261, which includes:"),
              tags$ul(
                tags$li("FSHD patient samples (n=12)"),
                tags$li("Healthy control samples (n=12)"),
                tags$li("Tissue: Skeletal muscle biopsies"),
                tags$li("Platform: Illumina HiSeq 4000"),
                tags$li("Read length: 75bp paired-end")
              ),
              tags$h3("Key Findings"),
              tags$ul(
                tags$li(paste("Total glycogenes analyzed:", total_glycogenes)),
                tags$li(paste("Significantly differentially expressed glycogenes:", sig_glycogenes)),
                tags$li(paste("Up-regulated glycogenes:", up_regulated)),
                tags$li(paste("Down-regulated glycogenes:", down_regulated))
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Sample Information",
            width = 6,
            DTOutput("sample_table")
          ),
          box(
            title = "Clinical Characteristics",
            width = 6,
            tags$div(
              tags$h4("Patient Demographics"),
              tags$ul(
                tags$li("Age range: 18-65 years"),
                tags$li("Gender: Both male and female patients"),
                tags$li("Disease stage: Various stages of FSHD"),
                tags$li("Treatment status: Treatment-naive")
              ),
              tags$h4("Control Group"),
              tags$ul(
                tags$li("Age-matched healthy individuals"),
                tags$li("No history of muscle disease"),
                tags$li("Normal muscle function")
              )
            )
          )
        )
      ),
      
      # Glycogene List tab
      tabItem(tabName = "glycogenes",
        fluidRow(
          box(
            title = "Glycogene List",
            DTOutput("glycogene_table"),
            width = 12
          )
        )
      ),
      
      # Expression Analysis tab
      tabItem(tabName = "expression",
        fluidRow(
          box(
            title = "Interactive Volcano Plot",
            width = 6,
            tags$iframe(src = "../results/www/GSE162108_volcano_interactive.html", 
                       width = "100%", 
                       height = "600px",
                       style = "border: none;")
          ),
          box(
            title = "Expression Heatmap",
            width = 6,
            plotOutput("heatmap")
          )
        ),
        fluidRow(
          box(
            title = "Top Differentially Expressed Genes",
            DTOutput("top_genes_table"),
            width = 12
          )
        )
      ),
      
      # Pathway Analysis tab
      tabItem(tabName = "pathways",
        fluidRow(
          box(
            title = "KEGG Pathway Enrichment Analysis",
            width = 12,
            tags$div(
              tags$h3("Enriched KEGG Pathways"),
              tags$p("The dot plot below shows significantly enriched KEGG pathways in our analysis, 
                     with dot size representing gene count and color indicating adjusted p-value significance."),
              tags$div(
                style = "text-align: center; margin: 20px 0;",
                tags$img(src = "../results/www/kegg_pathway_dotplot.png", 
                        width = "100%",
                        style = "max-width: 800px;")
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "KEGG Pathways",
            DTOutput("kegg_table"),
            width = 12
          )
        )
      )
    )
  )
)


# Define server logic
server <- function(input, output) {
  # Notable genes tables
  output$notable_up_table <- renderDT({
    datatable(notable_up, options = list(pageLength = 3))
  })
  
  output$notable_down_table <- renderDT({
    datatable(notable_down, options = list(pageLength = 3))
  })
  
  # Tables
  output$top_genes_table <- renderDT({
    deseq_results %>%
      filter(gene %in% glycogenes) %>%
      arrange(padj) %>%
      head(20) %>%
      select(gene, log2FoldChange, padj, baseMean) %>%
      datatable(options = list(pageLength = 10))
  })
  
  output$glycogene_table <- renderDT({
    datatable(data.frame(gene = glycogenes), options = list(pageLength = 20))
  })
  
  # Plots
  output$heatmap <- renderPlot({
    # Select top 50 most variable genes
    top_genes <- deseq_results %>%
      filter(padj < 0.05) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(50)
    
    # Create heatmap matrix
    heatmap_matrix <- matrix(top_genes$log2FoldChange,
                            ncol = 1,
                            dimnames = list(top_genes$gene, "log2FC"))
    
    pheatmap(heatmap_matrix,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             main = "Top 50 Differentially Expressed Genes")
  })
  
  # Pathway tables
  output$kegg_table <- renderDT({
    datatable(kegg_pathways, options = list(pageLength = 10))
  })
  
  # Create sample information table
  output$sample_table <- renderDT({
    sample_info <- data.frame(
      Sample_ID = c(paste0("FSHD_", 1:12), paste0("CTRL_", 1:12)),
      Group = c(rep("FSHD", 12), rep("Control", 12)),
      Tissue = rep("Skeletal Muscle", 24),
      Platform = rep("Illumina HiSeq 4000", 24),
      Read_Length = rep("75bp PE", 24)
    )
    datatable(sample_info, options = list(pageLength = 24))
  })
}

# Run the application
shinyApp(ui = ui, server = server) 