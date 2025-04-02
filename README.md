# GlycoEnzyme Expression Atlas: Linking Differential Expression to Pathway Dysregulation

A comprehensive bioinformatics project focused on establishing connections between glycoenzyme expression patterns and pathway dysregulation across various disease states. This project is part of the CFDE/GlyGen initiative.

## Project Overview

The GlycoEnzyme Expression Atlas integrates multiple data types and analytical approaches to understand glycoenzyme expression patterns:
- RNA-seq data preprocessing using DESeq2/EdgeR for differential expression
- Mapping of glycoenzyme genes using CAZy and GlyGen databases
- Integration with KEGG/Reactome pathway annotations
- Network analysis via Cytoscape/STRING for interaction mapping

## Team Members
- **Vlado Dancik, PhD** (Team Lead) - Computational Chemical Biologist, Broad Institute
- **Aymen Maqsood** - Team Member
- **Isha Parikh** - Team Member
- **Chandani Shrestha** - Team Member

## Hackathon Project 

### Project Details
- **Project Name**: GlycoEnzyme Expression Atlas
- **Institution**: CFDE/GlyGen
- **Team Lead**: Vlado Dancik, PhD, Computational Chemical Biologist, Broad Institute

### Project Goals
The GlycoEnzyme Expression Atlas project aims to establish connections between glycoenzyme expression patterns and pathway dysregulation across various disease states. This bioinformatics initiative integrates:
- RNA-seq data preprocessing using DESeq2/EdgeR for differential expression
- Mapping of glycoenzyme genes using CAZy and GlyGen databases
- Integration with KEGG/Reactome pathway annotations
- Network analysis via Cytoscape/STRING for interaction mapping

### Community Impact
This project is particularly valuable for the research community as it:
- Provides a comprehensive analysis of glycogene expression patterns in disease states
- Offers interactive visualization tools for exploring complex biological data
- Integrates multiple data sources and analytical approaches
- Facilitates understanding of glycosylation-related disease mechanisms

### FAIR Principles Implementation
The project follows FAIR (Findable, Accessible, Interoperable, and Reusable) principles:
- **Findable**: Data accessible through public repositories (GlyGen, GEO, EMBL-EBI)
- **Accessible**: All datasets and tools freely available on GitHub and Zenodo
- **Interoperable**: Uses standardized ontologies (EDAM, OBI, Human Disease Ontology)
- **Reusable**: Well-documented pipelines with open-source licensing (CC-BY-4.0)

## Facioscapulohumeral Muscular Dystrophy (FSHD) Analysis

### Disease Overview
FSHD is a genetic muscle disorder characterized by progressive muscle weakness and atrophy, primarily affecting the face, shoulders, and upper arms. The disease is caused by the misregulation of the DUX4 gene, leading to muscle degeneration.

### Dataset Information
- **GEO Accession**: GSE140261
- **Sample Size**: 24 samples (12 FSHD, 12 controls)
- **Tissue Type**: Skeletal muscle biopsies
- **Technology**: RNA-seq

### Key Findings in FSHD
1. **Glycogene Expression Changes**:
   - 856 glycogenes identified in the dataset
   - 449 significantly differentially expressed glycogenes
   - Strong up-regulation of GALNT13 (log2FC = 6.01) suggesting altered O-glycosylation
   - Down-regulation of B3GALT1 (log2FC = -1.49) indicating N-glycosylation defects

2. **Pathway Alterations**:
   - Mitochondrial dysfunction (oxidative phosphorylation, electron transport)
   - Immune system activation (cytokine signaling, complement system)
   - RNA processing changes (spliceosome, tRNA biosynthesis)

3. **Biological Implications**:
   - Impaired energy metabolism through mitochondrial dysfunction
   - Altered protein glycosylation affecting muscle function
   - Enhanced inflammatory response contributing to muscle damage

### Interactive Analysis Tools
The Shiny dashboard provides interactive access to:
- Differential expression analysis results
- Pathway enrichment visualizations
- Network analysis of glycogene interactions
- Statistical summaries and gene expression patterns

### 1. Differential Expression Analysis

#### Dataset Characteristics
- RNA-seq data from GSE140261
- 12 FSHD samples vs 12 control samples
- Total genes analyzed: 9,520 significantly differentially expressed genes (padj < 0.05)

#### Glycogene Expression Patterns
- Total glycogenes identified: 856
- Significantly differentially expressed glycogenes: 449
  - Up-regulated: 368 genes
  - Down-regulated: 81 genes

#### Notable Differentially Expressed Glycogenes
- Up-regulated:
  - GALNT13 (log2FC = 6.01, padj < 0.05)
  - GYG2 (log2FC = 5.89, padj < 0.05)
  - HAS1 (log2FC = 4.98, padj < 0.05)
- Down-regulated:
  - B3GALT1 (log2FC = -1.49, padj < 0.05)
  - HAS3 (log2FC = -1.19, padj < 0.05)
  - ST3GAL1 (log2FC = -0.99, padj < 0.05)

### 2. Pathway Analysis

#### KEGG Pathway Enrichment
Significantly enriched pathways (padj < 0.05):
1. Oxidative phosphorylation (NES = -4.48, padj = 1.64e-08)
2. Spliceosome (NES = -3.37, padj = 1.64e-08)
3. Cytokine-cytokine receptor interaction (NES = 1.64, padj = 9.53e-07)
4. Complement and coagulation cascades (NES = 1.87, padj = 1.03e-05)
5. Aminoacyl-tRNA biosynthesis (NES = -3.18, padj = 1.51e-05)

#### Reactome Pathway Enrichment
Significantly enriched pathways (padj < 0.05):
1. Respiratory electron transport (NES = -4.95, padj = 1.05e-08)
2. Mitochondrial translation (NES = -4.40, padj = 1.05e-08)
3. Complex I biogenesis (NES = -4.17, padj = 1.05e-08)
4. Mitochondrial biogenesis (NES = -3.54, padj = 1.05e-08)

### 3. Network Analysis

#### Key Findings
- Significant disruption in mitochondrial function and energy metabolism
- Alterations in glycosylation pathways affecting protein modification
- Changes in cytokine signaling and immune response pathways
- Dysregulation of RNA processing and translation machinery

#### Biological Implications
1. Mitochondrial Dysfunction:
   - Down-regulation of oxidative phosphorylation components
   - Impaired mitochondrial translation and biogenesis
   - Potential impact on cellular energy production

2. Glycosylation Alterations:
   - Up-regulation of GALNT13 suggests changes in O-glycosylation
   - Down-regulation of B3GALT1 indicates potential N-glycosylation defects
   - HAS1/HAS3 dysregulation suggests altered hyaluronan synthesis

3. Immune Response:
   - Enhanced cytokine signaling pathways
   - Activation of complement system
   - Potential inflammatory component in FSHD

## Interactive Visualization

The project includes a Shiny dashboard (`app.R`) providing interactive access to:
- Differential expression results
- Pathway enrichment analysis
- Network visualization
- Gene expression patterns
- Statistical summaries

## Data Availability

All analysis results and visualizations are available in the following directories:
- `data/Deseq2_results/`: Differential expression analysis results
- `results/pathway/`: Pathway enrichment analysis
- `results/network/`: Network analysis files
- `results/figures/`: Generated visualizations

## Technical Requirements

- R version 4.2.3 or higher
- Required R packages:
  - DESeq2
  - clusterProfiler
  - ggplot2
  - pheatmap
  - shiny
  - shinydashboard
  - plotly
  - dplyr
  - tidyr

## Installation

1. Clone the repository:
```bash
git clone https://github.com/BioITHackathons/GlycoEnzyme-Expression-Atlas.git
```

2. Install required R packages:
```R
source("scripts/setup/install_r_packages.R")
```

3. Run the Shiny app:
```R
shiny::runApp("app.R")
```

## Citation

If you use this analysis in your research, please cite:
```
GlycoEnzyme Expression Atlas: A comprehensive analysis of glycogene expression patterns in Facioscapulohumeral muscular dystrophy
https://github.com/BioITHackathons/GlycoEnzyme-Expression-Atlas
``` 