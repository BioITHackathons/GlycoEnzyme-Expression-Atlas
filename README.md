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
- **Isha Parikh** - Team Member
- **Aymen Maqsood** - Team Member
- **Chandani Shrestha** - Team Member

## Methodology

### 1. Disease-Glycosylation Association Analysis
Through systematic analysis of glycogene-disease associations and glycosylation sites, we identified muscular dystrophy as the most promising disease model for RNA-seq analysis. This selection was based on:

![Disease-Glycosylation Network](results/figures/disease_glycosylation_network.pdf)
*Network visualization of disease-glycosylation associations*

![Disease-Glycosylation Barplot](results/figures/disease_glycosylation_barplot.pdf)
*Distribution of glycosylation-related genes across disease categories*

### 2. RNA-seq Analysis Pipeline
1. **Data Collection and Preprocessing**
   - RNA-seq count matrix preparation
   - Quality control and normalization
   - Metadata organization and validation

2. **Differential Expression Analysis**
   - DESeq2 implementation for RNA-seq analysis
   - Statistical testing and multiple testing correction
   - Fold change and p-value calculations

3. **Pathway Analysis**
   - KEGG pathway enrichment analysis
   - Reactome pathway mapping
   - Gene set enrichment analysis (GSEA)

4. **Network Analysis**
   - STRING protein-protein interaction network construction
   - Cytoscape visualization and analysis
   - Module identification and pathway mapping

### 3. Interactive Visualization
The project includes a Shiny dashboard (`app.R`) providing interactive access to:
- Differential expression results
- Pathway enrichment analysis
- Network visualization
- Gene expression patterns
- Statistical summaries

## Project Structure

```
GlycoEnzyme-Expression-Atlas/
├── data/
│   ├── raw/           # Original, immutable data
│   ├── processed/     # Cleaned and processed data
│   └── interim/       # Intermediate data that has been transformed
├── scripts/
│   ├── scraping/     # Scripts for data collection and web scraping
│   ├── analysis/     # Scripts for data analysis and visualization
│   └── utils/        # Utility functions and helper scripts
├── docs/            # Documentation files
├── notebooks/       # Jupyter notebooks for exploration and analysis
└── results/
    ├── figures/     # Generated graphics and figures
    └── tables/      # Generated tables and CSV files
```

## Prerequisites

- Python 3.x
- R version ≥ 4.0.0
- Bioconductor 3.12 or higher
- Cytoscape 3.8 or higher (optional, for network visualization)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/GlycoEnzyme-Atlas.git
cd GlycoEnzyme-Atlas
```

2. Install required dependencies:
```bash
pip install pandas numpy
```

## FAIR Principles

This project follows FAIR (Findable, Accessible, Interoperable, and Reusable) principles:

- **Findable**: 
  - Data findable through public repositories (GlyGen, GEO, EMBL-EBI Expression Atlas)
  - Standardized metadata and persistent identifiers
  - Integration with existing databases

- **Accessible**: 
  - All datasets, analysis pipelines, and visualization tools freely available
  - Open-access policies on GitHub and Zenodo
  - Clear documentation and usage guidelines

- **Interoperable**: 
  - Adherence to standardized ontologies (EDAM, OBI, Human Disease Ontology)
  - Standard data formats (FASTA, CSV, JSON)
  - Integration with existing bioinformatics tools

- **Reusable**: 
  - Well-documented pipelines and methodologies
  - Open-source licensing (Creative Commons CC-BY-4.0)
  - Comprehensive documentation and examples

## License
This project is licensed under the Creative Commons CC-BY-4.0 license.

## Acknowledgments
- CFDE/GlyGen for project support
- GlyCosmos and CAZy databases
- Bioconductor project and all package maintainers
- Broad Institute for computational resources

## Analysis Scripts

The `scripts/analysis` directory contains:
- `find_common_genes.py`: Script to identify common genes across different datasets, handling alternative gene names

## Getting Started

1. Ensure you have Python 3.x installed
2. Install required dependencies:
   ```bash
   pip install pandas numpy
   ```
3. Run analysis scripts from the project root directory

## Results

Analysis results will be stored in the `results` directory, including:
- Common genes across datasets
- Unique genes for each dataset
- Comparison tables and visualizations

## Contributing

Please read the documentation in the `docs` directory for guidelines on contributing to this project. 

## Shiny App

The Shiny app can be run using the following R code:
```R
shiny::runApp("app.R")
```

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