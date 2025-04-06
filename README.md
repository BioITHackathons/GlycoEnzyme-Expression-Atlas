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

## Citation

If you use this analysis in your research, please cite:
```
GlycoEnzyme Expression Atlas: A comprehensive analysis of glycogene expression patterns in disease states
https://github.com/BioITHackathons/GlycoEnzyme-Expression-Atlas
``` 