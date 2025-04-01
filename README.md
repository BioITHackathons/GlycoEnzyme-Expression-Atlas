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

## Current Progress

### Pre-analysis Phase (Completed)
1. **Data Collection and Organization**
   - Established project structure following FAIR principles
   - Collected and organized glycoenzyme datasets:
     - `human_protein_glycohydrolase.csv`: Glycohydrolase proteins
     - `human_protein_glycotransferase.csv`: Glycotransferase proteins
     - `human_protein_glycogenes.csv`: Glycogenes and alternative names

2. **Gene Analysis**
   - Developed and implemented `find_common_genes.py` for cross-dataset analysis
   - Identified common genes across datasets:
     - 74 direct matches between hydrolase and glycogene lists
     - 205 direct matches between transferase and glycogene lists
     - Additional matches through alternative gene names
   - Generated detailed comparison reports and unique gene lists

### Next Steps
1. **Data Collection**
   - Download and process RNA-seq count matrices for identified genes
   - Prepare metadata files for sample information
   - Validate data quality and completeness

2. **Analysis Pipeline Development**
   - Implement DESeq2/EdgeR analysis pipeline
   - Develop pathway annotation integration
   - Create network analysis workflows

3. **Visualization and Reporting**
   - Generate publication-ready plots
   - Create interactive network visualizations
   - Develop comprehensive analysis reports

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