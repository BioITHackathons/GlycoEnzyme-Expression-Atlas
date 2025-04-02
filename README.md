# GlycoEnzyme Expression Atlas: Linking Differential Expression to Pathway Dysregulation

Despite the critical role glycosylation plays in health and disease, the expression patterns of glycoenzymes and their impact on biological pathways remain poorly understood and underexplored. Existing resources do not effectively link differential glycoenzyme expression to pathway dysregulation across disease states, making it difficult for researchers to uncover disease mechanisms or identify therapeutic targets related to glycosylation.

## Why This Needs to Be Solved ?
 - Glycosylation is involved in key cellular processes, yet is often overlooked in genomic and pathway-level analyses.
 - Glycoenzymes are potential biomarkers and drug targets, but are not widely studied due to lack of integrated resources.
 - Current databases are fragmented, lacking a centralized, open-source platform that connects expression data to functional impact.
 - Solving this will empower researchers and clinicians to uncover glycosylation-related disease mechanisms, aiding precision medicine efforts.

## Objective
To build an open-source, interactive atlas that integrates transcriptomic data, glycoenzyme gene sets, and pathway/network analysis — enabling researchers to explore glycosylation-related dysregulation in diseases, starting with muscular dystrophy

## Workflow
![GlycoEnzyme Expression Atlas Workflow](images/workflow.png)

## Methodology

The GlycoEnzyme Expression Atlas project was conducted in four key phases, each building toward identifying and visualizing glycoenzyme-driven pathway dysregulation in disease states.

### Phase 1: Disease Model Selection
- Identified 906 disease-associated glycosylation-related genes from databases like GlyGen. Data used:
- [GLY_000623 – Human Glycogenes](https://data.glygen.org/GLY_000623)
- [GLY_000004 – Human Glycosyltransferases](https://data.glygen.org/GLY_000004)
- [GLY_000025 – Human Glycoside Hydrolases](https://data.glygen.org/GLY_000025)

Through systematic analysis of glycogene-disease associations and glycosylation sites, we identified muscular dystrophy as the most promising disease model for RNA-seq analysis. This selection was based on:
- [GLY_000225 – Human Protein Glycosylation Diseases](https://data.glygen.org/GLY_000225)
- [GLY_000230 – Human Protein Glycosylation Disease Annotations](https://data.glygen.org/GLY_000230)
- [GLY_000308 – Human Glycogenes with Disease Associations](https://data.glygen.org/GLY_000308)

### Phase 2: RNA-Seq Differential Expression Analysis
**Disease-Glycosylation Association Analysis**
Selected two RNA-seq datasets related to muscular dystrophy from public repositories:
- [GSE140261 – RNA-seq data for Facioscapulohumeral Muscular Dystrophy (FSHD)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140261)
- [PMCID: PMC8756543 – Molecular Features of FSHD through Transcriptomic Analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8756543/)

1. **Data Collection and Preprocessing**
   - RNA-seq count matrix preparation
   - Quality control and normalization
   - Metadata organization and validation

2. **Differential Expression Analysis**
   - DESeq2/Edge R implementation for RNA-seq analysis
   - Statistical testing and multiple testing correction
   - Fold change and p-value calculations
   - These were filtered for glycoenzymes to focus on glycosylation-related dysregulation

### Phase 3: Pathway & Network Analysis
Mapped filtered glycoenzymes to biological pathways using KEGG.
KEGG enrichment analysis using clusterProfiler to identify significantly enriched pathways involving glycoenzymes.

### Phase 4: Interactive Visualization with R Shiny
To make our findings accessible and actionable, we developed an R Shiny application that allows researchers to:
- Input a gene of interest to explore its role in muscular dystrophy.
- View expression patterns and differential expression status.
- Visualize KEGG pathway associations and heatmaps across patients.
- Explore the gene’s network context via STRING-based interactions.
This tool bridges the gap between static analysis and interactive exploration, enabling users to dive deeper into glycosylation-related dysregulation in a user-friendly interface.

## Results

Analysis results will be stored in the `results` directory, including:
- Common genes across datasets
- Unique genes for each dataset
- Comparison tables and visualizations


## Repo Structure
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

## Getting Started

1. Ensure you have Python 3.x installed
2. Install required dependencies:
   ```bash
   pip install pandas numpy
   ```
3. Run analysis scripts from the project root directory
4. The Shiny app can be run using the following R code:
```R
shiny::runApp("app.R")
```
 
## Citation/References
## Acknowledgments
- CFDE/GlyGen for project support
- GlyCosmos and CAZy databases
- Bioconductor project and all package maintainers
- Broad Institute for computational resources
  
