# Glycogene Expression Atlas

## Study Overview
This interactive application provides a comprehensive analysis of glycogene expression patterns in Facioscapulohumeral muscular dystrophy (FSHD). The analysis focuses on understanding how glycosylation-related genes are differentially expressed in FSHD patients compared to healthy controls.

## Key Findings
- Analysis of 856 glycogenes in the FSHD dataset
- Identification of 449 significantly differentially expressed glycogenes
- 368 up-regulated and 81 down-regulated glycogenes
- Significant enrichment of glycosylation-related pathways

## Biological Implications
The analysis reveals several important insights into the role of glycosylation in FSHD:
1. Dysregulation of key glycosyltransferases and glycosidases
2. Alterations in glycan biosynthesis pathways
3. Changes in cell surface glycosylation patterns
4. Potential impact on muscle cell function and integrity

## Methodology
1. **Data Collection**
   - RNA-seq data from FSHD patients and controls
   - Comprehensive glycogene reference list from GlyGen database
   - Pathway information from KEGG and Reactome databases

2. **Analysis Pipeline**
   - Differential expression analysis using DESeq2
   - Pathway enrichment analysis
   - Network analysis of gene interactions
   - Visualization of expression patterns and pathways

3. **Quality Control**
   - Filtering for significant differential expression (padj < 0.05)
   - Cross-validation with reference glycogene list
   - Multiple testing correction
   - Pathway enrichment significance thresholds

## Interactive Features
This application provides:
- Interactive volcano plots for differential expression
- Heatmaps of glycogene expression patterns
- Pathway enrichment visualizations
- Network analysis tools for Cytoscape integration
- Detailed gene-level information and statistics 