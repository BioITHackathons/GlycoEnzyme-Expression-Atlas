#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNAseq DESeq2 Analysis for GlycoEnzyme-Expression-Atlas

This script:
1. Outlines how to select and download a relevant RNAseq dataset from the GlyGen database
2. Performs differential expression analysis using DESeq2 via rpy2
3. Filters results to focus on glycogenes
4. Creates visualizations for the differential expression results

Requirements:
- Python 3.7+
- R 4.0+ with DESeq2, ggplot2, and pheatmap packages installed
- rpy2 (Python interface to R)
- pandas, numpy, matplotlib, seaborn
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from pathlib import Path

# Set up directories
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = BASE_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

# Ensure directories exist
os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

def download_dataset_from_glygen():
    """
    Guidelines for selecting and downloading a dataset from GlyGen:
    
    1. Visit the GlyGen portal (https://www.glygen.org/)
    2. Navigate to the "Data" section
    3. Search for RNAseq datasets related to glycosylation or glycan processing
    4. Focus on datasets that compare:
       - Normal vs disease tissue where glycosylation is implicated
       - Wildtype vs glycosylation enzyme knockout
       - Glycosylation enzyme overexpression experiments
    
    Recommended dataset: GSE147592 - DESeq analysis of RNAseq from human lung cells with/without 
    coronavirus infection (which affects glycosylation pathways)
    
    For this script, we'll assume the data has already been downloaded to:
    - {RAW_DATA_DIR}/counts_matrix.csv (raw count matrix)
    - {RAW_DATA_DIR}/sample_metadata.csv (metadata about samples)
    
    Note: For actual implementation, use the GEOquery package in R or the 
    GEOparse package in Python to download the data directly
    """
    print("Dataset selection guidelines:")
    print("1. Visit the GlyGen portal (https://www.glygen.org/)")
    print("2. Navigate to the 'Data' section")
    print("3. Search for RNAseq datasets related to glycosylation or glycan processing")
    print("4. Select datasets comparing normal vs. disease states, wild-type vs. knockout, etc.")
    print("\nRecommended dataset: GSE147592 - Human lung cells with/without coronavirus infection")
    print("(Coronavirus infection affects glycosylation pathways)")
    print(f"\nData should be downloaded to:")
    print(f"- {RAW_DATA_DIR}/counts_matrix.csv (raw count matrix)")
    print(f"- {RAW_DATA_DIR}/sample_metadata.csv (metadata about samples)")
    
    # This is a placeholder for actual download code
    # For demonstration, we'll create a simulated dataset if real data isn't present
    if not os.path.exists(RAW_DATA_DIR / "counts_matrix.csv"):
        create_example_dataset()

def create_example_dataset():
    """Create example dataset for demonstration purposes"""
    print("Creating simulated dataset for demonstration...")
    
    # Create sample metadata
    samples = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6']
    conditions = ['control', 'control', 'control', 'treatment', 'treatment', 'treatment']
    metadata = pd.DataFrame({
        'sample': samples,
        'condition': conditions
    })
    
    # Create count matrix with some glycogenes
    # Using some known glycogenes
    genes = ['MGAT1', 'MGAT2', 'ST6GAL1', 'B4GALT1', 'MAN1A1', 'MAN2A1', 'FUT8', 
             'GALNT1', 'GALNT2', 'ALG1', 'ALG2', 'GANAB', 'NEU1', 'NEU2',
             'GENE' + str(i) for i in range(1, 986)]  # Add more genes to make 1000 total
    
    # Generate random counts with some differential expression
    np.random.seed(42)
    counts = np.random.negative_binomial(20, 0.3, size=(len(genes), len(samples)))
    
    # Make some glycogenes differentially expressed
    glycogene_indices = list(range(14))  # First 14 genes are glycogenes in our list
    for idx in glycogene_indices:
        if idx % 2 == 0:  # Even indices will be upregulated
            counts[idx, 3:] = np.random.negative_binomial(40, 0.3, size=3)
        else:  # Odd indices will be downregulated
            counts[idx, 3:] = np.random.negative_binomial(10, 0.3, size=3)
    
    # Create dataframe
    count_df = pd.DataFrame(counts, index=genes, columns=samples)
    
    # Save files
    os.makedirs(RAW_DATA_DIR, exist_ok=True)
    count_df.to_csv(RAW_DATA_DIR / "counts_matrix.csv")
    metadata.to_csv(RAW_DATA_DIR / "sample_metadata.csv", index=False)
    
    # Also save a list of glycogenes for filtering later
    glycogenes = genes[:14]  # First 14 genes in our example
    pd.DataFrame(glycogenes, columns=['gene_symbol']).to_csv(
        DATA_DIR / "isha_data" / "glycogenes_list.csv", index=False
    )
    
    print(f"Example dataset created and saved to {RAW_DATA_DIR}")

def load_glycogenes():
    """Load list of glycogenes for filtering results"""
    try:
        # Try to load from the all_glycogenes.txt file if it exists
        glycogenes_file = DATA_DIR / "isha_data" / "all_glycogenes.txt"
        if os.path.exists(glycogenes_file):
            with open(glycogenes_file, 'r') as f:
                glycogenes = [line.strip() for line in f.readlines()]
            return set(glycogenes)
        
        # Alternative: try to load from the CSV file we might have created
        glycogenes_file_csv = DATA_DIR / "isha_data" / "glycogenes_list.csv"
        if os.path.exists(glycogenes_file_csv):
            glycogenes_df = pd.read_csv(glycogenes_file_csv)
            return set(glycogenes_df['gene_symbol'].tolist())
        
        print("No glycogenes list found. Will analyze all genes.")
        return None
    except Exception as e:
        print(f"Error loading glycogenes: {e}")
        return None

def run_deseq2_analysis():
    """Run DESeq2 analysis on the dataset using rpy2"""
    print("Starting DESeq2 analysis...")
    
    # Check if the required files exist
    counts_file = RAW_DATA_DIR / "counts_matrix.csv"
    metadata_file = RAW_DATA_DIR / "sample_metadata.csv"
    
    if not os.path.exists(counts_file) or not os.path.exists(metadata_file):
        print("Required data files not found. Please download the dataset first.")
        return
    
    # Load data
    counts = pd.read_csv(counts_file, index_col=0)
    metadata = pd.read_csv(metadata_file)
    
    # Ensure counts are integers (required by DESeq2)
    counts = counts.astype(int)
    
    # Convert data to R format
    pandas2ri.activate()
    
    # Import R packages
    try:
        deseq = importr('DESeq2')
        base = importr('base')
        stats = importr('stats')
        grdevices = importr('grDevices')
        ggplot2 = importr('ggplot2')
        pheatmap = importr('pheatmap')
    except Exception as e:
        print(f"Error importing R packages: {e}")
        print("Please ensure R and the required packages (DESeq2, ggplot2, pheatmap) are installed.")
        return
    
    # Convert to R objects
    r_counts = pandas2ri.py2rpy(counts)
    r_metadata = pandas2ri.py2rpy(metadata)
    
    # Create DESeq2 dataset
    try:
        # Create formula
        formula = robjects.r('formula(~ condition)')
        
        # Create DESeq dataset
        dds = deseq.DESeqDataSetFromMatrix(
            countData=r_counts,
            colData=r_metadata,
            design=formula
        )
        
        # Run DESeq2 analysis
        dds = deseq.DESeq(dds)
        
        # Get results (treatment vs control)
        res = deseq.results(dds, contrast=robjects.StrVector(["condition", "treatment", "control"]))
        
        # Convert results to pandas DataFrame
        res_df = pandas2ri.rpy2py(base.as_data_frame(res))
        res_df.index = counts.index
        
        # Add gene names as a column
        res_df['gene'] = res_df.index
        
        # Save results
        res_df.to_csv(PROCESSED_DATA_DIR / "deseq2_results.csv")
        
        # Create MA plot in R
        grdevices.pdf(str(FIGURES_DIR / "ma_plot.pdf"))
        deseq.plotMA(res)
        grdevices.dev_off()
        
        # Generate PCA plot in R
        vsd = deseq.vst(dds, blind=robjects.BoolVector([False]))
        grdevices.pdf(str(FIGURES_DIR / "pca_plot.pdf"))
        deseq.plotPCA(vsd, intgroup=robjects.StrVector(["condition"]))
        grdevices.dev_off()
        
        # Generate heatmap of top differentially expressed genes
        rld = deseq.rlog(dds, blind=robjects.BoolVector([False]))
        
        # Get top 50 genes by adjusted p-value
        res_ordered = deseq.order(res_df.padj)
        top_genes = robjects.IntVector([i+1 for i in range(50)])  # Top 50 genes
        
        # Create matrix of normalized values for top genes
        mat = stats.model_matrix(robjects.Formula('~ condition'), r_metadata)
        mat = mat.rx(robjects.IntVector([i+1 for i in range(ncol)]), robjects.IntVector([2]))
        
        # Create heatmap
        grdevices.pdf(str(FIGURES_DIR / "heatmap_top50.pdf"))
        pheatmap.pheatmap(
            base.assay(rld).rx(top_genes, True),
            annotation_col=r_metadata.rx(True, robjects.StrVector(['condition'])),
            scale="row"
        )
        grdevices.dev_off()
        
        print("DESeq2 analysis completed successfully!")
        print(f"Results saved to {PROCESSED_DATA_DIR}/deseq2_results.csv")
        print(f"Figures saved to {FIGURES_DIR}")
        
        return res_df
        
    except Exception as e:
        print(f"Error during DESeq2 analysis: {e}")
        return None

def filter_results_for_glycogenes(results_df):
    """Filter DESeq2 results to focus on glycogenes"""
    if results_df is None:
        print("No results to filter.")
        return
    
    # Load glycogenes
    glycogenes = load_glycogenes()
    
    if glycogenes is None or len(glycogenes) == 0:
        print("No glycogenes list available. Showing all significant genes.")
        # Filter for significant genes (adjusted p-value < 0.05)
        sig_genes = results_df[results_df['padj'] < 0.05].sort_values('padj')
        print(f"Found {len(sig_genes)} significant genes (padj < 0.05)")
        return sig_genes
    
    # Filter for glycogenes
    glycogene_results = results_df[results_df.index.isin(glycogenes)]
    print(f"Found {len(glycogene_results)} glycogenes in the dataset")
    
    # Filter for significant glycogenes
    sig_glycogenes = glycogene_results[glycogene_results['padj'] < 0.05].sort_values('padj')
    print(f"Found {len(sig_glycogenes)} significantly differentially expressed glycogenes (padj < 0.05)")
    
    # Save filtered results
    sig_glycogenes.to_csv(PROCESSED_DATA_DIR / "significant_glycogenes.csv")
    print(f"Significant glycogenes saved to {PROCESSED_DATA_DIR}/significant_glycogenes.csv")
    
    return sig_glycogenes

def create_visualizations(results_df, sig_glycogenes=None):
    """Create additional visualizations for the results"""
    if results_df is None:
        print("No results to visualize.")
        return
    
    # Volcano plot for all genes
    plt.figure(figsize=(10, 8))
    plt.scatter(
        results_df['log2FoldChange'], 
        -np.log10(results_df['pvalue']),
        alpha=0.5, 
        s=5, 
        color='gray'
    )
    
    # Highlight significant glycogenes if available
    if sig_glycogenes is not None and len(sig_glycogenes) > 0:
        plt.scatter(
            sig_glycogenes['log2FoldChange'],
            -np.log10(sig_glycogenes['pvalue']),
            alpha=1,
            s=30,
            color='red',
            label='Significant Glycogenes'
        )
        
        # Label top 10 significant glycogenes
        top10 = sig_glycogenes.sort_values('padj').head(10)
        for idx, row in top10.iterrows():
            plt.annotate(
                idx,  # Gene name (index)
                (row['log2FoldChange'], -np.log10(row['pvalue'])),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8
            )
    
    # Add cutoff lines
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.3)
    plt.axvline(-1, color='blue', linestyle='--', alpha=0.3)
    plt.axvline(1, color='blue', linestyle='--', alpha=0.3)
    
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 p-value')
    plt.title('Volcano Plot of Differential Expression')
    
    if sig_glycogenes is not None and len(sig_glycogenes) > 0:
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "volcano_plot.png", dpi=300)
    print(f"Volcano plot saved to {FIGURES_DIR}/volcano_plot.png")
    
    # Barplot of top 20 significant genes (or glycogenes if available)
    if sig_glycogenes is not None and len(sig_glycogenes) > 0:
        data_to_plot = sig_glycogenes
        title_suffix = "Glycogenes"
    else:
        data_to_plot = results_df[results_df['padj'] < 0.05].sort_values('padj')
        title_suffix = "Genes"
    
    # Get top 20 (or fewer if less are available)
    top_n = min(20, len(data_to_plot))
    if top_n == 0:
        print("No significant genes to plot.")
        return
    
    top_genes = data_to_plot.sort_values('padj').head(top_n)
    
    plt.figure(figsize=(12, 8))
    bars = plt.barh(
        y=top_genes.index,
        width=top_genes['log2FoldChange'],
        color=[('red' if x < 0 else 'blue') for x in top_genes['log2FoldChange']]
    )
    
    # Add gene names and p-values
    for i, bar in enumerate(bars):
        gene = top_genes.index[i]
        padj = top_genes['padj'].iloc[i]
        plt.text(
            0,
            i,
            f"{gene} (p={padj:.2e})",
            ha='center',
            va='center',
            color='white',
            fontweight='bold',
            fontsize=8
        )
    
    plt.xlabel('log2 Fold Change')
    plt.title(f'Top {top_n} Differentially Expressed {title_suffix}')
    plt.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f"top_{top_n}_{title_suffix.lower()}.png", dpi=300)
    print(f"Barplot saved to {FIGURES_DIR}/top_{top_n}_{title_suffix.lower()}.png")

def main():
    """Main function to run the complete workflow"""
    print("Starting RNAseq analysis workflow...")
    
    # Download or simulate dataset
    download_dataset_from_glygen()
    
    # Run DESeq2 analysis
    results = run_deseq2_analysis()
    
    # Filter for glycogenes
    sig_glycogenes = filter_results_for_glycogenes(results)
    
    # Create visualizations
    create_visualizations(results, sig_glycogenes)
    
    print("RNAseq analysis workflow completed!")

if __name__ == "__main__":
    main() 