import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import os
from collections import Counter

def load_data():
    """Load the filtered disease data."""
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    data_file = os.path.join(base_dir, 'data/final_data/filtered_diseases.csv')
    return pd.read_csv(data_file)

def create_network_graph(data, output_dir):
    """Create a focused network graph of gene-disease relationships."""
    # Create a bipartite graph
    G = nx.Graph()
    
    # Add edges between genes and diseases
    for _, row in data.iterrows():
        G.add_edge(row['gene_symbol'], row['mondo_disease_name'])
    
    # Get top genes and diseases by degree
    gene_degrees = {node: G.degree(node) for node in G.nodes() if node in data['gene_symbol'].unique()}
    disease_degrees = {node: G.degree(node) for node in G.nodes() if node in data['mondo_disease_name'].unique()}
    
    # Select top 20 genes and their associated diseases
    top_genes = sorted(gene_degrees.items(), key=lambda x: x[1], reverse=True)[:20]
    top_genes = [gene for gene, _ in top_genes]
    
    # Create a subgraph with top genes and their neighbors
    subgraph_nodes = set(top_genes)
    for gene in top_genes:
        subgraph_nodes.update(G.neighbors(gene))
    
    H = G.subgraph(subgraph_nodes)
    
    # Create the plot
    plt.figure(figsize=(20, 20))
    
    # Use bipartite layout
    pos = nx.bipartite_layout(H, top_genes)
    
    # Draw nodes with different colors for genes and diseases
    gene_nodes = [node for node in H.nodes() if node in top_genes]
    disease_nodes = [node for node in H.nodes() if node not in top_genes]
    
    nx.draw_networkx_nodes(H, pos, nodelist=gene_nodes, node_color='lightblue', node_size=2000)
    nx.draw_networkx_nodes(H, pos, nodelist=disease_nodes, node_color='lightgreen', node_size=1500)
    
    # Draw edges
    nx.draw_networkx_edges(H, pos, edge_color='gray', alpha=0.5)
    
    # Draw labels
    nx.draw_networkx_labels(H, pos, font_size=8, font_weight='bold')
    
    plt.title('Top 20 Genes and Their Associated Diseases', pad=20, fontsize=16)
    plt.axis('off')
    plt.savefig(os.path.join(output_dir, 'network_graph.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_heatmap(data, output_dir):
    """Create a focused heatmap of disease-gene associations."""
    # Get top 20 genes and diseases by frequency
    top_genes = data['gene_symbol'].value_counts().head(20).index
    top_diseases = data['mondo_disease_name'].value_counts().head(20).index
    
    # Filter data for top genes and diseases
    filtered_data = data[
        (data['gene_symbol'].isin(top_genes)) & 
        (data['mondo_disease_name'].isin(top_diseases))
    ]
    
    # Create pivot table
    pivot_data = pd.crosstab(
        filtered_data['mondo_disease_name'], 
        filtered_data['gene_symbol']
    )
    
    # Create the heatmap
    plt.figure(figsize=(15, 15))
    sns.heatmap(pivot_data, 
                cmap='YlOrRd',
                cbar_kws={'label': 'Association'},
                square=True)
    
    plt.title('Top 20 Genes vs Top 20 Diseases Association Heatmap', pad=20, fontsize=16)
    plt.xlabel('Genes', fontsize=12)
    plt.ylabel('Diseases', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_distribution_plots(data, output_dir):
    """Create focused distribution plots."""
    # Diseases per gene (top 20)
    diseases_per_gene = data.groupby('gene_symbol')['mondo_disease_name'].nunique()
    top_20_genes = diseases_per_gene.nlargest(20)
    
    plt.figure(figsize=(15, 8))
    sns.barplot(x=top_20_genes.index, y=top_20_genes.values, palette='viridis')
    plt.title('Top 20 Genes by Number of Associated Diseases', pad=20, fontsize=16)
    plt.xlabel('Genes', fontsize=12)
    plt.ylabel('Number of Diseases', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'diseases_per_gene.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Genes per disease (top 20)
    genes_per_disease = data.groupby('mondo_disease_name')['gene_symbol'].nunique()
    top_20_diseases = genes_per_disease.nlargest(20)
    
    plt.figure(figsize=(15, 8))
    sns.barplot(x=top_20_diseases.index, y=top_20_diseases.values, palette='viridis')
    plt.title('Top 20 Diseases by Number of Associated Genes', pad=20, fontsize=16)
    plt.xlabel('Diseases', fontsize=12)
    plt.ylabel('Number of Genes', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'genes_per_disease.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_statistics(data, output_dir):
    """Create a clear summary statistics visualization."""
    # Calculate statistics
    stats = {
        'Total Genes': len(data['gene_symbol'].unique()),
        'Total Diseases': len(data['mondo_disease_name'].unique()),
        'Total Associations': len(data),
        'Avg Diseases per Gene': round(data.groupby('gene_symbol')['mondo_disease_name'].nunique().mean(), 2),
        'Avg Genes per Disease': round(data.groupby('mondo_disease_name')['gene_symbol'].nunique().mean(), 2)
    }
    
    # Create a more visually appealing bar plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(stats.keys(), stats.values(), color='skyblue')
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{height}',
                ha='center', va='bottom')
    
    plt.title('Summary Statistics of Gene-Disease Associations', pad=20, fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Count', fontsize=12)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'summary_statistics.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Create output directory if it doesn't exist
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    output_dir = os.path.join(base_dir, 'results/figures')
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    data = load_data()
    
    # Create visualizations
    print("Creating network graph...")
    create_network_graph(data, output_dir)
    
    print("Creating heatmap...")
    create_heatmap(data, output_dir)
    
    print("Creating distribution plots...")
    create_distribution_plots(data, output_dir)
    
    print("Creating summary statistics...")
    create_summary_statistics(data, output_dir)
    
    print(f"\nVisualizations have been saved to {output_dir}")

if __name__ == "__main__":
    main() 