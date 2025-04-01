import pandas as pd
import os

def read_gene_list(file_path):
    """Read genes from the text file."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def filter_diseases_by_genes(gene_list, disease_data, gene_col, disease_col):
    """Filter diseases based on the presence of genes in the gene list."""
    # Convert gene list to set for faster lookup
    gene_set = set(gene_list)
    
    # Filter the disease data to only include rows where genes are in our gene set
    filtered_data = disease_data[disease_data[gene_col].isin(gene_set)]
    
    # Select and rename columns to match our standard format
    filtered_data = filtered_data[[gene_col, disease_col]].copy()
    filtered_data.columns = ['gene_symbol', 'mondo_disease_name']
    
    return filtered_data

def process_disease_file(file_path, gene_set):
    """Process a single disease file and return filtered data."""
    disease_data = pd.read_csv(file_path)
    
    # Define column mappings for each file
    if 'human_protein_disease_glyco.csv' in file_path:
        gene_col = 'gene_symbol'
        disease_col = 'mondo_disease_name'
    elif 'human_protein_disease_cdg.csv' in file_path:
        gene_col = 'gene_symbol_human'
        disease_col = 'disease_name'
    elif 'human_protein_disease_genomics_england.csv' in file_path:
        gene_col = 'xref_id'
        disease_col = 'do_id'  # Using DO ID as disease identifier
    
    return filter_diseases_by_genes(gene_set, disease_data, gene_col, disease_col)

def main():
    # Define paths
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    gene_file = os.path.join(base_dir, 'data/isha_data/all_glycogenes.txt')
    output_file = os.path.join(base_dir, 'data/final_data/filtered_diseases.csv')
    
    # List of disease files to process
    disease_files = [
        'human_protein_disease_glyco.csv',
        'human_protein_disease_cdg.csv',
        'human_protein_disease_genomics_england.csv'
    ]
    
    # Read the gene list
    gene_list = read_gene_list(gene_file)
    gene_set = set(gene_list)
    
    # Process each disease file and combine results
    all_filtered_data = []
    for file in disease_files:
        file_path = os.path.join(base_dir, 'data/raw', file)
        print(f"\nProcessing {file}...")
        filtered_data = process_disease_file(file_path, gene_set)
        all_filtered_data.append(filtered_data)
        print(f"Found {len(filtered_data['mondo_disease_name'].unique())} unique diseases and {len(filtered_data['gene_symbol'].unique())} unique genes")
    
    # Combine all filtered data
    combined_data = pd.concat(all_filtered_data, ignore_index=True)
    
    # Remove duplicates based on gene_symbol and mondo_disease_name
    unique_data = combined_data.drop_duplicates(subset=['gene_symbol', 'mondo_disease_name'])
    
    # Save the filtered data
    unique_data.to_csv(output_file, index=False)
    print(f"\nFinal results:")
    print(f"Filtered data saved to {output_file}")
    print(f"Total number of unique diseases found: {len(unique_data['mondo_disease_name'].unique())}")
    print(f"Total number of unique genes matched: {len(unique_data['gene_symbol'].unique())}")

if __name__ == "__main__":
    main() 