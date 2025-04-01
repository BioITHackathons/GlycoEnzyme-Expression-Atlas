import pandas as pd
import os

def get_all_gene_names(df):
    """
    Get all possible gene names from a dataframe, including alternatives.
    Returns a dictionary mapping each gene name to its primary name.
    """
    gene_names = {}
    
    if 'gene_symbol' in df.columns:
        primary_col = 'gene_symbol'
    elif 'gene_name' in df.columns:
        primary_col = 'gene_name'
    else:
        print("No gene identifier column found")
        return gene_names

    # Add primary gene names
    for primary_name in df[primary_col].dropna().unique():
        gene_names[primary_name] = primary_name

    # Add alternative names if available
    if 'gene_names_alternative' in df.columns:
        for idx, row in df.iterrows():
            primary_name = row[primary_col]
            if pd.notna(row['gene_names_alternative']) and row['gene_names_alternative']:
                alt_names = row['gene_names_alternative'].split(';')
                for alt_name in alt_names:
                    if alt_name:  # Skip empty strings
                        gene_names[alt_name] = primary_name

    return gene_names

def read_csv_file(file_path):
    """
    Read a CSV file and return a set of unique gene symbols and a mapping of all names.
    """
    try:
        df = pd.read_csv(file_path)
        gene_names = get_all_gene_names(df)
        if not gene_names:
            print(f"Warning: No gene identifiers found in {file_path}")
            print("Available columns:", df.columns.tolist())
            return set(), {}
        return set(gene_names.values()), gene_names  # Return unique primary names and name mapping
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return set(), {}

def find_common_genes(set1, set2, name_map1, name_map2, output_file, set1_name="Set 1", set2_name="Set 2"):
    """
    Find common genes between two sets and save them to a CSV file.
    Also check for matches using alternative names.
    """
    # Find direct matches
    common_genes = set1.intersection(set2)
    
    # Find matches through alternative names
    all_names1 = set(name_map1.keys())
    all_names2 = set(name_map2.keys())
    
    # Find genes that match through alternative names
    alt_matches = set()
    for name in all_names1:
        if name in all_names2:
            primary1 = name_map1[name]
            primary2 = name_map2[name]
            if primary1 != primary2:  # Only add if they weren't already direct matches
                alt_matches.add((primary1, primary2))

    # Save results
    if common_genes or alt_matches:
        # Save direct matches
        results = []
        for gene in sorted(common_genes):
            results.append({'primary_name': gene, 'match_type': 'direct'})
        
        # Save alternative name matches
        for primary1, primary2 in sorted(alt_matches):
            results.append({
                'primary_name': primary1,
                'alternative_match': primary2,
                'match_type': 'alternative'
            })
        
        output_df = pd.DataFrame(results)
        output_df.to_csv(output_file, index=False)
        
        print(f"\nResults have been saved to {output_file}")
        print(f"Number of direct matches: {len(common_genes)}")
        if common_genes:
            print("First few direct matches:")
            print(', '.join(sorted(list(common_genes))[:5]))
        
        print(f"\nNumber of alternative name matches: {len(alt_matches)}")
        if alt_matches:
            print("First few alternative matches:")
            for primary1, primary2 in sorted(list(alt_matches))[:5]:
                print(f"{primary1} ↔ {primary2}")
    else:
        print("\nNo matches found.")
    
    # Report unique genes
    unique_to_set1 = set1 - set2
    unique_to_set2 = set2 - set1
    
    print(f"\nGenes unique to {set1_name}: {len(unique_to_set1)}")
    if unique_to_set1:
        print("First few unique genes:", ', '.join(sorted(list(unique_to_set1))[:5]))
        
        # Save unique genes to a separate file
        unique_file = f"unique_to_{set1_name.lower()}.csv"
        unique_df = pd.DataFrame({'gene_symbol': sorted(list(unique_to_set1))})
        unique_df.to_csv(unique_file, index=False)
        print(f"Unique genes saved to {unique_file}")
    
    print(f"\nGenes unique to {set2_name}: {len(unique_to_set2)}")
    if unique_to_set2:
        print("First few unique genes:", ', '.join(sorted(list(unique_to_set2))[:5]))
    
    return common_genes, alt_matches

def main():
    # Define the files to process
    files = {
        'hydrolase': 'human_protein_glycohydrolase.csv',
        'transferase': 'human_protein_glycosyltransferase.csv',
        'glycogene': 'human_protein_glycogenes.csv'
    }

    # Dictionary to store gene sets and name mappings
    gene_sets = {}
    name_mappings = {}

    # Read gene symbols from each file
    for file_type, file_path in files.items():
        if os.path.exists(file_path):
            gene_sets[file_type], name_mappings[file_type] = read_csv_file(file_path)
            print(f"\nFile: {file_path}")
            print(f"Number of unique genes: {len(gene_sets[file_type])}")
            print(f"Total gene names (including alternatives): {len(name_mappings[file_type])}")
        else:
            print(f"Warning: File not found - {file_path}")

    # Compare hydrolase genes with glycogene list
    if all(key in gene_sets for key in ['hydrolase', 'glycogene']):
        print("\n=== Comparing Hydrolase with Glycogene List ===")
        hydrolase_common, hydrolase_alt = find_common_genes(
            gene_sets['hydrolase'],
            gene_sets['glycogene'],
            name_mappings['hydrolase'],
            name_mappings['glycogene'],
            'common_hydrolase_glycogene_genes.csv',
            'Hydrolase',
            'Glycogene'
        )

    # Compare transferase genes with glycogene list
    if all(key in gene_sets for key in ['transferase', 'glycogene']):
        print("\n=== Comparing Transferase with Glycogene List ===")
        transferase_common, transferase_alt = find_common_genes(
            gene_sets['transferase'],
            gene_sets['glycogene'],
            name_mappings['transferase'],
            name_mappings['glycogene'],
            'common_transferase_glycogene_genes.csv',
            'Transferase',
            'Glycogene'
        )

    # Find genes that are common across all three files
    if all(key in gene_sets for key in ['hydrolase', 'transferase', 'glycogene']):
        print("\n=== Genes Common Across All Three Files ===")
        all_common = set.intersection(gene_sets['hydrolase'], gene_sets['transferase'], gene_sets['glycogene'])
        if all_common:
            output_df = pd.DataFrame({'gene_symbol': sorted(list(all_common))})
            output_df.to_csv('common_all_genes.csv', index=False)
            print(f"\nGenes common across all files have been saved to common_all_genes.csv")
            print(f"Number of genes common across all files: {len(all_common)}")
            print("\nFirst few common genes:")
            print(', '.join(sorted(list(all_common))[:5]))
        else:
            print("No genes are common across all three files.")

if __name__ == "__main__":
    main() 