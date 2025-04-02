import requests
import csv
import sys
import os

def control(gene):
    return {"value": gene, "name": "gene"}

def get_genes(filename, symbol_column='gene_symbol'):
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return None
    
    # Read the CSV file
    with open(filename, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        headers = reader.fieldnames
        if symbol_column not in headers:
            print(f"Column {symbol_column} does not exist in the file.")
            return None
        
        # Extract the gene symbols
        genes = set()
        for row in reader:
            if row[symbol_column] and row['padj'] and row['padj'] != "NA":
                genes.add(row[symbol_column].strip())

    controls = {
        "controls": [control(gene) for gene in genes]
    }
    url = "https://translator.broadinstitute.org/hgnc/genes/transform"
    response = requests.post(url, json=controls)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Error: {response.status_code} - {response.text}")


def run_gelinea(collection):
    url = 'https://translator.broadinstitute.org/gelinea/enrichment/transform'
    query = {
        'controls': [
            {
                "name": "maximum p-value",
                "value": "1e-03"
            },
            {
                "name": "network",
                "value": "STRING-human-700"
            },
            {
                "name": "gene-set collection",
                "value": "H - hallmark gene sets"
            },
            {
                "name": "gene-set collection",
                "value": "C2 - curated gene sets"
            }
        ],
        'collection': collection,
    }
    response = requests.post(url, json=query)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Error: {response.status_code} - {response.text}")


def get_gene(gene):
    return {
        "id": gene["id"],
        "identifiers": gene["identifiers"]
    }


def main():
    genes=get_genes(sys.argv[1])
    if genes:
        collection = [get_gene(gene) for gene in genes]
        print("analyzing", len(collection), "genes")
        gelinea = run_gelinea(collection)
        for pathway in gelinea:
            pvalue = 1.0
            for attribute in pathway["connections"][0]["attributes"]:
                if attribute["attribute_type_id"] == "biolink:p_value":
                    pvalue = attribute["value"]
            print(pathway["id"], pvalue, sep="\t")


if __name__ == "__main__":
    main()
