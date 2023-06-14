import pandas as pd

# Load data from HGNC file
with open("./data/HGNC.txt", "r") as f:
    data = [line.strip("\n").split("\t") for line in f if not line.startswith("Approved")]

# Split the data into individual lists
symbol, name, prev_symbol, alias_symbol, alias_name, prev_name = [[row[i] for row in data] for i in range(6)]

# Initialize the gene to symbol dictionary
gene_to_symbol = {s.upper(): s.upper() for s in symbol}

# Create a set to keep track of overlapping gene symbols
overlap = set()

# Function to parse and update gene to symbol dictionary
def update_dictionary(column, gene_to_symbol, overlap, entry_type):
    # Process each entry in the column
    for i, entry in enumerate(map(str.upper, column)):
        # Split the entry into terms based on entry_type
        terms = entry.split(",") if entry_type == "symbol" else entry.split("\",")
        terms = [term.strip().replace("\"", "") for term in terms]
        
        # If a term is already in the gene_to_symbol dictionary, add it to the overlap set. Otherwise, add it to gene_to_symbol.
        for term in terms:
            if term in gene_to_symbol:
                overlap.add(term)
                continue
            gene_to_symbol[term] = symbol[i].upper()
    return gene_to_symbol, overlap

# Lists of gene entries and corresponding types
gene_lists = [alias_symbol, name, alias_name, prev_name, prev_symbol]
entry_types = ["symbol", "name", "name", "name", "symbol"]

# Update gene_to_symbol dictionary
for gene_column, entry_type in zip(gene_lists, entry_types):
    gene_to_symbol, overlap = update_dictionary(gene_column, gene_to_symbol, overlap, entry_type)

# Remove overlapping entries from gene_to_symbol
gene_to_symbol = {k: v for k, v in gene_to_symbol.items() if k not in overlap}

# Add entries from 'symbol' list to gene_to_symbol
gene_to_symbol.update({s.upper(): s.upper() for s in symbol})

# If a gene symbol contains '-', add another entry without '-'
gene_to_symbol.update({g.replace("-", ""): v for g, v in gene_to_symbol.items() if "-" in g and g.replace("-", "") not in gene_to_symbol})

# Function to return gene_to_symbol dictionary
def get():
    return gene_to_symbol
