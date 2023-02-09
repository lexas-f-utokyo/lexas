#Gene_to_sym dictionary
import pandas as pd
with open("./data/HGNC.txt", "r") as f:
    data = [line.strip("\n").split("\t") for line in f if not line.startswith("Approved")]
    
sym, nam, psy, asy, ana, pna = [[row[i] for row in data] for i in range(6)]
gene_to_sym = {s.upper(): s.upper() for s in sym}
overlap=set()

def add_dic(col, gene_to_sym, overlap, typ):
    for i, a in enumerate(map(str.upper, col)):
        terms = a.split(",") if typ == "symbol" else a.split("\",")
        terms = [term.strip().replace("\"", "") for term in terms]
        for term in terms:
            if term in gene_to_sym:
                overlap.add(term)
                continue
            gene_to_sym[term] = sym[i].upper()
    return gene_to_sym, overlap

lists = [asy, nam, ana, pna, psy]
types = ["symbol", "name", "name", "name", "symbol"]

for col, typ in zip(lists, types):
    gene_to_sym, overlap = add_dic(col, gene_to_sym, overlap, typ)

gene_to_sym = {k: v for k, v in gene_to_sym.items() if k not in overlap}
gene_to_sym.update({s.upper(): s.upper() for s in sym})
gene_to_sym.update({g.replace("-", ""): v for g, v in gene_to_sym.items() if "-" in g and g.replace("-", "") not in gene_to_sym})


def get():
    return gene_to_sym