import pandas as pd
sym, asy, nam, ana, pna, psy, chr = [], [], [], [], [], [], []
# Downloaded from HGNC Bio-Mart at 2021-09-08
f = open("data/HGNC_0908.txt", "r")
for line in f:
    if line.startswith("Approved"):
        continue
    ls = line.strip("\n").split("\t")
    sym.append(ls[0])  # symbol
    nam.append(ls[1])  # name
    psy.append(ls[2])  # previous symbol
    asy.append(ls[3])  # alias symbol
    # chr.append(ls[4]) #chromosome
    ana.append(ls[5])  # alias name
    pna.append(ls[6])  # previous name
f.close()

gene_to_sym = {}
for n, s in enumerate(sym):
    gene_to_sym[s.upper()] = s.upper()
symbol_set = set([s.upper() for s in sym])
overlap = set()


def add_dic(col, gene_to_sym, overlap, type):
    for n in range(len(col)):
        a = str(col[n]).upper()
        if type == "symbol":
            terms = a.split(",")
            terms = [term.strip() for term in terms]
        if type == "name":
            terms = a.split("\",")
            terms = [term.replace("\"", "").strip() for term in terms]
        for term in terms:
            if term in symbol_set:
                continue
            if term in gene_to_sym:
                if gene_to_sym[term] != sym[n]:
                    overlap.add(term)
                    continue
            gene_to_sym[term] = sym[n]
    return gene_to_sym, overlap


gene_to_sym, overlap = add_dic(asy, gene_to_sym, overlap, "symbol")
gene_to_sym, overlap = add_dic(nam, gene_to_sym, overlap, "name")
gene_to_sym, overlap = add_dic(ana, gene_to_sym, overlap, "name")
gene_to_sym, overlap = add_dic(pna, gene_to_sym, overlap, "name")
gene_to_sym, overlap = add_dic(psy, gene_to_sym, overlap, "symbol")

for over in overlap:
    del gene_to_sym[over]

gene_to_sym2 = gene_to_sym.copy()
for g in gene_to_sym:
    if "-" in g and g.replace("-", "") not in gene_to_sym:
        gene_to_sym2[g.replace("-", "")] = gene_to_sym[g]

gene_to_sym = gene_to_sym2


def get():
    return gene_to_sym
