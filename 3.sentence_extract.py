import gene_to_sym
import pandas as pd
from ahocorapy.keywordtree import KeywordTree
dic_hgnc = KeywordTree(case_insensitive=True)
dic_expe = KeywordTree(case_insensitive=True)

gene_to_sym = gene_to_sym.get()

genes = set()
for key in gene_to_sym.keys():
    if len(key) < 3:
        continue
    else:
        genes.add(key.lower())

with open('data/terms_to_be_removed.csv') as remo:
    for word in remo:
        try:
            genes.remove(word.strip().lower())
        except KeyError:
            pass

for gene in genes:
    dic_hgnc.add(" " + gene + " ")
dic_hgnc.finalize()


expes = set()
with open('data/experiments_list.txt') as f:
    for line in f:
        line_splitted = line.strip().split('\t')
        expe = line_splitted[0]
        expes.add(expe.lower())
for e in expes:
    dic_expe.add(" " + e + " ")
dic_expe.finalize()


def replace(text, dicg, dice):
    text2 = " " + text + " "
    text2 = text2.replace(
        ",",
        "").replace(
        ".",
        "").replace(
            "(",
            " ").replace(
                ")",
        " ")
    ges = list(dicg.search_all(text2))
    if len(ges) == 0:
        return []
    tmp = []
    ret = []
    for g in ges:
        tmp.append((text2[:g[1]] + " [GENE] " +
                   text2[g[1] + len(g[0]):], g[0].strip()))
    for t, g in tmp:
        exs = list(dice.search_all(t))
        if len(exs) == 0:
            return []
        for e in exs:
            ret.append((t[:e[1]] + " [EXPE] " +
                       t[e[1] + len(e[0]):], g, e[0].strip()))
    return ret


def main(read, write):
    f2 = open(write, "w")
    with open(read, "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            if len(ls) == 3:
                ls2 = ls[2].strip().split("#####")  # textの集合体
                for sent in ls2:  # lllは文章
                    tmp = replace(sent, dic_hgnc, dic_expe)
                    if len(tmp) == 0:
                        continue
                    for t in tmp:
                        f2.write(
                            ls[0] +
                            "\t" +
                            ls[1] +
                            "\t" +
                            t[1] +
                            "\t" +
                            t[2] +
                            "\t" +
                            sent +
                            "\t" +
                            t[0] +
                            "\n")

    f.close()
    f2.close()


read = "data/result_sections_segmentation.txt"
write = "data/masked_sentences.txt"
main(read, write)
