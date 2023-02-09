
import xml.etree.ElementTree as ET
import scispacy
import en_core_sci_sm
import tqdm
from ahocorapy.keywordtree import KeywordTree

import lexas.gene_to_sym
gene_to_sym = lexas.gene_to_sym.get()

def join_and_remove_nt(sec):
    sent = " ".join([r for r in sec.itertext()]).strip()
    sent = sent.replace("\t"," ").replace("\n"," ").replace("\u2009","")
    sent = sent.replace("\xa0","")
    return sent

def parse(file):
    year,sent = 0,[]    
    try:
        tree = ET.parse(file)
        root = tree.getroot()
        year = root.find(".//article-meta").find("pub-date").find("year").text
        for sec in root.findall(".//sec"):
            try:
                sec_title = sec.find('title').text.lower()
                if "result" in sec_title:
                    for child in sec.findall(".//p"):
                        for fig in child.findall('fig'):
                            child.remove(fig)
                        sent.append(join_and_remove_nt(child))
            except:
                continue
    except:
        pass
    return year, " ".join(sent)


def segmentation(sent):
    en = en_core_sci_sm.load()
    doc = en(sent)
    ret = "#####".join([str(a) for a in doc.sents if len(str(a).split(" ")) < 500])
    return ret
 
def initialize_dictionaries():
    dic_hgnc = KeywordTree(case_insensitive=True)
    dic_expe = KeywordTree(case_insensitive=True)
    
    #Gene
    genes = {key.lower() for key in gene_to_sym if len(key) > 3}
    with open('./data/terms_to_be_removed.csv') as remo:
        to_remove = [word.strip().lower() for word in remo]
        genes = set([gene for gene in genes if gene not in to_remove])

    for gene in genes:
        dic_hgnc.add(" " + gene + " ")
    dic_hgnc.finalize()

    #Experiments
    expes = set()
    with open('./data/experiments_list.txt') as f:
        expes = [word.strip().lower() for word in f]
    
    for e in expes:
        dic_expe.add(" " + e + " ")
    dic_expe.finalize()
    
    return dic_hgnc,dic_expe

def mask(text, dicg, dice):
    text2 = " " + text + " "
    text2 = text2.replace(",","").replace(".","").replace("("," ").replace(")"," ")
    ges = list(dicg.search_all(text2))
    tmp = [(text2[:g[1]] + " [GENE] " + text2[g[1] + len(g[0]):], g[0].strip()) for g in ges]
    ret = []
    for t, g in tmp:
        exs = list(dice.search_all(t))
        ret += [[t[:e[1]] + " [EXPE] " +  t[e[1] + len(e[0]):], g, e[0].strip()] for e in exs]
    return ret
    