
import xml.etree.ElementTree as ET
import scispacy
import en_core_sci_sm
import tqdm
import re
import os
import concurrent.futures
from ahocorapy.keywordtree import KeywordTree

import lexas.gene_to_sym
gene_to_sym = lexas.gene_to_sym.get()

pattern = re.compile(r'\t|\n|\x0b|\x0c|\r|\x1c|\x1d|\x1e|\x1f|\x85|\xa0|\u1680|\u2000|\u2001|\u2002|\u2003|\u2004|\u2005|\u2006|\u2007|\u2008|\u2009|\u200a|\u2028|\u2029|\u202f|\u205f|\u3000')
def join_and_remove_nt(sec):
    sent = " ".join([r for r in sec.itertext()]).strip()
    sent = pattern.sub(' ', sent)
    return sent

def parse(file):
    year,sent,para = 0,[],[]   
    try:
        tree = ET.parse(file)
        root = tree.getroot()
        year = root.find(".//article-meta").find("pub-date").find("year").text
        for sec in root.findall(".//sec"):
            try:
                sec_title = sec.find('title').text.lower()
                if "result" in sec_title:
                    for n,child in enumerate(sec.findall(".//p")):
                        sent.append(join_and_remove_nt(child))
                        para.append(n)
            except:
                continue
    except:
        pass
    return year,sent,para


def result_extraction(article_dir="./articles/",output="./data/result_sections.txt"):
    with open(output, "w") as f:            
        for file in tqdm.tqdm(os.listdir(article_dir)):
            pmcid = file.split(".")[0]
            year,sent,para = lexas.sentence.parse(os.path.join(article_dir,file))
            if year != 0:
                for n in range(len(sent)):
                    sentences = sent[n]
                    paragraph = str(para[n])
                    segmented_sentences = lexas.sentence.segmentation(sentences)
                    f.write("\t".join([year, pmcid+"-"+paragraph, segmented_sentences]) + "\n")
                
en = en_core_sci_sm.load()        
def segmentation(sent):
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

dics = {}
def mask_gene_experiment(input="./data/result_sections.txt",output="./data/masked_sentences.txt"):
    if dics=={}:
        print("Initializing dictionaries...")
        dic_hgnc,dic_expe = lexas.sentence.initialize_dictionaries()
        dics["hgnc"] = dic_hgnc
        dics["expe"] = dic_expe
        print("Done")
    dic_hgnc,dic_expe = dics["hgnc"],dics["expe"]
    with open(output, "w") as f:
        with open(input, "r") as f2:    
            for line in tqdm.tqdm(f2):
                year,pmcid,sentences = line.strip("\n").split("\t")
                for sentence in sentences.split("#####"):
                    masked = lexas.sentence.mask(sentence,dic_hgnc,dic_expe)
                    for m in masked:
                        f.write("\t".join([year,pmcid] + m[1:3]+[sentence] + m[0:1])+"\n")
    