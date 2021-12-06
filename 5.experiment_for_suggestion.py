import tqdm
import gene_to_sym
gene_to_sym = gene_to_sym.get()
ches = ""
chep = ""
t = 0
tmp2 = ""
fw = open("data/experiments_for_xgboost.csv", "w")
f = open("data/masked_sentences_bert.txt", "r")
for line in tqdm.tqdm(f):
    ls = line.strip().split("\t")
    if chep != ls[2]:
        sentid = 0
        chep = ls[2]
    if ches != ls[-2]:
        ches = ls[-2]
        sentid += 1
    tmp = ",".join([ls[1], ls[2], str(sentid),
                   gene_to_sym[ls[3].upper()], ls[4]])
    if tmp != tmp2:
        fw.write(tmp + "\n")
        tmp2 = tmp
f.close()
fw.close()
