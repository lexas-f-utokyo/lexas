import numpy as np
import tqdm
import lexas.gene_to_sym
gene_to_sym = lexas.gene_to_sym.get()

with open("./data/symbols.txt","r") as f:
    symbols = [line.strip() for line in f]

categorical,numerical,string = {},{},{}

def feature_load():
  if categorical == {}:       
    #Loading dictionary
    print("Loading categorical features...")
    with open("./Repository/feature/categorical_feature.txt","r") as f:
        for line in f:
            ls = line.strip("\n").split("\t")
            if ls[0] not in categorical:
                categorical[ls[0]] = {}
            if ls[2] == "":
                categorical[ls[0]][ls[1]] = []
            else:
                categorical[ls[0]][ls[1]] = [s for s in ls[2].split(",")]
                
    print("Loading numerical features...")              
    with open("./Repository/feature/numerical_feature.txt","r") as f:
        for line in f:
            ls = line.strip("\n").split("\t")
            if ls[0] not in numerical:
                numerical[ls[0]] = {}
            if ls[2] == "":
                numerical[ls[0]][ls[1]] = []
            else:
                values = [float(s) for s in ls[2].split(",")]
                if np.std(values)==0:
                    numerical[ls[0]][ls[1]] = []
                else:
                    numerical[ls[0]][ls[1]] = values

    print("Loading STRING...")
    string = {}
    with open("./Repository/feature/string.txt","r") as f:
        for line in f:
            ls = line.strip("\n").split("\t")
            string[(ls[0],ls[1])] = int(ls[2])
  else:
    print("Features have been loaded!")

        
def numbering(feature_dic):
    all = set()
    for g in feature_dic:
        for f in feature_dic[g]:
            all.add(f)
    feature_list = sorted(list(all))
    ret = {feature:n for n,feature in enumerate(feature_list)}
    return ret

def choose_feature(cat_use,num_use):
    feature_list = []
    all_cat = {n:[] for n in categorical["GO"].keys()}
    for cat in cat_use:
        dic = categorical[cat]
        feature_list += numbering(dic)
        for symbol in dic:
            all_cat[symbol]+=list(dic[symbol])
    cat_num = {feature_list[n]:n for n in range(len(feature_list))}
    all_num = {f:numerical[f] for f in num_use}
    for f in all_num:
        for g in all_num[f]:
            all_num[f][g] = np.array(all_num[f][g]).astype("float")
    return feature_list,all_cat,cat_num,all_num

def experiment_context(input,output,threshold=0.5):
    sent, pmcid, cont, t = "","","",0
    f = open(input,"r")
    fw = open(output,"w")
    for line in tqdm.tqdm(f):
        ls = line.strip("\n").split("\t")
        if float(ls[0])>=threshold:
            if ls[3].upper() not in gene_to_sym:
                continue
            if pmcid != ls[2]:
                sentid = 0
                pmcid = ls[2]
            if sent != ls[-2]:
                sent = ls[-2]
                sentid += 1
            content = ",".join([ls[1], ls[2], str(sentid),
                   gene_to_sym[ls[3].upper()], ls[4]])
    
            #Avoid duplication
            if content != cont:
                fw.write(content + "\n")
                cont = content
    f.close()
    fw.close()
    print("Done!")


#sparse matrix
def make_sparse(exptup, col,row,data, k, C,fnum={},fdic={},mode="num"):
  for pre,nex in exptup:
    init_length = len(col)
    if mode=="string":
        if (pre,nex) in string:#stringはglobal
            col.append(C)
            data.append(string[(pre,nex)])
            row.append(k)
            
    elif mode=="cat":
        fpre,fnex = fdic[pre],fdic[nex]
        for fp in fpre:
          if fp in fnex:
            col.append(C+fnum[fp]*3)
          else:
            col.append(C+fnum[fp]*3+1)
        for fn in fnex:
          if fn not in fpre:
            col.append(C+fnum[fn]*3+2)
        row.extend([k] * (len(col)-init_length))
        data.extend([1] * (len(col)-init_length))
        
    elif mode=="num":
        fpre,fnex = fdic[pre],fdic[nex]
        if len(fnex)*len(fpre)!=0 and len(fnex)==len(fpre):
                 data.append(np.corrcoef(fpre,fnex)[0][1])
                 col.append(C)
                 row.append(k)
    k+=1        
  if mode=="cat":
      C += len(fnum)*3
  else:
      C += 1
  return col,row,data,C


def expespar(exptup, k, all_cat,cat_num,all_num,string=False):
    col,row,data = [],[],[]
    C = 0
    col,row,data,C = make_sparse(exptup, col,row,data, k, C,fdic=all_cat,fnum=cat_num,mode="cat")
    for f in all_num:
        dic = all_num[f]
        col,row,data,C = make_sparse(exptup, col,row,data, k, C,fdic=dic,fnum={},mode="num")
    if string:
        col,row,data,C = make_sparse(exptup, col,row,data, k, C,mode="string")
    return col,row,data, C

#generate tuples representing experimental contexts
from scipy.sparse import csr_matrix
import random
def make_tuple(path_to_csv,start,end,sampling=None):
    f = open(path_to_csv,"r")
    pmcid,sentid,genes = [],[],[]
    
    for line in f:
        ls = line.strip("\n").split(",")
        if ls[0].isnumeric():
            #Year Filter
            year = int(ls[0])
            if year > end or year < start:
                continue 

            #PMCID,sentid
            pmcid.append(ls[1])
            sentid.append(int(ls[2]))
            genes.append(ls[3])

    posi_tuple,nega_tuple=[],[]
    for n,p in enumerate(pmcid):
        n2 = 0
        while n+n2<len(pmcid) and sentid[n]+1 >= sentid[n+n2] and  p == pmcid[n+n2]:
            nn  = n+n2
            if sentid[n] +1 == sentid[nn]: 
                if genes[n]!=genes[nn]:
                     posi_tuple.append((genes[n],genes[nn]))
            n2 += 1
    f.close()
    
    if sampling != None:
        posi_tuple = random.sample(posi_tuple,sampling)

    #Negative sampling
    posi_tuple_all = set(posi_tuple)
    for g1,_ in posi_tuple:
            rand_g= random.sample(genes,20)
            for neg in rand_g:
                   if (g1,neg) not in posi_tuple_all and g1!=neg:
                        nega_tuple.append((g1,neg))
    nega_tuple = random.sample(nega_tuple,min(len(posi_tuple)*10,len(nega_tuple)))  
    return posi_tuple,nega_tuple       

def make_csr(posi_tuple,nega_tuple,all_cat,cat_num,all_num,string=False):
    print("Constructing CSR matrix...", end="  ")
    colp,rowp,datap, C = expespar(posi_tuple, 0,all_cat,cat_num,all_num,string)
    coln,rown,datan, C = expespar(nega_tuple, len(posi_tuple),all_cat,cat_num,all_num,string)
    row = rowp + rown
    col = colp + coln
    data = datap + datan
    y = [1 for _ in posi_tuple] + [0 for _ in nega_tuple]
    k = len(posi_tuple) + len(nega_tuple)
    x = csr_matrix((data, (row, col)), (k, C))
    print("Done")
    return x,y

def scoring(query,models,all_cat,cat_num,all_num,string=False):
    data,row,col,symbols,extup = [],[],[],[],[]
    
    with open("./data/symbols.txt","r") as f:
        for line in f:
            sym = line.strip()
            extup.append((query,sym))
            symbols.append(sym)
    
    extup = [(query,symbol) for symbol in symbols]
    extups = []
    extups.append(extup[:10000]) #memory
    extups.append(extup[10000:20000])
    extups.append(extup[20000:])
    scores={"Symbol":symbols}
    count = 0
    
    scores.update({m:[] for m in models})
    for n in range(len(extups)):
        col,row,data,C = expespar(extups[n],0,all_cat,cat_num,all_num,string)
        for m in models:
            model = models[m]
            score = model.predict_proba(csr_matrix((data, (row,col)),(len(extups[n]), C)))[:,1].tolist()
            scores[m] += score
        count += 1
        if count % 2000==0:
            print(count,"/",len(symbols))
    
    return scores



