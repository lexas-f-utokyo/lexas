import tqdm
import os
import numpy as np
import pandas as pd
import csv
from sklearn.metrics import roc_auc_score
from lexas import prediction
symbols = prediction.symbols

def generate_dic_for_eval(path_to_csv, start_year, end_year, mode):
    """
    Extract gene combinations from a CSV file.

    Parameters:
    - path_to_csv: path to the input CSV file
    - start_year, end_year: range of accepted published years
    - mode: 
    
    Returns:
    - a set of tuple combinations of genes
    """
    if mode=="just_next":
        gene_combinations, _ = prediction.generate_experiment_tuples(path_to_csv, start_year, end_year)
        return tuple_to_dic(gene_combinations)
    elif mode=="all_following":
        gene_combinations = generate_tuples_all_following(path_to_csv, start_year, end_year)
        return tuple_to_dic(gene_combinations)
    else:
        raise Exception("Mode should be \"just_next\" or \"all_following\"")
    
def generate_tuples_all_following(path_to_csv, start, end):
    # Resulting set to hold gene combinations
    gene_combinations = set()

    # Current list of genes for a specific PMCID
    current_genes = []

    # Initialize with an invalid PMCID
    current_pmcid = None

    with open(path_to_csv, "r") as file:
        # Using csv.reader for better handling of CSV files
        reader = csv.reader(file)
        
        for row in tqdm.tqdm(reader):
            year, _, pmcid, gene = int(row[0]), row[1], row[2], row[3]

            if not (start <= year <= end):
                continue
            
            # If we have a new PMCID, process the accumulated genes for the previous PMCID
            if pmcid != current_pmcid:
                for i, gene_a in enumerate(current_genes):
                    for gene_b in current_genes[i+1:]:
                        gene_combinations.add((gene_a, gene_b))
                current_genes = []
                current_pmcid = pmcid

            current_genes.append(gene)

    return gene_combinations

# Dictionary to hold results
def tuple_to_dic(expe_tuple):
    next_gene_dic = {s:set() for s in symbols}
    for first, second in expe_tuple:
        next_gene_dic[first].add(second)
    return next_gene_dic


import csv
def save_dic(train_dic,dev_dic,test_dic,mode):
    with open('./eval/dic_{}.csv'.format(mode), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['key', 'train_value', 'dev_value', 'test_value'])
        for key in train_dic.keys():
            train_value = ','.join(map(str, train_dic[key]))
            dev_value = ','.join(map(str, dev_dic[key]))
            test_value = ','.join(map(str, test_dic[key]))
            writer.writerow([key, train_value, test_value, dev_value])
            
def load_dic(mode):
    train_dic, dev_dic, test_dic = {}, {}, {}
    with open('./eval/dic_{}.csv'.format(mode), 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            key, train_value, dev_value, test_value = row
            train_dic[key] = train_value.split(',')
            dev_dic[key] = dev_value.split(',')
            test_dic[key] = test_value.split(',')
    return train_dic,dev_dic,test_dic

def prepare_df_for_eval(df, query_gene, prev_dic, answer_dic):
    df_empty = pd.DataFrame(index=symbols)
    df = pd.merge(df_empty, df, left_index=True, right_index=True,how="outer").reindex(index=symbols)
    df["prev"],df["answer"] = False,False
    
    for gene in prev_dic[query_gene]:
        if gene.strip():
            df.at[gene,"prev"] = True
    
    for gene in answer_dic[query_gene]:
        if gene.strip():
            df.at[gene,"answer"] = True
    
    filtered_df = df.drop(index=["NA", query_gene])
    return filtered_df

def auc_at_k(df, query_gene, prev_dic, answer_dic, top_k):
    df = prepare_df_for_eval(df, query_gene, prev_dic, answer_dic)
    
    # Genes examined in previous data are not considered false, so they are removed from evaluation.
    df = df[(df["prev"]==False) | (df["answer"]==True)]
    
    y_true = df["answer"].to_numpy()# True labels for AUC computation
  
    if sum(y_true):
        threshold_score = -1*np.sort(-1*df.iloc[:,0])[top_k]
        y_pred = [0 if np.isnan(score) or score < threshold_score else score for score in df.iloc[:,0]]
        auc_score = roc_auc_score(y_true,y_pred)
    else:
        return None
        
    return auc_score

def calculate_auc_for_many_genes(result_dir, model, prev_dic, answer_dic, top_k, genes=symbols):
    auc_results = []
    for query_gene in tqdm.tqdm(genes):
        path = os.path.join(result_dir, model, f"{query_gene}.csv")
        
        if answer_dic[query_gene]==[""]:
            auc_results.append(None)
            continue

        # Ensure that the file exists before trying to read it
        if os.path.exists(path):
            df = pd.read_csv(path, index_col=0)
            auc_score = auc_at_k(df, query_gene, prev_dic, answer_dic, top_k)
            auc_results.append(auc_score)
        else:
            print(f"Warning: File not found for gene {query_gene}. Skipping...")
            auc_results.append(None)
        
    df_results = pd.DataFrame({"Symbol":genes,f"AUC at {top_k} for {model}":auc_results})
    return df_results