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
    - mode: Determines which combinations to extract.
        * "just_next": Only consider experiments performed just after the current experiment.
        * "all_following": Consider all experiments performed after the current experiment in the article.
    
    Returns:
    - a set of tuple combinations of genes
    """
    modes = {
        "just_next": lambda: tuple_to_dic(prediction.generate_experiment_tuples(path_to_csv, start_year, end_year)[0]),
        "all_following": lambda: tuple_to_dic(generate_tuples_all_following(path_to_csv, start_year, end_year))
    }

    if mode not in modes:
        raise Exception("Mode should be \"just_next\" or \"all_following\"")
        
    return modes[mode]()

    
def generate_tuples_all_following(path_to_csv, start, end):
    """
    Extract all pairwise gene combinations from a CSV file for articles published within a specified year range.

    Parameters:
    - path_to_csv: Path to the input CSV file.
    - start: Start year (inclusive) for filtering articles.
    - end: End year (inclusive) for filtering articles.

    Returns:
    set: A set containing tuples of gene combinations.
    """
    
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
            year, pmcid, _, gene = int(row[0]), row[1], row[2], row[3]
            
            if "-" in pmcid:
                pmcid = pmcid.split("-")[0]

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
            
    # Last article
    for i, gene_a in enumerate(current_genes):
        for gene_b in current_genes[i+1:]:
            gene_combinations.add((gene_a, gene_b))

    return gene_combinations


def tuple_to_dic(expe_tuple):
    """
    Convert a list of gene tuple combinations into a dictionary representation.
    
    Parameters:
    - expe_tuple: List of gene combinations where each tuple is (gene_a, gene_b) indicating experiment on gene_a is followed by that on gene_b.
    
    Returns:
    dict: A dictionary where keys are individual genes and values are sets of genes that follow the key gene.
    """
    next_gene_dic = {s:set() for s in symbols}
    for first, second in expe_tuple:
        try:
            next_gene_dic[first].add(second)
        except KeyError:
            continue
    return next_gene_dic


import csv
def save_dic(train_dic,dev_dic,test_dic,mode,output_dir="./eval"):
    #Save dictionaries to a CSV file.
    with open(os.path.join(output_dir, f'dic_{mode}.csv'), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['key', 'train_value', 'dev_value', 'test_value'])
        for key in train_dic.keys():
            train_value = ','.join(map(str, train_dic[key]))
            dev_value = ','.join(map(str, dev_dic[key]))
            test_value = ','.join(map(str, test_dic[key]))
            writer.writerow([key, train_value, dev_value, test_value])
            
def load_dic(mode,input_dir="./eval"):
    #Load dictionaries to a CSV file.
    train_dic, dev_dic, test_dic = {}, {}, {}
    with open(os.path.join(input_dir, f'dic_{mode}.csv'), 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            key, train_value, dev_value, test_value = row
            train_dic[key] = train_value.split(',')
            dev_dic[key] = dev_value.split(',')
            test_dic[key] = test_value.split(',')
    return train_dic,dev_dic,test_dic


def prepare_df_for_eval(df, query_gene, prev_dic, answer_dic):
    """
    Prepare a DataFrame for evaluation by setting 'prev' and 'answer' columns.
    
    Parameters:
    - df (DataFrame): Input DataFrame containing results of prediction score.
    - query_gene (str): Gene of interest.
    - prev_dic (dict): Dictionary containing genes previously examined after the query_gene (neighter correct or incorrect).
    - answer_dic (dict): Dictionary containing correct genes examined after the query_gene.
    
    Returns:
    - DataFrame: A processed DataFrame, ready for evaluation.
    
    - 'prev' column indicates whether a gene was examined previously.
    - 'answer' column indicates the correct gene.
    """    
    
    df_empty = pd.DataFrame(index=symbols)
    df = pd.merge(df_empty, df, left_index=True, right_index=True,how="outer").reindex(index=symbols)
    df["prev"],df["answer"] = False,False
    
    for gene in prev_dic[query_gene]:
        if gene.strip():
            df.at[gene,"prev"] = True
    
    for gene in answer_dic[query_gene]:
        if gene.strip():
            df.at[gene,"answer"] = True
    
    filtered_df = df.drop(index=[query_gene]).dropna()
    return filtered_df


def auc_at_k(df, query_gene, prev_dic, answer_dic, top_k, remove_all_prev = False):
    """
    Calculate AUC score for a given query gene based on the top k predictions.
    
    Parameters:
    - df (DataFrame): Input DataFrame containing prediction scores.
    - query_gene (str): Gene of interest.
    - prev_dic (dict): Dictionary of previously examined genes.
    - answer_dic (dict): Dictionary of correct genes.
    - top_k (int): Number of top predictions to consider.
    
    Returns:
    - float: AUC score, or None if there's no true positive in y_true.
    """
    
    df = prepare_df_for_eval(df, query_gene, prev_dic, answer_dic)
    
    # Genes examined in previous data are not considered false, so they are removed from evaluation.
    df = df[(df["prev"]==False) | (df["answer"]==True)]
    
    if remove_all_prev:
        df =  df[df["prev"]==False]
    
    y_true = df["answer"].to_numpy().astype(int)# True labels for AUC computation

    # Add a small random value to scores to prevent identical scores.
    np.random.seed(100)
    df.iloc[:,0] = df.iloc[:,0] + np.random.uniform(0, 1e-10, len(df))
    
    if sum(y_true):
        try:
            threshold_score = -1*np.sort(-1*df.iloc[:,0])[top_k]
        except IndexError: #String_raw or Funoup_raw
            threshold_score = 0
        y_score = [0 if np.isnan(score) or score < threshold_score else score for score in df.iloc[:,0]]
        try:
            auc_score = roc_auc_score(y_true,y_score)
        except ValueError: #FunCoup_raw
            return None
    else:
        return None
        
    return auc_score

def calculate_auc_for_many_genes(result_dir, model, prev_dic, answer_dic, top_k, genes=symbols, remove_all_prev = False):
    """
    Calculate AUC scores for multiple query genes based on a model's predictions.
    
    Parameters:
    - result_dir (str): Directory containing the prediction results.
    - model (str): Model name used to find the corresponding prediction results.
    - prev_dic (dict): Dictionary of previously examined genes.
    - answer_dic (dict): Dictionary of correct genes.
    - top_k (int): Number of top predictions to consider.
    - genes (list): List of genes to calculate AUC scores for. 
    
    Returns:
    - DataFrame: A DataFrame containing the AUC scores for each gene.
    """
    
    auc_results = []
    for query_gene in tqdm.tqdm(genes):
        path = os.path.join(result_dir, f"{query_gene}.csv")
        
        if answer_dic[query_gene]==[""]:
            auc_results.append(None)
            continue

        # Ensure that the file exists before trying to read it
        if os.path.exists(path):
            df = pd.read_csv(path, index_col=0)
            df = df[~df.index.duplicated(keep='first')]
            auc_score = auc_at_k(df, query_gene, prev_dic, answer_dic, top_k, remove_all_prev=remove_all_prev)
            auc_results.append(auc_score)
        else:
            print(f"Warning: File not found for gene {query_gene}. Skipping...")
            auc_results.append(None)
        
    df_results = pd.DataFrame({"Symbol":genes,f"AUC at {top_k} for {model}":auc_results})
    return df_results