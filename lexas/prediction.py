import numpy as np
import tqdm
import lexas.gene_to_sym

# Initializing dictionaries at the module level
categorical, numerical, string, funcoup, semsim = {}, {}, {}, {}, {}

def load_categorical_features():
    """Loads categorical features into a dictionary"""
    print("Loading categorical features...")
    with open("../Repository/feature/categorical_feature.txt", "r") as f:
        for line in f:
            feature, subfeature, values = line.strip("\n").split("\t")
            categorical.setdefault(feature, {})
            categorical[feature][subfeature] = [] if not values else values.split(",")

def load_numerical_features():
    """Loads numerical features into a dictionary"""
    print("Loading numerical features...")
    with open("../Repository/feature/numerical_feature.txt", "r") as f:
        for line in f:
            feature, gene, values = line.strip("\n").split("\t")
            numerical.setdefault(feature, {})
            if values=="":
                numerical[feature][gene] = []
            else:
                values_list = [float(s) for s in values.split(",")]
                # Skip if all values are the same
                if len(set(values_list)) > 1:
                    numerical[feature][gene] = values_list
                else:
                    numerical[feature][gene] = []

def load_feature(filepath, feature_dict):
    """Loads features from a given filepath into a given dictionary"""
    print(f"Loading {filepath.split('/')[-1]}...")
    with open(filepath, "r") as f:
        for line in f:
            key1, key2, value = line.strip("\n").split("\t")
            feature_dict[tuple([key1,key2])] = float(value)

def feature_load():
    """Load features into dictionaries"""
    # Load categorical and numerical features
    if categorical == {}: 
        load_categorical_features()
        load_numerical_features()

        # Load STRING, FunCoup, and GoSemSim features
        load_feature("../Repository/feature/string11_rwr.txt", string)
        load_feature("../Repository/feature/funcoup5_rwr.txt", funcoup)
        load_feature("../Repository/feature/gosemsim.txt", semsim)

# Load gene_to_sym and symbols at the module level
gene_to_sym = lexas.gene_to_sym.get()

with open("./data/symbols.txt", "r") as f:
    symbols = [line.strip() for line in f]

def create_feature_mapping(feature_dict):
    """Generates a numerical mapping of the features"""
    all_features = set(feature for sub_dict in feature_dict.values() for feature in sub_dict)
    feature_list = sorted(all_features)
    return {feature: index for index, feature in enumerate(feature_list)}

def select_features(categorical_features, numerical_features):
    """Selects and formats the required features for further processing"""
    feature_mapping = []
    gene_to_categorical = {feature: [] for feature in categorical["GO"].keys()}
    for feature in categorical_features:
        sub_dict = categorical[feature]
        feature_mapping += create_feature_mapping(sub_dict)
        for symbol in sub_dict:
            gene_to_categorical[symbol] += list(sub_dict[symbol])

    gene_to_numerical = {feature: numerical[feature] for feature in numerical_features}
    for feature in gene_to_numerical:
        for symbol in gene_to_numerical[feature]:
            gene_to_numerical[feature][symbol] = np.array(gene_to_numerical[feature][symbol]).astype("float")

    return feature_mapping, gene_to_categorical, gene_to_numerical

def extract_context_from_experiments(input_file, output_file, threshold=0.5):
    """Extracts context from experiment data and writes to an output file"""
    sent, prev_pmcid, prev_sentence, prev_content, sentid = "", "", "", "", 0
    with open(input_file, "r") as input_f, open(output_file, "w") as output_f:
        for line in tqdm.tqdm(input_f):
            score, year, pmcid, gene, experiment, sentence, _  = line.strip("\n").split("\t")
            if float(score) >= threshold:
                if gene.upper() not in gene_to_sym:
                    continue
                if pmcid != prev_pmcid:
                    sentid = 0
                    prev_pmcid = pmcid
                if sentence != prev_sentence:
                    prev_sentence = sentence
                    sentid += 1
                content = ",".join([year, pmcid, str(sentid), gene_to_sym[gene.upper()], experiment])

                # Avoid duplication
                if content != prev_content:
                    output_f.write(content + "\n")
                    prev_content = content
    print("Done!")

# This function builds a sparse matrix for the given experiment tuple, depending on the mode.
def build_sparse_matrix(experiment_tuples, col, row, data, k, C, feature_num={}, feature_dict={}, mode="num"):
    for pre, nex in experiment_tuples:
        init_length = len(col)
        if mode == "string":
            if (pre, nex) in feature_dict:
                col.append(C)
                data.append(feature_dict[(pre, nex)])
                row.append(k)

        elif mode == "cat":
            fpre, fnex = feature_dict[pre], feature_dict[nex]
            for fp in fpre:
                col.append(C + feature_num[fp]*3 if fp in fnex else C + feature_num[fp]*3 + 1)
            for fn in fnex:
                if fn not in fpre:
                    col.append(C + feature_num[fn]*3 + 2)
            row.extend([k] * (len(col) - init_length))
            data.extend([1] * (len(col) - init_length))

        elif mode == "num":
            fpre, fnex = feature_dict[pre], feature_dict[nex]
            if len(fnex) * len(fpre) != 0 and len(fnex) == len(fpre):
                data.append(np.corrcoef(fpre, fnex)[0][1])
                col.append(C)
                row.append(k)
        k += 1
    C += len(feature_num)*3 if mode=="cat" else 1
    return col, row, data, C

# This function processes experiment data into a sparse matrix representation.
def generate_sparse_data(experiment_tuples, k, gene_categories, feature_list, gene_numbers, additional_features=[]):
    category_numbers = {feature_list[n]: n for n in range(len(feature_list))}

    col, row, data = [], [], []
    C = 0

    col, row, data, C = build_sparse_matrix(experiment_tuples, col, row, data, k, C, feature_dict=gene_categories, feature_num=category_numbers, mode="cat")

    for feature in gene_numbers:
        feature_dict = gene_numbers[feature]
        col, row, data, C = build_sparse_matrix(experiment_tuples, col, row, data, k, C, feature_dict=feature_dict, mode="num")

    if "String" in additional_features:
        col, row, data, C = build_sparse_matrix(experiment_tuples, col, row, data, k, C, mode="string", feature_dict=string)
    if "Funcoup" in additional_features:
        col, row, data, C = build_sparse_matrix(experiment_tuples, col, row, data, k, C, mode="string", feature_dict=funcoup)
    if "GOSemSim" in additional_features:
        col, row, data, C = build_sparse_matrix(experiment_tuples, col, row, data, k, C, mode="string", feature_dict=semsim)

    return col, row, data, C


from scipy.sparse import csr_matrix
import random

# Generate tuples representing experimental contexts
def generate_experiment_tuples(csv_path, start_year, end_year, sampling=None, negative_sampling=0):
    """
    Generate tuples representing experimental contexts. 
    Each tuple is a pair of gene symbols appearing in the same context.
    Negative sampling is also performed if 'negative_sampling' is not 0.

    :param csv_path: Path to the csv file.
    :param start_year: Start year for filtering the data.
    :param end_year: End year for filtering the data.
    :param sampling: If given, sample this amount of tuples from the entire dataset.
    :param negative_sampling: Ratio of negative to positive tuples.
    :return: A tuple of sets where the first element is the set of positive tuples, 
             and the second is the set of negative tuples.
    """
    with open(csv_path, "r") as file:
        pmcids, sentence_ids, gene_symbols = [], [], []

        for line in file:
            elements = line.strip("\n").split(",")
            if elements[0].isnumeric():
                # Filter by year
                year = int(elements[0])
                if start_year <= year <= end_year:
                    # PMCID, sentence id, and gene symbol
                    pmcids.append(elements[1])
                    sentence_ids.append(int(elements[2]))
                    gene_symbols.append(elements[3])

        positive_tuples, negative_tuples = set(), set()

        # Generate tuples
        for index, pmcid in enumerate(pmcids):
            offset = 0
            while index + offset < len(pmcids) and sentence_ids[index] + 1 >= sentence_ids[index + offset] and pmcid == pmcids[index + offset]:
                next_index = index + offset
                if sentence_ids[index] + 1 == sentence_ids[next_index] and gene_symbols[index] != gene_symbols[next_index]:
                    positive_tuples.add((gene_symbols[index], gene_symbols[next_index]))
                offset += 1

        # Negative sampling
        if negative_sampling != 0:
            for positive_tuple in positive_tuples:
                random_genes = random.sample(gene_symbols, negative_sampling * 3)
                for negative_gene in random_genes:
                    if (positive_tuple[0], negative_gene) not in positive_tuples and positive_tuple[0] != negative_gene:
                        negative_tuples.add((positive_tuple[0], negative_gene))

        # Sampling
        if sampling is not None:
            positive_tuples = set(random.sample(positive_tuples, sampling // (1 + negative_sampling)))
            negative_tuples = set(random.sample(negative_tuples, min(len(negative_tuples), sampling * negative_sampling // (1 + negative_sampling))))
        else:
            negative_tuples = set(random.sample(negative_tuples, min(len(negative_tuples), len(positive_tuples) * negative_sampling)))

        return positive_tuples, negative_tuples


from scipy.sparse import csr_matrix

def construct_csr_matrix(positive_tuples, negative_tuples, categorical_features, feature_list, numerical_features, additional_features=[]):
    """
    Constructs a CSR (Compressed Sparse Row) matrix using the provided tuples and features.

    :param positive_tuples: Set of tuples representing positive gene pairs.
    :param negative_tuples: Set of tuples representing negative gene pairs.
    :param categorical_features: Dictionary containing categorical features of genes.
    :param feature_list: List of all possible features.
    :param numerical_features: Dictionary containing numerical features of genes.
    :param additional_features: Additional features to be considered. Defaults to an empty list.
    :return: CSR matrix and labels.
    """
    print("Constructing CSR matrix...", end="  ")
    col_pos, row_pos, data_pos, C = generate_sparse_data(positive_tuples, 0, categorical_features, feature_list, numerical_features, additional_features=additional_features)
    col_neg, row_neg, data_neg, C = generate_sparse_data(negative_tuples, len(positive_tuples), categorical_features, feature_list, numerical_features, additional_features=additional_features)
    
    row = row_pos + row_neg
    col = col_pos + col_neg
    data = data_pos + data_neg
    labels = [1 for _ in positive_tuples] + [0 for _ in negative_tuples]
    num_rows = len(positive_tuples) + len(negative_tuples)
    
    matrix = csr_matrix((data, (row, col)), (num_rows, C))
    print("Done")
    return matrix, labels

def generate_scores(query, models, categorical_features, feature_list, numerical_features, additional_features=[]):
    """
    Generate prediction scores for gene pairs using the given models.

    :param query: Query gene symbol.
    :param models: Dictionary containing models to be used for prediction.
    :param categorical_features: Dictionary containing categorical features of genes.
    :param feature_list: List of all possible features.
    :param numerical_features: Dictionary containing numerical features of genes.
    :param additional_features: Additional features to be considered. Defaults to an empty list.
    :return: Dictionary containing prediction scores for each model.
    """
    data, row, col, symbols, experiment_tuples = [], [], [], [], []
    experiment_tuples_batched = []

    with open("./data/symbols.txt", "r") as file:
        for line in file:
            symbol = line.strip()
            experiment_tuples.append((query, symbol))
            symbols.append(symbol)

    # Divide experiment_tuples into batches due to memory limitations
    batch_size = 10000
    for i in range(0, len(experiment_tuples), batch_size):
        experiment_tuples_batched.append(experiment_tuples[i:i + batch_size])

    scores = {"Symbol": symbols}
    scores.update({model_name: [] for model_name in models})

    for batch in experiment_tuples_batched:
        col, row, data, C = generate_sparse_data(batch, 0, categorical_features, feature_list, numerical_features, additional_features=additional_features)
        for model_name, model in models.items():
            prediction_scores = model.predict_proba(csr_matrix((data, (row, col)), (len(batch), C)))[:, 1].tolist()
            scores[model_name] += prediction_scores

    return scores
