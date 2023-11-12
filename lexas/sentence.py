import xml.etree.ElementTree as ET
import os
import concurrent.futures
import re
from ahocorapy.keywordtree import KeywordTree

# Importing necessary modules for text processing
import scispacy
import en_core_sci_sm
import tqdm.notebook as tqdm

# Importing the module to convert gene names to symbols
import lexas.gene_to_sym
gene_to_sym = lexas.gene_to_sym.get()

# Creating a pattern for identifying and removing newline and tab characters
pattern = re.compile(r'\t|\n|\x0b|\x0c|\r|\x1c|\x1d|\x1e|\x1f|\x85|\xa0|\u1680|\u2000|\u2001|\u2002|\u2003|\u2004|\u2005|\u2006|\u2007|\u2008|\u2009|\u200a|\u2028|\u2029|\u202f|\u205f|\u3000')

# Function to join sections of text and remove newline and tab characters
def join_and_remove_nt(section):
    text = " ".join([r for r in section.itertext()]).strip()
    text = pattern.sub(' ', text)
    return text

# Function to parse the XML file and filter sections by title
def parse(xml_file, title_filter="result"):
    year = 0
    sentences = []
    paragraph_numbers = []
    try:
        # Parsing the XML file
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Extracting the publication year
        year = root.find(".//article-meta").find("pub-date").find("year").text

        # Extracting sections with the specified title
        for section in root.findall(".//sec"):
            try:
                section_title = section.find('title').text.lower()
                if title_filter in section_title:
                    for paragraph_number, child in enumerate(section.findall(".//p")):
                        sentences.append(join_and_remove_nt(child))
                        paragraph_numbers.append(paragraph_number)
            except:
                continue
    except:
        pass
    return year, sentences, paragraph_numbers

# Function to extract results from articles, optionally by paragraph
def extract_results(article_dir="./articles/", output_file="./data/result_sections.txt", title_filter="result", include_paragraph=True):
    with open(output_file, "w") as output:
        # Iterate over all files in the specified directory
        for filename in tqdm.tqdm(os.listdir(article_dir)):
            # Extracting the ID from the filename
            article_id = filename.split(".")[0]
            # Parsing the article
            year, sentences, paragraph_numbers = lexas.sentence.parse(os.path.join(article_dir,filename), title_filter=title_filter)
            
            # Checking if a valid year was found
            if year != 0:
                for idx in range(len(sentences)):
                    sentence = sentences[idx]
                    paragraph_number = str(paragraph_numbers[idx])
                    # Segmenting the sentences
                    segmented_sentence = lexas.sentence.segment_sentence(sentence)
                    # Writing the results to file, including the paragraph number if requested
                    if include_paragraph:
                        output.write("\t".join([year, article_id + "-" + paragraph_number, segmented_sentence]) + "\n")
                    else:
                        output.write("\t".join([year, article_id, segmented_sentence]) + "\n")


# Loading the ScispaCy model
en_model = en_core_sci_sm.load()

# Function to segment a sentence into separate sentences
def segment_sentence(sentence):
    doc = en_model(sentence)
    # Joining the sentences with '#####', ignoring sentences longer than 500 words
    segmented = "#####".join([str(sent) for sent in doc.sents if len(str(sent).split(" ")) < 500])
    return segmented


def initialize_dictionaries():
    """
    Initializes two KeywordTree objects for gene symbols and experiments.
    """
    dic_hgnc = KeywordTree(case_insensitive=True)
    dic_expe = KeywordTree(case_insensitive=True)

    # Gene symbols
    with open('./data/terms_to_be_removed.csv') as remo:
        to_remove = [word.strip().lower() for word in remo]
        
    for key in gene_to_sym:
        gene = key.lower()
        if len(gene) > 3 and gene not in to_remove:
            dic_hgnc.add(" " + gene + " ")
    dic_hgnc.finalize()

    # Experiments
    with open('./data/experiments_list.txt') as f:
        for experiment in f:
            dic_expe.add(" " + experiment.strip().lower() + " ")
    dic_expe.finalize()

    return dic_hgnc, dic_expe

def mask_text(text, dic_gene, dic_expe):
    """
    Masks gene symbols and experiment names in the given text.
    """
    cleaned_text = " " + text + " "
    cleaned_text = cleaned_text.replace(",","").replace(".","").replace("("," ").replace(")"," ")
    masked_texts = []

    for gene in dic_gene.search_all(cleaned_text):
        masked_text, gene_name = cleaned_text[:gene[1]] + " [GENE] " + cleaned_text[gene[1] + len(gene[0]):], gene[0].strip()
        for experiment in dic_expe.search_all(masked_text):
            result = [masked_text[:experiment[1]] + " [EXPE] " +  masked_text[experiment[1] + len(experiment[0]):], gene_name, experiment[0].strip()]
            masked_texts.append(result)
    
    return masked_texts

def load_dictionaries():
    """
    Loads the dictionaries if they aren't already loaded, and returns them.
    """
    if not hasattr(load_dictionaries, "dic_hgnc") or not hasattr(load_dictionaries, "dic_expe"):
        print("Initializing dictionaries...")
        load_dictionaries.dic_hgnc, load_dictionaries.dic_expe = initialize_dictionaries()
        print("Done")
    return load_dictionaries.dic_hgnc, load_dictionaries.dic_expe

def process_line(line):
    """
    Processes a line from the input file and returns its constituent parts.
    """
    try:
        year, pmcid, sentences = line.strip("\n").split("\t")
        return year, pmcid, sentences.split("#####")
    except:
        return None

def mask_gene_experiment(input_file_path="./data/result_sections.txt", output_file_path="./data/masked_sentences.txt"):
    dic_hgnc, dic_expe = load_dictionaries()
    
    with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
        for line in tqdm.tqdm(input_file):
            processed_line = process_line(line)
            
            if processed_line is not None:
                year, pmcid, sentences = processed_line
                for sentence in sentences:
                    masked_texts = mask_text(sentence, dic_hgnc, dic_expe)
                    for masked in masked_texts:
                        output_file.write("\t".join([year, pmcid] + masked[1:3] + [sentence] + masked[0:1]) + "\n")

