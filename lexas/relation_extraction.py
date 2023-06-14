import torch
import tqdm
from transformers import BertTokenizer, BertForSequenceClassification

def predict(device, input_filepath, output_filepath, threshold=0.5):
    """Predicts relations based on a pretrained BERT model.

    Args:
        device: The device (cpu or gpu) to use for prediction.
        input_filepath: Path to the file containing sentences to predict relations for.
        output_filepath: Path to the file to write the predictions to.
        threshold: Probability threshold for classifying a sentence as having a relation.
    """

    # Load the pretrained BERT model and its tokenizer
    tokenizer = BertTokenizer.from_pretrained('../Repository/biobert/tokenizer_biobert', do_lower_case=True)
    model = BertForSequenceClassification.from_pretrained(
        "../Repository/biobert/model_biobert-finetuning",
        num_labels = 2,
        output_attentions = False,
        output_hidden_states = False,  # Changed to False as hidden states are not used
    ).to(device)

    # Set model to evaluation mode
    model.eval()

    softmax = torch.nn.Softmax(dim=1)

    # Open input and output files
    with open(input_filepath, 'r') as input_file, open(output_filepath, 'w') as output_file:
        # Iterate over lines in the input file
        for line in tqdm.tqdm(input_file):
            sentence = line.strip().split("\t")[5]

            # Check if the number of tokens in the sentence is under the limit
            if len(tokenizer.tokenize(sentence)) > 500:
                continue

            # Tokenize the sentence and convert to Tensor
            input_ids = tokenizer.encode(sentence, add_special_tokens=True, return_tensors='pt').to(device)

            # Make prediction
            outputs = model(input_ids)
            probs = softmax(outputs.logits)

            # If the probability of the positive class exceeds the threshold, write to the output file
            if probs[0, 1] > threshold:
                output_file.write(f"{probs[0, 1].item()}\t{line}")
