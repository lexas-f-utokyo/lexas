import transformers
import torch
import tqdm
from transformers import BertTokenizer
from transformers import BertForSequenceClassification

def predict(device,input,output,threshold = 0.5):
    tokenizer = BertTokenizer.from_pretrained('Repository/biobert/tokenizer_biobert', do_lower_case=True)
    model = BertForSequenceClassification.from_pretrained(
         "Repository/biobert/model_biobert-finetuning",
         num_labels = 2, 
         output_attentions = False, # Whether the model returns attentions weights.
         output_hidden_states = True, # Whether the model returns all hidden-states.
    ).to(device)
    
    model.eval()
    m = torch.nn.Softmax()
    
    with torch.no_grad():
        f = open(input,"r")
        f2 = open(output,"w")
        for line in tqdm.tqdm(f):
            sent = line.strip().split("\t")[5]
            if len(sent.split(" ")) > 500:
                continue
            indexed_tokens = tokenizer.encode(sent,add_special_tokens = True)
            tokens_tensor = torch.tensor([indexed_tokens])
            all_encoder_layers = model(tokens_tensor.to(device))
            out= m(all_encoder_layers[0])
            if out[0][1] > threshold:
                 f2.write(str(out[0][1].item())+"\t"+line)
        f.close()
        f2.close() 