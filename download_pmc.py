from Bio import Entrez
Entrez.email = 'A.N.Other@example.com'#your email address

def download(pmc_id):
    article_xml = Entrez.efetch(db='pmc',
                     id=pmc_id, 
                     rettype='xml').read()
    file_name = f'{pmc_id}.nxml' 
    with open("articles/"+file_name, 'w') as file:
        file.write(article_xml.decode())

if __name__ == "__main__":
    pmc_ids = ["PMC3089914",
              "PMC6472344",
              "PMC4220463",
              "PMC6777370",
              "PMC4291483",
              "PMC6829667",
              "PMC5007451",
              "PMC7812875",
              "PMC5509424",
              "PMC7836272",
              "PMC6389942",
              "PMC7938805"]
    for pmc_id in pmc_ids:
        download(pmc_id)