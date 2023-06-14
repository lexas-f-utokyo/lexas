import gdown

def download_from_gdrive(url: str, output: str):
    gdown.download(url, output, quiet=False)

if __name__ == "__main__":
    url= 'https://drive.google.com/uc?id=1s8Na00l2GbxW112sfcb5q1g0jS-Ujp9h'
    output = '../Repository.tar.gz'
    download_from_gdrive(url, output)
