import os
import gzip
import shutil
import requests
from bs4 import BeautifulSoup

def download_assembly_summary():
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    r = requests.get(url)
    with open("assembly_summary_refseq.txt", "wb") as f:
        f.write(r.content)

def parse_assembly_summary(file_path="assembly_summary_refseq.txt"):
    acc_to_path = {}
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 20:
                continue
            acc, ftp_path = cols[0], cols[19]
            acc_to_path[acc] = ftp_path
    return acc_to_path

def read_accessions(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if line.strip()]

def list_ftp_files(ftp_url):
    page = requests.get(ftp_url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [a['href'] for a in soup.find_all('a') if a['href'].endswith('.gz')]

def download_and_unpack(url, outdir):
    local_filename = os.path.join(outdir, os.path.basename(url))
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    if local_filename.endswith(".gz"):
        with gzip.open(local_filename, 'rb') as f_in:
            with open(local_filename[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(local_filename)

def download_genomes(accessions, mapping, outdir="downloads", wanted_suffixes=(".fna.gz", ".gff.gz")):
    os.makedirs(outdir, exist_ok=True)
    for acc in accessions:
        if acc not in mapping:
            print(f"{acc} nicht in assembly_summary gefunden.")
            continue
        ftp_url = mapping[acc].replace("ftp://", "https://") + "/"
        try:
            files = list_ftp_files(ftp_url)
        except Exception as e:
            print(f"Fehler beim Auflisten des FTP-Verzeichnisses für {acc}: {e}")
            continue

        for file in files:
            if any(file.endswith(suffix) for suffix in wanted_suffixes):
                file_url = ftp_url + file
                try:
                    print(f"Downloading: {file_url}")
                    download_and_unpack(file_url, outdir)
                except Exception as e:
                    print(f"Fehler beim Download von {file_url}: {e}")

if __name__ == "__main__":
    from bs4 import BeautifulSoup  # Sicherstellen, dass bs4 installiert ist

    if not os.path.exists("assembly_summary_refseq.txt"):
        download_assembly_summary()

    mapping = parse_assembly_summary("assembly_summary_refseq.txt")
    accessions = read_accessions("accessions.txt")
    download_genomes(accessions, mapping)
