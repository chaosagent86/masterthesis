import os
import gzip
import shutil
import requests

def read_accessions(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if line.strip()]

def download_and_unpack(url, outdir):
    local_filename = os.path.join(outdir, os.path.basename(url))
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    # Entpacken falls .gz
    if local_filename.endswith(".gz"):
        with gzip.open(local_filename, 'rb') as f_in:
            with open(local_filename[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(local_filename)

def get_ftp_url(gcf):
    # GCF_043076235.1 → GCF/043/076/235
    parts = gcf.split('_')[1].split('.')
    digits = parts[0].zfill(9)
    path = f"{digits[:3]}/{digits[3:6]}/{digits[6:]}"
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{path}/{gcf}"

def download_genomes(accession_list, outdir="downloads"):
    os.makedirs(outdir, exist_ok=True)

    for acc in accession_list:
        ftp_base = get_ftp_url(acc)

        # Dateien: fna + gff
        files = [
            f"{acc}_genomic.fna.gz",
            f"{acc}_genomic.gff.gz"
        ]

        for file in files:
            url = f"{ftp_base}/{file}"
            try:
                print(f"Downloading: {url}")
                download_and_unpack(url, outdir)
            except Exception as e:
                print(f"Fehler beim Download von {acc}: {e}")

if __name__ == "__main__":
    accession_file = "accessions.txt"  # Eine Textdatei mit einer GCF-Nummer pro Zeile
    accessions = read_accessions(accession_file)
    download_genomes(accessions)
