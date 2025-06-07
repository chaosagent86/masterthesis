import os
import time
import pandas as pd
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "sebasbaier@gmail.com"

# Spezies aus Datei einlesen
with open("00_List_of_streptomyces_basis.txt") as f:
    species_list = [line.strip() for line in f if line.strip()]

# Qualitätsstufen in Präferenzreihenfolge
quality_order = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]

# Output-Ordner
os.makedirs("genomes", exist_ok=True)

overview = []
not_downloaded = []

def fetch_best_assembly(species_name):
    search_handle = Entrez.esearch(db="assembly", term=species_name + "[Organism]", retmax=50)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    ids = search_results['IdList']
    best = None

    for asm_id in ids:
        summary_handle = Entrez.esummary(db="assembly", id=asm_id, report="full")
        summary = Entrez.read(summary_handle)
        summary_handle.close()
        info = summary['DocumentSummarySet']['DocumentSummary'][0]
        quality = info['AssemblyStatus']
        if quality in quality_order:
            if best is None or quality_order.index(quality) < quality_order.index(best['quality']):
                ftp_path = info['FtpPath_RefSeq'] or info['FtpPath_GenBank']
                if ftp_path:
                    best = {
                        "species": species_name,
                        "accession": info['AssemblyAccession'],
                        "ftp_path": ftp_path,
                        "quality": quality
                    }
    return best

for species in tqdm(species_list, desc="Processing species"):
    try:
        best_assembly = fetch_best_assembly(species)
        if not best_assembly:
            not_downloaded.append(species)
            continue

        ftp = best_assembly["ftp_path"]
        file_basename = os.path.basename(ftp)
        fasta_url = f"{ftp}/{file_basename}_genomic.fna.gz"
        folder = f"genomes/{best_assembly['quality'].lower().replace(' ', '_')}"
        os.makedirs(folder, exist_ok=True)

        output_path = f"{folder}/{file_basename}_genomic.fna.gz"
        os.system(f"wget -q -O {output_path} {fasta_url}")
        overview.append(best_assembly)

        time.sleep(0.3)
    except Exception as e:
        not_downloaded.append(species)

# Speichern der Übersicht
df = pd.DataFrame(overview)
df.to_csv("genomes/genome_overview.csv", index=False)

# Speichern der Fehlerliste
with open("genomes/genomes_not_downloaded.txt", "w") as f:
    for name in not_downloaded:
        f.write(name + "\n")
