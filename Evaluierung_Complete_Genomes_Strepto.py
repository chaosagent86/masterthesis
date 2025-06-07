# Hier ist das überarbeitete Python-Skript, das für alle Streptomyces-Genome mit dem Status "Complete Genome" die folgenden Felder extrahiert:
# Species – wissenschaftlicher Name
# Accession – Assembly-Accession
# AssemblyName – interner Assembly-Name
# Status – Assembly-Status
# ReleaseDate – Veröffentlichungsdatum
# SeqLength – Genomlänge in Basen
# Submitter – Einreichende Institution
# Contigs – Anzahl Contigs (aus assembly_report.txt)
# TaxID – NCBI Taxonomie-ID
# FTP_Link – Downloadlink

import os
import pandas as pd
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "sebasbaier@gmail.com"

results = []

def search_complete_streptomyces():
    query = 'Streptomyces[Organism] AND "complete genome"[Assembly Level]'
    handle = Entrez.esearch(db="assembly", term=query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def get_assembly_metadata(asm_id):
    handle = Entrez.esummary(db="assembly", id=asm_id, report="full")
    summary = Entrez.read(handle)
    handle.close()
    doc = summary['DocumentSummarySet']['DocumentSummary'][0]
    return {
        "Species": doc["SpeciesName"],
        "Accession": doc["AssemblyAccession"],
        "AssemblyName": doc["AssemblyName"],
        "Status": doc["AssemblyStatus"],
        "ReleaseDate": doc["ReleaseDate"],
        "SeqLength": doc["SeqLength"],
        "Submitter": doc["SubmitterOrganization"],
        "TaxID": doc["Taxid"],
        "FTP_Link": doc["FtpPath_RefSeq"] or doc["FtpPath_GenBank"]
    }

def get_contig_count(ftp_url):
    try:
        file_base = os.path.basename(ftp_url)
        url = f"{ftp_url}/{file_base}_assembly_report.txt"
        local = "temp_assembly_report.txt"
        os.system(f"wget -q -O {local} {url}")
        contigs = 0
        with open(local) as f:
            for line in f:
                if not line.startswith("#"):
                    contigs += 1
        os.remove(local)
        return contigs if contigs > 0 else None
    except:
        return None

ids = search_complete_streptomyces()

for asm_id in tqdm(ids, desc="Verarbeite Assemblies"):
    try:
        meta = get_assembly_metadata(asm_id)
        contigs = get_contig_count(meta["FTP_Link"])
        meta["Contigs"] = contigs
        results.append(meta)
    except Exception as e:
        continue

df = pd.DataFrame(results)
df.to_csv("streptomyces_complete_genomes.csv", index=False)
