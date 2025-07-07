import os
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import logging

# Konfiguration
Entrez.email = "sebasbaier@gmail.com"  # <- Bitte setzen
MAX_GENOMES = 50  # Für Testzwecke z.B. auf 20 setzen, für vollständig: None

# Logging
logging.basicConfig(filename="failed_assemblies.log", level=logging.WARNING,
                    format="%(asctime)s - %(levelname)s - %(message)s")

results = []

def search_streptomyces_genomes():
    query = (
        'Streptomyces[Organism] AND latest[Property] AND '
        '("complete genome"[Assembly Level] OR "chromosome"[Assembly Level] OR "scaffold"[Assembly Level])'
    )
    handle = Entrez.esearch(db="assembly", term=query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def get_assembly_metadata(asm_id):
    try:
        handle = Entrez.esummary(db="assembly", id=asm_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        doc = summary['DocumentSummarySet']['DocumentSummary'][0]
        return {
            "Species": doc.get("SpeciesName"),
            "Accession": doc.get("AssemblyAccession"),
            "Status": doc.get("AssemblyStatus"),
            "TaxID": doc.get("Taxid"),
            "FTP_Link": doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank")
        }
    except Exception as e:
        logging.warning(f"Fehler beim Holen von Metadaten für {asm_id}: {e}")
        return None

def get_contig_count(ftp_url):
    if not ftp_url:
        return None
    try:
        file_base = os.path.basename(ftp_url)
        url = f"{ftp_url}/{file_base}_assembly_report.txt"
        local = "temp_assembly_report.txt"
        os.system(f"wget -q -O {local} {url}")
        contigs = 0
        with open(local) as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    contigs += 1
        os.remove(local)
        return contigs if contigs > 0 else None
    except Exception as e:
        logging.warning(f"Fehler beim Lesen von assembly_report für {ftp_url}: {e}")
        return None

# Hauptlogik
ids = search_streptomyces_genomes()

for i, asm_id in enumerate(tqdm(ids, desc="Verarbeite Assemblies")):
    if MAX_GENOMES and i >= MAX_GENOMES:
        break
    meta = get_assembly_metadata(asm_id)
    if not meta or not meta.get("FTP_Link"):
        logging.warning(f"Unvollständige Metadaten bei ID {asm_id}: {meta}")
        continue
    contigs = get_contig_count(meta["FTP_Link"])
    if contigs is None:
        logging.warning(f"Keine Contigs gezählt bei Assembly {meta['Accession']}")
        continue
    results.append({
        "Species": meta["Species"],
        "Contigs": contigs,
        "Accession": meta["Accession"],
        "Status": meta["Status"],
        "TaxID": meta["TaxID"]
    })

# DataFrame erstellen und speichern
df = pd.DataFrame(results)
if "Species" in df.columns:
    df = df.sort_values("Species")

df.to_csv("streptomyces_genomes.csv", index=False)
print(f"✅ {len(df)} Genome gespeichert in 'streptomyces_genomes.csv'")
