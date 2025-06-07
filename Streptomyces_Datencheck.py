# sucht nach synonymen

import requests
import time

INPUT_FILE = "00_List_of_streptomyces_basis.txt"
OUTPUT_FILE = "01_List_of_streptomyces_LookedUp.tsv"
EMAIL = "sebasbaier@gmail.com"  # Empfohlen für NCBI API-Nutzung
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def read_species_list(filename):
    with open(filename, encoding="utf-8") as f:
        return sorted(set(line.strip() for line in f if line.strip()))

def search_ncbi_taxonomy(species):
    params = {
        "db": "taxonomy",
        "term": species,
        "retmode": "json",
        "email": EMAIL
    }
    r = requests.get(NCBI_BASE + "esearch.fcgi", params=params)
    r.raise_for_status()
    result = r.json()
    id_list = result.get("esearchresult", {}).get("idlist", [])
    return id_list[0] if id_list else None

def fetch_summary(tax_id):
    params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "json",
        "email": EMAIL
    }
    r = requests.get(NCBI_BASE + "esummary.fcgi", params=params)
    r.raise_for_status()
    data = r.json()
    return data.get("result", {}).get(tax_id, {})

def main():
    species_list = read_species_list(INPUT_FILE)
    results = []

    for species in species_list:
        print(f"Checking: {species}")
        tax_id = search_ncbi_taxonomy(species)
        time.sleep(0.4)  # NCBI empfiehlt max. 3 Requests/Sek.

        if tax_id:
            summary = fetch_summary(tax_id)
            scientific_name = summary.get("scientificname", "")
            synonyms = "; ".join(summary.get("othernames", {}).get("synonym", []))
            results.append((species, "FOUND", tax_id, scientific_name, synonyms))
        else:
            results.append((species, "NOT_FOUND", "", "", ""))

    with open(OUTPUT_FILE, "w", encoding="utf-8") as out:
        out.write("Original_Name\tStatus\tTax_ID\tScientific_Name\tSynonyms\n")
        for row in results:
            out.write("\t".join(row) + "\n")

    print(f"\nFertig! Ergebnisse in: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
