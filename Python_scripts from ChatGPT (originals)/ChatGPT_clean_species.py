#Entfernt leerzeilen
#Entfernt duplikate
#sortiert alphabetisch

from pathlib import Path

# Ursprungsdatei
input_file = "streptomyces_species.txt"

# Zieldatei mit Suffix "_curated.txt"
output_file = Path(input_file).with_name(Path(input_file).stem + "_curated.txt")

# Zeilen lesen, filtern, deduplizieren, sortieren
with open(input_file, "r", encoding="utf-8") as infile:
    lines = {line.strip() for line in infile if line.strip()}  # Set für Duplikate

# Alphabetisch sortieren
sorted_lines = sorted(lines, key=lambda s: s.lower())

# Gesäuberte Datei schreiben
with open(output_file, "w", encoding="utf-8") as outfile:
    outfile.write("\n".join(sorted_lines) + "\n")

print(f"✅ Gesäubert, sortiert und gespeichert unter: {output_file}")
