import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# === INPUT ===
# Pfade zu AntiSMASH GenBank-Output und featureCounts-Read-Zähltabelle
gbk_file = "data/final_assembly.gbk"
counts_file = "data/38_read_counts.tsv"
outfile = "bgc_expression_summary.tsv"
outfile_annotations = "bgc_gene_annotations.tsv"

# === 1. Einlesen der featureCounts-Datei ===
print("[INFO] Lade Read Count Tabelle...")
counts = pd.read_csv(counts_file, sep='\t', comment='#', index_col=0)
counts.index.name = "GeneID"
print(f"[INFO] Anzahl Gene in Counts-Datei: {len(counts)}")

# === Spaltennamen automatisch ermitteln ===
count_columns = []
for col in counts.columns:
    if col.lower() not in ["chr", "start", "end", "strand", "length"]:
        count_columns.append(col)
if not count_columns:
    raise ValueError("Keine gültige Count-Spalte gefunden. Bitte überprüfe die featureCounts-Datei.")
count_col = count_columns[0]
print(f"[INFO] Verwende Count-Spalte: {count_col}")

# === 2. Parsen der GenBank-Datei zur Cluster/Gene-Zuordnung ===
print("[INFO] Parse GenBank Datei zur Cluster/Gene-Zuordnung...")
bgc_data = defaultdict(list)
bgc_coords = []
gene_annotations = {}
total_clusters_detected = 0
all_clusters = []
cluster_types = {}

# Typen, die als Cluster zählen sollen
target_types = ["region", "cand_cluster", "misc_feature"]

for record in SeqIO.parse(gbk_file, "genbank"):
    cluster_count = 0
    for feature in record.features:
        if feature.type in target_types:
            cluster_name = None
            if "product" in feature.qualifiers:
                cluster_name = feature.qualifiers["product"][0]
            elif "note" in feature.qualifiers:
                cluster_name = feature.qualifiers["note"][0]
            elif "region_number" in feature.qualifiers:
                cluster_name = f"Cluster_{feature.qualifiers['region_number'][0]}"

            if cluster_name:
                cluster_name = cluster_name.strip().replace(" ", "_")
                cluster_name_full = f"{record.id}_{feature.type}_{cluster_name}_{int(feature.location.start)}_{int(feature.location.end)}"
                start, end = int(feature.location.start), int(feature.location.end)
                bgc_coords.append((cluster_name_full, record.id, start, end))
                all_clusters.append(cluster_name_full)
                cluster_types[cluster_name_full] = feature.type
                cluster_count += 1

    print(f"[DEBUG] {record.id} hat {cluster_count} Cluster-Einträge")
    total_clusters_detected += cluster_count

    for cluster_name, contig, c_start, c_end in bgc_coords:
        for feature in record.features:
            if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                start, end = int(feature.location.start), int(feature.location.end)
                if record.id == contig and start >= c_start and end <= c_end:
                    gene_id = feature.qualifiers["locus_tag"][0]
                    bgc_data[cluster_name].append(gene_id)
                    if "product" in feature.qualifiers:
                        gene_annotations[gene_id] = feature.qualifiers["product"][0]

print(f"[INFO] Gesamtzahl erkannter Cluster (alle Typen): {total_clusters_detected}")
print(f"[INFO] Erkannte Cluster mit Genen: {len(bgc_data)}")

# === 3. Zusammenfassen der Expression pro BGC ===
summary = []

deduplicated_clusters = set(all_clusters)

for cluster in deduplicated_clusters:
    genes = bgc_data[cluster] if cluster in bgc_data else []
    found = []
    missing = []
    total_reads = 0

    for g in genes:
        if g in counts.index:
            read_count = counts.loc[g, count_col]
            if read_count > 0:
                found.append(g)
            else:
                missing.append(g)
            total_reads += read_count
        else:
            missing.append(g)

    expression_per_gene = total_reads / len(genes) if genes else 0
    completeness_ratio = len(found) / len(genes) if genes else 0
    completeness_label = "empty" if completeness_ratio == 0 else ("complete" if completeness_ratio == 1 else "partly")

    print(f"[DEBUG] {cluster}: {len(genes)} Gene, {len(found)} exprimiert, {len(missing)} nicht")

    summary.append({
        "Cluster": cluster,
        "Cluster_type": cluster_types.get(cluster, "unknown"),
        "Genes_total": len(genes),
        "Genes_expressed": len(found),
        "Reads_total": total_reads,
        "Mean_reads_per_gene": round(expression_per_gene, 2),
        "Completeness": f"{completeness_ratio:.1%}",
        "Status": completeness_label,
        "Missing_genes": ", ".join(missing) if missing else "-"
    })

# === 4. Export der Ergebnisse als TSV ===
summary_df = pd.DataFrame(summary)
summary_df.to_csv(outfile, sep='\t', index=False)
print(f"[INFO] Analyse abgeschlossen. Ergebnisse in: {outfile}")

# === 5. Export der Gen-Annotationen ===
annotations_df = pd.DataFrame(sorted(gene_annotations.items()), columns=["GeneID", "Product"])
annotations_df.to_csv(outfile_annotations, sep='\t', index=False)
print(f"[INFO] Annotationen gespeichert in: {outfile_annotations}")
