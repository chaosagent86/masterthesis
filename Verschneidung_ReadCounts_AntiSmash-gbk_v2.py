import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# === INPUT ===
# Pfade zu AntiSMASH GenBank-Output und featureCounts-Read-Zähltabelle
gbk_file = "data/final_assembly.gbk"
counts_file = "data/68_read_counts.tsv"
outfile = "output/68_crude_glycerol_pH4_bgc_expression_summary.tsv"
outfile_annotations = "output/68_crude_glycerol_pH4_bgc_gene_annotations.tsv"

# === 1. Einlesen der featureCounts-Datei ===
print("[INFO] Lade Read Count Tabelle...")
counts = pd.read_csv(counts_file, sep='\t', comment='#', index_col=0)
counts.index.name = "GeneID"
print(f"[INFO] Anzahl Gene in Counts-Datei: {len(counts)}")

# === Spaltennamen der featureCounts Tabelle festlegen ===
count_col = "mapped.sorted.bam"
if count_col not in counts.columns:
    raise ValueError(f"Spalte '{count_col}' nicht in der Counts-Datei gefunden. Verfügbare Spalten: {list(counts.columns)}")
print(f"[INFO] Verwende Count-Spalte: {count_col}")

# === 2. Parsen der GenBank-Datei zur Cluster/Gene-Zuordnung ===
print("[INFO] Parse GenBank Datei zur Cluster/Gene-Zuordnung...")
bgc_data = defaultdict(list)
bgc_coords = []  # speichert: (cluster_full, contig, start, end, name)
gene_annotations = {}
total_clusters_detected = 0
all_clusters = []
cluster_types = {}

target_types = ["region", "cand_cluster", "protocluster", "CDS", "misc_feature"]

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
                cname = cluster_name.strip().replace(" ", "_")
                start, end = int(feature.location.start), int(feature.location.end)
                full_id = f"{record.id}_{feature.type}_{cname}_{start}_{end}"
                bgc_coords.append((full_id, record.id, start, end, cname))
                all_clusters.append(full_id)
                cluster_types[full_id] = feature.type
                cluster_count += 1

    print(f"[DEBUG] {record.id} hat {cluster_count} Cluster-Einträge")
    total_clusters_detected += cluster_count

    for full_id, contig, c_start, c_end, cname in bgc_coords:
        for feature in record.features:
            if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                s, e = int(feature.location.start), int(feature.location.end)
                if record.id == contig and s >= c_start and e <= c_end:
                    gid = feature.qualifiers["locus_tag"][0]
                    bgc_data[full_id].append(gid)
                    if "product" in feature.qualifiers:
                        gene_annotations[gid] = feature.qualifiers["product"][0]

print(f"[INFO] Gesamtzahl erkannter Cluster (alle Typen): {total_clusters_detected}")
print(f"[INFO] Erkannte Cluster mit Genen: {len(bgc_data)}")

# Hilfsdict für Cluster-Metadaten
bgc_info = {full: (contig, cname, start, end)
            for full, contig, start, end, cname in bgc_coords}

# === 3. Zusammenfassen der Expression pro BGC ===
summary = []
deduplicated = set(all_clusters)
for full_id in deduplicated:
    contig, cname, start, end = bgc_info[full_id]
    ctype = cluster_types.get(full_id, "unknown")
    genes = bgc_data.get(full_id, [])
    found, missing, total_reads = [], [], 0

    for g in genes:
        if g in counts.index:
            rc = counts.loc[g, count_col]
            (found if rc > 0 else missing).append(g)
            total_reads += rc
        else:
            missing.append(g)

    expr = total_reads / len(genes) if genes else 0
    ratio = len(found) / len(genes) if genes else 0
    status = "empty" if ratio == 0 else ("complete" if ratio == 1 else "partly")

    summary.append({
        "Region": contig,
        "Cluster_type": ctype,
        "Cluster-Name": cname,
        "Start": start,
        "End": end,
        "Genes_total": len(genes),
        "Genes_expressed": len(found),
        "Reads_total": total_reads,
        "Mean_reads_per_gene": round(expr, 2),
        "Completeness": f"{ratio:.1%}",
        "Status": status,
        #"Missing_genes": ", ".join(missing) if missing else "-"
    })

# === 4. Export der Ergebnisse als TSV ===
summary_df = pd.DataFrame(summary)
# cols = ["Region", "Cluster_type", "Cluster-Name", "Start", "End",
#         "Genes_total", "Genes_expressed", "Reads_total", "Mean_reads_per_gene", "Completeness", "Status"]
# summary_df[cols].to_csv(outfile, sep='\t', index=False)
summary_df.to_csv(outfile, sep='\t', index=False)
print(f"[INFO] Analyse abgeschlossen. Ergebnisse in: {outfile}")

# === 5. Export der Gen-Annotationen ===
annotations_df = pd.DataFrame(sorted(gene_annotations.items()), columns=["GeneID", "Product"])
annotations_df.to_csv(outfile_annotations, sep='\t', index=False)
print(f"[INFO] Annotationen gespeichert in: {outfile_annotations}")

# Debug: alle Feature-Typen
types = set()
for rec in SeqIO.parse(gbk_file, "genbank"):
    types.update(f.type for f in rec.features)
print(sorted(types))
