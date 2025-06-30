import pandas as pd

# 1. Einlesen der CSV-Datei und Anzeige der ersten Zeilen
df = pd.read_csv('data/All_RNASeqs_Together.csv', sep=',')  # DataFrame mit allen Rohdaten
print("[INFO] CSV eingelesen: All_RNASeqs_Together.csv")
print("[DEBUG] Original Table Head:")
print(df.head(), "\n")  # Zeigt die Struktur und erste Zeilen der Daten an

# 2. Sortierung der Daten für konsistente Reihenfolge
print("[INFO] Sortiere Daten nach Region, Cluster_type, Cluster-Name, Start, End...")
df_sorted = df.sort_values(
    by=['Region', 'Cluster_type', 'Cluster-Name', 'Start', 'End'],
    ascending=[True, True, True, True, True]
)
print(f"[DEBUG] Erste Zeilen nach Sortierung:\n{df_sorted.head()}\n")

# 3. Pivotierung: Dreht die Tabelle um, sodass pro Kombi aus Region/Cluster nur noch eine Zeile existiert
print("[INFO] Pivot-Tabelle erstellen...")
df_pivot = (
    df_sorted
    .pivot_table(
        index=['Region', 'Cluster_type', 'Cluster-Name', 'Start', 'End', 'Genes_total'],
        columns='Sample',
        values=['Genes_expressed', 'Reads_total', 'Completeness'],
        aggfunc='first'  # Da jede Kombination nur einmal pro Sample vorkommt, reicht 'first'
    )
    .reset_index()  # Setzt den Index zurück, damit Region/Cluster-Spalten wieder normale Spalten sind
)
print(f"[DEBUG] Pivot-Tabelle Head:\n{df_pivot.head()}\n")

# 4. Umbenennung der mehrstufigen Spaltennamen in flache Strings
print("[INFO] Mehrstufige Spalten umbenennen...")
new_cols = []
for col_tuple in df_pivot.columns.to_flat_index():
    # col_tuple ist ein Tupel wie ('Genes_expressed', '38 Pure Glycerol pH7')
    if isinstance(col_tuple, tuple) and col_tuple[1] not in [None, 'Region','Cluster_type','Cluster-Name','Start','End','Genes_total']:
        new_cols.append(f"{col_tuple[0]}_{col_tuple[1]}")
    else:
        new_cols.append(col_tuple[0] if isinstance(col_tuple, tuple) else col_tuple)
# Weist die neuen Spaltennamen zu
df_pivot.columns = new_cols
print(f"[DEBUG] Neue Spaltennamen:\n{df_pivot.columns.tolist()}\n")

# 5. Speichern der vollständigen, umstrukturierten Tabelle mit Komma als Dezimaltrennzeichen
df_pivot.to_csv('output/RNASeq_reshaped_output_full.csv', index=False, decimal=',')
print("[INFO] Gespeichert: reshaped_output_full.csv (gesamt, Dezimaltrennzeichen ',')")

# 6. Identifikation der Completeness-Spalten für alle Samples
completeness_cols = [col for col in df_pivot.columns if col.startswith('Completeness_')]
print(f"[DEBUG] Gefundene Completeness-Spalten: {completeness_cols}")

# 7. Maske für identische Completeness-Werte über alle Samples
print("[INFO] Erstelle Maske für identische Completeness-Werte...")
mask_identical = df_pivot[completeness_cols].nunique(axis=1) == 1
print(f"[DEBUG] Anzahl identischer Zeilen: {mask_identical.sum()}")

# 8. Erstellen und Speichern der Tabelle mit nur unterschiedlichen Completeness-Werten
print("[INFO] Speichere nur Unterschiede...")
df_only_diff = df_pivot[~mask_identical]
df_only_diff.to_csv('output/RNASeq_reshaped_output_only_differences.csv', index=False, decimal=',')
print(f"[INFO] Gespeichert: only_differences.csv (nur Unterschiede, Dezimaltrennzeichen ',')")
print(df_only_diff.head(), "\n")

# 9. Erstellen und Speichern der Tabelle mit nur identischen Completeness-Werten
print("[INFO] Speichere nur identische Werte...")
df_ident = df_pivot[mask_identical]
df_ident.to_csv('output/RNASeq_reshaped_output_only_ident.csv', index=False, decimal=',')
print(f"[INFO] Gespeichert: ident.csv (nur identisch, Dezimaltrennzeichen ',')")
print(df_ident.head())
