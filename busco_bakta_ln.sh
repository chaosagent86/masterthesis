#!/usr/bin/bash
set -e
#https://ioflood.com/blog/bash-exit-on-error/#:~:text=The%20simplest%20way%20to%20make,beginning%20of%20your%20bash%20script.

# 1) In Busco-Umgebung wechseln
eval "$(conda shell.bash hook)"
conda activate busco_env

# 2) BUSCO-Analyse für jedes Genom
for fna in ./00_genomes/*.fna; do
  sample=$(basename "$fna" .fna)
  outdir="busco_${sample}"
  busco \
    -i "$fna" \
    -l streptomyces_odb12 \
    -m genome \
    --cpu 60 \
    --out "$outdir" 
done

# 3) In Bakta-Umgebung wechseln
eval "$(conda shell.bash hook)"
conda activate bakta_env

# 4) BAKTA-Annotation für jedes Genom
for fna in ./00_genomes/*.fna; do
  sample=$(basename "$fna" .fna)
  outdir="bakta_${sample}"
  bakta \
    --db ../baktadb/db/ \
    --output "$outdir" \
    "$fna"
done
