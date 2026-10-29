#!/bin/bash
# Download the official NCBI protein sets (protein.faa) for the 45 Lamiaceae assemblies.
# Needs: ncbi-datasets-cli  (conda install -c conda-forge ncbi-datasets-cli)  and internet access from the server.
# Usage: bash download_ncbi_proteomes.sh accessions_45.tsv [OUT_DIR]
set -uo pipefail
ACC="$1"; OUT="${2:-$HOME/project/proteomes_ncbi}"
mkdir -p "$OUT/raw"; cd "$OUT"
printf "accession\tspecies\tstatus\tproteins\n" > status.tsv
while IFS=$'\t' read -r acc sp; do
  tmp="raw/$acc"; mkdir -p "$tmp"
  datasets download genome accession "$acc" --include protein --filename "$tmp/$acc.zip" > "$tmp/log.txt" 2>&1
  if [ -s "$tmp/$acc.zip" ] && unzip -qo "$tmp/$acc.zip" -d "$tmp" 2>/dev/null && ls "$tmp"/ncbi_dataset/data/*/protein.faa >/dev/null 2>&1; then
    awk -v sp="$sp" '/^>/{id=substr($1,2); print ">" sp "|" id; next} {sub(/\*$/,""); gsub(/\*/,"X"); print}' "$tmp"/ncbi_dataset/data/*/protein.faa > "$sp.faa"
    n=$(grep -c '^>' "$sp.faa"); printf "%s\t%s\tOK\t%s\n" "$acc" "$sp" "$n" >> status.tsv
  else
    printf "%s\t%s\tNO_PROTEIN_SET\t0\n" "$acc" "$sp" >> status.tsv
  fi
done < "$ACC"
echo "OK: $(grep -c 'OK' status.tsv)   without protein set: $(grep -c NO_PROTEIN status.tsv)"
cat status.tsv
