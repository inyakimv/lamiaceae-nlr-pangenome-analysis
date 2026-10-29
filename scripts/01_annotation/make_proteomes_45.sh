#!/bin/bash
# Build the 45 clean Lamiaceae proteomes (ANNEVO) for functional annotation (e.g. FANTASIA).
# Usage: bash make_proteomes_45.sh  [IN_DIR]  [OUT_DIR]
#   IN_DIR  : folder with the AGAT outputs  *.annevo_proteins.fasta  (default ~/project/annevo_proteins)
#   OUT_DIR : output folder                                     (default ~/project/proteomes_45)
set -euo pipefail
IN_DIR="${1:-$HOME/project/annevo_proteins}"
OUT_DIR="${2:-$HOME/project/proteomes_45}"
EXCLUDE="Paraleonurus"          # excluded from all analyses (annotation problem)
mkdir -p "$OUT_DIR"; : > "$OUT_DIR/summary.tsv"
printf "file\tspecies\tproteins\tinternal_stops\tmin_len\tmax_len\n" > "$OUT_DIR/summary.tsv"
n=0
for f in "$IN_DIR"/*.annevo_proteins.fasta; do
  base=$(basename "$f" .annevo_proteins.fasta)
  [[ "$base" == *"$EXCLUDE"* ]] && { echo "skip $base"; continue; }
  # species name = genus_species taken from the file name (adjust the sed if your names differ)
  sp=$(echo "$base" | sed -E 's/^GC[AF]_[0-9.]+_//; s/^([A-Za-z]+_(x_)?[a-z]+).*/\1/')
  out="$OUT_DIR/${sp}.faa"
  awk -v sp="$sp" '
    /^>/ { if (seq!="") flush(); id=substr($1,2); sub(/[ \t].*/,"",id); hdr=">" sp "|" id; seq=""; next }
    { gsub(/[ \t\r]/,""); seq=seq $0 }
    END { if (seq!="") flush() }
    function flush(   s,internal) {
      s=seq; sub(/\*$/,"",s); internal=(index(s,"*")>0)
      gsub(/\*/,"X",s)                                   # internal stops -> X (counted in summary)
      print hdr; print s
      if (internal) print id > "/dev/stderr"
    }' "$f" > "$out" 2> "$OUT_DIR/${sp}.internal_stops.txt"
  np=$(grep -c '^>' "$out"); ni=$(wc -l < "$OUT_DIR/${sp}.internal_stops.txt")
  read mn mx < <(awk '!/^>/{l=length($0); if(min==""||l<min)min=l; if(l>max)max=l} END{print min,max}' "$out")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$base" "$sp" "$np" "$ni" "$mn" "$mx" >> "$OUT_DIR/summary.tsv"
  n=$((n+1))
done
echo "species written: $n (expected 45)"
# duplicate-ID check across species
cat "$OUT_DIR"/*.faa | grep '^>' | sort | uniq -d | head -3
chmod -R g+rX,o+rX "$OUT_DIR"
cat "$OUT_DIR/summary.tsv"
