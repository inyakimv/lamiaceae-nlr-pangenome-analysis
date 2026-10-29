#!/bin/bash
# Keep the longest isoform per gene and drop short peptides, for every proteome in IN_DIR.
# Usage: bash clean_proteomes.sh [IN_DIR] [OUT_DIR] [MIN_LEN]
set -euo pipefail
IN="${1:-$HOME/project/proteomes_45}"; OUT="${2:-$HOME/project/proteomes_45_clean}"; MINLEN="${3:-50}"
mkdir -p "$OUT"
printf "species\tinput_proteins\tgenes\tkept_longest_per_gene\tshorter_than_%s_removed\tfinal\n" "$MINLEN" > "$OUT/summary_clean.tsv"
for f in "$IN"/*.faa; do
  sp=$(basename "$f" .faa)
  awk -v min="$MINLEN" -v sp="$sp" -v OFS='\t' '
    /^>/ { if (id!="") store(); id=substr($0,2); seq=""; next }
    { seq=seq $0 }
    END { if (id!="") store()
          for (g in best) { nk++; if (length(bestseq[g])<min) { nshort++; continue }
                            print ">" bestid[g]; print bestseq[g]; nf++ }
          print sp, nin, nk, nk, nshort+0, nf+0 > "/dev/stderr" }
    function store(   g) { nin++; g=id; sub(/\.t[0-9]+$/,"",g)
        if (!(g in best) || length(seq)>length(bestseq[g])) { best[g]=1; bestid[g]=id; bestseq[g]=seq } }
  ' "$f" > "$OUT/$sp.faa" 2>> "$OUT/summary_clean.tsv"
done
chmod -R g+rX,o+rX "$OUT"
cat "$OUT/summary_clean.tsv"
