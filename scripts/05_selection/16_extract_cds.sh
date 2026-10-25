#!/bin/bash
# Extract CDS for one species from the ANNEVO GFF3 + genome, and check that the translation equals the ANNEVO protein.
# Usage: bash 16_extract_cds.sh GCA_xxxxxxxxx.1_Genus_species_genomic [OUTDIR]
set -uo pipefail
P=/data/tfg_users/mollai/project
b=$1; out=${2:-cds}; mkdir -p "$out"
gff=$P/00b_annevo_annotation_final/$b.annevo.gff3
prot=$P/annevo_proteins/$b.annevo_proteins.fasta
g=$P/00_genomes/$b.fna; tmpg=""
if [ ! -f "$g" ]; then
  [ -f "$g.gz" ] || { echo "no genome for $b"; exit 1; }
  tmpg="$out/$b.tmp.fna"; gunzip -c "$g.gz" > "$tmpg"; g=$tmpg
fi
[ -f "$gff" ] || { echo "no gff for $b"; exit 1; }
gffread "$gff" -g "$g" -x "$out/$b.cds.fna" -y "$out/$b.pep.faa" 2>"$out/$b.gffread.log"
[ -n "$tmpg" ] && rm -f "$tmpg" "$tmpg.fai"
python3 - "$out/$b.pep.faa" "$prot" "$b" <<'PY'
import sys
def rd(f):
    d, k = {}, None
    for l in open(f):
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        elif k: d[k].append(l.strip())
    return {a: "".join(v).replace("*", "").rstrip(".") for a, v in d.items()}
mine, ref = rd(sys.argv[1]), rd(sys.argv[2])
common = [k for k in ref if k in mine]
same = sum(mine[k] == ref[k] for k in common)
print(f"{sys.argv[3]}: annevo proteins {len(ref)}, translated CDS {len(mine)}, ids in common {len(common)}, identical {same}")
bad = [k for k in common if mine[k] != ref[k]][:3]
if bad: print("  example mismatches:", bad)
PY
