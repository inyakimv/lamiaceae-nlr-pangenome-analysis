#!/bin/bash
# Per-clade codon alignment for selection analysis.
# Usage: bash 17_codon_alignment.sh GENE_LIST.txt CDS_DIR OUTDIR [THREADS]
#   GENE_LIST.txt: one id per line, form Species|contig-gN.t1   (clade_gene_lists/*.txt)
#   CDS_DIR: folder with <GCx_..._Species_genomic>.cds.fna and .pep.faa (from 16_extract_cds.sh)
# Steps: extract CDS + protein of the listed genes -> MAFFT protein alignment -> pal2nal codon alignment.
set -uo pipefail
list=$1; cds=$2; out=$3; thr=${4:-8}
name=$(basename "$list" .txt); mkdir -p "$out"
python3 - "$list" "$cds" "$out/$name" <<'PY'
import sys, glob, os, re
lst, cdsdir, pref = sys.argv[1:4]
want = [l.strip() for l in open(lst) if l.strip()]
byfile = {}
for f in glob.glob(os.path.join(cdsdir, "*.cds.fna")):
    sp = re.sub(r"^GC[AF]_\d+\.\d+_", "", os.path.basename(f)[:-8]).replace("_genomic", "")
    byfile[sp] = f[:-8]
def rd(f):
    d, k = {}, None
    for l in open(f):
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        elif k: d[k].append(l.strip())
    return {a: "".join(v) for a, v in d.items()}
cds, pep, miss = {}, {}, []
for sp in sorted({w.split("|")[0] for w in want}):
    if sp not in byfile: miss += [w for w in want if w.startswith(sp + "|")]; continue
    c, p = rd(byfile[sp] + ".cds.fna"), rd(byfile[sp] + ".pep.faa")
    for w in want:
        s, g = w.split("|")
        if s != sp: continue
        if g in c and g in p: cds[w], pep[w] = c[g], p[g].rstrip("*").rstrip(".")
        else: miss.append(w)
bad = [w for w in cds if len(cds[w]) not in (3 * len(pep[w]), 3 * len(pep[w]) + 3)]
for w in bad: cds.pop(w); pep.pop(w)
open(pref + ".cds.fna", "w").write("".join(f">{w}\n{cds[w]}\n" for w in cds))
open(pref + ".pep.faa", "w").write("".join(f">{w}\n{pep[w]}\n" for w in cds))
print(f"{os.path.basename(pref)}: requested {len(want)}, written {len(cds)}, missing {len(miss)}, length-mismatch {len(bad)}")
if miss: print("  missing e.g.", miss[:3])
PY
mafft --auto --thread "$thr" "$out/$name.pep.faa" > "$out/$name.pep.aln" 2>"$out/$name.mafft.log"
if command -v pal2nal.pl >/dev/null; then
  pal2nal.pl "$out/$name.pep.aln" "$out/$name.cds.fna" -output fasta -codontable 1 > "$out/$name.codon.fasta" 2>"$out/$name.pal2nal.log"
  echo "$name: codon alignment $(grep -c '>' "$out/$name.codon.fasta") sequences, $(grep -v '>' "$out/$name.codon.fasta" | head -1 | tr -d '\n' | wc -c) columns"
else
  echo "pal2nal.pl not found: install it (conda install -c bioconda pal2nal) and rerun"
fi
