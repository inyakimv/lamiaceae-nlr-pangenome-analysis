#!/bin/bash
# Prepare a trimmed codon alignment for HyPhy and run FEL, MEME and BUSTED.
# Usage: bash 19_hyphy_run.sh CLADE   (expects codon_aln/CLADE.trim.fasta from 18_trim_and_saturation.py)
# Rules (same for every clade): '|' in names -> '_'; in-frame stop codons -> gaps; GTR+gamma FastTree; p<=0.1 for FEL/MEME.
set -uo pipefail
c=$1; d=codon_aln; o=hyphy_out; mkdir -p "$o"
python3 - "$d/$c.trim.fasta" "$o/$c.hy.fasta" <<'PY'
import sys
STOP = {"TAA", "TAG", "TGA"}
names, seqs = [], []
for l in open(sys.argv[1]):
    if l.startswith(">"): names.append(l[1:].split()[0].replace("|", "_")); seqs.append([])
    else: seqs[-1].append(l.strip().upper())
seqs = ["".join(s) for s in seqs]
n_stop = 0; out = []
for s in seqs:
    cod = [s[i:i+3] for i in range(0, len(s), 3)]
    for k, x in enumerate(cod):
        if x in STOP: cod[k] = "---"; n_stop += 1
    out.append("".join(cod))
with open(sys.argv[2], "w") as o:
    for nm, s in zip(names, out): o.write(f">{nm}\n{s}\n")
print(f"{sys.argv[2]}: {len(out)} seqs, {len(out[0])//3} codons, in-frame stop codons masked: {n_stop}")
PY
ft=$(command -v FastTree || command -v fasttree || command -v FastTreeMP)
[ -n "$ft" ] || { echo "FastTree not found"; exit 1; }
"$ft" -nt -gtr -gamma "$o/$c.hy.fasta" > "$o/$c.nwk" 2>"$o/$c.fasttree.log"
for m in fel meme busted; do
  case $m in busted) extra="";; *) extra="--pvalue 0.1";; esac
  "${HYPHY:-hyphy}" CPU=${CPU:-32} $m --alignment "$o/$c.hy.fasta" --tree "$o/$c.nwk" $extra --output "$o/$c.$m.json" > "$o/$c.$m.log" 2>&1
  if [ -s "$o/$c.$m.json" ]; then echo "$c $m OK ($(du -h "$o/$c.$m.json" | cut -f1))"; else echo "$c $m FAILED: no JSON, see $o/$c.$m.log"; fi
done
