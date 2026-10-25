#!/bin/bash
# Rarefaction control for the 45-species dataset: ONE replicate of one clade. Same rules as 21_subsample_hyphy.sh, but sequences of the
# excluded species are removed from the source alignment BEFORE subsampling (all-gap codons are dropped afterwards, as in script 21).
# Usage: SRC=hyphy_out OUT=hyphy_sub_45 EXCL=Paraleonurus HYPHY=$(which hyphy) CPU=10 bash 33_subsample_hyphy_excl.sh CLADE N REP
# Needs $SRC/CLADE.hy.fasta (from 19_hyphy_run.sh). Seed = REP. Output: $OUT/CLADE_nN_rK.{hy.fasta,nwk,fel.json,meme.json,busted.json}
set -uo pipefail
c=$1; n=$2; r=$3; src=${SRC:-hyphy_out}; o=${OUT:-hyphy_sub_45}; mkdir -p "$o"
ft=$(command -v FastTree || command -v fasttree || command -v FastTreeMP)
[ -n "$ft" ] || { echo "FastTree not found"; exit 1; }
id="${c}_n${n}_r${r}"
python3 - "$src/$c.hy.fasta" "$o/$id.hy.fasta" "$n" "$r" "${EXCL:-Paraleonurus}" <<'PY'
import sys, random
names, seqs = [], []
for l in open(sys.argv[1]):
    if l.startswith(">"): names.append(l[1:].strip()); seqs.append([])
    else: seqs[-1].append(l.strip())
seqs = ["".join(s) for s in seqs]
excl = [e for e in sys.argv[5].split(",") if e]
keep_i = [i for i, nm in enumerate(names) if not any(e in nm for e in excl)]
n, seed = int(sys.argv[3]), int(sys.argv[4])
pick = sorted(random.Random(seed).sample(keep_i, min(n, len(keep_i))))
sub = [seqs[i] for i in pick]
L = len(sub[0]); keep = [k for k in range(0, L, 3) if any(s[k:k+3] != "---" for s in sub)]
with open(sys.argv[2], "w") as out:
    for i, s in zip(pick, sub): out.write(">%s\n%s\n" % (names[i], "".join(s[k:k+3] for k in keep)))
print("%s: %d of %d seqs (%d removed as excluded), %d of %d codons kept" % (sys.argv[2], len(sub), len(keep_i), len(names) - len(keep_i), len(keep), L // 3))
PY
"$ft" -nt -gtr -gamma "$o/$id.hy.fasta" > "$o/$id.nwk" 2>"$o/$id.fasttree.log"
for m in fel meme busted; do
  case $m in busted) extra="";; *) extra="--pvalue 0.1";; esac
  "${HYPHY:-hyphy}" CPU=${CPU:-10} $m --alignment "$o/$id.hy.fasta" --tree "$o/$id.nwk" $extra --output "$o/$id.$m.json" > "$o/$id.$m.log" 2>&1
  if [ -s "$o/$id.$m.json" ]; then echo "$id $m OK ($(du -h "$o/$id.$m.json" | cut -f1))"; else echo "$id $m FAILED: no JSON, see $o/$id.$m.log"; fi
done
