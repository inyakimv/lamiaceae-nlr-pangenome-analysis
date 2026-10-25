#!/usr/bin/env python3
"""Saturation table: 3rd-codon-position p-distance per codon alignment, with Jukes-Cantor correction and threshold sensitivity.
Usage: python3 23_saturation_table.py OUT.tsv ALN1.fasta [ALN2.fasta ...]
Works on codon_aln/*.trim.fasta (full clades) and on hyphy_sub/*.hy.fasta (rarefaction replicates); gaps ('-') are skipped pairwise.
Columns: name, seqs, codons, median_p3, IQR, median JC distance (substitutions/site; 'inf' if p3>=0.75), % pairs p3>0.6,
         flag at thresholds 0.50 / 0.55 / 0.60 / 0.65 / 0.70 (S = saturated by the median rule, . = not)."""
import sys, os, numpy as np
out, files = sys.argv[1], sys.argv[2:]
TH = [0.50, 0.55, 0.60, 0.65, 0.70]
def jc(p): return float("inf") if p >= 0.75 else -0.75 * np.log(1 - p / 0.75)
rows = []
for f in files:
    names, seqs = [], []
    for l in open(f):
        if l.startswith(">"): names.append(l[1:].strip()); seqs.append([])
        else: seqs[-1].append(l.strip().upper())
    seqs = ["".join(s) for s in seqs]
    n, L = len(seqs), len(seqs[0]) // 3
    A = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(n, L, 3)
    gap = (A == ord("-")).any(axis=2); b3 = A[:, :, 2]
    p3 = []
    for i in range(n):
        for j in range(i + 1, n):
            ok = ~(gap[i] | gap[j])
            if ok.sum() >= 100: p3.append((b3[i, ok] != b3[j, ok]).mean())
    p3 = np.array(p3); med = float(np.median(p3))
    name = os.path.basename(f).replace(".trim.fasta", "").replace(".hy.fasta", "")
    rows.append([name, n, L, "%.3f" % med, "%.3f-%.3f" % tuple(np.percentile(p3, [25, 75])), "%.2f" % jc(med),
                 "%.1f" % (100 * np.mean(p3 > 0.6))] + ["S" if med > t else "." for t in TH])
keys = ["alignment", "seqs", "codons", "median_p3", "IQR_p3", "JC_dist_of_median", "pct_pairs_p3>0.6"] + ["sat@%.2f" % t for t in TH]
with open(out, "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print(open(out).read())
