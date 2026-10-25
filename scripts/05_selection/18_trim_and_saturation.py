#!/usr/bin/env python3
"""Uniform trimming of a codon alignment + saturation check (same rules for every clade, fixed beforehand).
Usage: python3 18_trim_and_saturation.py IN.codon.fasta OUT_PREFIX [MAX_COL_GAP=0.5] [MAX_SEQ_GAP=0.3]
Rules: (1) drop codon columns with a gap in more than MAX_COL_GAP of the sequences;
       (2) drop sequences with more than MAX_SEQ_GAP gap codons in the kept columns;
       (3) report median pairwise p-distance at 3rd codon positions and % of pairs above 0.6 (saturation proxy).
Writes OUT_PREFIX.trim.fasta and prints one summary line."""
import sys, numpy as np
fin, pref = sys.argv[1], sys.argv[2]
mc = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
ms = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
names, seqs, cur = [], [], None
for l in open(fin):
    if l.startswith(">"): names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.strip().upper())
seqs = ["".join(s) for s in seqs]
n, L = len(seqs), len(seqs[0]) // 3
A = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(n, L, 3)
gap = (A == ord("-")).any(axis=2)
keepc = gap.mean(axis=0) <= mc
A, gap = A[:, keepc], gap[:, keepc]
keeps = gap.mean(axis=1) <= ms
A, gap = A[keeps], gap[keeps]
names = [x for x, k in zip(names, keeps) if k]
n2, L2 = A.shape[:2]
# pairwise p-distance at 3rd positions and at all positions
p3 = []
for i in range(n2):
    for j in range(i + 1, n2):
        ok = ~(gap[i] | gap[j])
        if ok.sum() < 100: continue
        p3.append((A[i, ok, 2] != A[j, ok, 2]).mean())
p3 = np.array(p3) if p3 else np.array([np.nan])
with open(pref + ".trim.fasta", "w") as o:
    for k, nm in enumerate(names): o.write(f">{nm}\n{A[k].tobytes().decode()}\n")
flag = "SATURATED" if np.nanmedian(p3) > 0.6 else "ok"
print(f"{pref.split('/')[-1]}\tseqs {n}->{n2}\tcodons {L}->{L2}\tgaps_after {100*gap.mean():.1f}%\tmedian_p3 {np.nanmedian(p3):.3f}\tpairs_p3>0.6 {100*np.nanmean(p3>0.6):.1f}%\t{flag}")
