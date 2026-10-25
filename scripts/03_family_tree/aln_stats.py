import sys
import numpy as np

fa = "nbarc_domains_aligned_plus_outgroup.fasta"
names, seqs, cur = [], [], []
for line in open(fa):
    line = line.rstrip()
    if line.startswith(">"):
        if cur: seqs.append("".join(cur)); cur = []
        names.append(line[1:])
    else:
        cur.append(line)
seqs.append("".join(cur))

n = len(seqs); L = len(seqs[0])
assert all(len(s) == L for s in seqs), "sequences of different length"
A = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(n, L)
gap = (A == b"-") | (A == b"X") | (A == b"?")
print(f"Sequences: {n}")
print(f"Alignment length (columns): {L}")

gapfrac = gap.mean(axis=0)
for thr in (0.2, 0.5, 0.8):
    print(f"Columns with <{int(thr*100)}% gaps: {(gapfrac < thr).sum()}")
print(f"Columns with no gaps at all: {(gapfrac == 0).sum()}")
print(f"Mean residues per sequence (non-gap): {(~gap).sum(axis=1).mean():.1f}")

aas = [c.encode() for c in "ACDEFGHIKLMNPQRSTVWY"]
cnt = np.stack([(A == a).sum(axis=0) for a in aas])   # 20 x L
n_with2 = (cnt >= 2).sum(axis=0)
informative = n_with2 >= 2
variable = (cnt > 0).sum(axis=0) >= 2
conserved = ((cnt > 0).sum(axis=0) == 1) & (gapfrac < 0.5)

print(f"Parsimony-informative columns: {informative.sum()}")
print(f"Informative AND <50% gaps: {(informative & (gapfrac < 0.5)).sum()}")
print(f"Variable columns: {variable.sum()}")
print(f"Fully conserved columns (<50% gaps): {conserved.sum()}")
