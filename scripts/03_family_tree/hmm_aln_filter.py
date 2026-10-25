import sys
a2m = "nbarc_hmmalign.a2m"
M = 248
MIN_COV = 0.70

names, seqs, cur = [], [], []
for line in open(a2m):
    line = line.rstrip()
    if line.startswith(">"):
        if cur: seqs.append("".join(cur)); cur = []
        names.append(line[1:].split()[0])
    else:
        cur.append(line)
seqs.append("".join(cur))

# Keep only match columns: uppercase letters and '-' (drop lowercase inserts and '.')
match = ["".join(c for c in s if c == "-" or c.isupper()) for s in seqs]
bad = [n for n, s in zip(names, match) if len(s) != M]
print(f"Sequences: {len(names)}; with length != {M}: {len(bad)}")
if bad:
    print("Example of unexpected length:", bad[0], len(match[names.index(bad[0])]))
    sys.exit(1)

cov = [sum(1 for c in s if c != "-") / M for s in match]
srt = sorted(cov)
print("Coverage quantiles (5, 25, 50, 75, 95%):",
      [round(srt[int(len(srt) * q)], 2) for q in (0.05, 0.25, 0.5, 0.75, 0.95)])
for thr in (0.5, 0.6, 0.7, 0.8, 0.9):
    print(f"  coverage >= {thr}: {sum(c >= thr for c in cov)} sequences")

keep = [i for i, c in enumerate(cov) if c >= MIN_COV]
print(f"Kept with coverage >= {MIN_COV}: {len(keep)}; removed: {len(names) - len(keep)}")
removed_alga = [names[i] for i in range(len(names)) if i not in set(keep) and "GAQ" in names[i]]
print("Algal sequences removed:", removed_alga)

with open("nbarc_hmmaligned_cov70.fasta", "w") as out:
    for i in keep:
        out.write(f">{names[i]}\n{match[i]}\n")

# Alignment statistics
import numpy as np
A = np.frombuffer("".join(match[i] for i in keep).encode(), dtype="S1").reshape(len(keep), M)
gap = A == b"-"
gf = gap.mean(axis=0)
aas = [c.encode() for c in "ACDEFGHIKLMNPQRSTVWY"]
cnt = np.stack([(A == a).sum(axis=0) for a in aas])
inf = (cnt >= 2).sum(axis=0) >= 2
print(f"Alignment length: {M}; mean residues/sequence: {(~gap).sum(axis=1).mean():.1f}")
for t in (0.2, 0.5, 0.8):
    print(f"Columns with <{int(t*100)}% gaps: {(gf < t).sum()}")
print(f"Informative columns (all): {inf.sum()}; informative with <50% gaps: {(inf & (gf < 0.5)).sum()}")
