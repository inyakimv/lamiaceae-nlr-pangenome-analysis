#!/usr/bin/env python3
"""Compare the conserved core between the main partition and the re-rooted one.
Usage: compare_core.py MAIN_PREFIX REROOT_PREFIX
  (each prefix has _clade_members.tsv and _clade_conservation.tsv from 01/02)
Core = >=35 species and <=3 mean copies per species."""
import sys, pandas as pd, collections
def load(p):
    mem = collections.defaultdict(set)
    for l in open(p + "_clade_members.tsv"):
        c, leaf = l.rstrip("\n").split("\t"); mem[c].add(leaf)
    cons = pd.read_csv(p + "_clade_conservation.tsv", sep="\t")
    core = cons[(cons.n_species >= 35) & (cons.copies_per_species <= 3)]
    return mem, cons, core
A, cA, coreA = load(sys.argv[1]); B, cB, coreB = load(sys.argv[2])
print(f"main: {len(A)} clades, {len(coreA)} core | re-rooted: {len(B)} clades, {len(coreB)} core")
same = sum(1 for c, s in A.items() if any(s == t for t in B.values()))
print(f"clades with identical membership: {same} of {len(A)}")
print("\nMain core clades and their best match in the re-rooted partition:")
for c in coreA.clade:
    s = A[c]; best = max(B, key=lambda k: len(s & B[k]) / len(s | B[k]))
    j = len(s & B[best]) / len(s | B[best])
    row = cB[cB.clade == best].iloc[0] if (cB.clade == best).any() else None
    stat = f"{int(row.n_species)} species, {row.copies_per_species} copies/species" if row is not None else "no Lamiaceae genes"
    print(f"  {c} ({len(s)} leaves) -> {best} (Jaccard {j:.2f}; {stat}; core: {best in set(coreB.clade)})")
print("\nRe-rooted core clades:", ", ".join(coreB.clade))
unmatched = [c for c, s in A.items() if not any(s == t for t in B.values())]
print("Clades whose membership changed:", ", ".join(unmatched) or "none")
