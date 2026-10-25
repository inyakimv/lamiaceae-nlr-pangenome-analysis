#!/usr/bin/env python3
"""Deterministic choice of representative genes per clade for structure prediction.
Farthest-point sampling on patristic distance, at most one gene per species, fixed start
(the candidate with the smallest id), so the choice does not depend on judgement.
Usage: python3 12_select_reps.py TREE CANDIDATES.txt N OUT.tsv LABEL
 TREE: newick (leaf names must match the candidate ids after normalisation)
 CANDIDATES.txt: one leaf id per line."""
import sys, re
from Bio import Phylo
sys.setrecursionlimit(100000)
tree, cand, n, out, label = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
norm = lambda x: re.sub(r"^GC[AF]_\d+\.\d+_", "", x.strip()).replace("_genomic", "")
t = Phylo.read(tree, "newick")
leaf = {norm(l.name): l for l in t.get_terminals()}
want = [norm(l) for l in open(cand) if l.strip()]
C = sorted(w for w in want if w in leaf)
species = lambda x: x.split("|")[0]
chosen = [C[0]]
while len(chosen) < n:
    best = None
    for c in C:
        if species(c) in {species(x) for x in chosen}: continue
        d = min(t.distance(leaf[c], leaf[x]) for x in chosen)
        if best is None or d > best[0]: best = (d, c)
    if best is None: break
    chosen.append(best[1])
with open(out, "a") as f:
    for c in chosen: f.write(f"{label}\t{c}\n")
print(label, len(C), "candidates ->", chosen)
