#!/usr/bin/env python3
"""Summarise the Pfam-A hmmscan of the proteins added by the relaxed criterion.
Usage: python3 40_summarise_added_domains.py ADDED.domtblout ADDED.faa [EVALUE=1e-3]
Reports: how many have an LRR-family domain that the main profiles (LRR_1/3/4/8) missed, which other domains they carry, and how many have nothing but NB-ARC."""
import sys, collections
dt, fa = sys.argv[1:3]; ev = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
LRR_FAM = lambda n: n.startswith("LRR") or n in ("LRRNT", "LRRNT_2", "LRRCT", "LRR_adjacent")
lens = {}; name = None
for l in open(fa):
    if l.startswith(">"): name = l[1:].split()[0]; lens[name] = 0
    else: lens[name] += len(l.strip())
doms = collections.defaultdict(set)
for l in open(dt):
    if l.startswith("#"): continue
    c = l.split()
    if float(c[12]) <= ev: doms[c[3]].add(c[0])
N = len(lens); L = sorted(lens.values())
print("added proteins: %d | median length %d aa | shorter than 600 aa: %d (%.0f %%)" % (N, L[N // 2], sum(x < 600 for x in L), 100 * sum(x < 600 for x in L) / N))
lrr_any = [p for p in lens if any(LRR_FAM(x) for x in doms[p])]
lrr_main = [p for p in lens if any(x in ("LRR_1", "LRR_3", "LRR_4", "LRR_8") for x in doms[p])]
print("with any LRR-family domain (E<=%g): %d (%.1f %%) | of which found by the main profiles LRR_1/3/4/8: %d" % (ev, len(lrr_any), 100 * len(lrr_any) / N, len(lrr_main)))
fam = collections.Counter(x for p in lrr_any for x in doms[p] if LRR_FAM(x)); print("LRR families:", dict(fam.most_common(10)))
only_nb = [p for p in lens if doms[p] <= {"NB-ARC"}]
print("with nothing besides NB-ARC: %d (%.1f %%)" % (len(only_nb), 100 * len(only_nb) / N))
top = collections.Counter(x for p in lens for x in doms[p] if x != "NB-ARC" and not LRR_FAM(x))
print("top non-LRR domains (number of proteins):"); [print("  %-18s %d" % kv) for kv in top.most_common(15)]
