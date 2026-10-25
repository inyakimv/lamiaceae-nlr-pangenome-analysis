#!/usr/bin/env python3
"""Conserved-core statistics with and without one or more excluded species (same logic as 02_clade_stats.py).

Usage: 28_clade_stats_exclude.py CDHIT_CLSTR CLADE_MEMBERS_TSV CLADE_INFO_CSV OUTPREFIX EXCLUDE [MIN_SPECIES=35] [MAX_COPIES=3]
  EXCLUDE = comma-separated substrings of species names, e.g. Paraleonurus_japonicus
Writes OUTPREFIX_clade_conservation_excl.tsv (excluded run) and prints, for the all-species run and the excluded run,
the core clades (n_species >= MIN_SPECIES and copies_per_species <= MAX_COPIES) and, as a second line, the same with MIN_SPECIES-1.
The clade partition is NOT recomputed: genes of the excluded species are simply not counted (cd-hit clusters are still expanded)."""
import re, sys, statistics, collections
import pandas as pd
clstr, members_tsv, info_csv, prefix, excl = sys.argv[1:6]
ms = int(sys.argv[6]) if len(sys.argv) > 6 else 35
mc = float(sys.argv[7]) if len(sys.argv) > 7 else 3
excl = [e for e in excl.split(",") if e]
members, cur = {}, []
def flush():
    global cur
    rep = [m for m, star in cur if star]
    if rep: members[rep[0]] = [m for m, _ in cur]
    cur = []
for line in open(clstr):
    if line.startswith(">Cluster"): flush(); continue
    m = re.search(r">(.+?)\.\.\.", line)
    if m: cur.append((m.group(1), "*" in line.split("...")[1][:6]))
flush()
clade = {}
for line in open(members_tsv):
    c, leaf = line.rstrip("\n").split("\t"); clade[leaf] = c
info = pd.read_csv(info_csv).set_index("clade_id")
species = lambda name: name.split("|")[0].split("_genomic")[0]
def stats(skip):
    per = collections.defaultdict(collections.Counter)
    for rep, mem in members.items():
        c = clade.get(rep)
        if not c: continue
        for m in mem:
            if m.startswith(("Athaliana", "Knitens", "REF")): continue
            if skip and any(e in m for e in skip): continue
            per[c][species(m)] += 1
    rows = []
    for c, cnt in per.items():
        cp = list(cnt.values())
        rows.append(dict(clade=c, n_species=len(cnt), n_genes=sum(cp), median_copies=statistics.median(cp),
                         copies_per_species=round(sum(cp) / len(cp), 2)))
    return pd.DataFrame(rows).sort_values(["n_species", "n_genes"], ascending=[False, True])
a, b = stats([]), stats(excl)
b.to_csv(prefix + "_clade_conservation_excl.tsv", sep="\t", index=False)
def core(t, thr): return t[(t.n_species >= thr) & (t.copies_per_species <= mc)]
for label, t in (("all species", a), ("excluding " + ",".join(excl), b)):
    for thr in (ms, ms - 1):
        c = core(t, thr)
        print("%-34s min_species=%d  core (%d): %s" % (label, thr, len(c), ", ".join("%s[%d sp, %.2f cp]" % (r.clade, r.n_species, r.copies_per_species) for r in c.itertuples())))
m = a.merge(b, on="clade", suffixes=("_all", "_excl"))
ch = m[(m.n_species_all != m.n_species_excl) & (m.n_species_all >= ms - 3)]
print("\nClades near the threshold whose species count changes (%d):" % len(ch))
print(ch[["clade", "n_species_all", "n_species_excl", "copies_per_species_all", "copies_per_species_excl"]].to_string(index=False))
