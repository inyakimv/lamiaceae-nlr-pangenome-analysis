#!/usr/bin/env python3
"""Does the 'complete-like' filter (NB-ARC + >=1 LRR) change the between-species comparison?
Usage: python3 36_filter_robustness.py DOMTBLOUT_DIR OUT_PREFIX [PLOIDY_CSV] [EXCLUDE=substr,substr]
Per species: n_NBARC_any = proteins with an NB-ARC hit; n_complete_like = those that also have a LRR* hit (rule of extract_nbarc_coords.py).
Reports Spearman rho between both counts and, if PLOIDY_CSV (species,ploidy) is given, the polyploid vs diploid difference (Mann-Whitney) for each count.
Needs numpy/scipy (conda base)."""
import sys, os, glob
from scipy.stats import spearmanr, mannwhitneyu
d, outp = sys.argv[1:3]
pl = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "-") else None
excl = [e for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
rows = []
for f in sorted(glob.glob(os.path.join(d, "*.domtblout"))):
    sp = os.path.basename(f)[:-len(".domtblout")]
    if any(e in sp for e in excl): continue
    nb, lrr = set(), set()
    for l in open(f):
        if l.startswith("#"): continue
        c = l.split()
        if c[0] == "NB-ARC": nb.add(c[3])
        elif c[0].startswith("LRR"): lrr.add(c[3])
    rows.append((sp, len(nb), len(nb & lrr)))
def short(sp):  # GCA_xxx.1_Genus_species_genomic -> Genus species
    p = sp.replace("_genomic", "").split("_")[2:]; return " ".join(p)
ploidy = {}
if pl:
    import csv
    for r in csv.DictReader(open(pl)): ploidy[r["species"].strip()] = r["ploidy"].strip().lower()
with open(outp + "_counts.tsv", "w") as o:
    o.write("species\tn_NBARC_any\tn_complete_like\tfraction_complete_like\tploidy\n")
    for sp, a, b in rows: o.write("%s\t%d\t%d\t%.3f\t%s\n" % (short(sp), a, b, b / a if a else 0, ploidy.get(short(sp), "NA")))
a = [r[1] for r in rows]; b = [r[2] for r in rows]
rho, p = spearmanr(a, b)
print("species: %d | NB-ARC any: median %d (range %d-%d) | complete-like: median %d (range %d-%d)" % (
    len(rows), sorted(a)[len(a)//2], min(a), max(a), sorted(b)[len(b)//2], min(b), max(b)))
print("Spearman rho between the two counts: %.3f (p=%.2g)" % (rho, p))
print("fraction complete-like per species: min %.2f, max %.2f" % (min(y / x for x, y in zip(a, b)), max(y / x for x, y in zip(a, b))))
if ploidy:
    pol = [r for r in rows if ploidy.get(short(r[0])) == "polyploid"]; dip = [r for r in rows if ploidy.get(short(r[0])) == "diploid"]
    print("diploid n=%d, polyploid n=%d" % (len(dip), len(pol)))
    if not pol or not dip: sys.exit("ERROR: no species matched the ploidy table (check names such as 'Mentha x piperita')")
    unk = [short(r[0]) for r in rows if short(r[0]) not in ploidy]
    if unk: print("WARNING no ploidy for:", unk)
    for name, k in (("NB-ARC any", 1), ("complete-like", 2)):
        x = [r[k] for r in pol]; y = [r[k] for r in dip]
        u, pv = mannwhitneyu(x, y, alternative="two-sided")
        print("%-14s median polyploid %d vs diploid %d | Mann-Whitney p=%.3g" % (name, sorted(x)[len(x)//2], sorted(y)[len(y)//2], pv))
