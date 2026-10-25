#!/usr/bin/env python3
"""Sensitivity analysis (NOT the main criterion): relaxed 'complete-like' = NB-ARC AND (any LRR* hit OR >= MIN_CT aa after the NB-ARC).
Main criterion (extract_nbarc_coords.py) = NB-ARC AND any LRR* hit. NB-ARC = best hit per protein, as in the main pipeline (no extra E-value cut).
Usage: python3 37_relaxed_criterion.py DOMTBLOUT_DIR OUT_PREFIX [PLOIDY_CSV|-] [EXCLUDE=substr,substr] [MIN_CT=300]
Outputs: OUT_PREFIX_counts.tsv (per species: main, relaxed, added, fraction, ploidy) and a printed summary.
Needs scipy (conda base)."""
import sys, os, glob, csv
from scipy.stats import spearmanr, mannwhitneyu
d, outp = sys.argv[1:3]
pl = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "-") else None
excl = [e for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
min_ct = int(sys.argv[5]) if len(sys.argv) > 5 else 300
def short(sp): return " ".join(sp.replace("_genomic", "").split("_")[2:])
ploidy = {}
if pl:
    for r in csv.DictReader(open(pl)): ploidy[r["species"].strip()] = r["ploidy"].strip().lower()
rows = []
for f in sorted(glob.glob(os.path.join(d, "*.domtblout"))):
    sp = os.path.basename(f)[:-len(".domtblout")]
    if any(e in sp for e in excl): continue
    best = {}; lrr = set(); qlen = {}
    for l in open(f):
        if l.startswith("#"): continue
        c = l.split(); dom, p = c[0], c[3]; qlen[p] = int(c[5])
        if dom == "NB-ARC":
            ie = float(c[12])
            if p not in best or ie < best[p][0]: best[p] = (ie, int(c[18]))
        elif dom.startswith("LRR"): lrr.add(p)
    main = {p for p in best if p in lrr}
    relaxed = {p for p, (ie, end) in best.items() if p in lrr or qlen[p] - end >= min_ct}
    rows.append((short(sp), len(best), len(main), len(relaxed)))
with open(outp + "_counts.tsv", "w") as o:
    o.write("species\tn_NBARC_any\tn_main\tn_relaxed\tadded\tfraction_relaxed_over_any\tploidy\n")
    for s, a, m, r in rows: o.write("%s\t%d\t%d\t%d\t%d\t%.3f\t%s\n" % (s, a, m, r, r - m, r / a if a else 0, ploidy.get(s, "NA")))
M = [r[2] for r in rows]; R = [r[3] for r in rows]; A = [r[1] for r in rows]
med = lambda v: sorted(v)[len(v) // 2]
print("species: %d | NB-ARC any: median %d | main (NB-ARC+LRR): median %d (range %d-%d) | relaxed (>=%d aa C-term): median %d (range %d-%d)" % (
    len(rows), med(A), med(M), min(M), max(M), min_ct, med(R), min(R), max(R)))
print("total genes: main %d, relaxed %d (+%d, %.0f %%)" % (sum(M), sum(R), sum(R) - sum(M), 100 * (sum(R) - sum(M)) / sum(M)))
rho, p = spearmanr(M, R); print("Spearman rho main vs relaxed: %.3f (p=%.2g)" % (rho, p))
if ploidy:
    pol = [r for r in rows if ploidy.get(r[0]) == "polyploid"]; dip = [r for r in rows if ploidy.get(r[0]) == "diploid"]
    if not pol or not dip: sys.exit("ERROR: no species matched the ploidy table")
    print("diploid n=%d, polyploid n=%d" % (len(dip), len(pol)))
    for name, k in (("main", 2), ("relaxed", 3)):
        x = [r[k] for r in pol]; y = [r[k] for r in dip]
        print("%-8s median polyploid %d vs diploid %d | Mann-Whitney p=%.3g" % (name, med(x), med(y), mannwhitneyu(x, y, alternative="two-sided")[1]))
