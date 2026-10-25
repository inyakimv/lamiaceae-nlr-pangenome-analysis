#!/usr/bin/env python3
"""Tandem organisation of NB-ARC genes per species, from ANNEVO gene identifiers only.

Usage: python3 29_tandem_annevo.py CDHIT_CLSTR OUT.csv [GAPS=30,10] [EXCLUDE=species,substrings]
Genes are the members of the cd-hit (85%) clusters (all NB-ARC proteins predicted by ANNEVO, reference sequences excluded).
Gene index = the number in the identifier CONTIG-g<N>.t<k> (same definition as 04_inspect_clade.py).
Two NB-ARC genes are in tandem if they are on the same contig and their gene indices differ by <= GAP (default 30,
the rule used for the core clades; 10 is reported as a stricter alternative). A cluster = chain of consecutive NB-ARC genes
linked this way; a gene is 'in a tandem cluster' if its cluster has >= 2 genes.
OUT.csv: one row per species with n_nbarc_genes and, per gap G: pct_in_tandem_gG, n_clusters_gG, largest_cluster_gG.
The first gap in GAPS is the primary one; plot with plot_tandem_hist.py --tandem OUT.csv --col pct_in_tandem_g30."""
import re, sys, collections
import pandas as pd
clstr, out = sys.argv[1], sys.argv[2]
gaps = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "30,10").split(",")]
excl = [e for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
genes = []
for l in open(clstr):
    if l.startswith(">Cluster"): continue
    m = re.search(r">(.+?)\.\.\.", l)
    if m: genes.append(m.group(1))
rows, bad = [], 0
for g in genes:
    if g.startswith(("REF|", "Athaliana", "Knitens")) or any(e in g for e in excl): continue
    parts = g.split("|")
    m = re.match(r"(.+)-g(\d+)\.t\d+$", parts[1]) if len(parts) > 1 else None
    if not m: bad += 1; continue
    rows.append((parts[0].replace("_genomic", ""), m.group(1), int(m.group(2))))
print("NB-ARC genes parsed: %d (unparsed: %d)" % (len(rows), bad))
d = pd.DataFrame(rows, columns=["species", "contig", "gene_index"]).drop_duplicates()
res = []
for sp, s in d.groupby("species"):
    r = {"species": sp, "n_nbarc_genes": len(s)}
    for gap in gaps:
        sizes = []
        for ct, g in s.groupby("contig"):
            idx = sorted(g.gene_index); cur = 1
            for a, b in zip(idx, idx[1:]):
                if b - a <= gap: cur += 1
                else: sizes.append(cur); cur = 1
            sizes.append(cur)
        clustered = [x for x in sizes if x >= 2]
        r["pct_in_tandem_g%d" % gap] = round(100 * sum(clustered) / len(s), 1)
        r["n_clusters_g%d" % gap] = len(clustered)
        r["largest_cluster_g%d" % gap] = max(sizes)
    res.append(r)
t = pd.DataFrame(res).sort_values("species")
t.to_csv(out, index=False)
print("species:", len(t))
for gap in gaps:
    p = t["pct_in_tandem_g%d" % gap]
    print("gap <= %d genes: median %.1f%%, mean %.1f%%, range %.1f-%.1f%%; species < 80%%: %d" % (gap, p.median(), p.mean(), p.min(), p.max(), int((p < 80).sum())))
