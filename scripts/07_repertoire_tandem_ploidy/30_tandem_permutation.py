#!/usr/bin/env python3
"""Permutation control for the tandem organisation of NB-ARC genes (companion of 29_tandem_annevo.py).

Usage: python3 30_tandem_permutation.py CDHIT_CLSTR ANNEVO_GFF3_DIR OUT.csv [GAP=30] [NPERM=1000] [SEED=1] [EXCLUDE=species,substrings]
Null model, per species: on every contig keep the number k of NB-ARC genes and draw k distinct gene indices uniformly among
all genes of that contig (N = number of annotated genes on the contig in the ANNEVO GFF3, at least the largest observed index).
Same tandem rule as 29 (consecutive NB-ARC genes within GAP gene indices; a gene is clustered if its cluster has >= 2 genes).
Per species: observed %, null mean and 95% range, excess = observed - null mean, one-sided empirical p = (1 + #null >= observed) / (1 + NPERM).
Also prints the correlation of the excess with the number of NB-ARC genes (is the Spearman rho of 29 only a gene-number effect?)."""
import re, sys, os, collections
import numpy as np
import pandas as pd
from scipy import stats
clstr, gffdir, out = sys.argv[1:4]
gap = int(sys.argv[4]) if len(sys.argv) > 4 else 30
nperm = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
rng = np.random.default_rng(int(sys.argv[6]) if len(sys.argv) > 6 else 1)
excl = [e for e in (sys.argv[7] if len(sys.argv) > 7 else "").split(",") if e]
rows = []
for l in open(clstr):
    if l.startswith(">Cluster"): continue
    m = re.search(r">(.+?)\.\.\.", l)
    if not m: continue
    g = m.group(1)
    if g.startswith(("REF|", "Athaliana", "Knitens")) or any(e in g for e in excl): continue
    p = g.split("|"); mm = re.match(r"(.+)-g(\d+)\.t\d+$", p[1]) if len(p) > 1 else None
    if mm: rows.append((p[0].replace("_genomic", ""), mm.group(1), int(mm.group(2))))
d = pd.DataFrame(rows, columns=["species", "contig", "idx"]).drop_duplicates()
def genes_per_contig(sp):
    f = os.path.join(gffdir, sp + "_genomic.annevo.gff3")
    cnt = collections.Counter()
    if not os.path.exists(f): return None
    for l in open(f):
        if l.startswith("#"): continue
        c = l.split("\t")
        if len(c) > 8 and c[2] == "gene": cnt[c[0]] += 1
    return cnt
def frac_clustered(contigs):
    n_clu = 0; n = 0
    for idx in contigs:
        idx = np.sort(idx); n += len(idx)
        if len(idx) < 2: continue
        close = np.diff(idx) <= gap
        # a gene is clustered if it has a close neighbour on either side
        left = np.concatenate(([False], close)); right = np.concatenate((close, [False]))
        n_clu += int(np.sum(left | right))
    return 100.0 * n_clu / n
res = []
for sp, s in d.groupby("species"):
    gpc = genes_per_contig(sp)
    if gpc is None: print("WARNING: no GFF3 for", sp, "- using the largest observed index per contig as N")
    obs_c, size_k, size_N = [], [], []
    for ct, g in s.groupby("contig"):
        idx = np.array(sorted(g.idx)); obs_c.append(idx)
        N = max(int(idx.max()), (gpc or {}).get(ct, 0)); size_k.append(len(idx)); size_N.append(N)
    obs = frac_clustered(obs_c)
    null = np.empty(nperm)
    for i in range(nperm):
        null[i] = frac_clustered([rng.choice(N, size=k, replace=False) + 1 for k, N in zip(size_k, size_N)])
    res.append(dict(species=sp, n_nbarc_genes=len(s), observed_pct=round(obs, 1), null_mean_pct=round(null.mean(), 1),
                    null_lo95=round(np.percentile(null, 2.5), 1), null_hi95=round(np.percentile(null, 97.5), 1),
                    excess=round(obs - null.mean(), 1), p_one_sided=round((1 + np.sum(null >= obs)) / (1 + nperm), 4)))
t = pd.DataFrame(res).sort_values("species"); t.to_csv(out, index=False)
print("gap <= %d, %d permutations, %d species" % (gap, nperm, len(t)))
print("observed median %.1f%% | null median %.1f%% | median excess %.1f points" % (t.observed_pct.median(), t.null_mean_pct.median(), t.excess.median()))
print("species with observed > null (p < 0.05): %d of %d; with p < 0.001: %d" % ((t.p_one_sided < 0.05).sum(), len(t), (t.p_one_sided < 0.001).sum()))
r, pr = stats.spearmanr(t.n_nbarc_genes, t.excess)
print("Spearman rho (n genes vs excess over null) = %.2f, p = %.3g" % (r, pr))
