#!/usr/bin/env python3
"""Robustness of the saturation criterion: pairwise synonymous distance (dS, Nei-Gojobori 1986 with Jukes-Cantor correction)
computed on random pairs of every codon alignment, next to the 3rd-position p-distance rule used in the paper.

Usage: python3 47_ds_robustness.py OUT.tsv ALN1.fasta [ALN2.fasta ...]   (alignments in frame, '-' = gap)
Options by environment variables: PAIRS (default 3000 random pairs per alignment), SEED (default 1).
For every alignment: sequences, pairs used, median dS over pairs where dS is defined, % pairs with undefined dS
(synonymous sites saturated: the NG86 correction fails), % pairs with dS > 1, > 2 or undefined, and the median 3rd-position p-distance
for comparison. Pairs are scored on codons that are ungapped in both sequences."""
import sys, os, random, statistics, itertools
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
out, files = sys.argv[1], sys.argv[2:]
PAIRS = int(os.environ.get("PAIRS", 3000)); SEED = int(os.environ.get("SEED", 1))
STOPS = {"TAA", "TAG", "TGA"}
def read(f):
    d, k = {}, None
    for l in open(f):
        l = l.strip()
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        elif k: d[k].append(l.upper())
    return {a: "".join(b) for a, b in d.items()}
def pair_codons(a, b):
    ca = [a[i:i+3] for i in range(0, len(a) - 2, 3)]; cb = [b[i:i+3] for i in range(0, len(b) - 2, 3)]
    keep = [(x, y) for x, y in zip(ca, cb) if set(x + y) <= set("ACGT") and x not in STOPS and y not in STOPS]
    return keep
def p3(keep):
    return sum(x[2] != y[2] for x, y in keep) / len(keep) if keep else float("nan")
def ds_ng86(keep):
    if len(keep) < 30: return None
    s1 = "".join(x for x, _ in keep); s2 = "".join(y for _, y in keep)
    try:
        dn, ds = cal_dn_ds(CodonSeq(s1), CodonSeq(s2), method="NG86")
    except Exception:
        return None
    if ds is None or ds != ds or ds < 0 or ds > 50: return None
    return ds
rows = []
for f in files:
    s = read(f); names = list(s); n = len(names)
    if n < 3: continue
    random.seed(SEED)
    allp = n * (n - 1) // 2
    if allp <= PAIRS: pairs = list(itertools.combinations(names, 2))
    else:
        pairs = set()
        while len(pairs) < PAIRS:
            a, b = random.sample(names, 2); pairs.add((a, b) if a < b else (b, a))
        pairs = sorted(pairs)
    dss, undef, p3s = [], 0, []
    for a, b in pairs:
        keep = pair_codons(s[a], s[b]); p3s.append(p3(keep)); d = ds_ng86(keep)
        if d is None: undef += 1
        else: dss.append(d)
    m = len(pairs)
    over = lambda t: 100 * (sum(d > t for d in dss) + undef) / m
    rows.append([os.path.basename(f).split(".")[0], n, m, f"{statistics.median(dss):.2f}" if dss else "NA",
                 f"{100*undef/m:.1f}", f"{over(1):.1f}", f"{over(2):.1f}", f"{statistics.median([x for x in p3s if x == x]):.3f}"])
hdr = ["alignment", "seqs", "pairs", "median_dS_defined", "pct_dS_undefined", "pct_dS>1_or_undef", "pct_dS>2_or_undef", "median_p3"]
with open(out, "w") as o:
    o.write("\t".join(hdr) + "\n"); [o.write("\t".join(map(str, r)) + "\n") for r in rows]
print("\t".join(hdr)); [print("\t".join(map(str, r))) for r in rows]
