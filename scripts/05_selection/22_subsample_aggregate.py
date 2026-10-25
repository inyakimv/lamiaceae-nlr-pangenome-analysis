#!/usr/bin/env python3
"""Aggregate the rarefaction replicates (output of 20_hyphy_summary.py on hyphy_sub) per clade.
Usage: python3 22_subsample_aggregate.py SUB.tsv FULL.tsv OUT.tsv [REFERENCE_CLADE=Clade76]
SUB.tsv  : 20_hyphy_summary.py hyphy_sub SUB.tsv        (rows named CLADE_nN_rK)
FULL.tsv : 20_hyphy_summary.py hyphy_out FULL.tsv        (only REFERENCE_CLADE is taken, as one replicate: it already has N sequences)
Per clade: replicates, mean and range of % purifying (FEL), % positive sites (FEL), % episodic sites (MEME), and how many replicates have BUSTED p<0.05."""
import sys, csv, re
sub, full, out = sys.argv[1:4]; ref = sys.argv[4] if len(sys.argv) > 4 else "Clade76"
def rows(p): return list(csv.DictReader(open(p), delimiter="\t"))
recs = {}
def add(clade, r):
    nc = int(r["n_codons"])
    d = {"pur": 100 * int(r["fel_purifying"]) / nc, "pos": 100 * int(r["fel_positive"]) / nc}
    d["epi"] = 100 * int(r["meme_episodic"]) / nc if r["meme_episodic"] != "" else None
    d["bp"] = float(r["busted_p"]) if r["busted_p"] != "" else None
    d["n"] = int(r["n_seqs"]); recs.setdefault(clade, []).append(d)
for r in rows(sub):
    m = re.match(r"(.+)_n\d+_r\d+$", r["clade"])
    if m: add(m.group(1), r)
for r in rows(full):
    if r["clade"] == ref: add(ref, r)
def fmt(v):
    v = [x for x in v if x is not None]
    return "NA" if not v else "%.1f (%.1f-%.1f)" % (sum(v) / len(v), min(v), max(v))
keys = ["clade", "seqs_per_rep", "n_reps", "pct_purifying_mean(range)", "pct_positive_FEL", "pct_episodic_MEME", "BUSTED_p<0.05"]
with open(out, "w") as o:
    o.write("\t".join(keys) + "\n")
    for c, L in sorted(recs.items()):
        b = [x["bp"] for x in L if x["bp"] is not None]
        o.write("\t".join([c, str(L[0]["n"]), str(len(L)), fmt([x["pur"] for x in L]), fmt([x["pos"] for x in L]),
                           fmt([x["epi"] for x in L]), "%d/%d" % (sum(p < 0.05 for p in b), len(b))]) + "\n")
print(open(out).read())
