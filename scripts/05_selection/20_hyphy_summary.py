#!/usr/bin/env python3
"""Summarise HyPhy FEL / MEME / BUSTED JSON outputs per clade.
Usage: python3 20_hyphy_summary.py HYPHY_OUT_DIR OUT.tsv [P=0.1]
Per clade: n sequences, n codons, FEL sites under purifying (beta<alpha) and positive (beta>alpha) selection at p<=P,
% purifying, MEME episodic sites (p<=P), BUSTED p-value."""
import sys, os, json, glob
d, out = sys.argv[1], sys.argv[2]
P = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
def load(path):
    """Return parsed JSON, or None if the file is missing, empty (analysis still running) or unreadable."""
    try:
        if os.path.getsize(path) == 0: return None
        return json.load(open(path))
    except (OSError, ValueError):
        return None
def load(path):
    try:
        if os.path.getsize(path) == 0: return None
        return json.load(open(path))
    except (OSError, ValueError):
        return None
def col(js, name):
    heads = [h[0] for h in js["MLE"]["headers"]]
    return heads.index(name)
rows = []
for fel in sorted(glob.glob(os.path.join(d, "*.fel.json"))):
    c = os.path.basename(fel)[:-9]
    j = load(fel)
    if j is None:
        print("skip %s: FEL JSON empty or unreadable (still running?)" % c); continue
    rec = {"clade": c}
    rec["n_seqs"] = j["input"]["number of sequences"]; rec["n_codons"] = j["input"]["number of sites"]
    ia, ib, ip = col(j, "alpha"), col(j, "beta"), col(j, "p-value")
    cont = [r for part in j["MLE"]["content"].values() for r in part]
    rec["fel_purifying"] = sum(r[ip] <= P and r[ib] < r[ia] for r in cont)
    rec["fel_positive"] = sum(r[ip] <= P and r[ib] > r[ia] for r in cont)
    rec["pct_purifying"] = round(100 * rec["fel_purifying"] / len(cont), 1)
    mp = os.path.join(d, c + ".meme.json")
    m = load(mp)
    if m is not None:
        im = col(m, "p-value")
        rec["meme_episodic"] = sum(r[im] <= P for part in m["MLE"]["content"].values() for r in part)
    bp = os.path.join(d, c + ".busted.json")
    b = load(bp)
    if b is not None:
        tr = b["test results"]; rec["busted_p"] = tr.get("p-value", tr.get("p"))
    rows.append(rec)
keys = ["clade", "n_seqs", "n_codons", "fel_purifying", "fel_positive", "pct_purifying", "meme_episodic", "busted_p"]
with open(out, "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows: o.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
print(open(out).read())
