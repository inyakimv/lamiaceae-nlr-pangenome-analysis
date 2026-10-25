#!/usr/bin/env python3
"""Estimate how many of the proteins ADDED by the relaxed criterion (NB-ARC, no Pfam LRR, >= MIN_CT aa after NB-ARC) really carry a LRR region.
Method: count non-overlapping matches of the plant-NLR LRR motif L..L.L..[NC] in the C-terminal region (first WIN aa after the NB-ARC).
  positives  = proteins with NB-ARC + Pfam LRR (main criterion)       -> sensitivity of the motif count
  negatives  = proteins with NO hit at all in the domtblout (no NB-ARC, no LRR, no TIR); a random window of WIN aa -> false-positive rate
  test set   = proteins added by the relaxed criterion                  -> observed positive rate
  estimated true-LRR fraction among the added = (observed - FPR) / (sensitivity - FPR)   (Rogan-Gladen correction)
Usage: python3 38_lrr_motif_calibration.py DOMTBLOUT_DIR PROTEIN_DIR OUT_PREFIX [EXCLUDE=substr,substr] [MIN_CT=300] [WIN=300] [NNEG=3000] [SEED=1]
Protein files: PROTEIN_DIR/<species>.annevo_proteins.fasta (species = domtblout file name without .domtblout)."""
import sys, os, re, glob, random
d, pdir, outp = sys.argv[1:4]
excl = [e for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
min_ct = int(sys.argv[5]) if len(sys.argv) > 5 else 300
win = int(sys.argv[6]) if len(sys.argv) > 6 else 300
nneg = int(sys.argv[7]) if len(sys.argv) > 7 else 3000
rng = random.Random(int(sys.argv[8]) if len(sys.argv) > 8 else 1)
MOT = re.compile(r"L..L.L..[NC]")
def count(seq): return len(MOT.findall(seq))
def read_fa(f):
    out, name, buf = {}, None, []
    for l in open(f):
        if l.startswith(">"):
            if name: out[name] = "".join(buf)
            name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if name: out[name] = "".join(buf)
    return out
pos, add, neg = [], [], []; per_sp = []
for f in sorted(glob.glob(os.path.join(d, "*.domtblout"))):
    sp = os.path.basename(f)[:-len(".domtblout")]
    if any(e in sp for e in excl): continue
    pf = os.path.join(pdir, sp + ".annevo_proteins.fasta")
    if not os.path.exists(pf): print("WARNING no protein file for", sp); continue
    best = {}; lrr = set(); seen = set()
    for l in open(f):
        if l.startswith("#"): continue
        c = l.split(); dom, p = c[0], c[3]; seen.add(p)
        if dom == "NB-ARC":
            ie = float(c[12])
            if p not in best or ie < best[p][0]: best[p] = (ie, int(c[18]))
        elif dom.startswith("LRR"): lrr.add(p)
    prot = read_fa(pf)
    s_pos = s_add = 0
    for p, (ie, end) in best.items():
        if p not in prot: continue
        ct = prot[p][end:end + win]; n = count(ct)
        if p in lrr: pos.append(n); s_pos += 1
        elif len(prot[p]) - end >= min_ct: add.append(n); s_add += 1
    cand = [p for p, s in prot.items() if p not in seen and len(s) >= win + 100]
    for p in rng.sample(cand, min(len(cand), nneg // 45 + 1)):
        s = prot[p]; i = rng.randint(0, len(s) - win); neg.append(count(s[i:i + win]))
    per_sp.append((sp, s_pos, s_add))
print("proteins: positives (NB-ARC+Pfam LRR) %d | added by relaxed criterion %d | negatives (random windows, no domain hit) %d" % (len(pos), len(add), len(neg)))
print("k = min. motif matches in the %d aa after NB-ARC | sensitivity (positives) | false-positive rate (negatives) | observed rate in added | estimated true-LRR fraction among added | estimated real LRR genes added" % win)
rows = []
for k in range(1, 9):
    sens = sum(x >= k for x in pos) / len(pos); fpr = sum(x >= k for x in neg) / len(neg); obs = sum(x >= k for x in add) / len(add)
    est = max(0.0, min(1.0, (obs - fpr) / (sens - fpr))) if sens > fpr else float("nan")
    rows.append((k, sens, fpr, obs, est, est * len(add)))
    print("k>=%d | sens %.2f | FPR %.3f | observed %.2f | estimated true fraction %.2f | ~%d of %d" % (k, sens, fpr, obs, est, est * len(add), len(add)))
with open(outp + "_calibration.tsv", "w") as o:
    o.write("k\tsensitivity\tFPR\tobserved_added\test_true_fraction\test_true_genes\n")
    for r in rows: o.write("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.0f\n" % r)
