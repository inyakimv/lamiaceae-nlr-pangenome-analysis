#!/usr/bin/env python3
"""TNL census per species from the per-species hmmscan domtblout files (nlr_domain_search/per_species), independent of the family tree.
Usage: python3 35_tnl_census.py DOMTBLOUT_DIR OUT_PREFIX [CLSTR] [EXCLUDE=substr,substr]
A TNL = protein with an NB-ARC hit (i-E <= 1e-3) and a TIR or TIR_2 hit (i-E <= 1e-5) that STARTS BEFORE the first NB-ARC hit.
strict TNL = same, TIR i-E <= 1e-10 and >= 60 % of the TIR HMM covered.  If CLSTR is given, 'in_tree' tells whether the gene is in the NB-ARC family tree.
Outputs: OUT_PREFIX_per_species.tsv, OUT_PREFIX_genes.tsv"""
import sys, os, re, glob, collections
d, outp = sys.argv[1:3]
clstr = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "-") else None
excl = [e for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
intree = set()
if clstr:
    for l in open(clstr):
        m = re.search(r">(.+?)\.\.\.", l)
        if m and "|" in m.group(1): intree.add(m.group(1))
rows = []; tab = collections.defaultdict(collections.Counter)
for f in sorted(glob.glob(os.path.join(d, "*.domtblout"))):
    sp = os.path.basename(f)[:-len(".domtblout")]
    if any(e in sp for e in excl): continue
    nb = collections.defaultdict(list); tir = collections.defaultdict(list); lrr = set(); qlen = {}; nbend = {}
    for l in open(f):
        if l.startswith("#"): continue
        c = l.split(); dom, prot, ie = c[0], c[3], float(c[12])
        a0 = int(c[17]); cov = (int(c[16]) - int(c[15]) + 1) / int(c[2])
        if dom.startswith("LRR"): lrr.add(prot)   # same rule as extract_nbarc_coords.py: any LRR* hit in the file
        qlen[prot] = int(c[5])
        if dom == "NB-ARC" and ie <= 1e-3:
            nb[prot].append(a0)
            if prot not in nbend or a0 < nbend[prot][0]: nbend[prot] = (a0, int(c[18]))
        elif dom in ("TIR", "TIR_2") and ie <= 1e-5: tir[prot].append((a0, ie, cov))
    tab[sp]["NBARC_proteins"] = len(nb)
    for p in nb:
        start = min(nb[p]); t = [x for x in tir.get(p, []) if x[0] < start]
        if not t: continue
        strict = any(x[1] <= 1e-10 and x[2] >= 0.6 for x in t)
        it = (sp + "|" + p) in intree
        hl = p in lrr
        cterm = qlen[p] - nbend[p][1]
        rows.append((sp, p, int(strict), int(it) if clstr else "NA", int(hl), qlen[p], cterm))
        tab[sp]["TNL"] += 1; tab[sp]["TNL_strict"] += int(strict); tab[sp]["TNL_in_tree"] += int(it); tab[sp]["TNL_with_LRR"] += int(hl)
with open(outp + "_genes.tsv", "w") as o:
    o.write("species\tgene\tstrict\tin_tree\twith_LRR\tprotein_len\tcterm_len_after_NBARC\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
cols = ["NBARC_proteins", "TNL", "TNL_strict", "TNL_with_LRR", "TNL_in_tree"]
with open(outp + "_per_species.tsv", "w") as o:
    o.write("species\t" + "\t".join(cols) + "\n")
    for sp in sorted(tab): o.write(sp + "\t" + "\t".join(str(tab[sp][c]) for c in cols) + "\n")
    tot = collections.Counter()
    for sp in tab: tot.update(tab[sp])
    o.write("TOTAL\t" + "\t".join(str(tot[c]) for c in cols) + "\n")
print("species: %d | with TNL: %d | with strict TNL: %d" % (len(tab), sum(tab[s]["TNL"] > 0 for s in tab), sum(tab[s]["TNL_strict"] > 0 for s in tab)))
print(open(outp + "_per_species.tsv").read().splitlines()[-1], "(NBARC_proteins, TNL, TNL_strict, TNL_with_LRR, TNL_in_tree)")

# C-terminal length after the NB-ARC: truncated gene (short) or room for an LRR region that Pfam does not see (long)
def med(v): v = sorted(v); return v[len(v)//2] if v else "NA"
for lab, sel in (("TNL with LRR", [r for r in rows if r[4]]), ("TNL without LRR", [r for r in rows if not r[4]])):
    ct = [r[6] for r in sel]
    print("%-16s n=%3d | median C-terminal length after NB-ARC: %s aa | <150 aa (truncated-like): %d | >=300 aa (room for LRRs): %d" % (
        lab, len(ct), med(ct), sum(x < 150 for x in ct), sum(x >= 300 for x in ct)))
