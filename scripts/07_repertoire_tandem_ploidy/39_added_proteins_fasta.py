#!/usr/bin/env python3
"""Write a FASTA with the proteins ADDED by the relaxed criterion (NB-ARC, no Pfam LRR, >= MIN_CT aa after the NB-ARC), headers 'species|protein'.
Then scan them against Pfam-A (hmmscan) and summarise with 40_summarise_added_domains.py.
Usage: python3 39_added_proteins_fasta.py DOMTBLOUT_DIR PROTEIN_DIR OUT.faa [EXCLUDE=substr,substr] [MIN_CT=300]"""
import sys, os, glob
d, pdir, out = sys.argv[1:4]
excl = [e for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
min_ct = int(sys.argv[5]) if len(sys.argv) > 5 else 300
def read_fa(f):
    o, name, buf = {}, None, []
    for l in open(f):
        if l.startswith(">"):
            if name: o[name] = "".join(buf)
            name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if name: o[name] = "".join(buf)
    return o
n = 0
with open(out, "w") as fo:
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
        prot = read_fa(os.path.join(pdir, sp + ".annevo_proteins.fasta"))
        for p, (ie, end) in best.items():
            if p in lrr or p not in prot or qlen[p] - end < min_ct: continue
            fo.write(">%s|%s\n%s\n" % (sp, p, prot[p])); n += 1
print("added proteins written:", n, "->", out)
