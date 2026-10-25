#!/usr/bin/env python3
"""Mean pLDDT per AlphaFold3 model (B-factor column of the mmCIF) and % residues with pLDDT < 70.
Usage: python3 15_plddt_summary.py MODELS_DIR OUT.tsv
Uses CA atoms only (one per residue). Flag: mean < 70 = 'low'."""
import sys, os, glob
def read(path):
    cols, vals, inloop = [], [], False
    for l in open(path):
        if l.startswith("_atom_site."): cols.append(l.strip().split(".")[1]); continue
        if cols and l.startswith(("ATOM", "HETATM")): vals.append(l.split())
        elif cols and vals and (l.startswith("#") or l.startswith("loop_")): break
    i_atom, i_b = cols.index("label_atom_id"), cols.index("B_iso_or_equiv")
    return [float(v[i_b]) for v in vals if v[i_atom] == "CA"]
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.cif"))):
    b = read(f)
    rows.append((os.path.basename(f)[:-4], len(b), sum(b) / len(b), 100 * sum(x < 70 for x in b) / len(b)))
with open(sys.argv[2], "w") as o:
    o.write("model\tn_residues\tmean_plddt\tpct_below_70\tflag\n")
    for n, k, m, p in rows: o.write(f"{n}\t{k}\t{m:.1f}\t{p:.1f}\t{'low' if m < 70 else 'ok'}\n")
print(len(rows), "models; low:", sum(m < 70 for _, _, m, _ in rows))
