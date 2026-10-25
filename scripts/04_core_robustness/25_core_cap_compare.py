#!/usr/bin/env python3
"""Compare the conserved core clades obtained with different clade-size caps.
Usage: python3 25_core_cap_compare.py OUTDIR REF_CAP OUT.tsv [MIN_SPECIES=35] [MAX_COPIES=3]
Expects OUTDIR/capN/stats_clade_conservation.tsv and OUTDIR/capN/refs_clade_members.tsv (from 24_core_cap_sensitivity.sh).
Core rule (same as 03_fig_core_clades_circular.py): n_species >= MIN_SPECIES and copies_per_species <= MAX_COPIES.
For each core clade at each cap: size, and the reference-cap core clade that it overlaps most, with two fractions:
  ref_in_new  = share of that reference clade's leaves that fall inside the new clade (is the reference clade kept whole?)
  new_in_ref  = share of the new clade's leaves that come from that reference clade (did the new clade absorb other lineages?)"""
import sys, os, glob, collections
import pandas as pd
out, ref_cap, outtsv = sys.argv[1], sys.argv[2], sys.argv[3]
ms = int(sys.argv[4]) if len(sys.argv) > 4 else 35
mc = float(sys.argv[5]) if len(sys.argv) > 5 else 3
def load(cap):
    d = os.path.join(out, "cap%s" % cap)
    st = pd.read_csv(os.path.join(d, "stats_clade_conservation.tsv"), sep="\t")
    mem = collections.defaultdict(set)
    for l in open(os.path.join(d, "refs_clade_members.tsv")):
        c, leaf = l.rstrip("\n").split("\t"); mem[c].add(leaf)
    core = st[(st.n_species >= ms) & (st.copies_per_species <= mc)]
    return core, mem
rcore, rmem = load(ref_cap)
print("Reference cap %s core clades: %s" % (ref_cap, ", ".join(rcore.clade)))
rows = []
caps = sorted((int(os.path.basename(p)[3:]) for p in glob.glob(os.path.join(out, "cap*")) if os.path.exists(os.path.join(p, "stats_clade_conservation.tsv"))))
for cap in caps:
    core, mem = load(cap)
    for _, r in core.iterrows():
        leaves = mem[r.clade]
        best = max(rcore.clade, key=lambda rc: len(leaves & rmem[rc]))
        ov = len(leaves & rmem[best])
        rows.append(dict(cap=cap, clade=r.clade, leaves=len(leaves), n_species=r.n_species, genes=r.n_genes,
                         copies_per_species=r.copies_per_species, closest_ref_core_clade=best if ov else "none",
                         ref_in_new="%.2f" % (ov / len(rmem[best])) if ov else "0",
                         new_in_ref="%.2f" % (ov / len(leaves)) if ov else "0"))
t = pd.DataFrame(rows)
t.to_csv(outtsv, sep="\t", index=False)
print(t.to_string(index=False))
print("\nNumber of core clades per cap:")
print(t.groupby("cap").size().to_string())
print("\nReference core clades recovered (ref_in_new >= 0.8 and new_in_ref >= 0.8) at each cap:")
for rc in rcore.clade:
    s = t[(t.closest_ref_core_clade == rc) & (t.ref_in_new.astype(float) >= 0.8) & (t.new_in_ref.astype(float) >= 0.8)]
    print("  %s: %s" % (rc, ", ".join(str(c) for c in sorted(set(s.cap))) or "none"))
