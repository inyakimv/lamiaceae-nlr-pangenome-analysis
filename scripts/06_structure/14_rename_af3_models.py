#!/usr/bin/env python3
"""Rename AlphaFold Server model files to the gene ids used by 13_structure_pairs.py.
Usage: python3 14_rename_af3_models.py structure_reps.tsv IN_DIR OUT_DIR
IN_DIR holds fold_<job name>_model_0.cif files (job names are lowercased/sanitised by the server).
AtADR1 is added automatically. Prints files that could not be matched."""
import sys, os, re, shutil
reps, indir, outdir = sys.argv[1:4]
safe = lambda x: re.sub(r"[^A-Za-z0-9_.-]", "_", x)
key = lambda x: re.sub(r"[^a-z0-9]", "", x.lower())
genes = [l.rstrip("\n").split("\t")[1] for l in open(reps) if l.strip()] + ["AtADR1"]
table = {key(safe(g)): safe(g) for g in genes}
os.makedirs(outdir, exist_ok=True)
done, bad = set(), []
for f in sorted(os.listdir(indir)):
    m = re.match(r"fold_(.+)_model_0\.cif$", f)
    if not m: continue
    k = key(m.group(1))
    hit = [v for kk, v in table.items() if kk == k or kk.startswith(k) or k.startswith(kk)]
    if len(hit) != 1: bad.append((f, hit)); continue
    shutil.copy(os.path.join(indir, f), os.path.join(outdir, hit[0] + ".cif")); done.add(hit[0])
print("renamed", len(done), "of", len(genes))
print("not matched:", bad)
print("missing models:", sorted(set(table.values()) - done))
