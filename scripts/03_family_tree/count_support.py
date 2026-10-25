import re, sys
tree = open("nbarc_family_tree.treefile").read()
vals = [float(x) for x in re.findall(r"\)([0-9]*\.?[0-9]+):", tree)]
n = len(vals)
print(f"Internal nodes with support value: {n}")
for thr in (0.5, 0.7, 0.8, 0.9, 0.95):
    k = sum(v >= thr for v in vals)
    print(f">= {thr}: {k} ({100*k/n:.1f}%)")
