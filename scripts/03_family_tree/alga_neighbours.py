from ete3 import Tree
t = Tree("nbarc_family_tree_plus_outgroup.treefile", format=1)
leaf = [l for l in t.get_leaves() if "GAQ83946.1" in l.name][0]
print("Terminal branch length:", round(leaf.dist, 3))
node = leaf
for level in range(1, 8):
    node = node.up
    if node is None: break
    names = node.get_leaf_names()
    n_alga = sum("GAQ" in x for x in names)
    n_ath = sum(x.startswith("Athaliana") for x in names)
    n_lam = len(names) - n_alga - n_ath
    print(f"level {level}: leaves={len(names)} alga={n_alga} Arabidopsis={n_ath} Lamiaceae={n_lam} support={node.support}")
