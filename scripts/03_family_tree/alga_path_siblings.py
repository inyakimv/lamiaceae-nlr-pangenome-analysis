import re
from collections import Counter
from ete3 import Tree

t = Tree("nbarc_family_tree_plus_outgroup.treefile", format=1)
leaf = [l for l in t.get_leaves() if "GAQ83946.1" in l.name][0]
total = len(t.get_leaves())

def species(name):
    m = re.match(r"^(?:GC[AF]_\d+\.\d+_)?(.+?)_genomic\|", name)
    return m.group(1) if m else ("alga" if "GAQ" in name else name)

def comp(names):
    n_alga = sum("GAQ" in x for x in names)
    n_ath = sum(x.startswith("Athaliana") for x in names)
    return n_alga, n_ath, len(names) - n_alga - n_ath

# Complement (the other side of each branch on the path)
print(f"Total leaves: {total}")
print("\n== Other side of each branch (everything outside the clade) ==")
node = leaf
for level in range(1, 9):
    node = node.up
    if node is None: break
    inside = set(node.get_leaf_names())
    outside = [x for x in t.get_leaf_names() if x not in inside]
    a, r, l = comp(outside)
    print(f"outside level {level}: leaves={len(outside)} alga={a} Arabidopsis={r} Lamiaceae={l}")

# Sibling subtrees hanging off the path from the alga towards the root
print("\n== Sibling groups along the path (nearest first) ==")
node = leaf
for level in range(1, 13):
    parent = node.up
    if parent is None: break
    for sib in parent.children:
        if sib is node: continue
        names = sib.get_leaf_names()
        a, r, l = comp(names)
        top = Counter(species(x) for x in names).most_common(3)
        print(f"step {level}: sibling leaves={len(names)} alga={a} Arabidopsis={r} "
              f"Lamiaceae={l} support={sib.support} branch={sib.dist:.3f} top={top}")
    node = parent
