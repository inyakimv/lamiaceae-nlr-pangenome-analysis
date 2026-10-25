import sys
from ete3 import Tree

def report(path):
    t = Tree(path, format=0)
    t.set_outgroup(t.get_midpoint_outgroup())
    leaves = t.get_leaves()
    print(f"\n=== {path}: {len(leaves)} leaves ===")

    sup = [n.support for n in t.traverse() if not n.is_leaf() and not n.is_root()]
    n = len(sup)
    for thr in (0.5, 0.7, 0.8, 0.9, 0.95):
        k = sum(s >= thr for s in sup)
        print(f"  support >= {thr}: {k} of {n} ({100*k/n:.1f}%)")

    is_ath = lambda name: name.startswith("Athaliana")
    pure = {}
    for node in t.traverse("postorder"):
        if node.is_leaf():
            pure[node] = is_ath(node.name)
        else:
            pure[node] = all(pure[c] for c in node.children)
    blocks = [nd for nd in t.traverse("preorder")
              if pure[nd] and (nd.is_root() or not pure[nd.up])]
    sizes = sorted((len(nd.get_leaves()) for nd in blocks), reverse=True)
    print(f"  Arabidopsis: {sum(sizes)} leaves in {len(blocks)} pure monophyletic blocks")
    print(f"  Block sizes: {sizes}")

    big = max(blocks, key=lambda nd: len(nd.get_leaves()))
    print(f"  Largest block: {len(big.get_leaves())} leaves, stem support {big.support}")
    ath = [l for l in leaves if is_ath(l.name)]
    mrca = t.get_common_ancestor(ath)
    print(f"  MRCA of all Arabidopsis: {len(mrca.get_leaves())} leaves")

for p in sys.argv[1:]:
    report(p)
