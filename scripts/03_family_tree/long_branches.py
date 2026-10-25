import sys, statistics
from ete3 import Tree
t = Tree(sys.argv[1], format=0)
leaves = t.get_leaves()
d = sorted(((l.dist, l.name) for l in leaves), reverse=True)
med = statistics.median(x for x, _ in d)
print(f"{len(leaves)} leaves; median terminal branch {med:.3f}")
for thr in (2, 3, 5):
    print(f"  > {thr}x median: {sum(x > thr*med for x, _ in d)} leaves")
print("Top 15 longest terminal branches:")
for x, n in d[:15]:
    print(f"  {x:.3f}  {n}")
