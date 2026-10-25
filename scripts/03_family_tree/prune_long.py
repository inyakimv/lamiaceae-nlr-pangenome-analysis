import statistics
from ete3 import Tree
t = Tree("nbarc_hmm.treefile", format=0)
leaves = t.get_leaves()
med = statistics.median(l.dist for l in leaves)
drop = {l.name for l in leaves if l.dist > 5 * med}
open("ids_longbranch.txt", "w").write("\n".join(sorted(drop)) + "\n")
keep = True; n_in = n_out = 0
with open("nbarc_hmmaligned_cov70.fasta") as f, open("nbarc_hmm_noLB.fasta", "w") as o:
    for line in f:
        if line.startswith(">"):
            keep = line[1:].strip().split()[0] not in drop
            n_in += 1; n_out += keep
        if keep:
            o.write(line)
print(f"threshold {5*med:.3f}; removed {len(drop)}; kept {n_out} of {n_in}")
print("by species:"); 
from collections import Counter
for k, v in Counter(n.split('|')[0] for n in drop).most_common(): print(f"  {v}  {k}")
