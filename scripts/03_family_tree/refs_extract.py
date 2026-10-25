import re
seqs, order, cur = {}, [], None
for line in open("refplantnlr_AA.fasta"):
    line = line.rstrip()
    if line.startswith(">"):
        cur = line[1:].split()[0]; seqs[cur] = []; order.append(cur)
    else:
        seqs[cur].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
best = {}
for line in open("refs_nbarc.domtbl"):
    if line.startswith("#"): continue
    f = line.split()
    name, ie, es, ee = f[0], float(f[12]), int(f[19]), int(f[20])
    if name not in best or ie < best[name][0]:
        best[name] = (ie, es, ee)
n = 0
with open("refs_nbarc_domains.fasta", "w") as out:
    for name in order:
        if name in best:
            ie, s, e = best[name]
            ident = re.sub(r"[^A-Za-z0-9_.\-]", "_", name)
            out.write(f">REF|{ident}\n{seqs[name][s-1:e]}\n"); n += 1
print(f"Proteins: {len(order)}; with NB-ARC hit: {n}; without: {len(order)-n}")
