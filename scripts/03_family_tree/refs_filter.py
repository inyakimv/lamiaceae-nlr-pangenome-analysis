M, MIN_COV = 248, 0.70
names, seqs, cur = [], [], []
for line in open("refs_hmmalign.a2m"):
    line = line.rstrip()
    if line.startswith(">"):
        if cur: seqs.append("".join(cur)); cur = []
        names.append(line[1:].split()[0])
    else:
        cur.append(line)
seqs.append("".join(cur))
match = ["".join(c for c in s if c == "-" or c.isupper()) for s in seqs]
assert all(len(s) == M for s in match), "unexpected length"
cov = [sum(c != "-" for c in s) / M for s in match]
keep = [i for i, c in enumerate(cov) if c >= MIN_COV]
print(f"References aligned: {len(names)}; kept (coverage >= {MIN_COV}): {len(keep)}")
print("Removed:", [names[i] for i in range(len(names)) if i not in set(keep)])
with open("refs_hmmaligned_cov70.fasta", "w") as o:
    for i in keep: o.write(f">{names[i]}\n{match[i]}\n")
