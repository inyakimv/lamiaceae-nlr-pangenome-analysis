import re, subprocess
LABELS = ["Knitens", "Cbraunii", "Mpolymorpha", "Ppatens", "Smoellendorffii"]
M, MIN_COV = 248, 0.70
def read_fasta(p):
    d, k = {}, None
    for l in open(p):
        l = l.rstrip()
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l)
    return {a: "".join(b) for a, b in d.items()}
tot = {}
with open("outgroups_domains.fasta", "w") as out:
    for lab in LABELS:
        seqs = read_fasta(f"{lab}.faa"); best = {}
        for l in open(f"{lab}.nb.domtbl"):
            if l.startswith("#"): continue
            f = l.split(); ie = float(f[12])
            if f[0] not in best or ie < best[f[0]][0]: best[f[0]] = (ie, int(f[19]), int(f[20]))
        sel = []
        for l in open(f"{lab}_nbarc_hits.tsv"):
            p, e, tir = l.split("\t"); tir = tir.strip()
            if float(e) < 1e-10 and tir == "-": sel.append(p)
        for p in sel:
            _, s, e = best[p]
            out.write(f">OUT|{lab}|{re.sub(r'[^A-Za-z0-9_.-]', '_', p)}\n{seqs[p][s-1:e]}\n")
        tot[lab] = len(sel)
print("Selected (E<1e-10, no TIR):", tot)
subprocess.run(["hmmalign", "--amino", "--trim", "--outformat", "A2M", "-o", "outgroups.a2m",
                "../nbarc_only.hmm", "outgroups_domains.fasta"], check=True)
names, seqs, cur = [], [], []
for l in open("outgroups.a2m"):
    l = l.rstrip()
    if l.startswith(">"):
        if cur: seqs.append("".join(cur)); cur = []
        names.append(l[1:].split()[0])
    else: cur.append(l)
seqs.append("".join(cur))
match = ["".join(c for c in s if c == "-" or c.isupper()) for s in seqs]
assert all(len(s) == M for s in match)
cov = [sum(c != "-" for c in s) / M for s in match]
kept = {lab: 0 for lab in LABELS}
with open("outgroups_hmmaligned_cov70.fasta", "w") as o:
    for n, s, c in zip(names, match, cov):
        if c >= MIN_COV: o.write(f">{n}\n{s}\n"); kept[n.split("|")[1]] += 1
print("Kept (coverage >= 0.70):", kept)
