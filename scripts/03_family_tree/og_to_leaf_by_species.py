import os, subprocess, collections, glob
P = os.path.expanduser("~/project")
tsv = P + "/03_orthofinder/proteomes/OrthoFinder/Results_Jul25/Orthogroups/Orthogroups.tsv"
CORE10 = {"OG0000003","OG0000006","OG0000061","OG0000115","OG0000180","OG0000216","OG0000242","OG0000291","OG0000363","OG0000635"}
og_of = {}
with open(tsv) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    lam = [i for i, s in enumerate(hdr) if i > 0 and "Athaliana" not in s]
    for line in f:
        c = line.rstrip("\n").split("\t")
        nsp = sum(1 for i in lam if c[i].strip())
        if nsp >= 23 or c[0] in CORE10:
            for i in range(1, len(c)):
                for g in [x.strip() for x in c[i].split(",") if x.strip()]:
                    og_of[g] = c[0]
def fasta(path):
    d, name, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if name: d[name] = "".join(buf)
            name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if name: d[name] = "".join(buf)
    return d
q = fasta("core2/all_members.fasta")
# add any OG member missing from all_members
leaf = fasta("nbarc_domains_unaligned.fasta")
def lkey(h): return "ATH" if h.startswith("Athaliana") else h.split("|")[0][:15]
def qkey(h): return "ATH" if h.startswith("GCF_000001735.4") else h.split("|")[0][:15]
lsp, qsp = collections.defaultdict(dict), collections.defaultdict(dict)
for h, s in leaf.items(): lsp[lkey(h)][h] = s
for h, s in q.items():
    if h in og_of: qsp[qkey(h)][h] = s
os.makedirs("core2/sp", exist_ok=True)
print(f"queries with OG: {sum(len(v) for v in qsp.values())}; species: {len(qsp)}")
out = open("core2/og_to_leaf_species.tsv", "w")
out.write("OG\tgene\tleaf\tpident\taln_len\tqlen\tslen\tevalue\n")
for k, qs in qsp.items():
    if k not in lsp: print("no leaves for", k); continue
    with open(f"core2/sp/{k}.db.fa", "w") as o:
        for h, s in lsp[k].items(): o.write(f">{h}\n{s}\n")
    with open(f"core2/sp/{k}.q.fa", "w") as o:
        for h, s in qs.items(): o.write(f">{h}\n{s}\n")
    subprocess.run(f"makeblastdb -in core2/sp/{k}.db.fa -dbtype prot -out core2/sp/{k}.db > /dev/null", shell=True, check=True)
    r = subprocess.run(f"blastp -query core2/sp/{k}.q.fa -db core2/sp/{k}.db -outfmt '6 qseqid sseqid pident length qlen slen evalue' -max_target_seqs 5 -evalue 1e-10 2>/dev/null", shell=True, capture_output=True, text=True, check=True)
    for l in r.stdout.splitlines():
        c = l.split("\t"); out.write(f"{og_of[c[0]]}\t" + "\t".join(c) + "\n")
out.close()
print("done")
