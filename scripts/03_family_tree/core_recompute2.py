import os, glob, subprocess, collections, statistics, sys
import pandas as pd
P = os.path.expanduser("~/project")
tsv = P + "/03_orthofinder/proteomes/OrthoFinder/Results_Jul25/Orthogroups/Orthogroups.tsv"
OUT = "core2"; os.makedirs(OUT, exist_ok=True)
info, members = {}, {}
with open(tsv) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    lam = [i for i, s in enumerate(hdr) if i > 0 and "Athaliana" not in s]
    for line in f:
        c = line.rstrip("\n").split("\t")
        genes, nsp = [], 0
        for i in lam:
            g = [x.strip() for x in c[i].split(",") if x.strip()]
            if g: nsp += 1; genes += g
        if nsp >= 23:
            info[c[0]] = (nsp, len(genes))
            for g in genes: members[g] = c[0]
seqs = {}
for fa in glob.glob(f"{P}/03_orthofinder/proteomes/*.fa*"):
    name = None; buf = []
    for l in open(fa):
        if l.startswith(">"):
            if name in members: seqs[name] = "".join(buf)
            name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if name in members: seqs[name] = "".join(buf)
print(f"OGs (>=23 species): {len(info)}; genes: {len(members)}; sequences found: {len(seqs)}", file=sys.stderr)
with open(f"{OUT}/members.fasta", "w") as o:
    for g, s in seqs.items(): o.write(f">{g}\n{s}\n")
def run(c): subprocess.run(c, shell=True, check=True)
run(f"hmmscan --cpu 8 --noali --domtblout {OUT}/prof.dom {P}/nlr_domain_search/nlr_profiles.hmm {OUT}/members.fasta > /dev/null")
run(f"hmmsearch --cpu 8 --noali --domtblout {OUT}/rxn.dom eds1/rxn.hmm {OUT}/members.fasta > /dev/null")
run(f"hmmsearch --cpu 8 --noali --domtblout {OUT}/rpw8.dom eds1/rpw8.hmm {OUT}/members.fasta > /dev/null")
run(f"blastp -num_threads 8 -query {OUT}/members.fasta -db {P}/11_nbarc_family_tree/refplantnlr/refdb -outfmt '6 qseqid sseqid pident length evalue' -max_target_seqs 1 -evalue 1e-5 > {OUT}/blast.tsv")
ref = pd.read_excel("refplantnlr/TableS1_description.xlsx", sheet_name="RefNLRs")
struct = dict(zip(ref["ID"].astype(str), ref["Structure"].astype(str)))
def coarse(s):
    if s.startswith("T"): return "TNL"
    if s.startswith("R"): return "RNL"
    if s.startswith("C"): return "CNL"
    if s == "NL": return "NL(CNL-lineage)"
    return "other"
tir, rx, rp, bl = {}, {}, {}, {}
for fn, d, col in [("prof.dom", tir, 3), ("rxn.dom", rx, 0), ("rpw8.dom", rp, 0)]:
    for l in open(f"{OUT}/{fn}"):
        if l[0] == "#": continue
        c = l.split()
        if fn == "prof.dom" and c[0] != "TIR": continue
        d[c[col]] = min(d.get(c[col], 1.0), float(c[12]))
for l in open(f"{OUT}/blast.tsv"):
    c = l.split(); bl[c[0]] = (c[1], float(c[2]))
cnt = collections.defaultdict(lambda: collections.Counter())
tops = collections.defaultdict(collections.Counter)
with open(f"{OUT}/core2_members.tsv", "w") as o:
    o.write("OG\tgene\tdomain_class\tbest_RefPlantNLR\tpident\thomology_class\n")
    for g, og in members.items():
        if g not in seqs: continue
        dc = "TNL" if tir.get(g, 1) <= 1e-10 else ("RNL" if rp.get(g, 1) <= 1e-3 else ("CNL-like" if rx.get(g, 1) <= 1e-3 else "NBARC_only"))
        h = bl.get(g)
        hc = coarse(struct.get(h[0], "?")) if h else "no_hit"
        cnt[og][dc] += 1; cnt[og]["hom_" + hc] += 1
        if h: tops[og][h[0]] += 1
        o.write(f"{og}\t{g}\t{dc}\t{h[0] if h else 'NA'}\t{h[1] if h else 'NA'}\t{hc}\n")
cols = ["TNL", "RNL", "CNL-like", "NBARC_only", "hom_CNL", "hom_TNL", "hom_RNL", "hom_NL(CNL-lineage)", "hom_other", "hom_no_hit"]
print("OG\tn_species\tn_genes\t" + "\t".join(cols) + "\ttop_hits")
for og in sorted(info, key=lambda x: -info[x][0]):
    k = cnt[og]
    print(og, info[og][0], info[og][1], *[k[c] for c in cols], ";".join(f"{a}({b})" for a, b in tops[og].most_common(2)), sep="\t")
