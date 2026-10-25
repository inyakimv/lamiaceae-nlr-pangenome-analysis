#!/usr/bin/env python3
"""Classify the NB-ARC genes of the family tree (ANNEVO only) into TNL / RNL / CNL / other and count ADR1 and NRG1 per species.
Usage: python3 34_nlr_classes_tree.py CLSTR PROTEIN_DIR HMMS ADR1_IDS NRG1_IDS OUT_PREFIX [EXCLUDE=substr,substr] [CPU=32]
  CLSTR       cd-hit .clstr of the NB-ARC family tree (members 'SPECIES|GENE')
  PROTEIN_DIR annevo_proteins/ (files SPECIES.annevo_proteins.fasta)
  HMMS        comma-separated .hmm files, e.g. nlr_profiles.hmm,rpw8.hmm,rxn.hmm (TIR, TIR_2, NB-ARC, RPW8, Rx_N are taken from their NAME)
  ADR1_IDS, NRG1_IDS : text files with one 'SPECIES|GENE' per line (first token is used)
Rules (same for every gene, written in the table header):
  domain hit = hmmsearch independent E-value <= 1e-5 (domtblout i-Evalue); N-terminal = domain starts before the first NB-ARC hit.
  TNL (loose)  : TIR or TIR_2 hit, N-terminal.      TNL (strict): TIR (PF01582) or TIR_2 hit with i-E <= 1e-10 and >= 60 % of the HMM, N-terminal.
  RNL: RPW8 hit N-terminal.   CNL: Rx_N hit N-terminal (and no TIR/RPW8).   other: none of these (includes CC-less / truncated genes).
  Priority if several: TNL > RNL > CNL.
Outputs: OUT_PREFIX_per_species.tsv, OUT_PREFIX_gene_classes.tsv (one row per gene), OUT_PREFIX_hmmsearch.domtblout"""
import sys, re, os, subprocess, tempfile, collections
clstr, pdir, hmms, adr_f, nrg_f, outp = sys.argv[1:7]
excl = [e for e in (sys.argv[7] if len(sys.argv) > 7 else "").split(",") if e]
cpu = sys.argv[8] if len(sys.argv) > 8 else "32"
E_HIT, E_STRICT, COV_STRICT = 1e-5, 1e-10, 0.6

# 1. members of the tree
members = []
for l in open(clstr):
    m = re.search(r">(.+?)\.\.\.", l)
    if m and "|" in m.group(1): members.append(m.group(1))
members = [g for g in members if not any(e in g for e in excl)]
by_sp = collections.defaultdict(list)
for g in members: by_sp[g.split("|")[0]].append(g.split("|", 1)[1])
print("genes in tree:", len(members), "| species with ANNEVO ids:", len(by_sp))

# 2. protein file for hmmsearch
skipped = set(); tmp = tempfile.mkdtemp(); allfa = os.path.join(tmp, "tree_proteins.faa"); miss = 0
with open(allfa, "w") as out:
    for sp, genes in by_sp.items():
        f = os.path.join(pdir, sp + ".annevo_proteins.fasta")
        if not os.path.exists(f): print("skipping (no ANNEVO protein file, reference gene set):", sp, len(genes), "genes"); skipped.add(sp); continue
        want = set(genes); keep = False; seen = set()
        for l in open(f):
            if l.startswith(">"):
                gid = l[1:].split()[0]; keep = gid in want
                if keep: out.write(">%s|%s\n" % (sp, gid)); seen.add(gid)
            elif keep: out.write(l)
        miss += len(want - seen)
members = [g for g in members if g.split("|")[0] not in skipped]
print("genes analysed:", len(members), "| genes without sequence:", miss)

# 3. hmmsearch with all profiles
hmm_all = os.path.join(tmp, "all.hmm")
with open(hmm_all, "w") as o:
    for h in hmms.split(","): o.write(open(h).read().rstrip("\n") + "\n")
names = [l.split()[1] for l in open(hmm_all) if l.startswith("NAME")]
print("profiles:", names)
def role(n):
    nl = n.lower()
    if nl in ("tir", "tir_2"): return "TIR"
    if "nb-arc" in nl or nl == "nb_arc": return "NB"
    if "rpw8" in nl: return "RPW8"
    if re.search(r"rx_?n", nl): return "RXN"
    return None
roles = {n: role(n) for n in names}
for need in ("TIR", "NB", "RPW8", "RXN"):
    if need not in roles.values(): sys.exit("ERROR: no profile for %s among %s" % (need, names))
dom = outp + "_hmmsearch.domtblout"
subprocess.run(["hmmsearch", "--cpu", cpu, "-E", "1e-3", "--domE", "1e-3", "--domtblout", dom, hmm_all, allfa],
               stdout=subprocess.DEVNULL, check=True)

# 4. parse (hmmsearch orientation: col1 = protein, col4 = profile, col6 = profile length)
hits = collections.defaultdict(list)
for l in open(dom):
    if l.startswith("#"): continue
    c = l.split()
    prof, ie = c[3], float(c[12])
    if ie > E_HIT or roles.get(prof) is None: continue
    cov = (int(c[16]) - int(c[15]) + 1) / int(c[5])
    hits[c[0]].append((roles[prof], prof, ie, int(c[17]), cov))

def norm(x): return x.replace("_genomic|", "|", 1)   # ID files may lack the "_genomic" suffix of the species name
def ids(f): return {norm(l.split()[0]) for l in open(f) if l.strip()}
ADR1, NRG1 = ids(adr_f), ids(nrg_f)
rows = []; tab = collections.defaultdict(collections.Counter)
for g in members:
    sp, gid = g.split("|", 1); h = hits.get(g, [])
    nb = [x[3] for x in h if x[0] == "NB"]; start = min(nb) if nb else None
    def nterm(x): return start is not None and x[3] < start
    tir = [x for x in h if x[0] == "TIR" and nterm(x)]
    tnl_loose = bool(tir)
    tnl_strict = any(x[2] <= E_STRICT and x[4] >= COV_STRICT for x in tir)
    rnl = any(x[0] == "RPW8" and nterm(x) for x in h)
    cnl = any(x[0] == "RXN" and nterm(x) for x in h)
    cls = "TNL" if tnl_loose else "RNL" if rnl else "CNL" if cnl else "other"
    isadr, isnrg = norm(g) in ADR1, norm(g) in NRG1
    rows.append((sp, gid, cls, int(tnl_strict), int(isadr), int(isnrg), int(bool(nb))))
    t = tab[sp]; t["NBARC_in_tree"] += 1; t[cls] += 1; t["TNL_strict"] += int(tnl_strict); t["ADR1"] += int(isadr); t["NRG1"] += int(isnrg)
    t["no_NBARC_hit"] += int(not nb)
with open(outp + "_gene_classes.tsv", "w") as o:
    o.write("species\tgene\tclass\tTNL_strict\tADR1\tNRG1\thas_NBARC_hit\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
cols = ["NBARC_in_tree", "TNL", "TNL_strict", "RNL", "CNL", "other", "ADR1", "NRG1", "no_NBARC_hit"]
with open(outp + "_per_species.tsv", "w") as o:
    o.write("species\t" + "\t".join(cols) + "\n")
    for sp in sorted(tab): o.write(sp + "\t" + "\t".join(str(tab[sp][c]) for c in cols) + "\n")
    tot = collections.Counter()
    for sp in tab: tot.update(tab[sp])
    o.write("TOTAL\t" + "\t".join(str(tot[c]) for c in cols) + "\n")
n = len(tab)
print("ADR1 ids read: %d, NRG1 ids read: %d (genes matched: %d, %d)" % (len(ADR1), len(NRG1), sum(r[4] for r in rows), sum(r[5] for r in rows)))
print("species: %d | with TNL(loose): %d | with TNL(strict): %d | with ADR1: %d | with NRG1: %d | with RNL: %d" % (
    n, sum(tab[s]["TNL"] > 0 for s in tab), sum(tab[s]["TNL_strict"] > 0 for s in tab), sum(tab[s]["ADR1"] > 0 for s in tab),
    sum(tab[s]["NRG1"] > 0 for s in tab), sum(tab[s]["RNL"] > 0 for s in tab)))
print(open(outp + "_per_species.tsv").read().splitlines()[-1])
