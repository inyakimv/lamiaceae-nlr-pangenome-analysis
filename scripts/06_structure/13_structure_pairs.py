#!/usr/bin/env python3
"""Structure comparison of the representative genes (see structure_reps.tsv).
Two steps.
 1) seqs : write one FASTA per representative for the AlphaFold3 server
      python3 13_structure_pairs.py seqs structure_reps.tsv PROTEINS.faa OUTDIR [--extra NAME=FASTA ...]
    PROTEINS.faa: the ANNEVO proteins (headers 'GCA_..._genomic|contig-gN.t1'); --extra adds e.g. AtADR1=AtADR1.fa
 2) pairs: Foldseek TM-score + identity for all pairs inside each set, and for AtADR1 versus the ADR1 reps
      python3 13_structure_pairs.py pairs structure_reps.tsv MODELS_DIR OUT.tsv [ALL_VS_ALL.tsv]
    Runs foldseek createdb/search/convertalis exactly as run_foldseek_comparison_en.sh (alignment-type 1, -e 0.001, -a),
    columns fident and alntmscore as in the 210-pair table; pairs with no hit are listed, not invented.
    MODELS_DIR: AlphaFold3 models named <gene id with | replaced by _>.cif (best model, model_0)
Main TM metric = alntmscore (as in the 210-pair table); tm_min_qt = min(qtmscore, ttmscore) is a sensitivity column (normalised by full length).
Pre-registered criterion (before looking at results): a pair shows structure-sequence decoupling when
TM-score >= 0.5 and sequence identity <= 30%.  Foldseek parameters are those of the earlier 210-pair run."""
import sys, re, os, itertools, subprocess, csv, collections
def norm(x): return re.sub(r"^GC[AF]_\d+\.\d+_", "", x.strip()).replace("_genomic", "")
def safe(x): return re.sub(r"[^A-Za-z0-9_.-]", "_", x)
def read_reps(f):
    d = collections.OrderedDict()
    for l in open(f):
        s, g = l.rstrip("\n").split("\t"); d.setdefault(s, []).append(g)
    return d
def seqs(reps, faa, outdir, extra):
    os.makedirs(outdir, exist_ok=True)
    want = {g for v in read_reps(reps).values() for g in v}
    cur, buf, got = None, [], {}
    def flush():
        if cur and norm(cur) in want: got[norm(cur)] = "".join(buf)
    for l in open(faa):
        if l.startswith(">"):
            flush(); cur = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    flush()
    for e in extra:
        n, f = e.split("="); got[n] = "".join(x.strip() for x in open(f) if not x.startswith(">"))
    miss = [g for g in want if g not in got]
    for g, s in got.items():
        open(os.path.join(outdir, safe(g) + ".fa"), "w").write(f">{safe(g)}\n{s}\n")
    print(f"written {len(got)} fasta files; missing: {miss}")
def run_foldseek(models, workdir):
    """Same calls and parameters as run_foldseek_comparison_en.sh (all-vs-all, TM-align, -e 0.001, -a)."""
    os.makedirs(workdir, exist_ok=True)
    db, res, tmp, tsv = (os.path.join(workdir, x) for x in ("db", "results", "tmp", "all_vs_all.tsv"))
    subprocess.run(["foldseek", "createdb", models, db], check=True)
    os.makedirs(tmp, exist_ok=True)
    subprocess.run(["foldseek", "search", db, db, res, tmp, "--alignment-type", "1", "-e", "0.001", "-a"], check=True)
    subprocess.run(["foldseek", "convertalis", db, db, res, tsv, "--format-output",
                    "query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore"], check=True)
    return tsv
def clean(x):
    x = re.sub(r"\.cif(_[A-Za-z0-9]+)?$", "", x.strip())
    return x
def pairs(reps, models, out, tsv=None):
    sets = read_reps(reps)
    todo = []
    for s_, genes in sets.items(): todo += [(s_, a, b) for a, b in itertools.combinations(genes, 2)]
    adr1 = sets.get("core_ADR1_Clade18", [])
    if os.path.exists(os.path.join(models, "AtADR1.cif")):
        todo += [("ADR1_Ath_vs_Lam", "AtADR1", g) for g in adr1]
    tsv = tsv or run_foldseek(models, "foldseek_work")
    hit = {}
    for l in open(tsv):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        q, t = clean(f[0]), clean(f[1])
        if q == t: continue
        hit[(q, t)] = (float(f[2]), float(f[8]), min(float(f[6]), float(f[7])))
    rows, missing = [], []
    for s_, a, b in todo:
        a_, b_ = safe(a), safe(b)
        v = hit.get((a_, b_)) or hit.get((b_, a_))
        if v is None: missing.append((s_, a, b)); continue
        ident, tm, tmmin = v
        rows.append(dict(set=s_, a=a, b=b, identity=ident, tm_score=tm, tm_min_qt=tmmin,
                         decoupling=(tm >= 0.5 and ident <= 0.30)))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["set", "a", "b", "identity", "tm_score", "tm_min_qt", "decoupling"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("pairs written:", len(rows), "of", len(todo), "->", out)
    if missing: print("pairs NOT reported by Foldseek (no hit at e<=0.001):", missing)
if __name__ == "__main__":
    c = sys.argv[1]
    if c == "seqs":
        ex = [x for x in sys.argv[6:] if "=" in x]
        seqs(sys.argv[2], sys.argv[3], sys.argv[4], ex)
    elif c == "pairs": pairs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else None)
    else: sys.exit(__doc__)
