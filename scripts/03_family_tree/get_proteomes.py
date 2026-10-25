import gzip, json, re, shutil, subprocess, time, urllib.parse, urllib.request
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TARGETS = {  # label: assembly accession
    "Knitens": "GCA_000708835.1", "Cbraunii": "GCA_003427395.1", "Mendlicherianum": "GCA_009602735.1",
    "Mpolymorpha": "GCA_037833965.1", "Ppatens": "GCF_000002425.5",
    "Smoellendorffii": "GCF_000143415.4", "Atrichopoda": "GCF_000471905.2"}
def get(u):
    with urllib.request.urlopen(u, timeout=60) as r: return r.read()
res = {}
for lab, acc in TARGETS.items():
    ids = json.loads(get(f"{B}esearch.fcgi?db=assembly&term={urllib.parse.quote(acc + '[Assembly Accession]')}&retmode=json"))["esearchresult"]["idlist"]
    time.sleep(0.4)
    d = json.loads(get(f"{B}esummary.fcgi?db=assembly&id={ids[0]}&retmode=json"))["result"][ids[0]]
    time.sleep(0.4)
    ftp = (d["ftppath_refseq"] if acc.startswith("GCF") else d["ftppath_genbank"]).replace("ftp://", "https://")
    base = ftp.rstrip("/").split("/")[-1]
    url = f"{ftp}/{base}_protein.faa.gz"
    try:
        urllib.request.urlretrieve(url, f"{lab}.faa.gz")
        with gzip.open(f"{lab}.faa.gz") as i, open(f"{lab}.faa", "wb") as o: shutil.copyfileobj(i, o)
        n = sum(1 for l in open(f"{lab}.faa") if l.startswith(">"))
        print(f"{lab:16s} {acc}  proteins: {n}")
        res[lab] = n
    except Exception as e:
        print(f"{lab:16s} {acc}  NO PROTEIN FILE ({type(e).__name__}: {e})")

def best(dom, thr):
    b = {}
    for l in open(dom):
        if l.startswith("#"): continue
        f = l.split(); ie = float(f[12])
        if ie < thr and (f[0] not in b or ie < b[f[0]]): b[f[0]] = ie
    return b
print("\nlabel            proteins  NB-ARC(<1e-5)  NB-ARC(<1e-10)  with TIR(<1e-3)  NB-ARC without TIR")
for lab in res:
    for name, hmm in (("nb", "../nbarc_only.hmm"), ("tir", "tir_only.hmm"), ("tir2", "tir2_only.hmm")):
        subprocess.run(["hmmsearch", "--noali", "--domtblout", f"{lab}.{name}.domtbl", hmm, f"{lab}.faa"],
                       stdout=subprocess.DEVNULL, check=True)
    nb = best(f"{lab}.nb.domtbl", 1e-5); nb10 = best(f"{lab}.nb.domtbl", 1e-10)
    tir = set(best(f"{lab}.tir.domtbl", 1e-3)) | set(best(f"{lab}.tir2.domtbl", 1e-3))
    notir = [p for p in nb if p not in tir]
    open(f"{lab}_nbarc_hits.tsv", "w").write("".join(f"{p}\t{nb[p]:.1e}\t{'TIR' if p in tir else '-'}\n" for p in sorted(nb)))
    print(f"{lab:16s} {res[lab]:8d} {len(nb):14d} {len(nb10):15d} {len(set(nb)&tir):16d} {len(notir):18d}")
