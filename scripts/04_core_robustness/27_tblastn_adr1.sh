#!/usr/bin/env bash
# Search for ADR1 in a genome that has no ADR1 in the NB-ARC family tree (e.g. Paraleonurus japonicus, GCA_030762865.1).
# Usage: bash 27_tblastn_adr1.sh GENOME.fasta [QUERY.fasta] [OUT_PREFIX]
#   QUERY defaults to AtADR1_Q9FW44.fasta (same folder). Adding a Lamiaceae ADR1 (e.g. from a Lamioideae species) to the query file
#   makes the search more sensitive; keep AtADR1 in it as well.
# Output: OUT_PREFIX.tblastn.tsv (all hits, e-value <= 1e-5) and OUT_PREFIX.loci.tsv (hits merged into loci, with best query coverage).
set -euo pipefail
G="$1"; Q="${2:-$(dirname "$0")/AtADR1_Q9FW44.fasta}"; P="${3:-adr1_search}"
command -v tblastn >/dev/null || { echo "tblastn not found: activate an env with BLAST+ (conda install -c bioconda blast)"; exit 1; }
[ -s "$G" ] || { echo "genome not found: $G"; exit 1; }
DB="${P}_db"
makeblastdb -in "$G" -dbtype nucl -out "$DB" >/dev/null
tblastn -query "$Q" -db "$DB" -evalue 1e-5 -max_hsps 20 \
  -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send evalue bitscore" -out "${P}.tblastn.tsv"
echo "HSPs: $(wc -l < ${P}.tblastn.tsv)"
python3 - "${P}.tblastn.tsv" "${P}.loci.tsv" <<'PY'
import sys, collections
hits=[l.split("\t") for l in open(sys.argv[1]).read().splitlines() if l]
by=collections.defaultdict(list)
for h in hits:
    s,e=sorted((int(h[7]),int(h[8]))); by[h[1]].append((s,e,h))
out=[]
for c,v in by.items():
    v.sort(); cur=None
    for s,e,h in v:
        if cur and s<=cur[1]+20000:
            cur[1]=max(cur[1],e); cur[2].append(h)
        else:
            if cur: out.append((c,)+tuple(cur))
            cur=[s,e,[h]]
    if cur: out.append((c,)+tuple(cur))
rows=[]
for c,s,e,hs in out:
    cov=set()
    for h in hs: cov.update(range(int(h[4]),int(h[5])+1))
    qlen=int(hs[0][6]); best=max(float(h[10]) for h in hs); pid=max(float(h[2]) for h in hs)
    rows.append((best,c,s,e,len(hs),round(100*len(cov)/qlen,1),pid))
rows.sort(reverse=True)
with open(sys.argv[2],"w") as o:
    o.write("contig\tstart\tend\tn_hsps\tpct_query_covered\tmax_pident\tbest_bitscore\n")
    for b,c,s,e,n,cv,pid in rows: o.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\n"%(c,s,e,n,cv,pid,b))
print(open(sys.argv[2]).read()[:3000])
PY
