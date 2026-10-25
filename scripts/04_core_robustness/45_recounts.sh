#!/usr/bin/env bash
# Recounts for the 45-species dataset (Paraleonurus_japonicus excluded). Writes only to ~/structure_test/recounts_45/
# Usage: cd ~/structure_test && bash 45_recounts.sh
set -u
EXC=Paraleonurus; P=$HOME/project; S=$HOME/structure_test; O=$S/recounts_45; mkdir -p $O; cd $S
CLSTR=$P/11_nbarc_family_tree/nbarc_domains_cdhit85.fasta.clstr
FT=$P/11_nbarc_family_tree/family_tree_out

echo "### 1. Saturation (3rd-codon p-distance) without $EXC" | tee $O/01_saturation.log
mkdir -p $O/aln; rm -f $O/aln/*
SRC=codon_aln; ls $SRC/*.trim.fasta >/dev/null 2>&1 || SRC=hyphy_out
for f in $SRC/*.trim.fasta $SRC/*.hy.fasta; do [ -f "$f" ] || continue
  awk -v x="$EXC" '/^>/{keep=($0!~x)} keep' "$f" > $O/aln/$(basename $f); done
python3 23_saturation_table.py $O/saturation_45.tsv $O/aln/*.fasta 2>&1 | tee -a $O/01_saturation.log | tail -30

echo; echo "### 2. Core vs clade-size cap and thresholds (45 species)" | tee $O/02_core.log
for cap in 150 200 300 400 500; do d=out_cap_test/cap$cap; [ -f $d/refs_clade_members.tsv ] || continue
  echo "--- cap $cap (min 35 sp, max 3 copies)" | tee -a $O/02_core.log
  python3 28_clade_stats_exclude.py $CLSTR $d/refs_clade_members.tsv $d/refs_clade_info.csv $O/core_cap$cap $EXC 2>&1 | tee -a $O/02_core.log | tail -8; done
for ms in 30 33 35 38 40; do for mc in 2 3 4; do
  echo "--- cap 300, min_species $ms, max_copies $mc" | tee -a $O/02_core.log
  python3 28_clade_stats_exclude.py $CLSTR out_cap_test/cap300/refs_clade_members.tsv out_cap_test/cap300/refs_clade_info.csv $O/core_t${ms}_c${mc} $EXC $ms $mc 2>&1 | tail -3 | tee -a $O/02_core.log; done; done

echo; echo "### 3. Clade76 N-terminal (MADA-like) check without $EXC" | tee $O/03_mada.log
python3 $P/scripts/family_tree/05_mada_check.py $FT/Clade76_genes.tsv $P/annevo_proteins $O/Clade76_mada_all.tsv $P/refplantnlr/RefPlantNLR_proteins.fasta 2>&1 | tail -15 | tee -a $O/03_mada.log
head -1 $O/Clade76_mada_all.tsv > $O/Clade76_mada_45.tsv; grep -v "$EXC" $O/Clade76_mada_all.tsv | tail -n +2 >> $O/Clade76_mada_45.tsv
echo "rows all/45: $(($(wc -l < $O/Clade76_mada_all.tsv)-1)) / $(($(wc -l < $O/Clade76_mada_45.tsv)-1))" | tee -a $O/03_mada.log
head -3 $O/Clade76_mada_45.tsv | cut -c1-250 | tee -a $O/03_mada.log

echo; echo "### 4. TIR-NLRs in the family tree" | tee $O/04_tir_eds1.log
python3 - <<'PY' 2>&1 | tee -a $O/04_tir_eds1.log
import glob, pandas as pd
for f in sorted(glob.glob("tnl_census*gene*.tsv")):
    d = pd.read_csv(f, sep="\t"); print(f, d.shape, list(d.columns)[:12])
    c = [x for x in d.columns if "in_tree" in x]
    if c:
        t = d[d[c[0]].astype(str).str.lower().isin(["true","1","yes"])]
        print("in tree:", len(t)); print(t.head(10).to_string()[:1500])
PY
echo; echo "[ok] results in $O"
