#!/usr/bin/env bash
# Cross-validation ANNEVO+HMMER vs NLR-Annotator with 45 species (Paraleonurus_japonicus excluded).
# Works in a separate folder with symlinks; the original inputs and the TFM outputs are NOT touched.
# Usage: cd ~/structure_test && bash 46_xval_45.sh
set -u
EXC=Paraleonurus; P=$HOME/project; W=$HOME/structure_test/xval_45; M=$P/scripts/method_comparison
rm -rf $W; mkdir -p $W/02_nlr_annotation $W/00_genomes $W/nlr_domain_search/per_species; cd $W
for d in $P/02_nlr_annotation/*; do b=$(basename "$d"); [ "$b" = nlr_summary.csv ] && continue; case "$b" in *$EXC*) ;; *) ln -s "$d" 02_nlr_annotation/"$b";; esac; done
head -1 $P/02_nlr_annotation/nlr_summary.csv > 02_nlr_annotation/nlr_summary.csv
grep -v "$EXC" $P/02_nlr_annotation/nlr_summary.csv | tail -n +2 >> 02_nlr_annotation/nlr_summary.csv
for d in $P/00_genomes/*; do b=$(basename "$d"); case "$b" in *$EXC*) ;; *) ln -s "$d" 00_genomes/"$b";; esac; done
for f in $P/nlr_domain_search/per_species/*; do b=$(basename "$f"); case "$b" in *$EXC*) ;; *) ln -s "$f" nlr_domain_search/per_species/"$b";; esac; done
echo "species dirs: $(ls 02_nlr_annotation | grep -c '_genomic$')  domtblout: $(ls nlr_domain_search/per_species/*.domtblout | wc -l)  summary rows: $(($(wc -l < 02_nlr_annotation/nlr_summary.csv)-1))"
for s in compare_nlr_methods_refined gene_size_distribution; do
  sed -E 's/\b46 (Lamiaceae )?species/45 \1species/g' $M/$s.py > ${s}_45.py
done
echo "=== compare_nlr_methods_refined (45)"; python3 compare_nlr_methods_refined_45.py 2>&1 | tail -45 | tee compare_45.log
echo; echo "=== gene_size_distribution (45)"; python3 gene_size_distribution_45.py 2>&1 | tail -45 | tee gene_size_45.log
echo; echo "[ok] outputs in $W"; ls $W | grep -v "^02_\|^00_\|^nlr_domain"
