#!/bin/bash
cd ~/project/11_nbarc_family_tree/zar1
for f in ~/project/00_genomes/*_genomic.fna.gz; do
  g=$(basename $f _genomic.fna.gz)
  zcat $f | makeblastdb -in - -dbtype nucl -title $g -out db_$g > /dev/null 2>&1
  tblastn -query zar1.fasta -db db_$g -evalue 1e-20 -num_threads 8 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' > $g.tbl
  rm -f db_$g.*
  echo "$g done"
done
echo ALL_DONE
