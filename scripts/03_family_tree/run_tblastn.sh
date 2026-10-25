#!/bin/bash
f=$1
s=$(basename $f _genomic.fna.gz)
mkdir -p tbn_out tbn_db
if [ ! -s tbn_out/$s.tbl ]; then
  gunzip -c $f | makeblastdb -in - -dbtype nucl -title $s -out tbn_db/$s > /dev/null 2>&1
  tblastn -query tbn_queries.fasta -db tbn_db/$s -evalue 1e-5 -num_threads 8 \
    -outfmt "6 qseqid sseqid sstart send pident length bitscore evalue" > tbn_out/$s.tbl.tmp \
    && mv tbn_out/$s.tbl.tmp tbn_out/$s.tbl
  rm -f tbn_db/$s.n*
fi
echo "listo $s"
