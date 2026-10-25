#!/usr/bin/env bash
# ML tree of the ADR1/RNL clade (Clade18) of the NB-ARC family tree, rooted with ZAR1 and NRC0.
# This reproduces the commands that produced family_tree_out/Clade18_tree.treefile (17 Oct 2026):
#   1. sequences : 06_clade_fasta.py (ANNEVO proteins of the clade + RefPlantNLR ADR1/NRG1 references + ZAR1, NRC0 outgroups)
#   2. alignment : MAFFT --auto on full-length proteins (MAFFT v7.490, base environment)
#   3. tree      : IQ-TREE 3.1.4 (conda env "iqtree"), ModelFinder over LG,WAG,JTT,Q.plant,Q.pfam, 1000 UFBoot + 1000 SH-aLRT
#                  (Q.PLANT+R5 was selected by BIC)
# Usage: P=/path/to/project THREADS=8 bash 10b_clade18_adr1_tree.sh
# Needs family_tree_out/Clade18_genes.tsv (clade members, from 04_inspect_clade.py / the clade partition).
set -euo pipefail
P=${P:-$HOME/project}; THREADS=${THREADS:-8}
cd "$P/11_nbarc_family_tree"
python3 "$P/scripts/family_tree/06_clade_fasta.py" family_tree_out/Clade18_genes.tsv "$P/annevo_proteins" \
        "$P/refplantnlr/RefPlantNLR_proteins.fasta" 'ADR1|NRG1' family_tree_out/Clade18_tree.fasta '^ZAR1$|^NRC0$'
cd family_tree_out
mafft --auto --thread "$THREADS" Clade18_tree.fasta > Clade18_aln.fasta
iqtree3 -s Clade18_aln.fasta -m MFP -mset LG,WAG,JTT,Q.plant,Q.pfam -mrate G4 -B 1000 -alrt 1000 -T "$THREADS" \
        -o OUT_ZAR1,OUT_NRC0 --prefix Clade18_tree -redo
python3 "$P/scripts/family_tree/07_fig_clade_tree.py" Clade18_tree.treefile Fig_S_Clade18_tree_collapsed collapsed
python3 "$P/scripts/family_tree/07_fig_clade_tree.py" Clade18_tree.treefile Fig_S_Clade18_tree_full full
