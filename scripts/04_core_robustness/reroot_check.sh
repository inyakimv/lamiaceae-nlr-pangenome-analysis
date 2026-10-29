#!/bin/bash
# Re-rooting control for the conserved core.
# Re-partitions the NB-ARC tree after rooting on the stem of Clade 18 and compares the 85 clades and the core
# clades with the original (midpoint-rooted) partition.
#
# Usage: bash reroot_check.sh WORKDIR
#   WORKDIR must contain: nbarc_hmm_plus_refs.treefile, nbarc_domains_cdhit85.fasta.clstr,
#   refplantnlr/TableS1_description.xlsx and family_tree_out/refs_clade_members.tsv (original partition).
# Environment: SCRIPTS = folder with 01b_partition_refs_rerooted.py and compare_core.py (default: this repository's scripts);
#              02_clade_stats*.py (03_family_tree/02_clade_stats.py or the 45-species variant) must be in WORKDIR or on PATH.
# Output: REROOT_clade_members.tsv, REROOT_clade_info.csv, reroot_compare.txt (expected: 85/85 clades identical, same four core clades).
set -euo pipefail
WORKDIR=${1:?usage: bash reroot_check.sh WORKDIR}
SCRIPTS=${SCRIPTS:-$(cd "$(dirname "$0")/.." && pwd)}
cd "$WORKDIR"
TREE=nbarc_hmm_plus_refs.treefile
TABLES1=refplantnlr/TableS1_description.xlsx
CLSTR=nbarc_domains_cdhit85.fasta.clstr
ORIG=family_tree_out/refs
CLADE=Clade18
CAP=300
STATS=${STATS:-02_clade_stats_45.py}
[ -f "$STATS" ] || STATS="$SCRIPTS/03_family_tree/02_clade_stats.py"
awk -F'\t' -v c="$CLADE" '$1==c{print $2}' ${ORIG}_clade_members.tsv > clade18_leaves.txt
echo "leaves of $CLADE: $(wc -l < clade18_leaves.txt)"
REROOT_CLADE_MEMBERS=clade18_leaves.txt python3 "$SCRIPTS/03_family_tree/01b_partition_refs_rerooted.py" "$TREE" "$TABLES1" REROOT $CAP
python3 "$STATS" "$CLSTR" REROOT_clade_members.tsv REROOT_clade_info.csv REROOT
python3 "$SCRIPTS/04_core_robustness/compare_core.py" "$ORIG" REROOT | tee reroot_compare.txt
