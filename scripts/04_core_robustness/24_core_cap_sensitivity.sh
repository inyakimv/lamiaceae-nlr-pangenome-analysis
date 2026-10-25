#!/bin/bash
# Sensitivity of the conserved core to the clade-size cap (maximal subtrees of <= CAP leaves).
# Re-runs the partition (01) and the clade statistics (02) with several caps, then compares the core clades (25).
# Usage: TREE=nbarc_hmm_plus_refs.treefile S1=TableS1_description.xlsx CLSTR=nbarc_domains_cdhit85.fasta.clstr \
#        SCRIPTS=/path/to/family_tree_scripts bash 24_core_cap_sensitivity.sh OUTDIR [CAPS...]
# Default caps: 150 200 300 400 500. Cap 300 is the published partition and is the reference.
set -uo pipefail
out=$1; shift; caps=${@:-150 200 300 400 500}
: "${TREE:?set TREE}" "${S1:?set S1}" "${CLSTR:?set CLSTR}" "${SCRIPTS:?set SCRIPTS}"
mkdir -p "$out"
for cap in $caps; do
  d="$out/cap$cap"; mkdir -p "$d"
  python3 "$SCRIPTS/01_partition_refs.py" "$TREE" "$S1" "$d/refs" "$cap" > "$d/01.log" 2>&1 || { echo "cap $cap: 01 FAILED, see $d/01.log"; continue; }
  python3 "$SCRIPTS/02_clade_stats.py" "$CLSTR" "$d/refs_clade_members.tsv" "$d/refs_clade_info.csv" "$d/stats" > "$d/02.log" 2>&1 || { echo "cap $cap: 02 FAILED, see $d/02.log"; continue; }
  echo "cap $cap: $(tail -1 "$d/02.log")"
done
python3 "$(dirname "$0")/25_core_cap_compare.py" "$out" 300 "$out/core_cap_sensitivity.tsv"
