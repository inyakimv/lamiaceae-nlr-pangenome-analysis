#!/bin/bash
# ==========================================================
# reduce_nbarc_redundancy.sh
#
# Reduces sequence redundancy in the combined NB-ARC domain
# dataset (15,291 sequences across 47 organisms) using CD-HIT,
# collapsing near-identical paralogous copies (e.g., recent
# tandem duplicates within the same species) into a single
# representative sequence per cluster, at an 85% identity
# threshold.
#
# This step is necessary to make family-level phylogenetic
# reconstruction of the Lamiaceae NLR NB-ARC domain
# computationally tractable (alignment and tree-building on
# >15,000 near-identical sequences would be prohibitively slow
# and would add little additional phylogenetic signal, since
# near-identical duplicates barely diverge in sequence).
#
# Requirements: cd-hit (installed via conda, bioconda channel)
#
# Usage:
#   cd ~/project
#   bash scripts/reduce_nbarc_redundancy.sh
#
# Part of the Lamiaceae pan-NLRome family-level phylogeny
# pipeline (paper Figure — NB-ARC domain tree with Arabidopsis
# anchor points).
# ==========================================================

set -euo pipefail

INPUT_FASTA="nbarc_domains_combined.fasta"
OUTPUT_PREFIX="nbarc_domains_cdhit85"
IDENTITY_THRESHOLD=0.85

echo "Input sequences:"
grep -c "^>" "$INPUT_FASTA"

echo "Running CD-HIT at ${IDENTITY_THRESHOLD} identity threshold..."
cd-hit -i "$INPUT_FASTA" -o "$OUTPUT_PREFIX.fasta" \
    -c "$IDENTITY_THRESHOLD" \
    -n 5 \
    -M 8000 \
    -T 8 \
    -d 0

echo ""
echo "=============================================="
echo "Output sequences (after redundancy reduction):"
grep -c "^>" "$OUTPUT_PREFIX.fasta"
echo "=============================================="
echo "Representative sequences saved to: $OUTPUT_PREFIX.fasta"
echo "Cluster membership details saved to: $OUTPUT_PREFIX.fasta.clstr"
