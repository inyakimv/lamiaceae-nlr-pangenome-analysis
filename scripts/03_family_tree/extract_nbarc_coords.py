#!/usr/bin/env python3
"""
extract_nbarc_coords.py

Extracts NB-ARC domain coordinates from HMMER --domtblout output files
for 47 organisms (46 Lamiaceae species + Arabidopsis thaliana), keeping
only "complete-like" NLR proteins (those with both a detected NB-ARC
domain and at least one LRR domain), following the same criterion used
throughout this study for cross-validating NLR-Annotator against
ANNEVO+HMMER.

For each protein, if multiple NB-ARC domain hits are reported, only the
best-scoring hit (lowest independent e-value) is retained.

Input:
    - nlr_domain_search/per_species/*.domtblout  (46 Lamiaceae species)
    - nlr_domain_search/arabidopsis/Athaliana_genomic_ANNEVO.domtblout

Output:
    - nbarc_coords_summary.json: per-species dict of
      {protein_id: [ali_start, ali_end]} for complete-like NLR proteins
    - Console summary of protein counts per species

Usage:
    cd ~/project
    python3 scripts/extract_nbarc_coords.py

Part of the Lamiaceae pan-NLRome family-level phylogeny pipeline
(paper Figure — NB-ARC domain tree with Arabidopsis anchor points).
"""
import glob
import os
import json


def extract_nbarc_hits(domtblout_path):
    """Parses a single HMMER --domtblout file and returns a dict of
    {protein_id: (ali_start, ali_end)} for proteins classified as
    complete-like NLRs (NB-ARC domain present, with the best-scoring
    NB-ARC hit retained per protein, AND at least one LRR-family domain
    also detected in the same protein)."""
    nbarc_hits = {}  # protein_id -> (best_i_evalue, ali_start, ali_end)
    lrr_proteins = set()

    with open(domtblout_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 19:
                continue

            domain = fields[0]
            protein = fields[3]
            i_evalue = float(fields[12])
            ali_from = int(fields[17])
            ali_to = int(fields[18])

            if domain == 'NB-ARC':
                if protein not in nbarc_hits or i_evalue < nbarc_hits[protein][0]:
                    nbarc_hits[protein] = (i_evalue, ali_from, ali_to)
            elif domain.startswith('LRR'):
                lrr_proteins.add(protein)

    complete_like = {
        protein: (start, end)
        for protein, (_, start, end) in nbarc_hits.items()
        if protein in lrr_proteins
    }
    return complete_like


def main():
    domtblout_files = glob.glob('nlr_domain_search/per_species/*.domtblout')
    arabidopsis_file = 'nlr_domain_search/arabidopsis/Athaliana_genomic_ANNEVO.domtblout'
    domtblout_files.append(arabidopsis_file)

    print(f"Total domtblout files to process: {len(domtblout_files)}")

    records_to_extract = {}
    for domfile in domtblout_files:
        species_key = (
            os.path.basename(domfile)
            .replace('.domtblout', '')
            .replace('_ANNEVO', '')
        )
        complete_like = extract_nbarc_hits(domfile)
        records_to_extract[species_key] = complete_like
        print(f"{species_key}: {len(complete_like)} complete-like NB-ARC proteins")

    with open('nbarc_coords_summary.json', 'w') as f:
        json.dump(
            {k: {p: list(v) for p, v in d.items()} for k, d in records_to_extract.items()},
            f,
            indent=2,
        )

    total = sum(len(d) for d in records_to_extract.values())
    print(f"\nSummary saved to: nbarc_coords_summary.json")
    print(f"Total NB-ARC proteins to extract (47 organisms): {total}")


if __name__ == "__main__":
    main()
