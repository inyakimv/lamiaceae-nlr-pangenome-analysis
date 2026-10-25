#!/usr/bin/env python3
"""
extract_nbarc_sequences.py

Extracts the actual amino acid sequences of the NB-ARC domain for all
complete-like NLR proteins across 47 organisms (46 Lamiaceae species +
Arabidopsis thaliana), using the domain coordinates previously computed
by extract_nbarc_coords.py (saved in nbarc_coords_summary.json).

For each protein, the corresponding full-length sequence is read from
the ANNEVO protein FASTA file for that species, and only the NB-ARC
domain region (ali_start:ali_end, 1-based inclusive coordinates from
HMMER) is extracted. Sequence headers are standardized to
"<species_key>|<protein_id>" to keep species origin traceable in the
combined alignment.

Input:
    - nbarc_coords_summary.json (from extract_nbarc_coords.py)
    - annevo_proteins/<species>.annevo_proteins.fasta (46 Lamiaceae)
    - arabidopsis_outgroup/annevo/Athaliana_genomic.annevo_proteins.fasta

Output:
    - nbarc_domains_combined.fasta: single multi-FASTA file with all
      extracted NB-ARC domain sequences across the 47 organisms,
      ready for redundancy reduction (CD-HIT) and alignment (MAFFT).

Usage:
    cd ~/project
    python3 scripts/extract_nbarc_sequences.py

Part of the Lamiaceae pan-NLRome family-level phylogeny pipeline
(paper Figure — NB-ARC domain tree with Arabidopsis anchor points).
"""
import json
import os
from Bio import SeqIO


def get_proteome_path(species_key):
    """Returns the path to the ANNEVO protein FASTA file for a given
    species key, handling the special case of Arabidopsis, whose
    proteome lives in a different directory than the 46 Lamiaceae
    species."""
    if species_key == "Athaliana_genomic":
        return "arabidopsis_outgroup/annevo/Athaliana_genomic.annevo_proteins.fasta"
    return f"annevo_proteins/{species_key}.annevo_proteins.fasta"


def main():
    with open("nbarc_coords_summary.json") as f:
        coords_by_species = json.load(f)

    total_written = 0
    total_missing_proteome = 0
    total_protein_not_found = 0

    with open("nbarc_domains_combined.fasta", "w") as out_fasta:
        for species_key, protein_coords in coords_by_species.items():
            proteome_path = get_proteome_path(species_key)

            if not os.path.exists(proteome_path):
                print(f"[!] Proteome not found for {species_key}: {proteome_path}")
                total_missing_proteome += 1
                continue

            # Index the proteome once per species for fast lookup
            protein_index = SeqIO.index(proteome_path, "fasta")

            n_written_this_species = 0
            for protein_id, (start, end) in protein_coords.items():
                if protein_id not in protein_index:
                    total_protein_not_found += 1
                    continue

                full_seq = str(protein_index[protein_id].seq)
                # HMMER ali coordinates are 1-based inclusive
                domain_seq = full_seq[start - 1:end]

                header = f"{species_key}|{protein_id}"
                out_fasta.write(f">{header}\n{domain_seq}\n")
                n_written_this_species += 1
                total_written += 1

            protein_index.close()
            print(f"{species_key}: {n_written_this_species} NB-ARC domain sequences extracted")

    print(f"\nTotal NB-ARC domain sequences written: {total_written}")
    print(f"Species with missing proteome file: {total_missing_proteome}")
    print(f"Proteins not found in their proteome (ID mismatch): {total_protein_not_found}")
    print(f"Output: nbarc_domains_combined.fasta")


if __name__ == "__main__":
    main()
