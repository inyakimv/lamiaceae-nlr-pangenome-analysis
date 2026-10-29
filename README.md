# Pan-NLRome of Lamiaceae: analysis scripts

Scripts used for the manuscript on the NLR repertoire of 45 Lamiaceae genomes
(*Paraleonurus japonicus* was excluded after assembly checks, see below).
Only scripts that produced results or figures of the paper are included; raw genomes,
intermediate files and results are not part of this repository.
Software versions: see [`SOFTWARE_VERSIONS.md`](SOFTWARE_VERSIONS.md).

## Layout and order of execution

| Folder | Step |
|---|---|
| `scripts/01_annotation` | Genome QC (SeqKit, BUSCO), ANNEVO gene annotation, protein extraction, per-species proteome FASTAs (`make_proteomes_45.sh`, `clean_proteomes.sh`), check of official NCBI protein sets (`download_ncbi_proteomes.sh`, accessions in `data/accessions_45.tsv`), NLR-Annotator counts |
| `scripts/02_orthology` | OrthoFinder orthogroups (used for the orthogroup-based core, a robustness check) |
| `scripts/03_family_tree` | HMMER NB-ARC domain search, NB-ARC family tree (CD-HIT 85 %, FastTree), partition into clades with RefPlantNLR references; `01b_partition_refs_rerooted.py` repeats the partition after re-rooting on the stem of Clade 18 (re-rooting control) |
| `scripts/04_core_robustness` | Conserved core: clade-size cap sensitivity, exclusion of one species, ADR1 checks, `compare_core.py` (core clades under midpoint vs. re-rooted tree: same 85 clades, same four core clades; driver: `reroot_check.sh`) |
| `scripts/05_selection` | CDS extraction, codon alignments, saturation check, HyPhy (FEL, MEME, BUSTED), rarefaction |
| `scripts/06_structure` | Representative genes, AlphaFold 3 model handling, Foldseek/TM-score comparison, pLDDT |
| `scripts/07_repertoire_tandem_ploidy` | Repertoire per species, tandem organisation and permutation test, ploidy, TNL/RNL census, filter robustness |
| `scripts/08_validation_nlr_annotator` | Cross-validation of ANNEVO + HMMER with NLR-Annotator |

Scripts keep the numeric prefixes they had when they were run; within a folder, run them in numeric order.
Paths are given as command-line arguments or environment variables (see the header of each script).

## Definitions used throughout

* **Complete-like NLR**: best NB-ARC hit and at least one LRR domain (Pfam LRR_1/3/4/8) in the same protein.
* **Conserved core**: clade present in >= 35 of the 45 species with <= 3 mean copies per species,
  in subtrees of <= 300 leaves of the NB-ARC family tree (cap sensitivity: 150-500).
* **Excluded species**: *Paraleonurus japonicus*. The family tree and OrthoFinder were built with 46 genomes;
  the species is removed when counting (`04_core_robustness/28_clade_stats_exclude.py`), without recomputing the partition.

## Not included

Genome assemblies (public accessions listed in the manuscript), annotation outputs and large intermediates.

## License

The code in this repository is released under the MIT License (see `LICENSE`).

## Supplementary data

- `data/Table_S_threshold_sensitivity.tsv`: conserved-core clades under 12 combinations of species thresholds and mean-copy limits (Supplementary Methods S1).
- `scripts/03_family_tree/eds1/refs_extra.fasta` and `tbn_queries.fasta`: reference and query sequences used to identify EDS1, PAD4 and SAG101 (Supplementary Methods S3).
