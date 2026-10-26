# Software versions

All analyses were run on a Linux server (Ubuntu 22.04, kernel 5.15) with conda environments. Versions below were read from the environment exports in `environments/` and, where marked *log*, confirmed from the header of the output files of the step itself.

| Step | Tool | Version | Environment | Status |
|---|---|---|---|---|
| Gene annotation | ANNEVO | git commit `37bdd9aa62ddf24fa55941fb827061f7ed49ce53` (committed 7 Jun 2026) | `ANNEVO` (Python 3.10.20) | confirmed (`~/software/ANNEVO`) |
| Domain search (NB-ARC, LRR, TIR, RPW8, Rx_N, EDS1_EP) | HMMER (`hmmscan`, `hmmsearch`) | 3.3.2 (Nov 2020) | `base` (system) | confirm from `# Version` in the `.domtblout` files |
| Pfam profiles | Pfam-A (`nlr_profiles.hmm` subset: NB-ARC, LRR_1/3/4/8, TIR, TIR_2; full `Pfam-A.hmm` for the scan of the 9,934 added proteins) | `nlr_profiles.hmm`: Pfam 37.1 (stated in `run_nbarc_phylogeny_en.sh`); `/data/shared_dbs/pfam/Pfam-A.hmm`: file dated 2024-10-05, 23,794 profiles, no release file shipped (**ask the owner of the shared database which release it is**) | `nlr_domain_search/` | partly confirmed |
| Genome-motif NLR annotation (cross-validation) | NLR-Annotator | v2.1b (`NLR-Annotator-v2.1b.jar`), git commit `5d7f1ecfae8927f34baa1126998ccc9869995784` | Java: OpenJDK 25.0.2 (`base`) | confirmed |
| Redundancy reduction | CD-HIT | 4.8.1 | `base` | env export |
| Protein alignment | MAFFT | 7.490 (`base`); 7.525/7.526 in other environments | to confirm | confirm from `mafft.log` |
| Codon alignment | PAL2NAL | v14 | `base` | env |
| CDS extraction | gffread | 0.12.7 | `base` | env |
| Phylogenies (family tree, clade trees, HyPhy input trees) | FastTree | 2.1.11 (`base`; confirmed in `fasttree.log` and `nbarc_hmm.fasttree.log`: `FastTree -log … nbarc_hmmaligned_cov70.fasta`) | `base` | confirmed for the family tree; 2.2.0 exists only inside `busco_env`/`orthofinder_env` |
| Phylogeny (ADR1 tree, Q.PLANT+R5, UFBoot) | IQ-TREE | 3.1.4 (`iqtree` env, built 29 Sep 2026): ADR1/Clade18 tree, `iq_c50`, `iq_c60`; 3.1.2 (built 17 May 2026, `orthofinder_env`): the other 42 logs under `project/` | per-tree | confirmed from log headers |
| Selection (FEL, MEME, BUSTED) | HyPhy | **2.5.103** | `hyphy_env` | env (the 2.5.36 in `base` was not used) |
| Orthogroups | OrthoFinder | 3.1.5 (DIAMOND 2.2.2, MCL 22.282) | `orthofinder_env` | env; confirm from `Log.txt` |
| Genome completeness | BUSCO, lineage `eudicots_odb10` | 5.7.1 (`busco_env`, stated in `qc_seqkit_busco_flat.sh`); 6.1.0 also installed (`busco`) | `busco_env` | confirm from `short_summary*` |
| Genome-level search for ADR1 | BLAST+ (`tblastn`), miniprot | BLAST+ 2.12.0 (`base`) or 2.17.0 (`blast_env`); miniprot 0.18-r281 | to confirm | confirm which BLAST+ |
| Structure prediction | AlphaFold 3 | AlphaFold Server (web service, AlphaFold 3); models downloaded Oct 2026 (`14_rename_af3_models.py` renames the Server files) | not run on this server | record exact access dates in the README |
| Structure comparison | Foldseek (TM-score via `--alignment-type 1`) | Foldseek 10.941cd33 (`foldseek_env`); a standalone `/usr/bin/TMalign` also exists (version not recorded; `13_structure_pairs.py` does not call it) | `foldseek_env` | confirmed |
| Genome-level distance check (Mentha) | Mash | 2.3 | `mash_env` | env |
| Utilities | SeqKit 2.13.0; samtools 1.6 (`base`) / 1.24 (`ANNEVO`) | | `base` | env |
| Analysis scripts | Python 3.11.13; NumPy 2.4.6; pandas 3.0.3; SciPy 1.17.1; Matplotlib 3.11.1; Biopython 1.87; ete3 3.1.3 | | `base` | env |

Not used: DAMBE / Xia test (saturation was assessed with the third-position p-distance criterion; see Methods).

## Chain of the conserved core (for the Methods)
Tree: `nbarc_hmm_plus_refs.treefile` (FastTree 2.1.11 on the HMM-aligned NB-ARC domains; `FastTree -log … nbarc_hmmaligned_cov70.fasta`) with RefPlantNLR references (`TableS1_description.xlsx`) and CD-HIT 85 % clusters (`nbarc_domains_cdhit85.fasta.clstr`).
Partition into maximal subtrees of <= 300 leaves (`01_partition_refs.py`), clade statistics (`02_clade_stats.py`), cap sensitivity 150-500 (`24`, `25`), exclusion of *Paraleonurus japonicus* at counting time without recomputing the partition (`28_clade_stats_exclude.py`). Core rule: >= 35 of 45 species and <= 3 mean copies per species.
Note: the family tree and the OrthoFinder run (`Results_Jul25`) were built with the 46-genome set; the excluded species is removed only when counting.
