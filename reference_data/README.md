# reference_data/

Small reference files that the pipeline consumes but that are not part of
the Dryad bundle (they live elsewhere in the authors' analysis
environment on the cluster, so we ship them here).

| File | Size | What it is | Source on LiSC |
|------|------|------------|----------------|
| `BCnSSimakov2022.rbh` | 348 KB | BCnS ancestral linkage group database (ALG RBH file). Consumed by `egt decay-pairwise`, `egt alg-fusions`, `egt perspchrom-df-to-tree`. | BCnS ALG reference data used for this analysis. |
| `species_chrom_counts.tsv` | 243 KB | Haploid chromosome count per species (5,821 rows). Consumed by `egt newick-to-common-ancestors` and `egt chrom-number-vs-changes`. | Derived from the chromosome-count reference table used for this analysis. |

These are small enough to commit to git; everything else at scale ships
through the Dryad dataset (see `bin/download_data.sh`).
