# pipeline/

Nine ordered stages. Each `stepN_…/run.sh` sources the top-level
`config.yaml` (a copy of `config.template.yaml` with your paths filled in)
and invokes one or more `egt` subcommands.

Run one stage at a time:

```sh
cd pipeline/step1_generate_newick
bash run.sh
```

Or run all in order:

```sh
bash pipeline/run_all.sh
```

Clade-subset variants for steps 7 and 9 are in the same directory as the
primary `run.sh` — invoke them directly.

| Step | Purpose                                       | egt subcommand(s)                                        |
|------|-----------------------------------------------|----------------------------------------------------------|
| 1    | NCBI taxid → Newick topology                  | `egt taxids-to-newick`                                   |
| 2    | Newick + TimeTree → calibrated tree + node/edge TSVs | `egt newick-to-common-ancestors`                  |
| 3    | Pairwise ALG-decay + ALG-dispersion plots     | `egt decay-pairwise`, `egt alg-dispersion`               |
| 4    | Perspective-chromosome fusion/loss inference  | `egt alg-fusions`                                        |
| 5    | Chromosome-count vs rearrangement scatter     | `egt chrom-number-vs-changes`                            |
| 6    | Perspective-df → tree with Monte Carlo        | `egt perspchrom-df-to-tree`                              |
| 7    | Branch-wise rates against geologic time       | `egt branch-stats-vs-time` (+ clade variants)            |
| 8    | Annotated tree plot                           | `egt branch-stats-tree`                                  |
| 9    | Collapsed-tree + Fourier analyses             | `egt collapsed-tree` (+ clade variants)                  |
