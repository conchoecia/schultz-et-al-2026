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

### Step 3 is optional

`step3_dispersal_characterization` is a **leaf** stage — its outputs are
not consumed by steps 4–9. It also wants a `SYNTENY_DIR` (~58 GB of
odp-generated synteny PDFs) that is *not* part of the Dryad dataset, so a
typical reviewer environment won't have it. `run_all.sh` skips this step.

Run it manually if and only if you have the synteny directory available:

```sh
SYNTENY_DIR=/path/to/odp/step2-figures/synteny_nocolor  \
    bash pipeline/step3_dispersal_characterization/run.sh
```

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
