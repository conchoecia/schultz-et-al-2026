# overview

The `schultz-et-al-2026` repository is the glue that turns the published
analyses — described in the preprint at
[doi:10.1101/2024.07.29.605683](https://doi.org/10.1101/2024.07.29.605683)
— into an executable, reviewer-runnable pipeline.

## Who does what

- `chrombase` builds the chromosome-scale genome database from an NCBI
  species list. We vendor its snakefiles under `genome_database/` and
  ship the exact 5,821-assembly list used.
- `egt` ([`github.com/conchoecia/egt`](https://github.com/conchoecia/egt),
  on PyPI as `pip install egt`) provides every analysis CLI used by
  this repository.
- `schultz-et-al-2026` itself (this repo) contains **no analysis code**.
  It is a set of shell scripts driving `egt` and `chrombase` against
  paths declared in one top-level `config.yaml`.

## Pipeline shape

```
            genome_list.yaml
                   │
                   ▼
          ┌──────────────────┐
          │    chrombase     │   (genome_database/submit.sh)
          │   scrape + build │
          └────────┬─────────┘
                   │
                   ▼
    $GENOMES_DIR  (indexed chrom-scale genomes)
                   │
                   ▼
         ┌─────────────────────┐
         │   odp (not in repo) │   produce per-species RBH files against the
         │   reciprocal best   │   BCnSSimakov2022 ALG database
         └──────────┬──────────┘
                    │
                    ▼
              $RBH_DIR
                    │
                    ▼
  ┌───────────────────────────────────────────────────────────┐
  │ pipeline/ — nine ordered stages, each one `egt …` call(s) │
  │                                                           │
  │ step1 → step2 → step3 → step4 → step5 → step6 → step7 →   │
  │                                        step8 → step9     │
  └───────────────────────────────────────────────────────────┘
```

## The nine stages

See `pipeline/README.md` for the one-line summary of each stage. In
order they cover:

1. NCBI taxonomy → Newick topology
2. Calibration against TimeTree divergence times
3. Pairwise dispersal curves (figure: decay-over-time)
4. Perspective-chromosome inference (fusion/loss events on the tree)
5. Chromosome-count scatter against rearrangement-rate
6. Monte-Carlo-supported mapping of changes onto the tree
7. Branch-wise rates against geologic time (+ two clade-subset variants)
8. Annotated tree visualization
9. Collapsed-tree + Fourier analyses of rate periodicity (+ clade variants)
