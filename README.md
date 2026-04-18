# schultz-et-al-2026

Reproduction pipeline for the analyses in
[Schultz et al. (2024) bioRxiv](https://doi.org/10.1101/2024.07.29.605683).

Everything here is driven by the
[`egt`](https://github.com/conchoecia/egt) command-line interface — this
repository contains no analysis code of its own. It ships:

- `genome_database/` — the exact genome list used to build the comparative
  database, plus the vendored `chrombase` snakefiles that build it from
  the list.
- `pipeline/` — nine ordered stages (`step1_…` through `step9_…`), each a
  thin shell script that invokes one or more `egt` subcommands against a
  set of paths declared in `config.yaml`.
- `config.template.yaml` — every path and tunable parameter that the
  pipeline uses, with defaults matching the exact values the published
  analyses ran with.

## Quickstart

```sh
# 1. clone
git clone https://github.com/conchoecia/schultz-et-al-2026.git
cd schultz-et-al-2026

# 2. create the conda env (~5 min; uses mamba/micromamba if available)
#    On an HPC with lmod:  module load Conda/Miniforge3
bash bin/setup_env.sh
conda activate egt-repro

# 3. configure
cp config.template.yaml config.yaml
#  ... edit config.yaml: paths to RBH_DIR, ALG_RBH, GENOME_CONFIG_YAML, etc. ...

# 4. (optional) build the genome database from the 5,821-species list
#    Skip if you already have chromosome-scale genomes indexed.
bash genome_database/submit.sh

# 5. run the nine-stage pipeline
bash pipeline/run_all.sh
```

The setup script creates a conda env called `egt-repro` with Python 3.12
and installs `egt` from PyPI. `egt` transitively pulls in the scientific
Python stack the pipeline needs (numpy, pandas, scipy, scikit-learn,
matplotlib, networkx, umap-learn, bokeh, ete4, snakemake, …).

After `conda activate egt-repro`, `egt --help` should list all subcommands.
To capture exact dep versions for an archival run, `bash bin/freeze_env.sh`
writes `environment.lock.yml` + `requirements.lock.txt`.

## Layout

```
schultz-et-al-2026/
├── environment.yml           — conda env spec (python 3.12 + egt from PyPI)
├── config.template.yaml      — paths + tunables (copy → config.yaml)
├── bin/
│   ├── setup_env.sh          — create/update the conda env
│   └── freeze_env.sh         — capture exact dep versions (lockfile)
├── genome_database/          — genome list + vendored chrombase snakefiles
├── pipeline/                 — 9 ordered step directories
│   ├── run_all.sh
│   ├── step1_generate_newick/           → egt taxids-to-newick
│   ├── step2_divergence_times/          → egt newick-to-common-ancestors
│   ├── step3_dispersal_characterization/ → egt decay-pairwise, alg-dispersion
│   ├── step4_persp_chr/                 → egt alg-fusions
│   ├── step5_perspchangeplot/           → egt chrom-number-vs-changes
│   ├── step6_perspchrom_df_to_tree/     → egt perspchrom-df-to-tree
│   ├── step7_tree_analysis_branchstatsvtime/ → egt branch-stats-vs-time
│   ├── step8_plot_tree_analysis/        → egt branch-stats-tree
│   └── step9_fourier/                   → egt collapsed-tree, fourier-of-rates
├── data/                     — placeholder for processed intermediates (Dryad)
└── docs/
    ├── overview.md
    └── reproducing.md
```

## License

MIT — see [`LICENSE`](LICENSE).
