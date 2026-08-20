# schultz-et-al-2026

Reproduction pipeline for the analyses in
[Schultz et al. (2026) *Science Advances*](https://doi.org/10.1126/sciadv.adz5561).

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

# 2. conda env  (~5 min; uses mamba/micromamba if available)
#    On an HPC with lmod:  module load Conda/Miniforge3
bash bin/setup_env.sh
conda activate egt-repro
egt --help                       # sanity check

# 3. pull data from Dryad  (~1.2 GB download, 8 GB on disk after extract)
bash bin/download_data.sh

# 4. configure (the shipped defaults already point at dryad_data/)
cp config.template.yaml config.yaml
#  ... edit config.yaml only if your paths differ ...

# 5. (optional) build the genome database from the 5,821-species list
#    Skip if you only want to reproduce downstream analyses — the RBH
#    files pulled from Dryad already cover the compute-expensive part.
bash genome_database/submit.sh

# 6. run the nine-stage pipeline
bash pipeline/run_all.sh
```

### What gets downloaded

`bin/download_data.sh` fetches two tarballs from the Dryad dataset:

- `BCnSSimakov2022_current_rbh_202509.tar.gz` (568 MB → 3.6 GB) — the
  5,821-species reciprocal-best-hits database
- `newick_and_timetree_20251118.tar.gz` (643 MB → 4.7 GB) — the full
  published analysis workflow (trees, divergence times, per-clade stats,
  …). Also contains the user-supplied TimeTree newick and extinction-
  intensity TSV, so you don't have to supply those separately. Doubles as
  a reference for diffing your rerun against the published outputs.

The setup script creates a conda env called `egt-repro` with Python 3.12
and installs `egt` from PyPI. `egt` transitively pulls in the scientific
Python stack the pipeline needs (numpy, pandas, scipy, scikit-learn,
matplotlib, networkx, umap-learn, bokeh, ete4, snakemake, …).

To capture exact dep versions for an archival run, `bash bin/freeze_env.sh`
writes `environment.lock.yml` + `requirements.lock.txt`.

### Running from scratch vs. diffing against the published outputs

- **Run from scratch**: just `bash pipeline/run_all.sh`. Each stage
  regenerates its own outputs under `pipeline/stepN_*/`.
- **Diff against published**: after `run_all.sh` finishes, compare
  `pipeline/step*/` against `dryad_data/newick_and_timetree_20251118/step*/`.

### Step 3 is optional

`step3_dispersal_characterization` is a leaf stage (plots only; its
outputs feed no downstream stages) and it wants a 58 GB synteny-plot
directory that is not in the Dryad bundle. `run_all.sh` skips it. See
[`pipeline/README.md`](pipeline/README.md).

## Layout

```
schultz-et-al-2026/
├── environment.yml           — conda env spec (python 3.12 + egt from PyPI)
├── config.template.yaml      — paths + tunables (copy → config.yaml)
├── bin/
│   ├── setup_env.sh          — create/update the conda env
│   ├── download_data.sh      — pull the Dryad bundle into dryad_data/
│   └── freeze_env.sh         — capture exact dep versions (lockfile)
├── reference_data/           — small reference files shipped in-repo
│   ├── BCnSSimakov2022.rbh   — ALG RBH database (348 KB)
│   └── species_chrom_counts.tsv
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

## Citing this repository

If you use this pipeline in your work, please cite:

> Schultz, D.T., Blümel, A., Destanović, D., Sarigol, F., & Simakov, O. (2026).
> Topological mixing and irreversibility in animal chromosome evolution.
> *Science Advances*, **12**(34), eadz5561.
> [https://doi.org/10.1126/sciadv.adz5561](https://doi.org/10.1126/sciadv.adz5561)

See also [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
