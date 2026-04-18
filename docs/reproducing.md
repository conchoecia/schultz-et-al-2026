# reproducing.md

End-to-end walkthrough for reviewers. Assumes a linux HPC with
`Conda/Miniforge3` available (the pipeline was developed on LiSC at the
University of Vienna); adapt paths for your environment.

## 0. Prerequisites

- **conda / mamba / micromamba**. Every modern HPC ships one; on LiSC,
  `module load Conda/Miniforge3`.
- **`chrombase`** checked out locally if you need to rebuild the genome
  database: `git clone https://github.com/conchoecia/chrombase`.
- **`odp`** checked out if you need to regenerate the per-species RBH
  files: `git clone https://github.com/conchoecia/odp`. The ALG
  database used is `BCnSSimakov2022` under `odp/LG_db/`.
- **TimeTree export (Newick)** — download from
  [timetree.org](https://timetree.org) for the 4,454 species in the
  dataset. Place at `pipeline/step2_divergence_times/newick_timetree.nwk`.
- **Extinction intensity TSV** — a per-age extinction-rate table (e.g.
  derived from the Sepkoski dataset). Place at
  `pipeline/step7_tree_analysis_branchstatsvtime/extinction_intensity.tsv`.

## 1. Clone + environment + configure

```sh
git clone <repo-url>  schultz-et-al-2026
cd schultz-et-al-2026

# Create the conda env (Python 3.12 + egt + its scientific-Python stack).
bash bin/setup_env.sh
conda activate egt-repro
egt --help     # sanity check

# Configure.
cp config.template.yaml config.yaml
# Edit config.yaml so the paths match your environment. The defaults are
# the exact values used for the published run on LiSC.
```

(Optional archival step: after activating, `bash bin/freeze_env.sh` emits
`environment.lock.yml` + `requirements.lock.txt` capturing every resolved
version. Commit these alongside a paper submission for exact
reproducibility.)

Key fields to edit:

- `EGT_CONDA_ENV` — path to the conda env with `egt` installed.
- `GENOME_CONFIG_YAML` — points at `genome_database/genome_list.yaml`
  (after you've substituted `${GENOMES_DIR}` there) or at a pre-built
  odp config.
- `RBH_DIR` — where per-species `.rbh` files live.
- `ALG_RBH` — the BCnSSimakov2022 ALG RBH, shipped with
  `odp`'s `LG_db/`.
- `SYNTENY_DIR` — pre-generated odp synteny plots.
- `TIME_NEWICK` / `EXTINCTION_INTENSITY` — user-supplied inputs noted
  above.

## 2. Build the genome database (optional)

Skip this if you already have an indexed chromosome-scale genome tree
equivalent to chrombase's output.

```sh
export CHROMBASE_DIR=/path/to/chrombase/checkout
export GENOMES_DIR=/path/to/stage/genomes
cp genome_database/chrombase.config.yaml.template genome_database/chrombase.config.yaml
# edit genome_database/chrombase.config.yaml
bash genome_database/submit.sh
```

5,821 assemblies is a 1–2 TB stage; plan storage accordingly.

## 3. Generate per-species RBH files with odp

(Not orchestrated from this repo — run odp against `GENOMES_DIR` using
its own `CONFIG_odp.yaml`. See the odp README. The output directory is
what `RBH_DIR` in `config.yaml` should point at.)

## 4. Run the pipeline

```sh
bash pipeline/run_all.sh
```

Or one stage at a time (useful for re-running):

```sh
bash pipeline/step1_generate_newick/run.sh
bash pipeline/step2_divergence_times/run.sh
# ...
```

On a SLURM cluster, each `run.sh` has `#SBATCH` headers and can be
submitted directly:

```sh
sbatch pipeline/step4_persp_chr/run.sh
```

### Clade-subset variants

Steps 7 and 9 each ship a main `run.sh` and two clade-subset variants
(`run_protost_minus_clitellata.sh`,
`run_verte_minus_teleost.sh`). Run them after the main `run.sh` in the
same step has produced `branch_stats_output/`.

Step 9 additionally ships `run_verte_time_sweep.sh` (250 Mya → 550 Mya
max-time window sweep for the vertebrate subset).

## 5. Output filenames by stage

| Stage | Primary outputs                                                             |
|-------|------------------------------------------------------------------------------|
| 1     | `ncbi_tree.nwk`                                                              |
| 2     | `${TREE_PREFIX}.calibrated_tree.nwk`, `.divergence_times.txt`, `.node_information.tsv`, `.edge_information.tsv` |
| 3     | `odp_pairwise_decay/`, `cache/`                                              |
| 4     | `perspchrom.tsv`, `per_species_ALG_presence_fusions.tsv`, `changestring_checkpoints/`, `locdf.tsv`, `tree1.tsv.gz`, `unique_changes_summary.tsv` |
| 5     | `chrom_number_vs_changes.pdf`                                                |
| 6     | `statsdf.tsv`, `simulations/`                                                |
| 7     | `branch_stats_output/{modified_node_list.tsv, modified_edge_list.tsv, per_clade_analyses/, custom_clade_analyses/}` |
| 8     | `tree.pdf`, `tree_2d_bivariate*.pdf`, `tree_diagnostic_*.pdf`                |
| 9     | `collapsed_tree.pdf`, `per_clade_analyses/fourier_analysis/*`, `custom_clade_analyses/*/fourier_analysis/*` |

## Expected cost

- Step 4 (`alg-fusions`, 31 cores): ~12 hours wall-clock on the full dataset.
- Step 6 (`perspchrom-df-to-tree`, 25 cores, 1000 MC simulations): ~12–24
  hours.
- Step 8 (tree plots): ~1–4 hours depending on tree size.
- Other stages: minutes to low single-digit hours.

Total: one full pass on a mid-size cluster is ~2 days wall-clock plus
the DB build (1–3 days).
