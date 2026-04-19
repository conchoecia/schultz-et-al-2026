# Operational notes

Running log of measured resource usage for expensive pipeline steps,
so future runs can estimate realistic sbatch allocations instead of
guessing. Append new rows over time — don't overwrite.

## COO rebuild (`egt phylotreeumap combine-distances`)

Aggregates per-species `.gb.gz` distance files into the
`allsamples.coo.npz` sparse matrix used by defining-features, UMAP,
and downstream enrichment.

| date       | species (rows) | pairs (cols) | stored nnz    | output size | peak RAM | elapsed | alloc | job id  | notes                                   |
|------------|----------------|--------------|---------------|-------------|----------|---------|-------|---------|-----------------------------------------|
| 2026-04-19 | 5,821          | 2,785,980    | 920,682,774   | 5.9 GB      | 135.5 GB | 38 min  | 400 G | 2633259 | 202509 rebuild (post row-scramble fix). |

Sizing rule of thumb from this row: **~23 GB per 1,000 species** at this pair width; **~40 s per 100 M stored cells** on a single CPU. A 192 GB / 1 h sbatch is comfortable for the 5,821-species shape; scale linearly with species count.

## Grid-verify of COO (`grid_verify_new_coo.py`)

Samples a regular (n_rows × n_cols) grid across the COO and checks each
cell against the per-species `.gb.gz` — detects row-scramble, dropped
pairs, and value-level corruption anywhere in the matrix.

| date       | COO shape              | grid   | passed  | failed | elapsed | alloc | job id  | notes                                 |
|------------|------------------------|--------|---------|--------|---------|-------|---------|---------------------------------------|
| 2026-04-19 | 5,821 × 2,785,980      | 20×20  | 360     | 40     | 3 min   | 48 G  | 2631872 | Old scrambled COO (10 % fail). Evidence for rebuild. |
| 2026-04-19 | 5,821 × 2,785,980      | 50×50  | 2,500   | 0      | 3.7 min | 48 G  | 2634159 | Post-rebuild, clean.                  |

## COO density (202509 rebuild)

- Samples (rows): **5,821**
- Locus pairs (cols): **2,785,980**
- Stored nnz: **920,682,774**
- Total cells: 5,821 × 2,785,980 = 16,217,189,580
- Filled fraction: 920,682,774 / 16,217,189,580 = **5.68 %** (i.e. ~94.3 % sparse)
- On-disk npz size: 5.9 GB

A cell is "filled" iff the corresponding ortholog pair was observed on
the same scaffold in that species (and thus has a distance). Unfilled
cells mean either the pair wasn't in that species' RBH file or the two
families landed on different scaffolds.

## SupplementaryTable_16 selection rules

The publication table is a union-over-flags filter applied to each
per-clade `*_unique_pair_df.tsv.gz`. **Kept row** iff ≥ 1 of the 5 flags
below is set. Source of truth: `src/egt/legacy/defining_features_plot2.py`
(`add_ratio_columns`, `compute_z_scores`, `assign_flags`).

Pseudocount / log-ratio (with `pseudocount=1`):

```
mean_in_out_ratio_log = log10((mean_in + 1) / (mean_out + 1))
sd_in_out_ratio_log   = log10((sd_in   + 1) / (sd_out   + 1))
```

Z-score (computed *per clade* over rows with `occupancy_in >= 0.5`):

```
mean_sigma = (mean_in_out_ratio_log − mean_of_column) / std_of_column
sd_sigma   = (sd_in_out_ratio_log   − mean_of_column) / std_of_column
```

Default `sigma = 2` (the `--sigma` flag in `defining_features_plot2.py`;
SuppTable_16 used the default).

**Scope of the z-score reference population:** mean and std used for the
z-score are computed **only over pairs with `occupancy_in ≥ 0.5`**. The
resulting transform is then applied to every pair (high- or low-occupancy)
for plotting, but the reference population that defines "typical" is
strictly the high-occupancy subset. Source: `plot2.py::main` line 434
(`df = df[df["occupancy_in"] >= 0.5]`) before the sigma columns are
computed on lines 436–437. Why: low-occupancy pairs have noisy
out-clade stats (few in-clade species observed them) and their extreme
log-ratios would inflate the std and push the ±2σ cutoffs outward,
letting through high-occupancy pairs that aren't actually extreme. Since
the flags themselves also gate on `occupancy_in ≥ 0.5`, low-occupancy
points can never be flagged — they show up on the QC plots only for
visual context below the horizontal cutoff line.

| flag | rule |
|---|---|
| `close_in_clade`   | `mean_sigma < −σ`  AND `occupancy_in ≥ 0.5` |
| `distant_in_clade` | `mean_sigma > +σ`  AND `occupancy_in ≥ 0.5` |
| `stable_in_clade`  | `sd_sigma   < −σ`  AND `occupancy_in ≥ 0.5` |
| `unstable_in_clade`| `sd_sigma   > +σ`  AND `occupancy_in ≥ 0.5` |
| `unique_to_clade`  | `notna_out == 0` (no occupancy gate) |

Union filter (`defining_features_plot2.py:445`):

```python
df = df[df[columns_of_interest].sum(axis=1) > 0]
```

SuppTable_16 at `σ=2` had 142,815 rows across 20 clades (min 386 for Bilateria, max 22,251 for Scaphopoda).

## defining-features (28 clades) — post-rebuild

| date       | clades | peak RAM | elapsed                  | alloc | jobs                      | notes |
|------------|--------|----------|--------------------------|-------|---------------------------|-------|
| 2026-04-19 | 28     | 35.6 GB  | 10:04 (first 19 clades) + 8:56 (9 rerun) | 64 G | 2634374 + 2634907 | Sparse-native. First run had 9 clades with stale pre-rebuild outputs because of the "skip if output exists" short-circuit; quarantined to `out/.STALE/` and rerun in 2634907. |

## unique_pairs filter (post-rebuild)

Applied the SuppTable_16 selection rules (above) via `build_unique_pairs_tsv.py` with `--sigma 2`.

**First attempt had a filter-ordering bug** (`unique_pairs.pre_fix.tsv.gz` in `.STALE/`):
computed the z-score mean/std over the full per-clade df instead of the
`occupancy_in >= 0.5` subset. plot2.py's `main()` at line 434 filters BEFORE
z-scoring, so including low-occupancy rows (which have noisy out-clade stats)
was inflating the variance and making almost every pair pass the 2σ test for
large-in-clade groups. Fix: filter to `occupancy_in >= OCC_MIN` in
`build_unique_pairs_tsv.py` before calling `compute_z_scores`.

Post-fix counts (match the decay-branch filter math used for SuppTable_16):

| | old SuppTable_16 | new (post-rebuild) |
|---|---|---|
| total rows | 142,815 | **175,717** |
| clades | 20 | **28** (20 original + 8 new) |
| file size | n/a (xlsx) | 20.6 MB (tsv.gz) |

Per-clade sanity (20 overlapping clades, old → new):

| clade | old | new | note |
|---|---|---|---|
| Scaphopoda | 22,251 | 22,182 | match |
| Diptera | 15,575 | 14,487 | match |
| Mollusca | 11,651 | 11,502 | match |
| Bivalvia | 11,382 | 11,669 | match |
| Spiralia | 10,261 | 10,561 | match |
| Cnidaria | 10,007 | 9,840 | match |
| Gastropoda | 9,718 | 10,390 | match |
| Ctenophora | 9,201 | 5,567 | drop; different Ctenophora sampling |
| Echinodermata | 8,170 | 7,838 | match |
| Annelida | 7,885 | 7,910 | match |
| Nematoda | 7,488 | 6,444 | match |
| Coleoidea | 6,061 | 1,871 | drop |
| Porifera | 3,714 | 4,871 | up |
| Decapodiformes | 2,427 | 2,170 | match |
| Neoptera | 1,480 | 2,347 | up |
| Insecta | 1,466 | 2,276 | up |
| Deuterostomia | 1,398 | 1,416 | match |
| Vertebrata | 1,373 | 1,428 | match |
| Protostomia | 921 | 1,139 | match |
| Bilateria | 386 | 243 | match |

New 8 clades (no baseline): Arthropoda 1,913 · Cephalopoda 1,949 · Chordata 1,411 ·
Mammalia 3,356 · Teleostei 3,837 · Clitellata 3,840 · Hexapoda 2,212 ·
Platyhelminthes 21,048.

## GO enrichment sweep (post-rebuild)

Sweep of BH-FDR-corrected hypergeometric enrichment per clade, top-N
threshold ranged over (stability / closeness / intersection) ranking
axes. Source: `post_analyses/go_enrichment_sweep/run.sh`. Emits 8
artifacts: summary.tsv, significant_terms.tsv + _annotated, curves.pdf,
dotplots.pdf, heatmap.pdf, pair_distance.pdf, volcanos.pdf,
volcanos_fold3plus.pdf.

| date       | input unique_pairs | clades | elapsed | alloc | job id  | notes |
|------------|--------------------|--------|---------|-------|---------|-------|
| 2026-04-19 | post-rebuild (175,717 rows) | 28 | 1:27 sweep + ~2 min plots | 4 G | 2635758 | All 28 clades hit q ≤ 0.25 in ≥ 1 sweep config. Summary: 3,220 rows, 48,870 significant terms. |

## Downstream steps

(to be filled in as each step runs against the post-rebuild artifacts — UMAP regen)
