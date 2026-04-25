Inputs for Figure 4 and Figure 5 UMAP plotting scripts.

The compressed UMAP dataframe is the source table used by
`post_analyses/Figure_4_plots/run.sh` and
`post_analyses/Figure_5_plots/run.sh`. Pandas reads the `.gz` file directly.

The calibrated Newick is included with these inputs for figure provenance and
for linked-tree/Bokeh rerendering context. The static Figure 4/5 UMAP PDF
scripts currently assert that it is present, but do not draw from it directly.

Original source paths:

- `/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509_phyloSampling/features_plot/featuresplot_paperpalette_refresh2_20260423_g3jFJL/subsample_allsamples.neighbors_250.mind_1.0.missing_large.paper_palette.df`
- `/lisc/data/scratch/molevo/dts/manifold/newick_and_timetree_20251118/step2_download_newick_from_timetree/20251130Tree.calibrated_tree.nwk`

SHA-256 checksums:

```text
dc4756ded1b6acc8aa5b5418a21279a60c9e73c595b0fd5da0d506aa134a5379  subsample_allsamples.neighbors_250.mind_1.0.missing_large.paper_palette.df.gz
ca2baee0cf8b2476c67aa06cb9a77dd8a734d4241e7f77206be464ba5e75d94f  20251130Tree.calibrated_tree.nwk
```
