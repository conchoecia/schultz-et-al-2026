# stable_pair_examples

Builds the stable-pair-examples figure: a clade-wide distribution of stable-pair
distances plus locus tracks for top-ranked example pairs at three distance
scales (close <1 kb, intermediate 10–100 kb, long-range 100 kb–1 Mb), focused
on Vertebrata (human regulatory context) and Diptera (Drosophila context).

Driven by `egt`'s upstream output `unique_pairs.tsv.gz` from
`post_analyses/defining_features_pairs/out/`.

## Outputs

- `out/candidate_pairs.tsv` — all qualifying stable pairs with composite
  score, distance bin, GO/entanglement/fusion support columns.
- `out/selected_examples.tsv` — top-N per (clade, distance_bin).
- `out/stable_pair_distance_distribution.pdf` — Stage 2 distribution panel.
- `out/stable_pair_example_loci.pdf` — Stage 3 cross-species locus tracks.
- `out/epigenomics/` — cached ReMap + cCRE + REDfly BED files (Stage 0).

## Run

```
sbatch run.sh        # SLURM
bash run.sh          # login node (Stage 1+2 are light)
NO_FETCH=1 bash run.sh   # skip Stage 0 network
```

## Notes

- Stage 3 needs `GENOMES_DIR` exported so `genome_database/genome_list.yaml`
  resolves to actual `.chrFilt.chrom.gz` paths. If unset, Stage 3 writes a
  placeholder PDF and the rest of the pipeline still completes.
- Composite score weights are documented in `select_examples.py`.
- EIF2S1/ATP6V1D@Vertebrata is force-kept as the close exemplar.
