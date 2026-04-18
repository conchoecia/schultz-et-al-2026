# post_analyses/

Supplementary analyses that consume the outputs of the 9-stage pipeline.
Each subdirectory wraps one of the new `egt >= 0.2` subcommands (or a
short helper that produces its inputs) against the top-level
`config.yaml`.

These are not required to reproduce the main figure-generating pipeline.
Run any of them after `pipeline/run_all.sh` has finished.

| Directory | egt subcommand | Produces |
|-----------|----------------|----------|
| `palette_preview/` | `egt palette-preview` | Schematic-tree PDF + FigTree-ready Newick/NEXUS colored by the paper palette |
| `pigeonhole_check/` | `egt pigeonhole-check` | Null-model test of ALG-pair co-localization per clade |
| `build_family_naming_map/` | `egt build-family-naming-map` | BCnS family → human gene ID map (prerequisite for GO enrichment) |
| `divergence_vs_dispersal/` | `egt divergence-vs-dispersal` | Per-clade regression of protein divergence on ALG dispersal |
| `entanglement_browse/` | `egt entanglement-browse` | Clade-characteristic ALG fusion pairs |
| `entanglement_go_enrich/` | `egt entanglement-go-enrich` | GO enrichment on clade-characteristic ALG gene sets |

## Typical run order

```
palette_preview            # anytime, purely visual
pigeonhole_check           # anytime after step 4
divergence_vs_dispersal    # needs an upstream divergence.tsv
build_family_naming_map    # one-shot; feeds the next two
entanglement_browse
entanglement_go_enrich     # needs GOA download
```

Each `run.sh` sources `config.yaml` from the repo root and writes
outputs under `post_analyses/<name>/out/`.
