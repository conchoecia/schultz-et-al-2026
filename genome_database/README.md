# genome_database/

The exact genome list used to build the chromosome-scale genome database
for the analyses, plus a vendored snapshot of the
[`chrombase`](https://github.com/conchoecia/chrombase) snakefiles that
built it from that list.

## Contents

```
genome_list.yaml             — 5,821 species, paths templated on ${GENOMES_DIR}
genome_list.tsv              — same, flat: key / accession / taxid / genus / species / minscafsize
chrombase.config.yaml.template — chrombase-side config (copy to chrombase.config.yaml, edit)
snakefiles/                  — frozen copies of the chrombase snakefiles actually run
├── chrombase_scrape_genomes_NCBI.smk
├── chrombase_build_db_annotated_chr.smk
└── chrombase_build_db_unannotated_chr.smk
scripts/
└── plot_NCBI_genomes_history.py   — helper referenced by the scrape snakefile
submit.sh                    — wrapper invoking chrombase against genome_list.yaml
```

## Provenance

Vendored from `conchoecia/chrombase` at commit **`0b7557d`**
("updates to ignore list", 2025-10-18). The snakefiles depend on other
files (`bin/`, `src/`, `dependencies/`) that are part of chrombase proper
— to actually *run* them, install chrombase itself and point
`CHROMBASE_DIR` at that checkout. The `snakefiles/` copies here are for
reference and for auditing the exact logic that was executed.

## Usage

### Option A — Reproduce the database with chrombase

```sh
git clone https://github.com/conchoecia/chrombase
export CHROMBASE_DIR=$(pwd)/chrombase
export GENOMES_DIR=/path/to/where/genomes/will/live
cp chrombase.config.yaml.template chrombase.config.yaml
# edit chrombase.config.yaml for your env
bash submit.sh
```

This stages and indexes the 5,821 assemblies in `genome_list.yaml`. The
resulting directory is what the `RBH_DIR` / `GENOMES_DIR` paths in the
top-level `config.yaml` point at.

### Option B — Use your own database

If you already have chromosome-scale genomes staged in a different tree,
just point `GENOMES_DIR` at it (the `genome_list.yaml` paths will resolve
against that), or skip the DB build entirely and provide your own RBH
files to the pipeline.

### Existing annotated-genome tree

Either `export GENOMES_DIR=/path/to/annotated_genomes`, or symlink for
convenience:

```sh
ln -s /path/to/annotated_genomes \
      genome_database/annotated_genomes_link
export GENOMES_DIR="$(pwd)/genome_database/annotated_genomes_link"
```

The symlink is in `.gitignore`; the placeholder file `annotated_genomes_link.placeholder` documents this in a fresh clone. Never copy these files into the repo — read in place.

## Genome count

5,821 assemblies from 4,454 species — see `genome_list.tsv`.
