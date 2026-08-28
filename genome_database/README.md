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

## Citations for the assemblies

`genome_list.tsv` identifies each assembly by accession and NCBI submitter. That
credits institutions rather than the people who produced the data, and citation
indices do not read supplementary tables, so none of it counts as a citation.
Acknowledging only the largest contributors makes it worse: 27 groups
contributed 20 or more of these assemblies and account for 62.0% of them, while
the remaining 38.0% comes from 891 groups, 498 of which contributed exactly one
genome.

`citation_table.tsv` adds the originating publication for each assembly, one row
per accession, built with
[`chrombase`](https://github.com/conchoecia/chrombase)'s
`scripts/build_citation_table.py`:

```
citation_table.tsv                    5,821 assemblies -> DOI / PMID / authors / reference string
citation_table.bib                    the same references, deduplicated (3,459 entries)
ncbi_missing_publication_links.tsv    assemblies whose paper exists but NCBI does not link to it
```

Coverage:

| | assemblies |
| --- | --- |
| publication recovered | 4,800 (82.5%) |
| — `authoritative` (submitter linked it to the BioProject, dates consistent) | 798 |
| — `high` | 1,870 |
| — `medium` | 1,314 |
| — `low` (needs a human to confirm) | 833 |
| no publication recoverable | 1,006 |

Recovery is *better* for the groups the acknowledgements missed: **91.8%** of
assemblies from groups contributing fewer than 20 genomes resolve to a
publication.

Each row records how the publication was found (`evidence_route`), how much to
trust it (`confidence`), and the signals behind that (`notes`). Nothing is
invented — an assembly with no recoverable paper keeps its submitter so it can
still be credited at the group level. Anything below `medium` should be checked
by a human before being used as a citation.

`ncbi_missing_publication_links.tsv` is the actionable subset: **3,865
assemblies whose originating paper is public but which NCBI does not link to**,
because the BioProject record carries no `<Publication>` element. Only 815 of
the 5,171 BioProjects behind this dataset (15.8%) record one at all. Submitters
can fix this on their own BioProject records, which makes the link available to
everyone instead of only to whoever re-derives it.

To rebuild or refresh the table:

```sh
git clone https://github.com/conchoecia/chrombase
python chrombase/scripts/build_citation_table.py \
    --genome-list genome_database/genome_list.tsv \
    --out genome_database/citation_table.tsv \
    --bibtex genome_database/citation_table.bib \
    --unlinked-report genome_database/ncbi_missing_publication_links.tsv \
    --email you@example.org
```

## Genome count

5,821 assemblies from 4,454 species — see `genome_list.tsv`.
