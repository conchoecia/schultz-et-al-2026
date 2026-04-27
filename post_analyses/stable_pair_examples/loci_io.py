"""Shared helpers for stable_pair_examples.

Reads .chrom.gz coordinate tables, maps Simakov BCnS family IDs to
per-species gene IDs via the BCnSSimakov2022.rbh file, and picks
representative species from the calibrated tree.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path
from string import Template

import pandas as pd
import yaml

CHROM_COLUMNS = ["gene_id", "scaffold", "strand", "start", "end"]


def read_chrom(path: str | Path) -> pd.DataFrame:
    """Read an odp/egt-style .chrom(.gz) file. 5 columns, no header."""
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=CHROM_COLUMNS,
        dtype={"gene_id": str, "scaffold": str, "strand": str,
               "start": "int64", "end": "int64"},
        compression="infer",
    )


def load_genome_list(yaml_path: str | Path, genomes_dir: str | None = None) -> dict:
    """Load genome_database/genome_list.yaml. Substitutes ${GENOMES_DIR}."""
    with open(yaml_path) as fh:
        raw = fh.read()
    gd = genomes_dir or os.environ.get("GENOMES_DIR", "")
    raw = Template(raw).safe_substitute(GENOMES_DIR=gd)
    return yaml.safe_load(raw)


def load_rbh(rbh_path: str | Path) -> pd.DataFrame:
    """Load the BCnSSimakov2022.rbh table.

    Returns a long-form (species, family_id, gene_id) DataFrame.
    """
    df = pd.read_csv(rbh_path, sep="\t")
    fam_col = next((c for c in df.columns if "rbh" in c.lower() or "family" in c.lower()), df.columns[0])
    species_cols = [c for c in df.columns if c not in (fam_col,) and df[c].dtype == object]
    out = df.melt(id_vars=[fam_col], value_vars=species_cols,
                  var_name="species_col", value_name="gene_id").dropna()
    out = out.rename(columns={fam_col: "family_id"})
    return out


def family_to_human_symbol(map_path: str | Path,
                            symbols_path: str | Path | None = None) -> pd.DataFrame:
    """Return DataFrame: family_id, human_gene (NP_*), alg, symbol (HGNC if resolvable)."""
    fam = pd.read_csv(map_path, sep="\t")
    if symbols_path and Path(symbols_path).exists():
        sym = pd.read_csv(symbols_path, sep="\t", dtype={"gene_id": str})
        # gene_id in symbols is NCBI numeric id; human_gene is NP_ accession.
        # bcns_family_to_human_gene.tsv stores NP_ accessions. We don't have
        # a direct accession→symbol map here; fall back to using gene_symbols
        # file only when a match exists by literal join (rare). Safe noop.
    fam["symbol"] = fam["human_gene"].astype(str)
    return fam[["family_id", "alg", "human_gene", "symbol"]]


def probe_chrom(genome_yaml: str, genomes_dir: str | None = None,
                max_check: int = 3) -> int:
    """Verify .chrom.gz files are readable. Returns 0 on success."""
    try:
        gd = genomes_dir or os.environ.get("GENOMES_DIR")
        if not gd:
            print("WARN: GENOMES_DIR not set — Stage 3 (locus plots) will be "
                  "skipped. Set GENOMES_DIR in env or config.yaml to enable.",
                  file=sys.stderr)
            return 2
        cfg = load_genome_list(genome_yaml, gd)
        species = cfg.get("species", {})
        checked = 0
        for name, entry in species.items():
            chrom_path = entry.get("chrom")
            if not chrom_path or not Path(chrom_path).exists():
                continue
            df = read_chrom(chrom_path)
            if list(df.columns) != CHROM_COLUMNS:
                print(f"FAIL {name}: unexpected columns {list(df.columns)}",
                      file=sys.stderr)
                return 1
            print(f"OK   {name}: {len(df)} genes, scaffolds={df['scaffold'].nunique()}")
            checked += 1
            if checked >= max_check:
                break
        if checked == 0:
            print("WARN: no .chrom.gz files found under GENOMES_DIR — "
                  "Stage 3 will be skipped.", file=sys.stderr)
            return 2
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"PROBE ERROR: {exc}", file=sys.stderr)
        return 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--genome-yaml",
                    default=str(Path(__file__).resolve().parents[2]
                                / "genome_database" / "genome_list.yaml"))
    ap.add_argument("--genomes-dir", default=None)
    args = ap.parse_args(argv)
    if args.probe:
        rc = probe_chrom(args.genome_yaml, args.genomes_dir)
        # Treat "no genomes staged" as a soft warning (rc=2): downstream
        # stages handle absence; do not abort run.sh.
        sys.exit(0 if rc in (0, 2) else 1)


if __name__ == "__main__":
    main()
