#!/usr/bin/env python3
"""End-to-end smoke-test of `egt phylotreeumap combine-distances`.

Builds a tiny 5-species x 6-pair fixture in a tmp dir, runs the CLI as
a subprocess, loads the output npz, and checks shape + hand-verifiable
values. Confirms the CLI wrapper works before the 6-hour production
rebuild finishes (so a wrapper bug wouldn't waste compute).
"""
from __future__ import annotations
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from scipy.sparse import load_npz

# 4 BCnS families, C(4,2) = 6 pairs, 5 species.
FAMS = ["fam_A", "fam_B", "fam_C", "fam_D"]
PAIRS = [(FAMS[i], FAMS[j]) for i in range(len(FAMS)) for j in range(i + 1, len(FAMS))]
COMBO_TO_IX = {p: k for k, p in enumerate(PAIRS)}

SPECIES = [f"species_{i:03d}" for i in range(5)]

# distance[s][p] encodes (species, pair) uniquely as 1000*(s+1)+p
EXPECTED = {(s, p): 1000 * (s + 1) + p
            for s in range(len(SPECIES))
            for p in range(len(PAIRS))}


def build_fixture(root: Path) -> tuple[Path, Path, Path]:
    gb_dir = root / "gbgz"; gb_dir.mkdir()
    sampledf_rows = []
    for s, sample in enumerate(SPECIES):
        gbgz = gb_dir / f"{sample}.gb.gz"
        with gzip.open(gbgz, "wt") as fh:
            fh.write("rbh1\trbh2\tdistance\n")
            for p, (a, b) in enumerate(PAIRS):
                fh.write(f"{a}\t{b}\t{EXPECTED[(s, p)]}\n")
        sampledf_rows.append({"sample": sample,
                              "dis_filepath_abs": str(gbgz)})
    sampledf_path = root / "sampledf.tsv"
    pd.DataFrame(sampledf_rows).to_csv(sampledf_path, sep="\t")
    algcombo_path = root / "combo_to_index.txt"
    with open(algcombo_path, "w") as fh:
        for p, ix in COMBO_TO_IX.items():
            fh.write(f"{p}\t{ix}\n")
    return sampledf_path, algcombo_path, gb_dir


def main():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        sampledf, algcombo, _ = build_fixture(root)
        out = root / "test.coo.npz"
        egt = shutil.which("egt") or os.path.join(os.path.dirname(sys.executable), "egt")
        cmd = [egt, "phylotreeumap", "combine-distances",
               "--sampledf", str(sampledf),
               "--algcomboix", str(algcombo),
               "--output", str(out)]
        print(f"[smoke] $ {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print("STDERR:\n" + result.stderr, file=sys.stderr)
            sys.exit(f"CLI failed with exit {result.returncode}")
        assert out.exists(), f"expected output {out} missing"
        coo = load_npz(out)
        assert coo.shape == (len(SPECIES), len(PAIRS)), \
            f"shape mismatch: {coo.shape} != ({len(SPECIES)}, {len(PAIRS)})"
        csr = coo.tocsr()
        mismatches = []
        for (s, p), want in EXPECTED.items():
            got = float(csr[s, p])
            if got != want:
                mismatches.append((s, p, got, want))
        if mismatches:
            for s, p, got, want in mismatches[:5]:
                print(f"  MISMATCH row={s} col={p} got={got} want={want}")
            sys.exit(f"{len(mismatches)}/{len(EXPECTED)} cells mismatched")
        print(f"\n[smoke] PASS: shape={coo.shape}, all {len(EXPECTED)} cells match expected values.")


if __name__ == "__main__":
    main()
