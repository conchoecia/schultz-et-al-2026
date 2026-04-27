"""Stage 0: cache regulatory BED tracks for human (hg38) and Drosophila (dm6).

Primary source is ReMap 2022 (https://remap.univ-amu.fr/), a meta-catalog of
ChIP-seq peaks aggregated across hundreds of GEO/ENCODE/modENCODE studies.
Complementary tracks: ENCODE cCRE (human), phastCons27way (dm6), REDfly CRMs
(dm6, optional).

Idempotent: skips files that already exist and are non-empty. Writes
out/epigenomics/MANIFEST.tsv with source URL, fetch date, sha256.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import sys
import urllib.request
from pathlib import Path

DEFAULT_SOURCES = {
    "hg38/ReMap2022_nr_macs2_hg38.bed.gz":
        "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz",
    "hg38/ENCODE_cCRE_hg38.bed":
        "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
    "dm6/ReMap2022_nr_macs2_dm6.bed.gz":
        "https://remap.univ-amu.fr/storage/remap2022/dm6/MACS2/remap2022_nr_macs2_dm6_v1_0.bed.gz",
    # phastCons bigWigs (UCSC multiz) — read by pyBigWig in plot_examples.
    # hg38 is 5.9 GB; dm6 is 234 MB. Range-requested at plot time so the
    # full files only need to live on scratch, not in memory.
    "hg38/phastCons100way_hg38.bw":
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw",
    "dm6/phastCons27way_dm6.bw":
        "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phastCons27way/dm6.27way.phastCons.bw",
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(out_dir: Path, sources: dict[str, str], no_fetch: bool = False) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for relpath, url in sources.items():
        target = out_dir / relpath
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.stat().st_size > 0:
            print(f"SKIP {relpath} (cached, {target.stat().st_size} bytes)")
            rows.append({"path": relpath, "url": url, "status": "cached",
                         "bytes": target.stat().st_size,
                         "sha256": _sha256(target),
                         "date": _dt.date.today().isoformat()})
            continue
        if no_fetch:
            print(f"MISS {relpath} (--no-fetch set, not downloading)",
                  file=sys.stderr)
            rows.append({"path": relpath, "url": url, "status": "missing",
                         "bytes": 0, "sha256": "", "date": ""})
            continue
        print(f"GET  {relpath} <- {url}")
        try:
            urllib.request.urlretrieve(url, target)
            rows.append({"path": relpath, "url": url, "status": "fetched",
                         "bytes": target.stat().st_size,
                         "sha256": _sha256(target),
                         "date": _dt.date.today().isoformat()})
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL {relpath}: {exc}", file=sys.stderr)
            rows.append({"path": relpath, "url": url, "status": f"fail: {exc}",
                         "bytes": 0, "sha256": "", "date": ""})
    return rows


def write_manifest(out_dir: Path, rows: list[dict]) -> None:
    manifest = out_dir / "MANIFEST.tsv"
    cols = ["path", "url", "status", "bytes", "sha256", "date"]
    with open(manifest, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"manifest -> {manifest}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True,
                    help="e.g. .../stable_pair_examples/out/epigenomics")
    ap.add_argument("--no-fetch", action="store_true",
                    help="skip network; assume prior cache")
    args = ap.parse_args(argv)
    out_dir = Path(args.out_dir)
    rows = fetch(out_dir, DEFAULT_SOURCES, args.no_fetch)
    write_manifest(out_dir, rows)


if __name__ == "__main__":
    main()
