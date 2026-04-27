"""Score every candidate pair by how many of our 10+10 anchor species
have BOTH anchor genes resolved on the same scaffold.

This is what makes a "good" Stage 3 example: a pair that conserves
deeply across the focal clade's anchor lineages, not one that's
microsynteny in a single subclade.

Inputs:
    out/candidate_pairs.tsv     all qualifying stable pairs
    BCnSSimakov2022_*.rbh files (per-species RBH dir)
    genome_list.yaml + GENOMES_DIR (in-place .chrom.gz reads)

Output:
    out/candidate_pairs_anchor_coverage.tsv
    Adds columns:
        n_anchors_resolved_<clade>   how many of the clade's anchors have both
                                      ortholog1 + ortholog2 on the same scaffold
        n_subclades_resolved_<clade> independent sublineages covered
                                      (e.g. for Arthropoda: insects vs crustacea
                                      vs chelicerates)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import loci_io
import recompute_distances

REPO_ROOT = Path(__file__).resolve().parents[2]

# Sublineage tags for diversity scoring. For each anchor species, which broad
# subclade does it belong to? Used to compute n_subclades_resolved.
SUBCLADE = {
    # Vertebrata
    "Homosapiens-9606-GCF000001405.40":        "Mammalia",
    "Musmusculus-10090-GCF000001635.27":       "Mammalia",
    "Rattusnorvegicus-10116-GCF000001895.5":   "Mammalia",
    "Pantroglodytes-9598-GCF000001515.7":      "Mammalia",
    "Gallusgallus-9031-GCF000002315.6":        "Aves",
    "Anoliscarolinensis-28377-GCF000090745.2": "Reptilia",
    "Xenopustropicalis-8364-GCF000004195.4":   "Amphibia",
    "Daniorerio-7955-GCF000002035.6":          "Teleostei",
    "Lepisosteusoculatus-7918-GCF000242695.1": "Holostei",
    "Petromyzonmarinus-7757-GCF010993605.1":   "Cyclostomata",
    # Arthropoda
    "Drosophilamelanogaster-7227-GCF000001215.4":   "Diptera",
    "Anophelesfunestus-62324-GCF943734845.2":       "Diptera",
    "Aedesaegypti-7159-GCF002204515.2":             "Diptera",
    "Apismellifera-7460-GCF000002195.4":            "Hymenoptera",
    "Bombusterrestris-30195-GCF000214255.1":        "Hymenoptera",
    "Triboliumcastaneum-7070-GCF000002335.3":       "Coleoptera",
    "Bombyxmori-7091-GCF014905235.1":               "Lepidoptera",
    "Daphniapulex-6669-GCF021134715.1":             "Branchiopoda",
    "Penaeusmonodon-6687-GCF015228065.2":           "Malacostraca",
    "Parasteatodatepidariorum-114398-GCF043381705.1": "Arachnida",
}


def build_lookup(species_keys: list[str], rbh_dir: Path,
                  genome_cfg: dict) -> dict[str, dict[str, str]]:
    """For each species, return {family_id: scaffold} dict by joining its
    RBH (family→gene_id) with its chrom (gene_id→scaffold)."""
    lookup: dict[str, dict[str, str]] = {}
    for sk in species_keys:
        ent = genome_cfg["species"].get(sk)
        if not ent or not Path(ent.get("chrom", "")).exists():
            print(f"WARN: {sk} chrom.gz not found", file=sys.stderr)
            lookup[sk] = {}
            continue
        rbh_path = recompute_distances.per_species_rbh_path(rbh_dir, sk)
        fam2gene = recompute_distances.load_species_family_to_gene(rbh_path, sk)
        if not fam2gene:
            print(f"WARN: {sk} RBH empty", file=sys.stderr)
            lookup[sk] = {}
            continue
        chrom_df = loci_io.read_chrom(ent["chrom"])
        gene2scaf = dict(zip(chrom_df["gene_id"], chrom_df["scaffold"]))
        fam2scaf = {}
        for fam, gene in fam2gene.items():
            scaf = gene2scaf.get(gene)
            if scaf is not None:
                fam2scaf[fam] = scaf
        lookup[sk] = fam2scaf
        print(f"  {sk}: {len(fam2gene)} fam→gene, {len(fam2scaf)} fam→scaffold")
    return lookup


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--rbh-dir", default=str(recompute_distances.DEFAULT_RBH_DIR))
    ap.add_argument("--genome-yaml",
                    default=str(REPO_ROOT / "genome_database/genome_list.yaml"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    cands = pd.read_csv(args.candidates, sep="\t")
    cfg = loci_io.load_genome_list(args.genome_yaml)
    rbh_dir = Path(args.rbh_dir)
    if not rbh_dir.exists():
        sys.exit(f"rbh dir not found: {rbh_dir}")

    # Build lookup once per anchor species (the slow part).
    all_anchors = sorted({sk for spp in recompute_distances.ANCHORS.values()
                          for sk in spp})
    print(f"building lookup for {len(all_anchors)} anchor species ...")
    fam2scaf = build_lookup(all_anchors, rbh_dir, cfg)

    # For each candidate pair, count anchors resolved per clade and the
    # number of distinct subclades covered.
    for clade, anchors in recompute_distances.ANCHORS.items():
        n_resolved = []
        n_subclades = []
        for fam1, fam2 in zip(cands["ortholog1"], cands["ortholog2"]):
            count = 0
            sclades = set()
            for sk in anchors:
                s = fam2scaf.get(sk, {})
                a1, a2 = s.get(fam1), s.get(fam2)
                if a1 is not None and a1 == a2:
                    count += 1
                    sub = SUBCLADE.get(sk)
                    if sub:
                        sclades.add(sub)
            n_resolved.append(count)
            n_subclades.append(len(sclades))
        cands[f"n_anchors_resolved_{clade}"] = n_resolved
        cands[f"n_subclades_resolved_{clade}"] = n_subclades
        print(f"  {clade}: median anchors_resolved="
              f"{pd.Series(n_resolved).median():.0f} "
              f"max={max(n_resolved)}")

    cands.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out} rows={len(cands)}")


if __name__ == "__main__":
    main()
