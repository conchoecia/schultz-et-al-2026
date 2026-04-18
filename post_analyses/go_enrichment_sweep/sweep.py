#!/usr/bin/env python3
"""
Per-clade GO enrichment sweep over locus-pair rank thresholds.

Reads Dryad's SupplementaryTable_16 (per-clade locus pairs with z-scored
stability/closeness statistics), ranks each clade's pairs along two axes
(stability, closeness), and sweeps a top-N threshold to show how GO
enrichment varies as the pair set grows from most-constrained to
all-flagged. The full egt entanglement-go-enrich command operates at
ALG-letter granularity, which dilutes each clade's foreground to near
the background — this script works at the BCnS-family level.

Annotation join uses NCBI's gene2accession + gene2go (Entrez GeneID
keyed), not GOA GAF, because the family map's RefSeq NP_/XP_ protein
accessions map cleanly to GeneIDs via gene2accession.

Emits per-clade TSVs, cross-clade summary, every term reaching q <= 0.25,
and curves.pdf of -log10(top q) vs N.
"""
import argparse
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Statistics — native hypergeometric (one-tailed upper) + BH-FDR.
# ---------------------------------------------------------------------
def log_binom(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N_total, K, n):
    """P(X >= k) with X ~ Hypergeometric(N_total, K, n).

    Exact sum in log-space to dodge overflow.
    """
    if k <= 0:
        return 1.0
    if k > min(K, n):
        return 0.0
    denom = log_binom(N_total, n)
    logs = [log_binom(K, i) + log_binom(N_total - K, n - i) - denom
            for i in range(k, min(K, n) + 1)]
    m = max(logs)
    return math.exp(m) * sum(math.exp(x - m) for x in logs)


def bh_qvalues(pvals):
    """Benjamini-Hochberg. Returns ndarray of q-values aligned with input."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q_sorted = np.minimum(q_sorted, 1.0)
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return q


# ---------------------------------------------------------------------
# NCBI annotation parsers
# ---------------------------------------------------------------------
NCBI_CATEGORY_MAP = {
    "Process": "BP",
    "Function": "MF",
    "Component": "CC",
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

REFSEQ_PROTEIN_RE = re.compile(r"^(NP|XP|YP)_[0-9]+\.[0-9]+$")
UNIPROT_SP_RE = re.compile(r"^sp\|([^|]+)\|([^|]+)$")


def parse_gene2accession(path):
    """Return {protein_accession.version: GeneID:str} for all rows with a
    non-empty protein_accession.version in the gene2accession table.

    Expected columns (tab):
      0 tax_id, 1 GeneID, 2 status, 3 RNA_nucleotide_accession.version,
      4 RNA_nucleotide_gi, 5 protein_accession.version, 6 protein_gi, ...
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    prot_to_gene = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            prot, gene = f[5], f[1]
            if prot and prot != "-" and gene and gene != "-":
                prot_to_gene[prot] = gene
    return prot_to_gene


def parse_gene2go(path):
    """Return (gene_to_terms, term_namespace) from gene2go.

    Expected columns (tab):
      0 tax_id, 1 GeneID, 2 GO_ID, 3 Evidence, 4 Qualifier,
      5 GO_term, 6 Pubmed, 7 Category
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    gene_to_terms = defaultdict(set)
    term_ns = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            gene, go_id, qualifier, category = f[1], f[2], f[4], f[7]
            if "NOT" in qualifier.upper().split("|"):
                continue
            gene_to_terms[gene].add(go_id)
            ns = NCBI_CATEGORY_MAP.get(category)
            if ns and go_id not in term_ns:
                term_ns[go_id] = ns
    return gene_to_terms, term_ns


# ---------------------------------------------------------------------
# Family → GeneID mapping (from bcns_family_to_human_gene.tsv)
# ---------------------------------------------------------------------
def parse_family_map(path, prot_to_gene):
    """Return (fam_to_genes, stats) with stats reporting join coverage.

    For each family row, the `human_gene` cell is expected to be an
    NP_/XP_/YP_ accession (from the GCF_000001405 RefSeq human proteome).
    Maps to a GeneID via prot_to_gene; drops entries that don't map.
    """
    df = pd.read_csv(path, sep="\t")
    family_col = df.columns[0]
    gene_col = "human_gene" if "human_gene" in df.columns else df.columns[2]
    fam_to_genes = defaultdict(set)
    n_rows = 0
    n_empty = 0
    n_refseq_hit = 0
    n_refseq_miss = 0
    n_sp = 0
    n_other = 0
    for _, row in df.iterrows():
        n_rows += 1
        fam = row[family_col]
        val = row[gene_col]
        if not isinstance(val, str) or not val.strip():
            n_empty += 1
            continue
        v = val.strip()
        if REFSEQ_PROTEIN_RE.match(v):
            g = prot_to_gene.get(v)
            if g:
                fam_to_genes[fam].add(g)
                n_refseq_hit += 1
            else:
                # Version mismatch: try stripping the .N
                base = v.rsplit(".", 1)[0]
                # Search prot_to_gene for any key with that base
                # (rare; skip for now — keeps this O(1))
                n_refseq_miss += 1
        elif UNIPROT_SP_RE.match(v):
            n_sp += 1  # pass-2 entries; not mappable via gene2accession alone
        else:
            n_other += 1
    stats = dict(
        n_rows=n_rows,
        n_empty=n_empty,
        n_refseq_hit=n_refseq_hit,
        n_refseq_miss=n_refseq_miss,
        n_sp_uniprot=n_sp,
        n_other=n_other,
        n_families_mapped=len(fam_to_genes),
    )
    return fam_to_genes, stats


# ---------------------------------------------------------------------
# Enrichment core
# ---------------------------------------------------------------------
def enrich_for_foreground(foreground, background_to_terms, term_namespace,
                          namespaces=("all", "BP", "MF", "CC"),
                          min_term_hits=2):
    """One-shot hypergeometric + BH-FDR for a single foreground.

    background_to_terms: {gene_key: set(GO_ID)} restricted to the background.
    Returns {namespace: [dict(go_id, k, K, n, N, fold, p, q, namespace)]}.
    """
    fg = {g for g in foreground if g in background_to_terms}
    n = len(fg)
    N_total = len(background_to_terms)
    if n == 0 or N_total == 0:
        return {ns: [] for ns in namespaces}

    term_K = defaultdict(int)
    for g, terms in background_to_terms.items():
        for t in terms:
            term_K[t] += 1
    term_k = defaultdict(int)
    for g in fg:
        for t in background_to_terms[g]:
            term_k[t] += 1

    rows = []
    for t, k in term_k.items():
        if k < min_term_hits:
            continue
        K = term_K[t]
        p = hypergeom_sf(k, N_total, K, n)
        fold = (k / n) / (K / N_total) if K > 0 else float("inf")
        rows.append((t, k, K, n, N_total, fold, p))

    out_by_ns = {ns: [] for ns in namespaces}
    if not rows:
        return out_by_ns

    # Per-namespace BH correction.
    for ns in namespaces:
        if ns == "all":
            sub = rows
        else:
            sub = [r for r in rows if term_namespace.get(r[0]) == ns]
        if not sub:
            continue
        pvals = np.array([r[6] for r in sub])
        qvals = bh_qvalues(pvals)
        enriched = [
            dict(go_id=t, k=k, K=K, n=ntot_, N=N_, fold=fold, p=p, q=q,
                 term_namespace=term_namespace.get(t, "?"))
            for (t, k, K, ntot_, N_, fold, p), q in zip(sub, qvals)
        ]
        enriched.sort(key=lambda d: d["q"])
        out_by_ns[ns] = enriched
    return out_by_ns


# ---------------------------------------------------------------------
# Per-clade sweep
# ---------------------------------------------------------------------
OCCUPANCY_MIN = 0.5  # matches Dryad's close_in_clade / stable_in_clade threshold
N_GRID = [5, 10, 20, 50, 100, 200, 500, 1000, 2000]


def sweep_clade(clade_rows, fam_to_genes, background_to_terms, term_namespace):
    """Run the (N × axis × namespace) sweep for one clade."""
    df = clade_rows.copy()
    df = df[df["occupancy_in"].fillna(0) >= OCCUPANCY_MIN]
    df = df.dropna(subset=["sd_in_out_ratio_log_sigma",
                           "mean_in_out_ratio_log_sigma"])
    df = df.reset_index(drop=True)
    if df.empty:
        return [], {}

    stability_order = df.sort_values("sd_in_out_ratio_log_sigma").index.to_numpy()
    closeness_order = df.sort_values("mean_in_out_ratio_log_sigma").index.to_numpy()
    n_rows = len(df)
    n_grid = sorted(set([x for x in N_GRID if x <= n_rows] + [n_rows]))

    records = []
    curve_data = defaultdict(list)
    for N in n_grid:
        for axis in ("stability", "closeness", "intersection"):
            if axis == "stability":
                idxs = set(stability_order[:N].tolist())
            elif axis == "closeness":
                idxs = set(closeness_order[:N].tolist())
            else:
                idxs = (set(stability_order[:N].tolist())
                        & set(closeness_order[:N].tolist()))
            if not idxs:
                continue
            sub = df.loc[list(idxs)]
            families = pd.concat([sub["ortholog1"], sub["ortholog2"]]).dropna().unique()
            foreground = set()
            for fam in families:
                foreground |= fam_to_genes.get(fam, set())
            if not foreground:
                continue
            by_ns = enrich_for_foreground(foreground, background_to_terms,
                                           term_namespace)
            for ns in ("all", "BP", "MF", "CC"):
                res = by_ns.get(ns, [])
                n_q05 = sum(1 for r in res if r["q"] <= 0.05)
                n_q25 = sum(1 for r in res if r["q"] <= 0.25)
                top = res[0] if res else None
                top_q = top["q"] if top else float("nan")
                records.append(dict(
                    axis=axis, N_threshold=N, pairs_used=len(idxs),
                    namespace=ns,
                    foreground_size=len(foreground),
                    n_families=len(families),
                    n_terms_tested=len(res),
                    n_hits_q05=n_q05,
                    n_hits_q25=n_q25,
                    top_term=top["go_id"] if top else "",
                    top_term_fold=top["fold"] if top else float("nan"),
                    top_term_k=top["k"] if top else 0,
                    top_term_K=top["K"] if top else 0,
                    top_q=top_q,
                ))
                if not math.isnan(top_q) and top_q > 0:
                    curve_data[(axis, ns)].append((N, -math.log10(top_q)))
                elif top and top_q == 0:
                    curve_data[(axis, ns)].append((N, 300.0))
    return records, curve_data


def harvest_significant_terms(clade, clade_rows, records,
                               fam_to_genes, background_to_terms,
                               term_namespace):
    """For configs that produced q<=0.25 hits, re-run enrichment at the
    smallest such N to collect the per-term rows for the sig-terms TSV.
    """
    out = []
    df = clade_rows.copy()
    df = df[df["occupancy_in"].fillna(0) >= OCCUPANCY_MIN]
    df = df.dropna(subset=["sd_in_out_ratio_log_sigma",
                           "mean_in_out_ratio_log_sigma"])
    df = df.reset_index(drop=True)
    if df.empty:
        return out
    stability_order = df.sort_values("sd_in_out_ratio_log_sigma").index.to_numpy()
    closeness_order = df.sort_values("mean_in_out_ratio_log_sigma").index.to_numpy()

    for axis in ("stability", "closeness", "intersection"):
        cand = [r for r in records
                if r["axis"] == axis and r["namespace"] == "all"
                and r["n_hits_q25"] > 0]
        if not cand:
            continue
        cand.sort(key=lambda r: r["N_threshold"])
        r = cand[0]
        N = r["N_threshold"]
        if axis == "stability":
            idxs = stability_order[:N]
        elif axis == "closeness":
            idxs = closeness_order[:N]
        else:
            idxs = list(set(stability_order[:N].tolist())
                        & set(closeness_order[:N].tolist()))
        s2 = df.loc[list(idxs)]
        families = pd.concat([s2["ortholog1"], s2["ortholog2"]]).dropna().unique()
        fg = set()
        for f in families:
            fg |= fam_to_genes.get(f, set())
        ns_res = enrich_for_foreground(fg, background_to_terms, term_namespace)
        for ns, rows in ns_res.items():
            for rr in rows:
                if rr["q"] > 0.25:
                    continue
                out.append(dict(
                    clade=clade, axis=axis, N_threshold=N,
                    sweep_namespace=ns,
                    go_id=rr["go_id"], go_namespace=rr["term_namespace"],
                    k=rr["k"], K=rr["K"], n=rr["n"], N=rr["N"],
                    fold=rr["fold"], p=rr["p"], q=rr["q"],
                ))
    return out


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--supp-table", required=True,
                    help="Dryad SupplementaryTable_16.xlsx")
    ap.add_argument("--family-map", required=True,
                    help="bcns_family_to_human_gene.tsv")
    ap.add_argument("--gene2accession", required=True,
                    help="human-filtered gene2accession.gz")
    ap.add_argument("--gene2go", required=True,
                    help="human-filtered gene2go.gz")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    (out_dir / "per_clade").mkdir(parents=True, exist_ok=True)

    print(f"[load] SupplementaryTable_16: {args.supp_table}")
    df = pd.read_excel(args.supp_table, engine="openpyxl")
    print(f"  shape={df.shape}  clades={df['nodename'].nunique()}")
    required_cols = ("nodename", "ortholog1", "ortholog2", "occupancy_in",
                     "sd_in_out_ratio_log_sigma", "mean_in_out_ratio_log_sigma")
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: xlsx missing required columns: {missing}")

    print(f"[load] gene2accession: {args.gene2accession}")
    prot_to_gene = parse_gene2accession(args.gene2accession)
    print(f"  protein-accessions -> GeneID rows={len(prot_to_gene)}")

    print(f"[load] family map: {args.family_map}")
    fam_to_genes, stats = parse_family_map(args.family_map, prot_to_gene)
    print(f"  rows={stats['n_rows']}  "
          f"refseq_hit={stats['n_refseq_hit']}  "
          f"refseq_miss={stats['n_refseq_miss']}  "
          f"sp={stats['n_sp_uniprot']}  "
          f"other={stats['n_other']}  "
          f"empty={stats['n_empty']}  "
          f"mapped_families={stats['n_families_mapped']}")

    print(f"[load] gene2go: {args.gene2go}")
    gene_to_terms_all, term_namespace = parse_gene2go(args.gene2go)
    print(f"  GeneIDs-with-GO-annotations={len(gene_to_terms_all)}  "
          f"terms-in-namespace-map={len(term_namespace)}")

    # Background = GeneIDs that (a) appear in the family map AND (b) have
    # any GO annotation. This matches the hypergeom semantics: the
    # population from which foreground is drawn is the annotatable subset
    # of the mapped universe.
    background_gene_ids = set().union(*fam_to_genes.values())
    background_to_terms = {}
    for g in background_gene_ids:
        t = gene_to_terms_all.get(g)
        if t:
            background_to_terms[g] = t
    print(f"  background_family-map-GeneIDs={len(background_gene_ids)}  "
          f"background_with_any_term={len(background_to_terms)}")

    if len(background_to_terms) < 100:
        print("WARN: background is tiny; enrichment power will be low.",
              file=sys.stderr)

    # Per-clade sweep
    all_records = []
    all_curves = {}
    all_significant = []
    for clade in sorted(df["nodename"].dropna().unique()):
        sub = df[df["nodename"] == clade]
        records, curves = sweep_clade(sub, fam_to_genes,
                                       background_to_terms, term_namespace)
        any_hit = any(r["n_hits_q25"] > 0 for r in records)
        print(f"[clade] {clade}  rows={len(sub)}  sweep_configs={len(records)}  "
              f"any_q25_hit={any_hit}")
        if records:
            cdf = pd.DataFrame(records)
            cdf.insert(0, "clade", clade)
            cdf.to_csv(out_dir / "per_clade" / f"{clade}.tsv",
                       sep="\t", index=False)
            all_records.extend(cdf.to_dict("records"))
        if curves:
            all_curves[clade] = curves
        all_significant.extend(harvest_significant_terms(
            clade, sub, records, fam_to_genes,
            background_to_terms, term_namespace))

    # Summary
    if all_records:
        sdf = pd.DataFrame(all_records).sort_values(
            ["n_hits_q25", "n_hits_q05"], ascending=False)
        sdf.to_csv(out_dir / "summary.tsv", sep="\t", index=False)
        any_q25 = int((sdf["n_hits_q25"] > 0).sum())
        any_q05 = int((sdf["n_hits_q05"] > 0).sum())
        print(f"[write] summary.tsv  rows={len(sdf)}  "
              f"configs_with_q25_hits={any_q25}  "
              f"configs_with_q05_hits={any_q05}")

    if all_significant:
        sig_df = pd.DataFrame(all_significant).drop_duplicates(
            subset=["clade", "axis", "N_threshold", "sweep_namespace", "go_id"]
        ).sort_values(["clade", "q"])
        sig_df.to_csv(out_dir / "significant_terms.tsv", sep="\t", index=False)
        print(f"[write] significant_terms.tsv  rows={len(sig_df)}")
    else:
        pd.DataFrame(columns=["clade", "axis", "N_threshold", "sweep_namespace",
                              "go_id", "go_namespace", "k", "K", "n", "N",
                              "fold", "p", "q"]
                     ).to_csv(out_dir / "significant_terms.tsv",
                              sep="\t", index=False)
        print("[write] significant_terms.tsv  (empty — no q<=0.25 anywhere)")

    # Curves
    if all_curves:
        n_clades = len(all_curves)
        fig, axes = plt.subplots(n_clades, 3, figsize=(15, 3 * n_clades),
                                  squeeze=False)
        colors = {"all": "black", "BP": "C0", "MF": "C1", "CC": "C2"}
        for row, (clade, cdata) in enumerate(sorted(all_curves.items())):
            for col, axis in enumerate(("stability", "closeness", "intersection")):
                ax = axes[row][col]
                for ns in ("all", "BP", "MF", "CC"):
                    pts = cdata.get((axis, ns), [])
                    if not pts:
                        continue
                    xs, ys = zip(*pts)
                    ax.plot(xs, ys, marker="o", markersize=3, label=ns,
                            color=colors[ns])
                ax.axhline(-math.log10(0.05), ls="--", color="red", lw=0.5,
                           label="q=0.05" if row == 0 and col == 0 else None)
                ax.axhline(-math.log10(0.25), ls=":", color="orange", lw=0.5,
                           label="q=0.25" if row == 0 and col == 0 else None)
                ax.set_xscale("log")
                ax.set_title(f"{clade} — {axis}")
                if col == 0:
                    ax.set_ylabel("-log10(top q)")
                if row == n_clades - 1:
                    ax.set_xlabel("top-N threshold")
                if row == 0 and col == 0:
                    ax.legend(fontsize=7, loc="best")
        fig.tight_layout()
        fig.savefig(out_dir / "curves.pdf")
        plt.close(fig)
        print(f"[write] curves.pdf  clades={n_clades}")

    print("[done]")


if __name__ == "__main__":
    main()
