#!/bin/bash
# One-shot: fetch NCBI's gene2accession + gene2go, filtered to Homo
# sapiens (tax_id 9606). Streams through awk to avoid staging the full
# multi-GB globals.
#
# Produces:
#   ncbi_ref/human_gene2accession.tsv.gz  (~5 MB) — protein acc → GeneID
#   ncbi_ref/human_gene2go.tsv.gz         (~3 MB) — GeneID → GO term
#
# Run on a login node with outbound HTTPS.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$SCRIPT_DIR/ncbi_ref"
mkdir -p "$OUT"

BASE=https://ftp.ncbi.nlm.nih.gov/gene/DATA

fetch_filter() {
    local url="$1" target="$2"
    if [ -s "$target" ]; then
        echo "  [present] $(basename "$target")"
        return
    fi
    echo "  [fetch  ] $(basename "$target")  from $url"
    curl -L --fail --silent --show-error "$url" | zcat | \
        awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1=="9606" {print}' | \
        gzip > "$target.part"
    mv "$target.part" "$target"
}

fetch_filter "$BASE/gene2accession.gz" "$OUT/human_gene2accession.tsv.gz"
fetch_filter "$BASE/gene2go.gz"         "$OUT/human_gene2go.tsv.gz"

echo "Done."
ls -la "$OUT"
