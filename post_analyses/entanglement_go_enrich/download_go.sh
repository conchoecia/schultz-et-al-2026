#!/bin/bash
# One-shot helper: download the human GO annotation file and the GO ontology.
# Runs on a login node.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"

if [ ! -s "$OUT/goa_human.gaf.gz" ]; then
    curl -L -o "$OUT/goa_human.gaf.gz" \
        https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz
fi

if [ ! -s "$OUT/go-basic.obo" ]; then
    curl -L -o "$OUT/go-basic.obo" \
        https://current.geneontology.org/ontology/go-basic.obo
fi

echo "Downloaded:"
ls -la "$OUT"/goa_human.gaf.gz "$OUT"/go-basic.obo
