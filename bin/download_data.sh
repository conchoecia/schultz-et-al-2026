#!/bin/bash
# Fetch the Dryad-hosted data bundle for this pipeline and unpack it
# into $DATA_DIR (default: $REPO_ROOT/dryad_data/).
#
# Two tarballs + the authors' extract script + a checksums file come
# down. `extract_data.sh` (from Dryad) verifies the checksums and
# unpacks both tarballs; we defer to it so the integrity check stays
# exactly the same one the data producer uses.
#
# Override any file-stream ID via env var if the Dryad URLs migrate
# after publication. The current IDs point at the private review copy.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR="${DATA_DIR:-$REPO_ROOT/dryad_data}"
DRYAD_BASE="${DRYAD_BASE:-https://datadryad.org/downloads/file_stream}"

# File-stream IDs (override via env). Current values are the private-review copy.
RBH_ID="${RBH_ID:-4623678}"             # BCnSSimakov2022_current_rbh_202509.tar.gz    568 MB
WORKFLOW_ID="${WORKFLOW_ID:-4623679}"   # newick_and_timetree_20251118.tar.gz          643 MB
EXTRACT_ID="${EXTRACT_ID:-4623680}"     # extract_data.sh                                3 KB
CHECKSUMS_ID="${CHECKSUMS_ID:-4623681}" # checksums.md5                                146 B

# Filename on disk → Dryad ID.
declare -A FILES=(
    ["BCnSSimakov2022_current_rbh_202509.tar.gz"]="$RBH_ID"
    ["newick_and_timetree_20251118.tar.gz"]="$WORKFLOW_ID"
    ["extract_data.sh"]="$EXTRACT_ID"
    ["checksums.md5"]="$CHECKSUMS_ID"
)

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "==> Destination: $DATA_DIR"
echo

for name in "${!FILES[@]}"; do
    id="${FILES[$name]}"
    url="$DRYAD_BASE/$id"
    if [ -s "$name" ]; then
        echo "    [present] $name"
    else
        echo "    [fetch  ] $name  ($url)"
        # -L follow redirects, --fail on HTTP error, --continue-at - for resumable,
        # --retry for transient network flakes.
        curl -L --fail --retry 5 --continue-at - -o "$name" "$url"
    fi
done

echo
echo "==> Extracting (delegated to Dryad's extract_data.sh)"
chmod +x extract_data.sh
bash extract_data.sh

echo
echo "==> Done."
echo "   Per-species RBH dir : $DATA_DIR/BCnSSimakov2022_current_rbh_202509/"
echo "   Workflow outputs    : $DATA_DIR/newick_and_timetree_20251118/"
echo
echo "   Config defaults in config.template.yaml already point at these paths."
