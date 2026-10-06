#!/bin/bash
#
# unicon-sdf.sh -- Perceive bond orders and formal charges for every molecule
#                  that sdf2pdbs.py wrote into the scratch directory.
#
# Unicon (Schroedinger) is called one molecule at a time, preserving the
# original crystal coordinates (-g 1 = keep the input geometry, single
# conformer). The resulting SDFs are what assemble_cell.py stitches back
# together into the labelled cell.
#
# Usage:
#   unicon-sdf.sh [TMP_DIR] [UNICON_BIN]
#
# Defaults: TMP_DIR=$HOME/2507/tmp, UNICON_BIN=unicon

set -uo pipefail

TMP_DIR="${1:-$HOME/2507/tmp}"
UNICON_BIN="${2:-unicon}"

if [ ! -d "$TMP_DIR" ]; then
    echo "error: scratch directory $TMP_DIR does not exist" >&2
    exit 1
fi

if ! command -v "$UNICON_BIN" > /dev/null 2>&1; then
    echo "error: '$UNICON_BIN' not found on PATH." >&2
    echo "Unicon ships with Schroedinger; source its env.sh first, e.g." >&2
    echo "  source /path/to/schrodinger/env.sh" >&2
    exit 1
fi

count=0
find "$TMP_DIR" -maxdepth 1 -type f -name '*.pdb' | LC_COLLATE=C sort | while read -r pdb_file; do
    pdb_file="$(readlink -f "$pdb_file")"
    pdb_dir="$(dirname "$pdb_file")"
    sdf_file="${pdb_file%.pdb}.sdf"

    cd "$pdb_dir" || exit 1

    # -t single : one molecule  -p single : one protonation state
    # -g 1      : keep the input coordinates (do NOT re-embed)
    if "$UNICON_BIN" -i "$(basename "$pdb_file")" -o "$(basename "$sdf_file")" \
                    -t single -p single -g 1 > /dev/null 2>&1; then
        echo "ok:   $(basename "$pdb_file")"
    else
        echo "fail: $(basename "$pdb_file")" >&2
    fi
done

echo "Generated unicon-sdf files in $TMP_DIR"
