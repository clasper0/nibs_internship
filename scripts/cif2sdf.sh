#!/bin/bash
#
# cif2sdf.sh -- CIF -> SDF conversion using Jmol, as a licence-free fallback
#               for the Schroedinger ``structconvert`` step.
#
# Jmol expands the unit cell ({3 3 3} = the symmetry operators of the CIF) and
# writes a V3000 SDF. This path was used before the Schroedinger licence was
# available; ``structconvert`` gives better bond perception, so prefer it when
# you have it.
#
# Usage:
#   cif2sdf.sh [CIF_ROOT_DIR] [OUTPUT_DIR] [COD_ID_LIST]
#
# Requires: a java binary and jmol.jar. Set JMOL_JAR / JAVA_BIN below or in
# the environment.

set -uo pipefail

JMOL_JAR="${JMOL_JAR:-$HOME/softwares/jmol/jmol-16.3.31/Jmol.jar}"
JAVA_BIN="${JAVA_BIN:-java}"

CIF_ROOT_DIR="${1:-$HOME/2507/cif}"
OUTPUT_DIR="${2:-$HOME/2507/sdf}"
COD_ID_LIST="${3:-$CIF_ROOT_DIR/cod_ids.txt}"

if [ ! -f "$JMOL_JAR" ]; then
    echo "error: Jmol.jar not found at $JMOL_JAR" >&2
    echo "Set JMOL_JAR to its real location." >&2
    exit 1
fi

if [ ! -f "$COD_ID_LIST" ]; then
    echo "error: COD id list $COD_ID_LIST not found" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

LOG_SUCCESS="$OUTPUT_DIR/success.log"
LOG_FAILED="$OUTPUT_DIR/failed.log"
LOG_SKIPPED="$OUTPUT_DIR/skipped.log"
: > "$LOG_SUCCESS"
: > "$LOG_FAILED"
: > "$LOG_SKIPPED"

echo "Converting CIF -> SDF"
echo "  CIF dir : $CIF_ROOT_DIR"
echo "  SDF dir : $OUTPUT_DIR"
printf '%.0s=' {1..50}; echo

while read -r cod_id; do
    [ -z "$cod_id" ] && continue

    # COD stores files sharded as <d1>/<d2d3>/<d4d5>/<id>.cif
    cif_file="$CIF_ROOT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}/$cod_id.cif"
    sdf_dir="$OUTPUT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}"
    new_sdf_file="$sdf_dir/$cod_id.sdf"
    jmol_file="$(mktemp --suffix=.jmol)"

    mkdir -p "$sdf_dir"

    if [ -f "$new_sdf_file" ]; then
        echo "skip $cod_id"
        echo "$cod_id" >> "$LOG_SKIPPED"
        rm -f "$jmol_file"
        continue
    fi

    {
        echo 'set defaultLoadScript "unitcell"'
        echo "load '$cif_file' {3 3 3}"
        echo "write V3000 '$new_sdf_file'"
        echo "quit"
    } > "$jmol_file"

    if timeout 30 "$JAVA_BIN" -jar "$JMOL_JAR" -n -s "$jmol_file" > /dev/null 2>&1; then
        echo "success: $cod_id"
        echo "$cod_id" >> "$LOG_SUCCESS"
    else
        echo "failed: $cod_id"
        echo "$cod_id" >> "$LOG_FAILED"
        rm -f "$new_sdf_file"
    fi

    rm -f "$jmol_file"
    sleep 0.1
done < "$COD_ID_LIST"

printf '%.0s=' {1..50}; echo
success_count=$(wc -l < "$LOG_SUCCESS")
failed_count=$(wc -l < "$LOG_FAILED")
skipped_count=$(wc -l < "$LOG_SKIPPED")
echo "Conversion complete."
echo "  success=$success_count skipped=$skipped_count failed=$failed_count"
echo "  logs in $OUTPUT_DIR"
