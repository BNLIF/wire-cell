#!/bin/bash
# Run wire-cell-display-convert for one or more nue root files.
#
# Usage:
#   ./magnify_convert.sh -e <event>          - process all files with <event>
#   ./magnify_convert.sh -r <run>            - process all events in <run>
#   ./magnify_convert.sh -r <run> -e <event> - process one specific event
#
# Input files are expected in the same directory as this script:
#   nue_<run>_<subrun>_<event>.root
#
# Output files are written to the same directory:
#   track_com_<run>_<subrun>_<event>.root

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONVERT="$SCRIPT_DIR/build/mcs/wire-cell-display-convert"

# ── parse flags ───────────────────────────────────────────────────────────────

usage() {
    echo "Usage:"
    echo "  $0 -e <event>          - process all files with <event>"
    echo "  $0 -r <run>            - process all events in <run>"
    echo "  $0 -r <run> -e <event> - process one specific event"
    exit 0
}

RUN=""
EVENT=""

while getopts ":r:e:h" flag; do
    case "$flag" in
        r) RUN="$OPTARG" ;;
        e) EVENT="$OPTARG" ;;
        h) usage ;;
        :) echo "[magnify] ERROR: -$OPTARG requires an argument" >&2; exit 1 ;;
       \?) echo "[magnify] ERROR: unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done

if [[ -z "$RUN" && -z "$EVENT" ]]; then
    echo "[magnify] ERROR: at least one of -r <run> or -e <event> is required" >&2
    usage
fi

if [[ ! -x "$CONVERT" ]]; then
    echo "[magnify] ERROR: convert binary not found or not executable: $CONVERT" >&2
    exit 1
fi

# ── helper ────────────────────────────────────────────────────────────────────

process_file() {
    local filepath="$1"
    local basename="${filepath##*/}"               # nue_5384_130_6505.root

    # Parse run, subrun, event directly from filename: nue_<run>_<subrun>_<event>.root
    local inner="${basename#nue_}"; inner="${inner%.root}"  # 5384_130_6505
    local file_run="${inner%%_*}"                  # 5384
    local rest="${inner#*_}"                       # 130_6505
    local subrun="${rest%%_*}"                     # 130
    local ev="${rest#*_}"                          # 6505

    local outfile="$SCRIPT_DIR/track_com_${file_run}_${subrun}_${ev}.root"

    echo "[magnify] Processing: $basename"
    echo "[magnify]   -> output: ${outfile##*/}"
    "$CONVERT" -b"$filepath" -o"$outfile" -f1
    echo "[magnify]   done."
}

# ── main ──────────────────────────────────────────────────────────────────────

# Build glob pattern from whichever of -r / -e were supplied
if [[ -n "$RUN" && -n "$EVENT" ]]; then
    PATTERN="nue_${RUN}_*_${EVENT}.root"
elif [[ -n "$RUN" ]]; then
    PATTERN="nue_${RUN}_*.root"
else
    PATTERN="nue_*_${EVENT}.root"
fi

mapfile -t MATCHES < <(find "$SCRIPT_DIR" -maxdepth 1 -name "$PATTERN" | sort)

if [[ ${#MATCHES[@]} -eq 0 ]]; then
    echo "[magnify] ERROR: no $PATTERN found in $SCRIPT_DIR" >&2
    exit 1
fi

echo "[magnify] Found ${#MATCHES[@]} file(s) matching: $PATTERN"

if [[ -n "$RUN" && -n "$EVENT" ]]; then
    # ── mode: run + event (expect exactly one) ────────────────────────────────
    for f in "${MATCHES[@]}"; do
        process_file "$f"
    done

else
    # ── mode: run-only or event-only → possibly multiple files ───────────────
    FAILED=()
    for f in "${MATCHES[@]}"; do
        process_file "$f" || { echo "[magnify] WARNING: failed on ${f##*/}" >&2; FAILED+=("$f"); }
    done

    if [[ ${#FAILED[@]} -gt 0 ]]; then
        echo "[magnify] WARNING: ${#FAILED[@]} file(s) failed."
        exit 1
    fi
fi

echo "[magnify] All done."
