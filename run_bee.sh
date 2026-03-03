#!/bin/bash
# Wrapper around dump_json_parallel.py + upload_to_bee.sh
#
# Usage:
#   ./run_bee.sh [OPTIONS] [event1 event2 ...]
#
#   No arguments        - run on ../nue_01.root (from bee/'s perspective) with
#                         default dump options
#   event1 event2 ...   - collect ALL nue*_<event>.root files matching every
#                         given event number into ONE temporary directory, run
#                         dump_json_parallel.py once, upload, then clean up
#
# Options:
#   -o "opt1 opt2 ..."  - override dump options (default: simple charge deblob mc flash cluster)
#   -h                  - show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BEE_DIR="$SCRIPT_DIR/bee"
OPTIONS="simple charge deblob mc flash cluster"   # default dump options

# ── parse flags ──────────────────────────────────────────────────────────────

usage() {
    sed -n '3,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 0
}

while getopts ":o:h" flag; do
    case "$flag" in
        o) OPTIONS="$OPTARG" ;;
        h) usage ;;
        :) echo "[run_bee] ERROR: -$OPTARG requires an argument" >&2; exit 1 ;;
       \?) echo "[run_bee] ERROR: unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done
shift $(( OPTIND - 1 ))
EVENTS=("$@")

# ── helpers ──────────────────────────────────────────────────────────────────

# rootfile path must be relative to bee/ (or absolute)
run_dump() {
    local rootfile="$1"
    echo "[run_bee] Running dump_json_parallel.py on: $rootfile"
    echo "[run_bee] Dump options: $OPTIONS"
    cd "$BEE_DIR"
    python dump_json_parallel.py "$rootfile" $OPTIONS
}

upload() {
    echo "[run_bee] Uploading to_upload.zip ..."
    cd "$BEE_DIR"
    bash upload-to-bee.sh to_upload.zip
}

# ── main ─────────────────────────────────────────────────────────────────────

if [ ${#EVENTS[@]} -eq 0 ]; then
    # ── default mode: single file, from bee/'s perspective ../nue_01.root ────
    run_dump "../nue_01.root"
    upload

else
    # ── event-number mode ────────────────────────────────────────────────────
    # Gather ALL files matching every requested event into one temp directory,
    # then run dump_json_parallel.py once (its glob picks up all files in the dir).

    TMPDIR_PATH="$(mktemp -d "$BEE_DIR/tmp_bee_XXXXXX")"
    echo "[run_bee] Temporary directory: $TMPDIR_PATH"

    FIRST_FILE=""
    FAILED=()

    for EVENT in "${EVENTS[@]}"; do
        mapfile -t MATCHES < <(find "$SCRIPT_DIR" -maxdepth 1 -name "nue*_${EVENT}.root" | sort)

        if [ ${#MATCHES[@]} -eq 0 ]; then
            echo "[run_bee] WARNING: no nue*_${EVENT}.root found in $SCRIPT_DIR -- skipping" >&2
            FAILED+=("$EVENT")
            continue
        fi

        echo "[run_bee] Event ${EVENT}: found ${#MATCHES[@]} file(s):"
        for f in "${MATCHES[@]}"; do
            echo "[run_bee]   Copying: $f"
            cp "$f" "$TMPDIR_PATH/"
            # Remember the very first copied file to pass to dump_json_parallel.py
            [[ -z "$FIRST_FILE" ]] && FIRST_FILE="${f##*/}"
        done
    done

    if [[ -z "$FIRST_FILE" ]]; then
        echo "[run_bee] ERROR: no matching files found for any event -- aborting" >&2
        rm -rf "$TMPDIR_PATH"
        exit 1
    fi

    if [ ${#FAILED[@]} -gt 0 ]; then
        echo ""
        echo "[run_bee] WARNING: no files found for event(s): ${FAILED[*]}"
    fi

    # Path to the first file, relative to bee/
    TMPDIR_REL="${TMPDIR_PATH##$BEE_DIR/}"
    run_dump "${TMPDIR_REL}/${FIRST_FILE}"
    upload

    # Clean up
    echo "[run_bee] Removing temporary directory: $TMPDIR_PATH"
    rm -rf "$TMPDIR_PATH"

    # Exit with error if any events were missing
    [ ${#FAILED[@]} -gt 0 ] && exit 1
fi
