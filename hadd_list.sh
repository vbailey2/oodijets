#!/bin/bash
#
# hadd_list.sh - merge all ROOT files listed in a text file using hadd
#
# Usage:
#   ./hadd_list.sh <list_file> <output_file> [extra hadd args...]
#
# Example:
#   ./hadd_list.sh filelists/82703.txt merged/82703.root
#   ./hadd_list.sh filelists/82703.txt merged/82703.root -f -j 4

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <list_file> <output_file> [extra hadd args...]" >&2
    exit 1
fi

list_file="$1"
output_file="$2"
shift 2
extra_args=("$@")

if [ ! -f "$list_file" ]; then
    echo "Error: list file '$list_file' not found" >&2
    exit 1
fi

# Read files from list, skipping blank lines and comments (#)
mapfile -t files < <(grep -v -e '^[[:space:]]*$' -e '^[[:space:]]*#' "$list_file")

if [ "${#files[@]}" -eq 0 ]; then
    echo "Error: no files found in '$list_file'" >&2
    exit 1
fi

# Check all input files exist before running hadd
missing=0
for f in "${files[@]}"; do
    if [ ! -f "$f" ]; then
        echo "Warning: file not found: $f" >&2
        missing=1
    fi
done

if [ "$missing" -eq 1 ]; then
    echo "One or more input files are missing. Aborting." >&2
    exit 1
fi

mkdir -p "$(dirname "$output_file")"

echo "Merging ${#files[@]} files into $output_file ..."
hadd "${extra_args[@]}" "$output_file" "${files[@]}"
echo "Done."
