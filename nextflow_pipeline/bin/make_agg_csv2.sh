#!/bin/bash
set -euo pipefail

# Usage: make_agg_csv_by_suffix.sh counts_paths.txt output.csv
PATHS_FILE=$1
OUTPUT_FILE=$2

tmp=$(mktemp)

# Build "suffix<TAB>sample_id,molecule_h5" rows
tr ',' '\n' < "$PATHS_FILE" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' | sed '/^$/d' | \
while IFS= read -r path; do
  # Strip stray quotes if present
  path=${path%\"}; path=${path#\"}

  sample_id=$(basename "$path")
  molecule_h5="$path/outs/molecule_info.h5"

  # numeric suffix = digits after final underscore; fallback to huge number if missing
  suffix="${sample_id##*_}"
  [[ "$suffix" =~ ^[0-9]+$ ]] || suffix=999999

  printf '%s\t%s,%s\n' "$suffix" "$sample_id" "$molecule_h5" >> "$tmp"
done

# Sort by suffix (numeric) and emit CSV
{
  echo "sample_id,molecule_h5"
  sort -t $'\t' -k1,1n "$tmp" | cut -f2-
} > "$OUTPUT_FILE"

rm -f "$tmp"

# Show result
cat "$OUTPUT_FILE"
