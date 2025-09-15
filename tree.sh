#!/usr/bin/env bash
# Print concise TEPEAK layout without spamming giant dirs.
set -euo pipefail

ROOT="${1:-.}"

echo "=== TOP-LEVEL ==="
if command -v tree >/dev/null 2>&1; then
  # Hide common bulky/noisy dirs
  tree -a -L 2 -I '.git|.snakemake|__pycache__|.DS_Store|prefetch_tmp|picard/build|output|data|env|.venv|.gradle' "$ROOT"
else
  # fallback if 'tree' isn't available
  find "$ROOT" -maxdepth 2 -mindepth 1 \
    -not -path '*/.git*' \
    -not -path '*/.snakemake*' \
    -not -path '*/__pycache__*' \
    -not -path '*/.gradle*' \
    -not -name '.DS_Store' \
    -printf '%P\n' | sort
fi

echo
echo "=== data/ (first 3 levels; pruning fastq + big extracted ncbi_dataset) ==="
if [ -d "$ROOT/data" ]; then
  find "$ROOT/data" -maxdepth 3 \
    -not -path '*/fastq/*' \
    -not -path '*/ncbi_dataset/*' \
    -not -name '.DS_Store' \
    -printf '%P\n' | sort
else
  echo "data/ (missing)"
fi

echo
echo "=== output/ (first 3 levels) ==="
if [ -d "$ROOT/output" ]; then
  find "$ROOT/output" -maxdepth 3 -printf '%P\n' | sort
else
  echo "output/ (missing)"
fi

echo
echo "=== species dirs (indexes & key files) ==="
for spdir in "$ROOT"/data/*/ ; do
  [ -d "$spdir" ] || continue
  spname=$(basename "$spdir")
  echo "-- $spname --"
  ls -lh "$spdir" 2>/dev/null | grep -E '(\.fa(\..+)?$|\.gtf$|\.zip$)' || true
done

