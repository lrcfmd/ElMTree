#!/usr/bin/env bash
set -euo pipefail

# Nearest-neighbour search for every composition in a CSV with a 'composition'
# column (no S3 required).
# Usage: ./scripts/run_dataset.sh [csv_file] [k]
#
# Example:
#   ./scripts/run_dataset.sh data/input.csv 5
# The results are written to ./data/output/neighbours.csv

CSV_FILE="${1:-data/input.csv}"
K="${2:-10}"

uv run model-run --local \
  --parameters "$(jq -n \
    --arg file "file://${CSV_FILE}" --argjson k "${K}" \
    '{mode:"dataset", inputs:{file:{uri:$file, mime_type:"text/csv"}, k:{value:$k}}, parameters:{}}')"
