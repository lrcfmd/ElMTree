#!/usr/bin/env bash
set -euo pipefail

# Find the nearest neighbours of one composition (no S3 required).
# Usage: ./scripts/run_instance.sh [composition] [k]
#
# Example:
#   ./scripts/run_instance.sh LiFePO4 10
# Needs the index in src/model/artifacts/ (see README: uv run build-index ...).

COMPOSITION="${1:-LiFePO4}"
K="${2:-10}"

uv run model-run --local \
  --parameters "$(jq -n \
    --arg comp "${COMPOSITION}" --argjson k "${K}" \
    '{mode:"instance", inputs:{composition:{value:$comp}, k:{value:$k}}, parameters:{}}')"
