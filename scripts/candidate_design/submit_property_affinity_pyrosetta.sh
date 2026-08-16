#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/../.."
RUN_KIND="${1:-}"
if [[ "${RUN_KIND}" != "pilot" && "${RUN_KIND}" != "full_scan" ]]; then
  echo "Usage: bash scripts/candidate_design/submit_property_affinity_pyrosetta.sh pilot|full_scan" >&2
  exit 2
fi

mkdir -p logs/property_affinity_pyrosetta
if [[ "${RUN_KIND}" == "full_scan" ]]; then
  python -c 'import json; p="docs/result_artifacts/candidate_design/property_affinity_pyrosetta_pilot_20260816/property_affinity_scoring_gate.json"; d=json.load(open(p)); assert d["status"] == "pass" and d["release"] == "ready_for_full_property_affinity_scan"'
fi

sbatch --export="ALL,RUN_KIND=${RUN_KIND}" scripts/candidate_design/submit_property_affinity_pyrosetta.slurm
