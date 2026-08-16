#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/../.."

python -c 'import json; d=json.load(open("docs/result_artifacts/candidate_design/targeted_structure_review_plan_20260816/targeted_structure_review_contract.json")); assert d["status"] == "pass" and d["release"] == "ready_for_remote_targeted_structure_review" and d["candidate_count"] == 9 and d["hard_exclusion_count"] == 7'

for path in \
  results/candidate_design/targeted_structure_review_runtime_20260816 \
  docs/result_artifacts/candidate_design/targeted_structure_review_result_20260816 \
  docs/run_summaries/candidate_design/targeted_structure_review_runtime_20260816.json \
  docs/run_summaries/candidate_design/targeted_structure_review_result_20260816.json
do
  if [[ -e "${path}" ]]; then
    echo "Output already exists: ${path}" >&2
    exit 2
  fi
done

mkdir -p logs/targeted_structure_review
sbatch scripts/candidate_design/submit_targeted_structure_review.slurm
