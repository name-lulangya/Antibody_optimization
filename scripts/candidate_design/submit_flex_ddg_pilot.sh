#!/bin/bash

set -euo pipefail

PROJECT_ROOT="$(pwd -P)"
PLAN_DIR="docs/result_artifacts/candidate_design/flex_ddg_pilot_plan_20260812"

if [[ ! -f "${PLAN_DIR}/flex_ddg_pilot_plan.json" ]]; then
  echo "Run the local Flex ddG pilot plan builder before submission." >&2
  exit 1
fi

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/multi_ligand

python scripts/candidate_design/run_flex_ddg_task_pyrosetta.py \
  --plan-dir "${PLAN_DIR}" \
  --task-index 0 \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir results/candidate_design/flex_ddg_pilot_20260812/precheck-not-written \
  --check_only

cd "${PROJECT_ROOT}"
ARRAY_JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_flex_ddg_pilot_array.slurm)
SUMMARY_JOB_ID=$(sbatch --parsable \
  --dependency="afterok:${ARRAY_JOB_ID}" \
  scripts/candidate_design/submit_flex_ddg_pilot_summary.slurm)

echo "Flex ddG pilot array job: ${ARRAY_JOB_ID}"
echo "Flex ddG pilot summary job: ${SUMMARY_JOB_ID}"
