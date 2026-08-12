#!/usr/bin/env bash
#SBATCH --job-name=nb252-flex-prod
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/flex_ddg_production/production-%A_%a.log
#SBATCH --error=logs/flex_ddg_production/production-%A_%a.log

set -euo pipefail

if [[ -z "${FLEX_DDG_INDEX_FILE:-}" || ! -f "${FLEX_DDG_INDEX_FILE}" ]]; then
  echo "FLEX_DDG_INDEX_FILE is missing or invalid" >&2
  exit 2
fi
GLOBAL_TASK_INDEX=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${FLEX_DDG_INDEX_FILE}")
if [[ ! "${GLOBAL_TASK_INDEX}" =~ ^[0-9]+$ ]]; then
  echo "Cannot map local array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
TASK_ID=$(printf 'task_%04d' "${GLOBAL_TASK_INDEX}")
OUTPUT_DIR="results/candidate_design/flex_ddg_production_20260812/tasks/${TASK_ID}"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/multi_ligand

python scripts/candidate_design/run_flex_ddg_task_pyrosetta.py \
  --run-kind production \
  --plan-dir docs/result_artifacts/candidate_design/flex_ddg_production_plan_20260812 \
  --task-index "${GLOBAL_TASK_INDEX}" \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir "${OUTPUT_DIR}" \
  --backrub-trials 35000 \
  --backrub-temperature 1.2 \
  --backrub-neighborhood-angstrom 8.0
