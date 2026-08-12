#!/bin/bash
#SBATCH --job-name=nb252_flex_pilot
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --array=0-7%8
#SBATCH --time=06:00:00
#SBATCH --output=nb252-flex-pilot-%A_%a.log
#SBATCH --error=nb252-flex-pilot-%A_%a.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/multi_ligand

printf -v TASK_ID "task_%02d" "${SLURM_ARRAY_TASK_ID}"
OUTPUT_DIR="results/candidate_design/flex_ddg_pilot_20260812/tasks/${TASK_ID}"
mkdir -p "${OUTPUT_DIR}"

python scripts/candidate_design/run_flex_ddg_task_pyrosetta.py \
  --plan-dir docs/result_artifacts/candidate_design/flex_ddg_pilot_plan_20260812 \
  --task-index "${SLURM_ARRAY_TASK_ID}" \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir "${OUTPUT_DIR}" \
  --backrub-trials 35000 \
  --backrub-temperature 1.2 \
  --backrub-neighborhood-angstrom 8.0
