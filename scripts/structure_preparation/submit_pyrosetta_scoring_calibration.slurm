#!/bin/bash
#SBATCH --job-name=nb252_score_cal
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=pyrosetta-score-calibration-%j.log
#SBATCH --error=pyrosetta-score-calibration-%j.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/multi_ligand

python scripts/structure_preparation/calibrate_pyrosetta_scoring.py \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --import-gate-dir docs/result_artifacts/structure_preparation/pyrosetta_wt_import_20260810 \
  --experimental-structure data/structures/cxs_exports/NK2R-252__native.cif \
  --output-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_20260810 \
  --run-summary docs/run_summaries/structure_preparation/pyrosetta_scoring_calibration_20260810.json \
  --replicates 8
