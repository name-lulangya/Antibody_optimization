#!/usr/bin/env bash
#SBATCH --job-name=nb252-double-rose
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=logs/double_mutant_scan/pyrosetta-%j.log
#SBATCH --error=logs/double_mutant_scan/pyrosetta-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source /data/software/anaconda3/etc/profile.d/conda.sh
set +u
conda activate /data/software/env/luly25/multi_ligand
set -u
python scripts/candidate_design/score_double_mutants_pyrosetta.py --plan-dir docs/result_artifacts/candidate_design/double_mutant_plan_20260816 --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 --output-dir results/candidate_design/double_mutant_scan_20260816/pyrosetta --run-summary docs/run_summaries/candidate_design/double_mutant_pyrosetta_20260816.json --replicates 3 --base-seed 8163000
