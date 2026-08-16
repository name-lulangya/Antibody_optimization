#!/usr/bin/env bash
#SBATCH --job-name=nb252-double-join
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/double_mutant_scan/analysis-%j.log
#SBATCH --error=logs/double_mutant_scan/analysis-%j.log
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim
set -u
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-double-analysis-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_double_mutant_scan.py \
  --plan-dir docs/result_artifacts/candidate_design/double_mutant_plan_20260816 \
  --score-root results/candidate_design/double_mutant_scan_20260816 \
  --calibration-gate docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811/pyrosetta_scoring_calibration_gate.json \
  --output-dir docs/result_artifacts/candidate_design/double_mutant_scan_review_v2_1_20260816 \
  --run-summary docs/run_summaries/candidate_design/double_mutant_scan_review_v2_1_20260816.json
