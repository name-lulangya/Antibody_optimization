#!/usr/bin/env bash
#SBATCH --job-name=nb252-double-melt
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/double_mutant_scan/nanomelt-%j.log
#SBATCH --error=logs/double_mutant_scan/nanomelt-%j.log
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/nanomelt
set -u
python scripts/candidate_design/score_double_mutant_properties.py --tool nanomelt --plan-dir docs/result_artifacts/candidate_design/double_mutant_plan_20260816 --executable /data/software/env/luly25/nanomelt/bin/nanomelt --output-dir results/candidate_design/double_mutant_scan_20260816/nanomelt
