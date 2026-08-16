#!/usr/bin/env bash
#SBATCH --job-name=nb252-double-net
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/double_mutant_scan/netsolp-%j.log
#SBATCH --error=logs/double_mutant_scan/netsolp-%j.log
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/netsolp
set -u
python scripts/candidate_design/score_double_mutant_properties.py --tool netsolp --plan-dir docs/result_artifacts/candidate_design/double_mutant_plan_20260816 --tool-root /homes/Tianlab/luly25/software/netsolp --output-dir results/candidate_design/double_mutant_scan_20260816/netsolp
