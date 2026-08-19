#!/usr/bin/env bash
#SBATCH --job-name=nb252-prop-v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=logs/expression_property_completion_v2/run-%j.log
#SBATCH --error=logs/expression_property_completion_v2/run-%j.log

set -eo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/expression_property_completion_plan_v2_20260819"
ANTIFOLD_PLAN="docs/result_artifacts/candidate_design/antifold_validation_plan_20260815"
SCORE_ROOT="results/candidate_design/expression_property_completion_v2_20260819"
VALIDATION_RESULT="docs/result_artifacts/candidate_design/expression_property_reuse_validation_v2_20260819"
FINAL_RESULT="docs/result_artifacts/candidate_design/expression_property_complete_matrix_v2_20260819"

source /data/software/anaconda3/etc/profile.d/conda.sh
export PYTHONUNBUFFERED=1

echo "[1/7] Repeat NetSolP on WT plus 12 validation candidates"
conda activate /data/software/env/luly25/netsolp
set -u
python scripts/candidate_design/score_expression_property_netsolp.py \
  --sample-table "${PLAN_DIR}/reuse_validation_samples.csv" \
  --fasta "${PLAN_DIR}/reuse_validation_sequences.fasta" \
  --contract "${PLAN_DIR}/expression_property_completion_contract.json" \
  --netsolp-workdir /homes/Tianlab/luly25/software/netsolp \
  --output-dir "${SCORE_ROOT}/reuse_validation/netsolp" \
  --num-threads 12

echo "[2/7] Repeat NanoMelt on WT plus 12 validation candidates"
set +u
conda activate /data/software/env/luly25/nanomelt
set -u
python scripts/candidate_design/score_expression_property_nanomelt.py \
  --sample-table "${PLAN_DIR}/reuse_validation_samples.csv" \
  --fasta "${PLAN_DIR}/reuse_validation_sequences.fasta" \
  --contract "${PLAN_DIR}/expression_property_completion_contract.json" \
  --nanomelt-executable /data/software/env/luly25/nanomelt/bin/nanomelt \
  --immune-builder-refine /data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py \
  --output-dir "${SCORE_ROOT}/reuse_validation/nanomelt"

echo "[3/7] Repeat AntiFold once on the three frozen WT structural views"
set +u
conda activate /data/software/env/luly25/antifold
set -u
python scripts/candidate_design/score_antifold_validation.py \
  --plan-dir "${ANTIFOLD_PLAN}" \
  --output-dir "${SCORE_ROOT}/reuse_validation/antifold"

echo "[4/7] Enforce the repeat-validation gate before reuse"
set +u
conda activate /data/software/env/luly25/ab_optim
set -u
python scripts/candidate_design/validate_expression_property_reuse.py \
  --plan-dir "${PLAN_DIR}" \
  --antifold-plan-dir "${ANTIFOLD_PLAN}" \
  --netsolp-score-dir "${SCORE_ROOT}/reuse_validation/netsolp" \
  --nanomelt-score-dir "${SCORE_ROOT}/reuse_validation/nanomelt" \
  --antifold-score-dir "${SCORE_ROOT}/reuse_validation/antifold" \
  --output-dir "${VALIDATION_RESULT}" \
  --run-summary docs/run_summaries/candidate_design/expression_property_reuse_validation_v2_20260819.json

echo "[5/7] Score WT plus 126 missing candidates with NetSolP"
set +u
conda activate /data/software/env/luly25/netsolp
set -u
python scripts/candidate_design/score_expression_property_netsolp.py \
  --sample-table "${PLAN_DIR}/missing_property_samples.csv" \
  --fasta "${PLAN_DIR}/missing_property_sequences.fasta" \
  --contract "${PLAN_DIR}/expression_property_completion_contract.json" \
  --netsolp-workdir /homes/Tianlab/luly25/software/netsolp \
  --output-dir "${SCORE_ROOT}/missing/netsolp" \
  --num-threads 12

echo "[6/7] Score WT plus 126 missing candidates with NanoMelt"
set +u
conda activate /data/software/env/luly25/nanomelt
set -u
python scripts/candidate_design/score_expression_property_nanomelt.py \
  --sample-table "${PLAN_DIR}/missing_property_samples.csv" \
  --fasta "${PLAN_DIR}/missing_property_sequences.fasta" \
  --contract "${PLAN_DIR}/expression_property_completion_contract.json" \
  --nanomelt-executable /data/software/env/luly25/nanomelt/bin/nanomelt \
  --immune-builder-refine /data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py \
  --output-dir "${SCORE_ROOT}/missing/nanomelt"

echo "[7/7] Merge the complete 847-row raw score matrix"
set +u
conda activate /data/software/env/luly25/ab_optim
set -u
python scripts/candidate_design/finalize_expression_property_matrix.py \
  --v2-contract-dir docs/result_artifacts/input_baseline/vhh_conservation_consensus_v2_20260819 \
  --plan-dir "${PLAN_DIR}" \
  --validation-dir "${VALIDATION_RESULT}" \
  --legacy-property-dir docs/result_artifacts/candidate_design/unified_property_scoring_result_20260815 \
  --legacy-antifold-dir docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815 \
  --netsolp-score-dir "${SCORE_ROOT}/missing/netsolp" \
  --nanomelt-score-dir "${SCORE_ROOT}/missing/nanomelt" \
  --output-dir "${FINAL_RESULT}" \
  --run-summary docs/run_summaries/candidate_design/expression_property_complete_matrix_v2_20260819.json

echo "Expression-property completion pipeline finished"
