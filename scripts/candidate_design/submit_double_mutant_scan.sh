#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
PLAN=docs/result_artifacts/candidate_design/double_mutant_plan_20260816
ROOT=results/candidate_design/double_mutant_scan_20260816
test -s "$PLAN/double_mutant_plan_gate.json"
if [[ -e "$ROOT" || -e docs/run_summaries/candidate_design/double_mutant_pyrosetta_20260816.json || -e docs/result_artifacts/candidate_design/double_mutant_scan_result_20260816 || -e docs/run_summaries/candidate_design/double_mutant_scan_result_20260816.json ]]; then echo "Double-mutant scan output already exists; refusing to overwrite." >&2;exit 1;fi
mkdir -p logs/double_mutant_scan
NET=$(sbatch --parsable scripts/candidate_design/submit_double_mutant_netsolp.slurm);MELT=$(sbatch --parsable scripts/candidate_design/submit_double_mutant_nanomelt.slurm);TNP=$(sbatch --parsable scripts/candidate_design/submit_double_mutant_tnp.slurm);ROSE=$(sbatch --parsable scripts/candidate_design/submit_double_mutant_pyrosetta.slurm)
NET_ID=${NET%%;*};MELT_ID=${MELT%%;*};TNP_ID=${TNP%%;*};ROSE_ID=${ROSE%%;*};JOIN=$(sbatch --parsable --dependency="afterok:${NET_ID}:${MELT_ID}:${TNP_ID}:${ROSE_ID}" scripts/candidate_design/submit_double_mutant_analysis.slurm)
echo "NetSolP job: ${NET_ID}";echo "NanoMelt job: ${MELT_ID}";echo "TNP job: ${TNP_ID}";echo "PyRosetta job: ${ROSE_ID}";echo "Joint analysis job: ${JOIN%%;*} (after all four scans pass)"
