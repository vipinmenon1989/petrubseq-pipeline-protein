#!/bin/bash
#SBATCH --job-name=scp1064_analysis_smoke
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=../results/analysis/SCP1064_smoke/logs/slurm_%j.out
#SBATCH --error=../results/analysis/SCP1064_smoke/logs/slurm_%j.err
# Smoke test of the complete analysis (incl. the perturbation-distance stages) on the 4,500-cell SCP1064 subset (docs/HPC.md).
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)  commit: $(git rev-parse HEAD)"
petrubseq-protein --version
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
petrubseq-protein validate --config config/analysis_scp1064_smoke.yaml
/usr/bin/time -v petrubseq-protein run --config config/analysis_scp1064_smoke.yaml
python scripts/validate_processed.py ../results/analysis/SCP1064_smoke/processed/scp1064_analysis_smoke_processed.h5ad
echo "end: $(date)"
