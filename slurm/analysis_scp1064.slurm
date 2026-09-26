#!/bin/bash
#SBATCH --job-name=analysis_scp1064
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=16:00:00
#SBATCH --output=../results/analysis/SCP1064/logs/slurm_%j.out
#SBATCH --error=../results/analysis/SCP1064/logs/slurm_%j.err
# Scientific analysis run (docs/HPC.md); output separate from the regression/demo directories.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)  commit: $(git rev-parse HEAD)"
petrubseq-protein --version
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
petrubseq-protein validate --config config/analysis_scp1064.yaml
/usr/bin/time -v petrubseq-protein run --config config/analysis_scp1064.yaml
python scripts/validate_processed.py ../results/analysis/SCP1064/processed/scp1064_processed.h5ad
echo "end: $(date)"
