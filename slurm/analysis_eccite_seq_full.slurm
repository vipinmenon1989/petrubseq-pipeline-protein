#!/bin/bash
#SBATCH --job-name=analysis_eccite_seq_full
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=../results/analysis/ECCITE-seq/logs/slurm_%j.out
#SBATCH --error=../results/analysis/ECCITE-seq/logs/slurm_%j.err
# Production analysis of the full pooled ECCITE-seq screen (docs/HPC.md); output separate from the regression/demo directories.
#   mkdir -p ../results/analysis/ECCITE-seq/logs && sbatch slurm/analysis_eccite_seq_full.slurm
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)  commit: $(git rev-parse HEAD)"
petrubseq-protein --version
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
petrubseq-protein validate --config config/analysis_eccite_seq_full.yaml
/usr/bin/time -v petrubseq-protein run --config config/analysis_eccite_seq_full.yaml
python scripts/validate_processed.py ../results/analysis/ECCITE-seq/processed/eccite_seq_full_processed.h5ad
echo "end: $(date)"
