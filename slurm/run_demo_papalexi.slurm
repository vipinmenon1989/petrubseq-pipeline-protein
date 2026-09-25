#!/bin/bash
#SBATCH --job-name=papalexi_demo
# Cluster policy: docs/HPC.md
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=../results/demo/papalexi_eccite/logs/slurm_%j.out
#SBATCH --error=../results/demo/papalexi_eccite/logs/slurm_%j.err
#
# Papalexi ECCITE-seq demo through the production CLI (~1,800 cells; the
# 4,500-cell SCP1064 smoke needs ~1 GB and 80 s, so 16 GB / 1 h is ample).
#   mkdir -p ../results/demo/papalexi_eccite/logs && sbatch slurm/run_demo_papalexi.slurm
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)"
petrubseq-protein --version
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
petrubseq-protein validate --config config/demo_papalexi.yaml
/usr/bin/time -v petrubseq-protein run --config config/demo_papalexi.yaml
python scripts/validate_processed.py ../results/demo/papalexi_eccite/processed/papalexi_eccite_demo_processed.h5ad
echo "end: $(date)"
