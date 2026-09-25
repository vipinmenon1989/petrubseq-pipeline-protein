#!/bin/bash
#SBATCH --job-name=petrubseq_scp1064
# Cluster policy: docs/HPC.md (account ihc, partition ihc, CPU node ihc-grid-1-1-1).
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=../results/SCP1064/logs/slurm_%j.out
#SBATCH --error=../results/SCP1064/logs/slurm_%j.err
#
# Full SCP1064 processing run (218,331 cells x 23,712 genes x 20 proteins).
#
# Resource basis (docs/datasets/SCP1064_AUDIT.md section 9): the first full run (job
# 20838373) peaked at 23.9 GB RSS and took 28 min on 8 CPUs, so 48 GB / 2 h
# is ~2x / ~4x headroom. Six RNA chunk files are parsed in parallel
# (compute.n_jobs: 6), hence 8 CPUs.
#
# Submit from the repository root. Results go to ../results/SCP1064 (i.e.
# ${RESULTS_ROOT}/SCP1064, the CURRENT full output; override with
# PETRUBSEQ_RESULTS_ROOT and --output). The frozen regression reference lives in
# ../results/baselines/v0.1/SCP1064 and is never written (docs/REGRESSION.md).
#   mkdir -p ../results/SCP1064/logs
#   sbatch slurm/run_scp1064.slurm
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p ../results/SCP1064/logs

echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)"
echo "python: $(which python)  ($(python --version 2>&1))"
petrubseq-protein --version

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

petrubseq-protein validate --config config/scp1064.yaml
/usr/bin/time -v petrubseq-protein run --config config/scp1064.yaml

echo "end: $(date)"
