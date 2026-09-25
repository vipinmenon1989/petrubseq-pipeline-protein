#!/bin/bash
#SBATCH --job-name=papalexi_prepare
# Cluster policy: docs/HPC.md
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=../results/demo/papalexi_eccite/logs/prepare_%j.out
#SBATCH --error=../results/demo/papalexi_eccite/logs/prepare_%j.err
#
# Builds the deterministic Papalexi ECCITE-seq demo subset twice (into the
# canonical location and into a scratch copy) and proves they are identical.
# Resource basis: the streamed read of the 64 MB gzip / 1.5 GB dense cDNA
# matrix keeps ~2 GB resident; 16 GB is generous.
#   mkdir -p ../results/demo/papalexi_eccite/logs && sbatch slurm/prepare_papalexi_demo.slurm
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
DATA_ROOT="${PETRUBSEQ_DATA_ROOT:-$(pwd)/../data}"
SCRATCH="${PETRUBSEQ_SCRATCH:-$(pwd)/../results/demo/papalexi_eccite/prepare_check}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)"
/usr/bin/time -v python demo/prepare_papalexi_demo.py --source "$DATA_ROOT/ECCITE-seq" --output "$DATA_ROOT/demo/papalexi_eccite"
python demo/prepare_papalexi_demo.py --source "$DATA_ROOT/ECCITE-seq" --output "$SCRATCH/run2"
python demo/prepare_papalexi_demo.py --source "$DATA_ROOT/ECCITE-seq" --output "$SCRATCH/run3"
python demo/check_demo_reproducibility.py "$SCRATCH/run2" "$SCRATCH/run3"
python demo/check_demo_reproducibility.py "$DATA_ROOT/demo/papalexi_eccite" "$SCRATCH/run2"
echo "end: $(date)"
