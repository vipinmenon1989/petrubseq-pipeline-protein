#!/bin/bash
#SBATCH --job-name=papalexi_full_prep
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=../results/analysis/papalexi_eccite_full/logs/prepare_%j.out
#SBATCH --error=../results/analysis/papalexi_eccite_full/logs/prepare_%j.err
# Full pooled Papalexi screen: every HTO singlet, no downsampling (n-cells = eligible singlets).
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
DATA_ROOT="${PETRUBSEQ_DATA_ROOT:-$(pwd)/../data}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)"
/usr/bin/time -v python demo/prepare_papalexi_demo.py --source "$DATA_ROOT/ECCITE-seq" --output "$DATA_ROOT/analysis/papalexi_eccite_full" --n-cells 17473 --min-per-stratum 1 --seed 0
echo "end: $(date)"
