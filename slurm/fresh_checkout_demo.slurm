#!/bin/bash
#SBATCH --job-name=fresh_checkout_demo
# Cluster policy: docs/HPC.md
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Fresh-checkout test of the bundled Papalexi demo: runs the production CLI
# from a clean copy of the repository (FRESH_REPO) whose parent directory
# holds no data/ or results/, using only the files inside the checkout.
# The package is imported from the checkout's src/ (PYTHONPATH), never from
# the development install, and nothing is symlinked or copied in.
#   FRESH_REPO=/path/to/clean/petrubseq-pipeline-protein sbatch slurm/fresh_checkout_demo.slurm
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate petrubseq-protein
: "${FRESH_REPO:?set FRESH_REPO to the clean checkout}"
cd "$FRESH_REPO"
export PYTHONPATH="$FRESH_REPO/src"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  checkout: $FRESH_REPO  start: $(date)"
test ! -e ../data && echo "no ../data next to the checkout: OK"
test ! -e ../results && echo "no ../results next to the checkout: OK"
python -c "import petrubseq_protein, sys; print('package from', petrubseq_protein.__file__); assert petrubseq_protein.__file__.startswith('$FRESH_REPO')"
python -m petrubseq_protein.cli validate --config config/demo_papalexi.yaml
/usr/bin/time -v python -m petrubseq_protein.cli run --config config/demo_papalexi.yaml
python scripts/validate_processed.py ../results/demo/papalexi_eccite/processed/papalexi_eccite_demo_processed.h5ad
test -s ../results/demo/papalexi_eccite/report.html && echo "report.html created"
echo "end: $(date)"
