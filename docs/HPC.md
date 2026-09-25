# HPC execution policy

**Use the login node for development and lightweight inspection; use SLURM
for substantial real-data computation.** Cluster settings live only in
`slurm/*.slurm` and this file; the Python package and CLI stay portable.

## Login node (allowed)

File names / sizes / headers, small metadata reads, grep / find, config
inspection, Git, tiny synthetic tests and `pytest`, the 4,500-cell SCP1064
smoke run (~80 s, ~1 GB), small code-development tasks.

Not allowed: full real-data matrix loading, ECCITE-seq / SCP1064 full
processing, expensive PCA / UMAP, complete production runs, large report
generation, high-memory conversions, ML / GPU work. If in doubt, use SLURM.

## SLURM resources

| workload | account | partition | node | extra |
|---|---|---|---|---|
| CPU real-data jobs | `ihc` | `ihc` | `ihc-grid-1-1-1` | `--cpus-per-task`, `--mem`, `--time` sized from the workload |
| GPU (future) H200 | `ihc` | `ihc` | `ihc-h200-1` | `--gres=gpu:1` |
| GPU (future) L40S | `ihc` | `ihc` | `ihc-l40s-1` | `--gres=gpu:1` |

Size CPU / RAM / walltime from measurements: the full SCP1064 run uses
~24 GB peak RSS and ~24 min on 8 CPUs (`slurm/run_scp1064.slurm` asks for
48 GB / 2 h). Validating or comparing the 3.4 GB processed object needs
~20 GB (`srun ... --mem=80G` in `docs/REGRESSION.md`). Stage D (Papalexi
demo) is CPU-only; request a GPU only for GPU work and verify it is visible
inside the job (`nvidia-smi`).

## Script rules (`slurm/`)

Every substantial job has a reproducible script that sets account / partition
/ node, activates `petrubseq-protein`, `cd`s to the repository root
(`SLURM_SUBMIT_DIR`), uses the production CLI, runs with `set -euo pipefail`,
writes stdout / stderr under the run's `logs/`, and never writes into `data/`.
Submit with `sbatch`, record the job ID, monitor with `squeue -u $USER` /
`sacct -j <id> --format=State,Elapsed,MaxRSS,NodeList`, inspect only log tails,
and record job ID, node, exit status, walltime, peak memory and CPU allocation
in the stage report.

Interactive one-off commands on the CPU node:

```bash
srun -A ihc -p ihc --nodelist=ihc-grid-1-1-1 --mem=32G -c 4 --time=60:00 <command>
```
