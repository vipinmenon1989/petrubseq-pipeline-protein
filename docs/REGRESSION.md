# Regression protocol and results layout

```
results/
├── baselines/v0.1/SCP1064/          # immutable v0.1 full reference (job 20838373, 2026-09-24), read-only
├── baselines/v0.1/SCP1064_smoke/    # immutable v0.1 smoke reference, read-only
├── SCP1064/                         # CURRENT validated full output  (config/scp1064.yaml, slurm/run_scp1064.slurm)
├── SCP1064_smoke/                   # CURRENT validated smoke output (config/scp1064_smoke.yaml)
├── demo/<name>/                     # demo runs (config/demo_*.yaml)
└── archive/                         # superseded stage outputs kept for reference (may be deleted)
```

* `results/baselines/` is never a pipeline destination: `Config.validate` rejects
  an `output.dir` under a `baselines` directory, and the directory is `chmod a-w`.
* Never recompute or overwrite v0.1. A new baseline version gets a new
  `results/baselines/<version>/` directory.
* `data/` is input only; no run writes there.

## Protocol (after any change to the pipeline)

```bash
conda activate petrubseq-protein
pytest -q 2>&1 | tail -3                                     # login node
petrubseq-protein run --config config/scp1064_smoke.yaml 2>&1 | tail -3    # ~80 s, ~1 GB: login node acceptable
python scripts/validate_processed.py ../results/SCP1064_smoke/processed/scp1064_smoke_processed.h5ad \
    --expect-cells 4500 --expect-genes 23712 --expect-proteins 20 | tail -2
python scripts/compare_processed.py ../results/baselines/v0.1/SCP1064_smoke/processed/scp1064_smoke_processed.h5ad \
    ../results/SCP1064_smoke/processed/scp1064_smoke_processed.h5ad | tail -3
# full run (SLURM, docs/HPC.md): ~24 min, ~24 GB peak on 8 CPUs
mkdir -p ../results/SCP1064/logs && sbatch slurm/run_scp1064.slurm
# validate / compare the 3.4 GB object on the CPU node
srun -A ihc -p ihc --nodelist=ihc-grid-1-1-1 --mem=80G -c 4 --time=60:00 \
    python scripts/validate_processed.py ../results/SCP1064/processed/scp1064_processed.h5ad --expect-cells 218331 | tail -2
srun -A ihc -p ihc --nodelist=ihc-grid-1-1-1 --mem=80G -c 4 --time=60:00 \
    python scripts/compare_processed.py ../results/baselines/v0.1/SCP1064/processed/scp1064_processed.h5ad \
    ../results/SCP1064/processed/scp1064_processed.h5ad | tail -3
```

`compare_processed.py` statuses: `SAME`, `ADD` (slot present only in the new
object and on the additive whitelist: guide-count slots, protein annotation
columns, lane columns), `NOTE` (whitelisted provenance keys), `DIFF`
(anything else; non-zero exit). The acceptance criterion is **0 DIFF**; every
ADD/NOTE must be explainable by the change being tested.

## Reference values (Stage C, 2026-09-25)

| check | smoke | full |
|---|---|---|
| validator | 60 checks, 0 FAIL, 1 WARN (protein PCA variance 78 %) | 60 checks, 0 FAIL, 0 WARN |
| compare vs v0.1 | 67 slots, 0 DIFF, 1 ADD, 1 NOTE | 67 slots, 0 DIFF, 1 ADD, 1 NOTE |
| run | 78 s, 6 warnings | job 20846621: 24 min, 23.5 GB peak, 7 warnings |
