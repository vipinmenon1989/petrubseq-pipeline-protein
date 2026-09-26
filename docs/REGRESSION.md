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

## Reference-parity corrections (2026-09-26) and the v0.1 baseline

The corrections that make the RNA analyses reproduce `weili-lab/perturbseq-pipeline`
(docs/reference/PARITY_RESULTS.md) change three things that the v0.1 baseline comparison
reports as DIFF or ADD and that are **intentional**:

* `obsm['X_pca']` / `obsm['X_umap_rna']`: the RNA PCA scales the sparse HVG block on the
  reference code path (scanpy densifies it in float64 before `sc.tl.pca`) instead of a
  float32 dense copy. PCs move by ~1e-3 (v0.1 smoke: max|diff| 7e-4), which is exactly the
  difference that made the Leiden partition diverge from the reference; the reference PCs
  are now reproduced to 1e-8. The UMAP moves accordingly.
* `obsm['X_pca_protein']` / `obsm['X_umap_protein']`: the protein PCA now follows the same
  numerical principle. The CLR block (float32 in `obsm['protein']`, unchanged) is cast to an
  explicit float64 working matrix, scaled and decomposed by arpack in float64, and the
  coordinates / loadings are stored as float32 (a deterministic cast; the h5ad size is
  unchanged). Against v0.1 the protein PCs differ at the ~1e-5 level (float32 arpack in
  v0.1); repeat runs of the new implementation are bit-identical, including the protein
  UMAP. Protein counts, the CLR formula and axis, isotype handling, feature exclusion and
  every protein-effect / concordance statistic are untouched.
* `obs['pct_counts_hb']`, `obs['total_counts_hb']`, `var['hb']`: the reference
  haemoglobin gene class (whitelisted as additive slots in `compare_processed.py`).

Everything else in the regression objects is unchanged (`X`, counts, QC metrics, protein
CLR, guide slots). The v0.1 baselines stay frozen; a new baseline version
(`results/baselines/v0.2/`) can be frozen from the next validated SCP1064 run.

## Stage record

| stage | content | validation |
|---|---|---|
| A | software / data / results separation | smoke identical to v0.1 |
| B | lifecycle, QC filtering, reports | smoke + full: 0 DIFF vs v0.1 |
| C | input adapters, guide calling | smoke + full: 0 DIFF vs v0.1; cross-format equivalence |
| D | Papalexi demo, bundled | deterministic preparation; validator 0 FAIL |
| E | PS, lochNESS, programs / modules, protein effects, concordance (optional) | reference agreement; smoke 0 DIFF with analyses off |
| F | Leiden cell states, perturbation × cluster enrichment (optional) | smoke 0 DIFF with clustering off; see the Stage F completion report |
