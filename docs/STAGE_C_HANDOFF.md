# Stage C handoff (generic input adapters)

> **Status (2026-09-25, end of the Stage C session): Stage C is implemented.**
> `io/adapters.py` (canonical input), `io/tenx.py` (MTX, 10x H5),
> `io/h5ad.py`, `validation.py`, `preprocessing/guides.py` (guide calling),
> `config/examples/*.yaml`, `tests/test_cross_format.py` (25 tests) exist;
> docs/DEFAULTS.md ("Inputs", "Guide assignment"), docs/PROCESSED_OBJECT.md,
> docs/OUTPUTS.md and README.md describe the result. Sections 7-11 below are
> the plan that was implemented and are kept for context; the decisions taken
> on the open questions (assignment `source` semantics, single `multiplexed`
> block, protein availability rule, `<barcode>__<lane>` IDs, additive slots in
> `compare_processed.py`) are recorded in docs/DEFAULTS.md.

Written 2026-09-25 at the end of Stage B, for a fresh session. **The repository
on disk is authoritative**: inspect the existing implementation before
modifying anything; this document points at it rather than duplicating it.

Read first, in this order: `README.md`, this file, `docs/DEFAULTS.md`,
`docs/PROCESSED_OBJECT.md`, `docs/OUTPUTS.md`, `config/scp1064.yaml`,
`config/scp1064_smoke.yaml`.

**Context/token rule for the next session**: do not `cat` large files; use
`grep` and targeted reads with offsets/limits; never print matrices, h5ad
contents, CSVs, logs, figure lists or full test output; redirect long command
output to a file and inspect only head/tail/error lines; keep every command's
output short.

## 1. Filesystem organization (after Stage A)

```
<project>/      (one mount, two mount-point names)
├── petrubseq-pipeline-protein/   Git repository root (this repo)
├── data/SCP1064/                 SCP download in its original layout (cluster/ documentation/ expression/
│                                 metadata/ other/ file_supplemental_info.tsv) + smoke/ (4,500-cell subset)
├── results/                      generated outputs only (see section 5)
└── SCP1064/                      EMPTY leftover directory (old repo/data root); safe to rmdir
```

Repository layout: `src/petrubseq_protein/{cli,config,pipeline,audit}.py`,
subpackages `io/` (readers.py), `preprocessing/` (align, normalize,
embeddings), `qc/` (metrics, filtering, perturbation), `reporting/` (figures,
plots, html + templates/report.html, markdown, context, provenance, archive);
`config/`, `docs/`, `scripts/`, `slurm/`, `tests/`. Conda env
`petrubseq-protein` (Python 3.11.16, scanpy 1.11.5, anndata 0.12.19, jinja2
3.1), editable install points at this repo.

## 2. Stage A state (approved)

* Software / data / results are separated as above; nothing under the repo is
  data or output. `.gitignore` is a whitelist (py, yaml, md, html, slurm, toml).
* Git: root `petrubseq-pipeline-protein/`, branch `master` with **zero
  commits**, all 41 files untracked; `origin` = github vipinmenon1989/
  petrubseq-pipeline-protein, `origin/main` holds only an "Initial commit"
  README. Never commit/push/tag without explicit instruction.
* Path placeholders (implemented in `config.py`: `path_roots`, `expand_roots`):
  `${DATA_ROOT}` = `$PETRUBSEQ_DATA_ROOT` or `<repo>/../data`;
  `${RESULTS_ROOT}` = `$PETRUBSEQ_RESULTS_ROOT` or `<repo>/../results`;
  `${REPO_ROOT}`. Both SCP1064 YAMLs use them; no absolute HPC path in the repo.
* Regression after the move: 28 tests, smoke object identical to v0.1 (66/66
  slots), tables byte-identical.

## 3. Stage B state (implemented, validated, awaiting the user's review)

* **Lifecycle**: 23 stages in `pipeline.py::run_pipeline` (`_Stages` context
  manager gives `=== Stage k/23 ===` banners and timings; one `warnings`
  list; `notes`). Order: validate config, output tree, logging/resolved
  config, input audit, load, align, harmonize perturbations, RNA state/count
  layer (+ protein representations), permissive prefilter, pre-filter QC,
  BEFORE figures, strict filtering, audit table, post-filter summaries +
  perturbation QC, AFTER figures, representations (raw-count normalization
  is applied here, after filtering: `normalize.log_normalize_rna`),
  cross-modality diagnostics, h5ad, tables, figure manifest, report, run
  manifest, optional archive. Listed in `docs/OUTPUTS.md`.
* **QC semantics** (`qc/filtering.py`, `docs/DEFAULTS.md` "QC" section):
  prefilter (`qc.prefilter`, default 200 genes / 3 cells per gene, steps
  `prefilter_*`), strict filter (`qc.filter`, default on: rna min_genes 500,
  max_pct_mt 20, protein thresholds null, `perturbation.cells: all |
  assigned | single_guide`), flags (`qc.flags`, never remove). `FilterAudit`
  records every step with consistent arithmetic (`check()`); `strict_filter`
  returns `removed_by` (first removing step). Zero cells left -> ValueError.
* **Before/after architecture**: `reporting/plots.py::rna_qc_figures` and
  `protein_qc_figures` are called twice with `stage="before_filtering"` /
  `"after_filtering"` (`pipeline._qc_figures`); perturbation QC and
  representations are `post_filter` / `representation` stages.
* **Figure registry**: `reporting/figures.py::FigureRegistry.save(fig, name,
  section, stage, title, caption)` -> `figures/<section>[/<stage>]/<name>.png`;
  `manifest()` -> `tables/figure_manifest.csv`; `by_section(section, stage)`
  feeds the report. Sections: qc_rna, qc_protein, qc_perturbation,
  representations_rna, representations_protein, multimodal.
* **Reports**: one `reporting/context.py::ReportContext` is built in
  `pipeline.py` (stage 21) and rendered by `reporting/html.py::write_html`
  (jinja2 template `reporting/templates/report.html`, figures embedded as
  base64 when `report.embed_figures`) and `reporting/markdown.py::write_markdown`
  (`report.write_markdown`). Sections: summary cards, warnings, 1 inputs and
  modality audit, 2 cell/RNA QC, 3 protein QC, 4 perturbation QC,
  5 representations, 6 outputs and provenance.
* **Output tree / tables / figures**: `docs/OUTPUTS.md` (tables are CSV,
  `cell_qc*.csv.gz` gzipped; `qc_filtering_steps.csv`, `cell_qc_prefilter.csv.gz`
  with `qc_retained` and `removed_by`; optional `processed/<stem>_prefilter.h5ad`
  via `output.write_prefilter_h5ad`; optional `<run>_results.tar.gz` via
  `output.archive`, `reporting/archive.py`).
* **Processed-object additions** (`docs/PROCESSED_OBJECT.md`):
  `uns['petrubseq_protein']` now also holds `schema_version` ("0.2"),
  `run_name`, `inputs` (modality audit DataFrame `modalities`, `files`,
  `protein_embedding_representation`, `isotype_map`), `warnings`,
  `qc.prefilter`, `qc.filter` (cells before/after/removed),
  `qc.filtering_steps` (DataFrame), `perturbations.guide_counts_available`
  (False for SCP1064). Note: lists of dicts cannot be written to h5ad by
  anndata; store tables as DataFrames.
* **Config schema** (`config.py`; every key documented there, unknown keys
  fail): `run.name`, `dataset`, `inputs.{rna,protein,protein_counts,metadata,
  guide_assignments,embedding}`, `columns`, `perturbation`, `alignment`,
  `rna`, `protein`, `neighbors`, `umap`, `multimodal`,
  `qc.{prefilter,filter,flags}`, `report`, `compute`, `output`.
  `petrubseq-protein init-config x.yaml` writes the defaults. SCP1064 configs
  keep `qc.prefilter.enabled: false` and `qc.filter.enabled: false`
  (documented regression settings; the prefilter would drop 1 gene on the
  full data and 5,402 genes on the smoke subset).

## 4. Validation state at the end of Stage B

| check | result |
|---|---|
| `pytest` | 44 passed (28 v0.1 + 16 Stage B in `tests/test_stage_b.py`) |
| SCP1064 smoke (`results/SCP1064_smoke_stageB`) | 56 s, 1.0 GB; validator 54 checks 0 FAIL / 1 known WARN (protein PCA variance 78 %) |
| full SCP1064 (`results/SCP1064_stageB`, SLURM job 20840649) | 25 min 45 s, 24.0 GB peak, exit 0, 7 expected warnings; validator 56 checks 0 FAIL / 0 WARN |
| `compare_processed.py` vs v0.1 (smoke and full) | 66 slots: 65 identical/close, 0 DIFF, 1 NOTE (uns only) |
| tables vs v0.1 TSV | all renamed tables identical by content |

Whitelisted (expected) `uns['petrubseq_protein']` differences vs v0.1:
`config/*` (schema keys), `inputs/*`, `provenance/*` (paths, timestamps,
version, git root), `qc/filter/{cells_removed,cells_flagged}`,
`qc/filtering_steps`, `qc/prefilter/*`, `qc/protein/thresholds`,
`perturbations/guide_counts_available`, `run_name`, `schema_version`,
`version`, `warnings`, `schema`. Anything else differing in X,
reconstructed_counts, protein*, obs labels, QC columns, PCA, UMAP is a
regression.

## 5. Baselines and current outputs (under `../results/`)

> **Superseded on 2026-09-25 (Task A after Stage C).** The tree was normalized:
> `results/baselines/v0.1/{SCP1064,SCP1064_smoke}` (frozen v0.1, read-only),
> `results/{SCP1064,SCP1064_smoke}` (current validated Stage C outputs),
> `results/archive/` (superseded Stage A/B outputs). See `docs/REGRESSION.md`.
> The table below describes the historical layout at the end of Stage B.

| directory | role |
|---|---|
| `SCP1064` | frozen v0.1 full run (job 20838373): **regression baseline** |
| `SCP1064_smoke` | frozen v0.1 smoke run: **smoke regression baseline** |
| `SCP1064_smoke_stageA` | post-migration check (identical to baseline) |
| `SCP1064_smoke_stageB` | Stage B smoke output (report.html etc.) |
| `SCP1064_stageB` | Stage B full run output |

The repo's `slurm/run_scp1064.slurm` still writes to `../results/SCP1064`
(the baseline). Do not overwrite the baseline; use `--outdir` (the Stage B
run used a scratch copy of the script with `--outdir ../results/SCP1064_stageB`).

## 6. Commands (run from the repo root, `conda activate petrubseq-protein`)

```bash
pytest -q 2>&1 | tail -3
petrubseq-protein run --config config/scp1064_smoke.yaml --outdir ../results/SCP1064_smoke_stageC 2>&1 | tail -3
python scripts/validate_processed.py ../results/SCP1064_smoke_stageC/processed/scp1064_smoke_processed.h5ad \
       --expect-cells 4500 --expect-genes 23712 --expect-proteins 20 | tail -3
python scripts/compare_processed.py ../results/SCP1064_smoke/processed/scp1064_smoke_processed.h5ad \
       ../results/SCP1064_smoke_stageC/processed/scp1064_smoke_processed.h5ad | tail -5
# full run: needs a compute node (login node: 4 CPUs / 15 GB). sbatch requires --account=ihc on this cluster.
sed 's#../results/SCP1064/logs#../results/SCP1064_stageC/logs#g; s#--config config/scp1064.yaml#--config config/scp1064.yaml --outdir ../results/SCP1064_stageC#' \
    slurm/run_scp1064.slurm > /tmp/run_scp1064_stageC.slurm && mkdir -p ../results/SCP1064_stageC/logs && sbatch /tmp/run_scp1064_stageC.slurm
# validate/compare the 3.4 GB full object on a compute node:
srun -A ihc -p all --mem=80G -c 4 --time=60:00 python scripts/validate_processed.py ../results/SCP1064_stageC/processed/scp1064_processed.h5ad --expect-cells 218331 | tail -2
srun -A ihc -p all --mem=80G -c 4 --time=60:00 python scripts/compare_processed.py ../results/SCP1064/processed/scp1064_processed.h5ad ../results/SCP1064_stageC/processed/scp1064_processed.h5ad | tail -3
```
`compare_processed.py --allow-uns PREFIX ...` extends the whitelist.

## 7. Modules relevant to Stage C (what each does now)

| module | role now | Stage C relevance |
|---|---|---|
| `io/readers.py` (360 l.) | `Matrix` (cells, features, CSR, sources), streamed dense CSV/TSV reader (row chunks, multi-file, parallel), SCP tables (TYPE row), guide-list reader, embedding reader, `detect_value_state`, `discover_files` | becomes the dense_csv adapter; `Matrix` gains feature metadata |
| `config.py` (636 l.) | dataclass schema + validation + placeholders | `inputs.*` generalization, `perturbation.assignment`, `protein.feature_table` |
| `pipeline.py` (597 l.) | 23-stage lifecycle; stage 5 "load data" calls `pio.read_dense_matrix` per modality directly | replace stage 5 with `adapters.load_inputs(cfg)`; store guide counts |
| `preprocessing/align.py` | `align_cells` (exact ID intersection, duplicates abort, overlap threshold), `harmonize_perturbations` (guide lists -> classes), `build_obs` | consume guide lists derived from counts |
| `preprocessing/normalize.py` | RNA state detection / reconstruction / deferred normalization; protein primary + extra representations, isotype inference (`verify_isotype_ratio`), `protein_features` table | merge structured antibody metadata |
| `qc/metrics.py`, `qc/perturbation.py` | cell QC + flags; guide/target/condition coverage tables and `perturbation_qc` | add guide-count diagnostics columns (`guide_total_counts`, `n_guides_detected`, `guide_top_count`, `guide_second_count`) |
| `reporting/plots.py` | figures; `perturbation_qc_figures` already draws guide-count panels when `guide_total_counts` exists | mostly unchanged |
| `audit.py` | `audit --data/--config` inventory and role-aware audit (dense CSV only) | make format-aware |
| `scripts/validate_processed.py`, `scripts/compare_processed.py` | contract validator (56 checks), regression gate | add guide_counts / protein feature checks |
| `tests/make_synthetic.py`, `tests/test_pipeline.py`, `tests/test_stage_b.py` | synthetic fixture (dense CSV serialization only), 44 tests | fixture must also emit guide counts and MTX / h5ad serializations |

Reference implementation to study (read-only, clone into the scratchpad if
needed): weili-lab/perturbseq-pipeline `src/perturbseq_pipeline/io.py`
(`_load_mtx`, `_read_guide_mtx`, `_load_h5ad`, `split_features`,
`read_guide_table`, lane concatenation) and `guides.py` (`top_two_guides`,
dominant-guide rule `min_umi 3`, `dominance_ratio 2`, `max_second_umi -1`,
`ambiguous`/`unassigned`).

## 8. Stage C goals (approved plan)

* One canonical input model (`CanonicalInput`: `rna`, `protein_counts`,
  `protein_normalized`, `guide_counts` as `Matrix` with `features_meta`,
  `guide_assignments`, `metadata`, `embedding`, `provenance`) produced by an
  adapter layer (`io/adapters.py` or similar) so nothing after stage 5 knows
  the source format.
* Formats: dense CSV/TSV (existing), combined 10x MTX split by feature type
  (Gene Expression / Antibody Capture / CRISPR Guide Capture, configurable
  type names), separate RNA/ADT/GUIDE MTX directories, H5AD with
  slot-addressed modalities (X/layer/obsm/var feature type/obs column),
  10x H5 if it fits the same abstraction, multi-lane inputs with lane
  metadata joined on `lane_id`.
* Guide-count ingestion (sparse `obsm['guide_counts']` + `uns['guide_features']`)
  and guide calling (`perturbation.assignment`: dominant rule from the
  reference and a threshold rule producing multi-guide lists) that feeds
  the existing `harmonize_perturbations`; `ambiguous` recorded explicitly.
* Structured protein feature metadata in `uns['protein_features']`
  (`feature_id`, `antibody_name`, `protein_name`, `gene_symbol`, `clone`,
  `feature_type`, `isotype`, `isotype_control`, `annotation_source`; NA when
  unavailable, never invented); isotypes still excluded from embeddings.
* Centralized input validation (`validation.py`): duplicate barcodes,
  incompatible cell sets, duplicate feature IDs, missing metadata, negative or
  non-integer "counts", unknown RNA state when normalization requested,
  protein/guide name mismatches, malformed guide-target maps, impossible
  modality declarations, missing h5ad slots. Errors when the output would be
  ambiguous, warnings otherwise.
* Input provenance in `uns['petrubseq_protein']['inputs']` (format, files,
  splits, guide source, `guide_counts_available`).
* Synthetic cross-format equivalence test: one experiment serialized as
  dense CSV, combined MTX, separate MTX, h5ad (two layouts), run through the
  pipeline, compared with the same tolerances as `compare_processed.py`.
* Preserve Stage B reporting and the SCP1064 numerical behaviour (dense_csv
  path unchanged; SCP1064 configs remain valid without edits beyond additive keys).

## 9. Stage C non-goals

No Papalexi/GSE153056 demo yet (Stage D), no PS score, no lochNESS, no
clustering, no WNN, no DSB, no differential perturbation biology, no
RNA->protein ML, no caching layer, no new dependencies unless clearly needed
(scanpy/h5py/scipy.io already cover MTX, 10x H5 and h5ad).

## 10. Known Stage-B limitations

* Guide-count figures/columns are dormant until guide-count inputs exist.
* `protein_antibody_summary` uses a symlog axis, so antibodies with median 0
  show no bar (isotypes, background-dominated); cosmetic.
* With `qc.filter.enabled: false` the after-filter figures duplicate the
  before-filter set (by design; a report note says so).
* `output.archive` is tested on synthetic data only, not at full scale.
* Full-run `report.html` embeds 48 PNGs (5 MB); fine for email.

## 11. Files most likely touched in Stage C

`src/petrubseq_protein/io/{readers,adapters,mtx,h5ad}.py` (new adapters),
`validation.py` (new), `preprocessing/guides.py` (new), `config.py`,
`pipeline.py` (stage 5 and guide-count storage), `preprocessing/{align,
normalize}.py`, `qc/{metrics,perturbation}.py`, `audit.py`,
`scripts/{validate_processed,compare_processed}.py`, `tests/make_synthetic.py`,
`tests/test_cross_format.py` (new), `docs/{PROCESSED_OBJECT,DEFAULTS,OUTPUTS}.md`,
`README.md`, `config/examples/*.yaml` (new), `environment.yml`/`pyproject.toml`
only if a dependency is truly needed.

## 12. Working rules carried over

The user approves one stage at a time; STOP and report after Stage C with
the same regression protocol (pytest, smoke, validator, compare, full SLURM
run). Do not commit, push, tag or rewrite history. Do not modify
`../data/`, the frozen baselines, or the v0.1 numerical behaviour. Keep
command output small.
