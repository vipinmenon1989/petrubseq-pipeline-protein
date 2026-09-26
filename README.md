# petrubseq-pipeline-protein

A reproducible Perturb-CITE-seq pipeline for RNA, CRISPR guide and
optional ADT/protein measurements: input harmonization, QC, guide assignment,
RNA/protein preprocessing, Leiden clustering with perturbation × cluster
enrichment, perturbation-response analysis, gene programs, protein effects and a
self-contained report, in one command.

```
input → QC → guide assignment → RNA (+ protein) preprocessing → PCA / neighbours / UMAP
      → Leiden clustering → perturbation strength → perturbation × cluster enrichment
      → gene programs / perturbation modules → PS → lochNESS
      → protein effects → RNA–protein concordance   (both only when protein exists)
      → report.html + processed .h5ad
```

```bash
petrubseq-protein run --config config/demo_papalexi.yaml
```

Every run produces:

1. a processed **`.h5ad`** (aligned modalities, QC, perturbation labels, embeddings, analysis results),
2. a self-contained **`report.html`** (every figure embedded; `report.md` mirrors it),
3. **CSV tables and figures** on disk for everything in the report, plus the resolved config and a run manifest.

Protein is optional: RNA + guides alone runs the full pipeline; the protein
analyses skip themselves when no protein is present.

---

## What it does

**1 · Input harmonization.** Dense CSV/TSV, 10x MTX, 10x H5 or `.h5ad`; one combined
feature matrix (Gene Expression / Antibody Capture / CRISPR Guide Capture) or one
matrix per modality; several lanes (`<barcode>__<lane>` IDs, per-lane sample sheet).
Everything becomes one canonical model and is validated before processing: duplicate
barcodes/features, negative or non-integer "counts", unclaimed feature types, missing
h5ad slots, metadata sharing no cell with the RNA matrix.

**2 · QC.** Reference Perturb-seq QC: prefilter (200 genes, genes in ≥ 3 cells),
QC metrics (mt, ribo, hb) on the raw counts, strict filter (≥ 1,000 genes, < 20 %
mito, optional counts / % hb, then genes in ≥ 3 retained cells) recorded step by step;
protein and perturbation criteria optional; flags that never remove cells. RNA,
protein (ADT depth, antibodies detected, isotype background) and guide QC (guide UMIs,
guides per cell, top-vs-second guide, cells per target / per guide) figures.

**3 · Guide assignment.** Provided per-cell guide lists / metadata columns, or calls
from guide counts: `dominant` (top guide ≥ 3 UMIs and > 2× the runner-up, else
`ambiguous` / `unassigned`) or `threshold` (every guide ≥ 3 UMIs, multi-guide).
Classes `single_targeting`, `single_control`, `multi_targeting`, `multi_control`,
`mixed_control_targeting`, `ambiguous`, `unassigned`. When both a provided list and
guide counts exist the run asks you to choose and reports their agreement.

**4 · RNA and protein preprocessing.** Raw counts normalized once (median library
size, log1p, as the reference); log-normalized input preserved and integer counts
reconstructed when a library size is given. 3,000 HVGs, scaled HVG block, 50 PCs,
15-neighbour graph, UMAP, Leiden (resolution 1.0, igraph, seeded), bit-identical to
the reference pipeline. Protein: raw ADT counts kept, per-protein CLR across cells,
isotypes used for QC and excluded from the protein PCA/UMAP; supplied normalized
values preserved, never treated as counts.

**5 · Perturbation strength.** For every target whose gene is measured, the gene's
*own* expression in perturbed vs control cells (non-targeting and other-target arms):
log2FC on de-logged means, Kolmogorov–Smirnov and Mann–Whitney tests, BH-FDR per
arm, effective knockdown = KS FDR < 0.05 and log2FC < 0; volcano, waterfall,
control comparison and one ECDF / violin / UMAP figure per target.

**6 · Perturbation × cluster enrichment.** Fisher exact test (or Cochran–Mantel–
Haenszel across a sample/lane column) of each target's single-guide cells per
Leiden cluster against both reference arms (`other` drives the calls), Haldane odds
ratios, BH within arm, guide concordance for significant pairs, omnibus permutation
test, composition / phenocopy / volcano / effect-magnitude figures and one per target.

**7 · Co-functional modules and gene programs.** Perturbation × gene log2FC matrix
over the union of Leiden cluster markers, vs non-targeting controls; genes clustered
into programs (`P1..`, 1 − Pearson), perturbations into modules (`M1..`, 1 − Spearman)
after Zhou et al. 2023; module × program strength, alluvial, TF hubs / edges /
module connectivity, program activity per cell (`score_genes`) and per cluster.

**8 · Per-cell perturbation response (PS)** (Song et al. 2025, PS_python definition
re-implemented): per-cell score in [0, 1], knockdown / escaper / non-responder /
low-signal quadrants against the target's own expression, agreement with the
group-level test, and the supervised LDA embedding with per-target maps.

**9 · lochNESS.** Cluster-free local enrichment of each perturbation among its 300
nearest neighbours in PCA space (Huang et al. 2023): self-enrichment ranking,
distributions, mean per cluster, self score on the UMAP and one map per
perturbation. Cluster enrichment and lochNESS are shown side by side, never combined.

**10 · Perturbation distance, distance space and phenotype modules** *(reference
stages 10–11; off by default as in the reference, enabled in the demo and analysis
configs)*. Energy distance between each target's cells and the controls in PCA space
with a seeded permutation DistanceTest and BH-FDR; the target × target distance matrix,
its PCoA phenotype space, nearest phenotypic neighbours and average-linkage **phenotype
modules** (cell-state similarity, compared with the gene-effect modules by ARI / NMI);
ranking, atlas, PS-vs-distance map and phenotype-space figures. The **master
perturbation table** (`tables/master_perturbation_table.csv`) consolidates every
target-level result without a composite score. docs/PERTURBATION_DISTANCE.md.

**11 · Protein extension** *(needs protein)*. Per target × protein: CLR mean
difference, Cohen's d, Mann–Whitney p, BH-FDR, sign consistency per sample and per
guide. Concordance: PS ↔ protein within target, lochNESS ↔ protein across targets,
program activity ↔ protein per cell and per target, distance ↔ protein across targets
(signed and absolute) and phenotype space ↔ protein (Mantel, module-wise), the
integrated table and the target × protein master table with every statistic labelled
CELL_LEVEL or TARGET_LEVEL. Associations only.

Stages 5–10 reproduce the reference `weili-lab/perturbseq-pipeline` numerically
(docs/reference/PARITY_RESULTS.md); stage 11 is the multimodal addition.

---

## Install

```bash
git clone https://github.com/vipinmenon1989/petrubseq-pipeline-protein.git
cd petrubseq-pipeline-protein
conda env create -f environment.yml      # creates 'petrubseq-protein' (Python 3.11)
conda activate petrubseq-protein
pip install -e .
petrubseq-protein --help
```

`pip install -e .` also works in any Python ≥ 3.10 environment (`pyproject.toml`
lists the runtime dependencies; `igraph` is needed for Leiden).

---

## Quick start: the demo

A deterministic 1,800-cell subset of the Papalexi et al. 2021 ECCITE-seq screen
(GEO GSE153056, pooled screen in IFN-γ stimulated THP-1 cells) is bundled in
`demo/data/papalexi_eccite/` (~20 MB): 18,649 genes, 4 ADTs (CD86, PD-L1, PD-L2,
TIM-3), 111 guides (25 targets, 9 non-targeting, 1 eGFP), 3 hashed samples, 8 lanes.

```bash
petrubseq-protein run --config config/demo_papalexi.yaml
python scripts/validate_processed.py ../results/demo/papalexi_eccite/processed/papalexi_eccite_demo_processed.h5ad
```

About a minute, under 1 GB. Output goes to `../results/demo/papalexi_eccite/`
(`PETRUBSEQ_RESULTS_ROOT` overrides); open `report.html`. The notebook
[`notebooks/demo_papalexi_eccite.ipynb`](notebooks/demo_papalexi_eccite.ipynb)
walks through the outputs.

On this subset: the IFN-γ receptor/JAK/STAT targets (JAK2, IFNGR1, IFNGR2, STAT1)
occupy one Leiden cluster without control cells, have the strongest lochNESS
scores, lower the IFN-response program P1 (CXCL9, CXCL10, SOCS1, WARS, HLA-DRA)
and reduce PD-L1 protein, consistently across samples and guides. It is a
demonstration subset (1,800 cells, 21–102 cells per target, four proteins, every
sample stimulated), not a reanalysis of the paper.

To rebuild the subset from GEO (optional): `demo/fetch_papalexi_data.py`,
`demo/prepare_papalexi_demo.py`, `demo/check_demo_reproducibility.py`.

---

## Inputs

| data | needed | config |
|---|---|---|
| RNA counts or log-normalized RNA | always | `inputs.rna` or `inputs.multiplexed` |
| guide counts, or a per-cell guide list / metadata column | for perturbation analyses | `inputs.guide_counts`, `inputs.guide_assignments`, `columns.guide` |
| ADT raw counts and/or normalized protein | optional | `inputs.protein_counts`, `inputs.protein` |
| per-cell metadata (condition, sample, library size …) | required by default (`inputs.metadata.required: false` to drop) | `inputs.metadata`, `columns.*` |
| per-lane sample sheet, antibody table, precomputed embedding | optional | `inputs.lane_metadata`, `protein.feature_table`, `inputs.embedding` |

Formats: `dense_csv`, `mtx`, `10x_h5`, `h5ad`. One config per layout in
`config/examples/`. Start from the documented defaults:

```bash
petrubseq-protein init-config my_run.yaml
petrubseq-protein validate --config my_run.yaml     # config + every input exists
petrubseq-protein audit --config my_run.yaml        # cell-ID overlap, value states, protein panel
```

### Option 1 — 10x / matrix input

One combined feature-barcode matrix, read once and split by feature type
(`config/examples/mtx_combined.yaml`; separate MTX/H5 directories per modality:
`mtx_separate.yaml`; several lanes: `multilane.yaml`):

```yaml
dataset:
  name: my_screen
  input_dir: ${DATA_ROOT}/my_screen
inputs:
  multiplexed:
    format: mtx
    path: filtered_feature_bc_matrix
    feature_types: {rna: [Gene Expression], protein: [Antibody Capture], guide: [CRISPR Guide Capture]}
    var_names: name
  metadata: {file: cells.csv, format: auto}
columns: {cell_id: barcode, condition: condition}
perturbation:
  guide_target_regex: "^(?P<target>.+?)_(?P<index>\\d+)$"
  control_classes: {non_targeting: ["^NTC$"]}
analysis:
  clustering: {enabled: true}
  perturbation_effects: {enabled: true}
output:
  dir: ${RESULTS_ROOT}/my_screen
```

### Option 2 — existing `.h5ad`

Each modality is read from a slot of the same AnnData; cell metadata and an
upstream guide call come from `obs`. Perturb-seq, RNA + guide counts, no protein
(`config/examples/h5ad_perturbseq.yaml`):

```yaml
inputs:
  rna:          {format: h5ad, file: screen.h5ad, slot: layers, key: counts, state: raw_counts}
  guide_counts: {format: h5ad, file: screen.h5ad, slot: obsm, key: guide_counts}
  protein:      {required: false}
  metadata:     {file: screen.h5ad, format: h5ad}        # obs is the cell table
columns: {cell_id: obs_names, condition: condition, sample: sample}
alignment: {required: [rna, metadata]}
perturbation:
  guide_target_regex: "^(?P<target>.+?)_(?P<index>\\d+)$"
  control_classes: {non_targeting: ["^NTC$"]}
```

Perturb-CITE-seq, RNA + guide counts + ADT counts, with the upstream guide call
in `obs['guide_call']` driving the assignment
(`config/examples/h5ad_perturb_cite_seq.yaml`):

```yaml
inputs:
  rna:               {format: h5ad, file: screen.h5ad, slot: X, state: normalized}
  protein_counts:    {format: h5ad, file: screen.h5ad, slot: obsm, key: protein_counts}
  guide_counts:      {format: h5ad, file: screen.h5ad, slot: obsm, key: guide_counts}
  guide_assignments: {file: screen.h5ad, format: h5ad, guides_column: guide_call, list_separator: ";"}
  metadata:          {file: screen.h5ad, format: h5ad}
columns: {cell_id: obs_names, condition: condition, sample: sample}
perturbation:
  assignment: {source: provided}    # counts are still loaded, compared and stored
```

| purpose | AnnData location | config |
|---|---|---|
| RNA expression | `X` or `layers[<key>]` | `inputs.rna: {slot: X}` or `{slot: layers, key: <key>}` |
| RNA value state | – | `inputs.rna.state: raw_counts \| normalized \| auto` (auto = detected from the values) |
| raw ADT counts / normalized protein | `obsm[<key>]` | `inputs.protein_counts` / `inputs.protein` with `{slot: obsm, key: <key>}` |
| guide UMI counts | `obsm[<key>]` | `inputs.guide_counts: {slot: obsm, key: <key>}` |
| feature names of an `obsm` matrix | DataFrame columns, or `uns['<key>_features']` for a plain array | automatic |
| cell metadata | `obs` | `inputs.metadata: {file: <h5ad>, format: h5ad}`; `columns.cell_id: obs_names` or an obs column |
| upstream guide / perturbation call | `obs[<column>]` | `inputs.guide_assignments: {format: h5ad, guides_column: <column>}` |
| gene metadata | `var` | carried into the processed object |
| precomputed embedding | – | `inputs.embedding` (CSV only) |

One combined h5ad with a `var['feature_types']` column (as written by
`scanpy.read_10x_mtx(gex_only=False)`) can instead be given as
`inputs.multiplexed: {format: h5ad, path: …, feature_types: …}`, or with explicit
`slots` (`config/examples/h5ad_slots.yaml`).

RNA + guides only, no protein, no metadata: `inputs.protein.required: false`,
`inputs.metadata.required: false`, `alignment.required: [rna]`. `${DATA_ROOT}` and
`${RESULTS_ROOT}` default to `../data` and `../results` next to the repository.
Every key, default and rationale: `docs/DEFAULTS.md`; clustering, PS, lochNESS,
programs and the protein analyses can each be switched off under `analysis:`.

---

## Outputs

```
results/<run>/
├── report.html                      # self-contained; report.md alongside
├── processed/<name>.h5ad            # docs/PROCESSED_OBJECT.md
├── figures/<section>/<stage>/       # every figure, listed in tables/figure_manifest.csv
├── tables/                          # QC, filtering steps, perturbation coverage, protein QC ...
│   ├── cell_states/                 # cluster summary/composition, perturbation_cluster_enrichment
│   └── perturbation_effects/        # ps_targets, lochness_targets, gene_programs, perturbation_modules,
│                                    # protein_effects, *_association, perturbation_summary
└── logs/                            # run.log, resolved_config.yaml, run_manifest.json
```

Analysis results in the object: `obs['leiden']`, `obs['ps_score']`,
`obs['ps_quadrant']`, `obs['lochness_self']`; `obsm['ps_scores']`,
`obsm['lochness']`, `obsm['program_activity']`; the tables in `uns`; parameters and
provenance in `uns['petrubseq_protein']`. Full list: `docs/OUTPUTS.md`.

---

## Testing and validation

```bash
pytest -q            # 114 tests on synthetic data
```

The synthetic data have a known ground truth (targets with strong, weak and no
effects; enriched, depleted and null clusters; protein linked to some
perturbations) and are serialized in every supported layout; the tests assert that
each layout gives the same processed object and that each analysis recovers the
designed answer. In addition:

* `scripts/validate_processed.py` checks the processed-object contract (demo: 73 checks, 0 FAIL, 1 expected warning: no isotype controls);
* `scripts/compare_reference_perturbation.py` reproduces PS, lochNESS and the program/module effect matrix of the group's reference implementation exactly on shared inputs;
* `scripts/compare_processed.py` compares runs slot by slot against a frozen SCP1064 (Frangieh et al. 2021, 218,331 cells) reference; changes are accepted only with zero unexplained numerical differences (`docs/REGRESSION.md`).

---

## Repository layout

```
src/petrubseq_protein/   io (adapters, readers, tenx, h5ad) · validation · preprocessing (align, guides, normalize, embeddings)
                         · qc · analysis (clustering, cluster_enrichment, ps_score, lochness, modules, protein_effects, concordance)
                         · reporting (figures, plots, html, markdown) · pipeline · cli · config
config/                  demo_papalexi.yaml · scp1064*.yaml · examples/ (mtx, multilane, h5ad_perturbseq, h5ad_perturb_cite_seq, ...)
demo/                    data/papalexi_eccite (bundled subset) · fetch/prepare/check scripts
docs/                    DEFAULTS · OUTPUTS · PROCESSED_OBJECT · CELL_STATES · PERTURBATION_EFFECTS · PERTURBATION_DISTANCE · RNA_PROTEIN_LEVELS · REGRESSION · HPC · datasets/ · reference/
scripts/                 validate_processed · compare_processed · compare_reference_perturbation
slurm/                   job scripts (docs/HPC.md)
tests/                   synthetic generator + tests
notebooks/               demo_papalexi_eccite.ipynb
```

---

## License

MIT (`LICENSE`). The bundled demo data derive from GEO GSE153056; cite the original publication when using them.
