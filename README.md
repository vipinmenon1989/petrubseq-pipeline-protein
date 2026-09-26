# petrubseq-pipeline-protein

A reproducible Perturb-CITE-seq pipeline for RNA, CRISPR guide and
optional ADT/protein measurements: input harmonization, QC, guide assignment,
RNA/protein preprocessing, Leiden clustering, perturbation-response analysis,
gene programs, protein effects and a self-contained report, in one command.

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

**2 · QC.** Prefilter (200 genes, genes in ≥ 3 cells), strict filter (≥ 500 genes,
≤ 20 % mito; protein and perturbation criteria optional) recorded step by step,
flags that never remove cells. RNA, protein (ADT depth, antibodies detected,
isotype background) and guide QC (guide UMIs, guides per cell, top-vs-second guide)
figures, before and after filtering.

**3 · Guide assignment.** Provided per-cell guide lists / metadata columns, or calls
from guide counts: `dominant` (top guide ≥ 3 UMIs and > 2× the runner-up, else
`ambiguous` / `unassigned`) or `threshold` (every guide ≥ 3 UMIs, multi-guide).
Classes `single_targeting`, `single_control`, `multi_targeting`, `multi_control`,
`mixed_control_targeting`, `ambiguous`, `unassigned`. When both a provided list and
guide counts exist the run asks you to choose and reports their agreement.

**4 · RNA and protein preprocessing.** Raw counts normalized once (10,000, log1p);
log-normalized input preserved and integer counts reconstructed when a library size
is given. Up to 3,000 HVGs, 50 PCs, 15-neighbour graph, UMAP. Protein: raw ADT
counts kept, per-protein CLR across cells, isotypes used for QC and excluded from the
protein PCA/UMAP; supplied normalized values preserved, never treated as counts.

**5 · Cell states** *(optional, `analysis.clustering`)*. Leiden on the existing RNA
graph (igraph backend, resolution 1.0, seeded); cluster sizes, composition by
perturbation class / sample / lane, depth; flags for clusters dominated by one
sample, by ambiguous cells or by library depth. **Perturbation × cluster
enrichment**: Fisher exact test of each target's single-guide cells vs non-targeting
controls per cluster, BH-FDR over all pairs, guide agreement per pair, optional
Cochran–Mantel–Haenszel across a sample/lane column.

**6 · Perturbation response** *(optional, `analysis.perturbation_effects`)*.
**PS**: per-cell perturbation-response score (Song et al. 2025), max-normalized
within target, with knockdown / escaper / non-responder quadrants and a
cross-target separation (AUC, Cohen's d). **lochNESS**: cluster-free local
enrichment of each perturbation in PCA space (Huang et al. 2023), k = 300 capped
at 10 % of cells, permutation null, control-referenced delta. Cluster enrichment
and lochNESS are reported side by side, never combined.

**7 · Gene programs and perturbation modules** *(optional)*. Perturbation × response-
gene log2FC matrix vs non-targeting controls; genes clustered into programs
(`P1..`, 1 − Pearson), perturbations into modules (`M1..`, 1 − Spearman), after
Zhou et al. 2023; per-cell program activity. Numbered, never named.

**8 · Protein effects and RNA–protein concordance** *(optional, needs protein)*.
Per target × protein: CLR mean difference, Cohen's d, Mann–Whitney p, BH-FDR,
sign consistency per sample and per guide. Concordance: PS ↔ protein within target
(Spearman, FDR), lochNESS ↔ protein across targets, program activity ↔ protein
per cell and per target, plus one integrated target-level table. Associations only.

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

Formats: `dense_csv`, `mtx`, `10x_h5`, `h5ad` (slot-addressed). One config per
layout in `config/examples/`. Start from the documented defaults:

```bash
petrubseq-protein init-config my_run.yaml
petrubseq-protein validate --config my_run.yaml     # config + every input exists
petrubseq-protein audit --config my_run.yaml        # cell-ID overlap, value states, protein panel
```

Minimal combined-10x config:

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

RNA + guides only, no protein, no metadata: `inputs.protein.required: false`,
`inputs.metadata.required: false`, `alignment.required: [rna]`. `${DATA_ROOT}` and
`${RESULTS_ROOT}` default to `../data` and `../results` next to the repository.
Every key, default and rationale: `docs/DEFAULTS.md`.

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
config/                  demo_papalexi.yaml · scp1064*.yaml · examples/
demo/                    data/papalexi_eccite (bundled subset) · fetch/prepare/check scripts
docs/                    DEFAULTS · OUTPUTS · PROCESSED_OBJECT · CELL_STATES · PERTURBATION_EFFECTS · REGRESSION · HPC · datasets/ · reference/
scripts/                 validate_processed · compare_processed · compare_reference_perturbation
slurm/                   job scripts (docs/HPC.md)
tests/                   synthetic generator + tests
notebooks/               demo_papalexi_eccite.ipynb
```

---

## Citation

* Song B. et al. Decoding heterogeneous single-cell perturbation responses. *Nat Cell Biol* 27, 493–504 (2025) — PS.
* Huang X. et al. Single-cell, whole-embryo phenotyping of mammalian developmental disorders. *Nature* 623, 772–781 (2023) — lochNESS.
* Zhou P. et al. Single-cell CRISPR screens in vivo map T cell fate regulomes in cancer. *Nature* 624, 154–163 (2023) — gene programs / perturbation modules.
* Traag V.A., Waltman L., van Eck N.J. From Louvain to Leiden. *Sci Rep* 9, 5233 (2019) — clustering.
* Papalexi E. et al. Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens. *Nat Genet* 53, 322–331 (2021) — demo data, GEO GSE153056.
* Frangieh C.J. et al. Multimodal pooled Perturb-CITE-seq screens in patient models define mechanisms of cancer immune evasion. *Nat Genet* 53, 332–341 (2021) — SCP1064 regression data.

Engineering conventions follow [weili-lab/perturbseq-pipeline](https://github.com/weili-lab/perturbseq-pipeline).

## License

MIT (`LICENSE`). The bundled demo data derive from GEO GSE153056; cite the original publication when using them.
