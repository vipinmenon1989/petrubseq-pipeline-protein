# petrubseq-pipeline-protein

A reproducible pipeline for **Perturb-CITE-seq** and other multimodal Perturb-seq
screens: CRISPR perturbations read out by single-cell RNA, surface-protein (ADT)
and guide measurements. One command turns the raw matrices into a validated,
cell-aligned, QC-annotated AnnData object, per-perturbation effect analyses on
RNA and protein, and a self-contained HTML report.

```bash
petrubseq-protein run --config config/demo_papalexi.yaml
```

## What the pipeline does

```
INPUT            dense text / 10x MTX / 10x H5 / h5ad; combined or per-modality; one or many lanes
  ↓
RNA + protein + guides        one canonical input model, validated (cell barcodes, feature IDs, value states)
  ↓
QC / filtering                prefilter, strict filter with an audit table, flags; before/after figures
  ↓
representations               RNA normalize + HVG + PCA + UMAP; protein CLR + PCA + UMAP; cross-modality checks
  ↓
perturbation assignment       provided guide lists or guide calling from guide counts (dominant / threshold rule)
  ↓
PS / lochNESS                 per-cell response strength; local cell-state enrichment of each perturbation
  ↓
gene programs                 perturbation × response-gene effects → gene programs and perturbation modules
  ↓
protein perturbation effects  each perturbation's effect on each measured protein
  ↓
RNA–protein concordance       PS ↔ protein, lochNESS ↔ protein, gene program ↔ protein
  ↓
self-contained report         report.html (+ report.md), processed .h5ad, CSV tables, figures, provenance
```

Everything is driven by one YAML file; unknown keys are errors, every default is
documented (`docs/DEFAULTS.md`), and dataset-specific facts (file names, column
names, control-guide patterns, target aliases) live only in the YAML.

## Quick start

```bash
git clone https://github.com/vipinmenon1989/petrubseq-pipeline-protein.git
cd petrubseq-pipeline-protein

conda env create -f environment.yml
conda activate petrubseq-protein
pip install -e .

petrubseq-protein run --config config/demo_papalexi.yaml
```

The repository bundles a ~20 MB deterministic 1,800-cell subset of the Papalexi
et al. ECCITE-seq screen (`demo/data/papalexi_eccite/`), so the demo runs from a
fresh clone without downloading anything. It takes about a minute and under 1 GB
of memory and writes `../results/demo/papalexi_eccite/` (set
`PETRUBSEQ_RESULTS_ROOT` to write elsewhere). Open `report.html` there; it embeds
every figure. Check the processed object with

```bash
python scripts/validate_processed.py ../results/demo/papalexi_eccite/processed/papalexi_eccite_demo_processed.h5ad
```

## Demo dataset

Papalexi E., Mimitou E.P., Butler A.W., et al. *Characterizing the molecular
regulation of inhibitory immune checkpoints with multimodal single-cell screens.*
Nat Genet 53, 322–331 (2021). GEO **GSE153056**; the demo uses the pooled
ECCITE-seq screen (GSM4633614–GSM4633618) in IFN-γ stimulated THP-1 cells.

| | bundled subset |
|---|---|
| cells | 1,800 (of 20,729; HTO singlets, stratified by hashed sample × guide target, seed 0) |
| RNA | 18,649 genes, raw UMI counts |
| protein | 4 ADTs: CD86, PD-L1 (`PDL1`), PD-L2 (`PDL2`), TIM-3 (`CD366`); no isotype controls |
| guides | 111 sgRNAs: 25 gene targets × 3–4 guides, 9 non-targeting (`NT`), 1 eGFP |
| design | 3 hashed samples, 8 lanes; every sample IFN-γ stimulated |

The subset is a software demonstration dataset, not a reanalysis of the
publication: with 1,800 cells and four proteins its statistics are illustrative.
Provenance: `demo/data/papalexi_eccite/README.md`; audit of the source files:
`docs/datasets/PAPALEXI_ECCITE_AUDIT.md`.

## Supported inputs

| layout | config |
|---|---|
| dense CSV / TSV, either orientation, one file or several that partition the cells | `inputs.<modality>.format: dense_csv` |
| 10x MTX directory, 10x feature-barcode H5 | `format: mtx` / `10x_h5` |
| AnnData `.h5ad` with a feature-type column or explicit slots (`X`, `layers`, `obsm`, `obs`) | `format: h5ad` |
| one combined matrix (Gene Expression / Antibody Capture / CRISPR Guide Capture) read once and split | `inputs.multiplexed` |
| separate RNA / ADT / guide matrices | `inputs.rna`, `inputs.protein_counts`, `inputs.protein`, `inputs.guide_counts` |
| several lanes (cell IDs `<barcode>__<lane>`, per-lane sample sheet) | `lanes: [{id, path}]`, `inputs.lane_metadata` |
| provided guide assignments, per-cell metadata, antibody annotation, precomputed embedding | `inputs.guide_assignments`, `inputs.metadata`, `protein.feature_table`, `inputs.embedding` |

Examples for every layout: `config/examples/*.yaml`.

## Pipeline stages

* **Input harmonization.** Every layout becomes one canonical model (`io/adapters.py`) and is
  validated (`validation.py`): duplicate barcodes or feature IDs, negative or non-integer "counts",
  unclaimed feature types, missing h5ad slots, metadata that shares no cell with the RNA matrix.
* **QC.** A permissive prefilter (200 genes, genes in ≥ 3 cells), a strict filter (≥ 500 genes,
  ≤ 20 % mitochondrial; protein and perturbation criteria optional) recorded step by step in
  `qc_filtering_steps.csv`, and flags that never remove cells. RNA and protein QC figures are drawn
  before and after filtering.
* **RNA.** The value state is detected; raw counts are normalized once (`normalize_total` + `log1p`)
  after filtering, log-normalized input is kept and integer counts reconstructed when a library
  size is provided. 3,000 HVGs, scaled copy for PCA, 50 PCs, 15-neighbour graph, UMAP.
* **Protein.** Raw ADT counts are kept; per-protein CLR across cells is the default representation;
  provided normalized values are preserved and never treated as counts; isotype controls are
  excluded from embeddings and used for QC; structured antibody annotation.
* **Guides.** Provided assignments or guide calling from guide counts (dominant rule: top guide
  ≥ 3 UMIs and > 2 × the runner-up, else `ambiguous` / `unassigned`; or a threshold rule for
  multi-guide lists). Both sources are kept and compared when both exist.
* **Perturbation effects** (optional, `analysis.perturbation_effects`), **gene programs** and
  **multimodal analysis** — below.
* **Reporting.** `report.html` with summary cards, warnings, inputs, RNA / protein / perturbation QC,
  representations, perturbation effects, outputs and provenance.

Perturbed cells in every effect analysis are single-guide targeting cells of one target; controls
are single-guide non-targeting cells; ambiguous, multi-guide and unassigned cells enter neither group.

## Perturbation-response score (PS)

PS measures, per cell, how strongly the cell shows the transcriptional response of its
perturbation. For each target the top 100 response genes (t-test vs non-targeting cells) define
a signature; the signature is the OLS effect of the perturbation on those genes, and each cell is
projected onto it, centred on the controls and scaled to [0, 1]. Combined with the target gene's own
expression it separates successful knockdowns from escapers and non-responders. Scores are
comparable within a target; across targets the report gives the separation of perturbed from control
scores. Song B., Liu D., Dai W., et al. *Decoding heterogeneous single-cell perturbation responses.*
Nat Cell Biol 27, 493–504 (2025); implementation follows PS_python
(https://github.com/weili-lab/PS_python).

## lochNESS

lochNESS asks, for every cell and perturbation, whether the perturbation is over- or
under-represented among the cell's nearest neighbours in PCA space:
`local fraction / overall fraction − 1` (0 = chance). Averaged over a target's own cells it says
whether that perturbation occupies a distinct transcriptional state; a label-permutation null and a
control-referenced difference accompany it. Huang X., Henck J., Qiu C., et al. *Single-cell,
whole-embryo phenotyping of mammalian developmental disorders.* Nature 623, 772–781 (2023);
implementation ported from pertTF (https://github.com/davidliwei/pertTF) as in the group's
perturbseq-pipeline.

## Gene programs and perturbation modules

Both come from one perturbation × response-gene matrix of log2 fold changes against the
non-targeting controls, clustered along its two axes:

* a **gene program** (`P1..`) is a group of downstream genes whose responses are correlated across
  perturbations (genes clustered on 1 − Pearson);
* a **perturbation module** (`M1..`) is a group of perturbations with similar response profiles
  (perturbations clustered on 1 − Spearman).

The *perturbation target* is the gene a guide acts on; *response genes* are the genes that change
downstream. Programs and modules are numbered clusters, never given biological names. The design
follows Zhou P., Shi H., Huang H., et al. *Single-cell CRISPR screens in vivo map T cell fate
regulomes in cancer.* Nature 624, 154–163 (2023).

## Protein effects

For each perturbation and each measured protein, the CLR-normalized values of perturbed cells are
compared with non-targeting controls: difference of means (a log fold change on the CLR scale),
Cohen's d, Mann–Whitney U p-value, Benjamini–Hochberg FDR over all tests, and whether the direction
holds in each hashed sample and each guide.

## RNA–protein concordance

* **PS ↔ protein**: within one target, do cells with a stronger transcriptional response carry more
  or less of a protein? (Spearman over the target's cells; pooling targets is reported only as
  descriptive.)
* **lochNESS ↔ protein**: across targets, do perturbations that shift cell state more also shift a
  protein more? lochNESS is a population quantity, so this is a target-level comparison.
* **gene program ↔ protein**: does a program's per-cell activity track a protein (cell level), and do
  perturbations that move a program also move the protein (target level)?

These are correlations. They do not show that a program mediates a protein change, or that one
causes the other. `perturbation_summary.csv` puts every measure for one perturbation on one row.
Details and limitations: `docs/PERTURBATION_EFFECTS.md`.

## Outputs

```
results/<run>/
├── report.html                 self-contained report (figures embedded)
├── report.md                   the same content with linked figures
├── processed/<name>.h5ad       the processed object (docs/PROCESSED_OBJECT.md)
├── figures/<section>/<stage>/  every figure, listed in tables/figure_manifest.csv
├── tables/                     CSV tables; tables/perturbation_effects/ for the effect analyses
└── logs/                       run.log, resolved_config.yaml, run_manifest.json (provenance)
```

`docs/OUTPUTS.md` lists every table and figure.

## Configuration

```yaml
dataset: {name: my_screen, input_dir: ${DATA_ROOT}/my_screen}
inputs:
  multiplexed:
    format: mtx
    path: filtered_feature_bc_matrix
    feature_types: {rna: [Gene Expression], protein: [Antibody Capture], guide: [CRISPR Guide Capture]}
perturbation:
  guide_target_regex: "^(?P<target>.+?)_(?P<index>\\d+)$"
  control_classes: {non_targeting: ["^NTC$"]}
  assignment: {source: auto, method: dominant}
analysis:
  perturbation_effects:
    enabled: true
    target_gene_map: {PDL1: CD274}   # only where a target label differs from its gene symbol
output: {dir: ${RESULTS_ROOT}/my_screen}
```

`${DATA_ROOT}` and `${RESULTS_ROOT}` default to `../data` and `../results` next to the repository
(`PETRUBSEQ_DATA_ROOT`, `PETRUBSEQ_RESULTS_ROOT` override them); `${REPO_ROOT}` is the repository.
`petrubseq-protein init-config my.yaml` writes every key with its default. Other commands:
`petrubseq-protein audit --config CFG` (role-aware input audit) and `validate --config CFG`.

## Rebuilding the demo

Normal use does not need this. To rebuild the bundled subset from the public accession:

```bash
python demo/fetch_papalexi_data.py --dest ${DATA_ROOT}/ECCITE-seq            # accession-based GEO FTP; skips complete files
python demo/prepare_papalexi_demo.py --source ${DATA_ROOT}/ECCITE-seq --output /tmp/papalexi_eccite
python demo/check_demo_reproducibility.py /tmp/papalexi_eccite demo/data/papalexi_eccite   # byte-identical
```

## HPC

On a cluster, use the login node for development and lightweight inspection and SLURM for
substantial real-data runs. `docs/HPC.md` describes the policy used on the development cluster;
`slurm/` holds the job scripts. The package itself has no cluster-specific settings.

## Reproducibility

* A frozen reference object for the Frangieh et al. SCP1064 Perturb-CITE-seq dataset (218,331 cells;
  `config/scp1064.yaml`, `docs/datasets/SCP1064_AUDIT.md`) and its 4,500-cell smoke subset lives
  outside the repository; `scripts/compare_processed.py` compares any new run against it slot by
  slot, and a change is accepted only with zero unexplained numerical differences (`docs/REGRESSION.md`).
* The Papalexi demo is prepared deterministically (byte-identical on repeat) and runs through the
  same production command as real datasets.
* PS, lochNESS and the gene-program effect matrix reproduce the group's reference implementation
  exactly on shared inputs (`scripts/compare_reference_perturbation.py`).
* `pytest` runs 97 tests on synthetic data (inputs in every format, QC, reports, guide calling,
  every effect analysis).

## Citation

Methods: Song et al., Nat Cell Biol 2025 (PS); Huang et al., Nature 2023 (lochNESS);
Zhou et al., Nature 2023 (gene programs / perturbation modules). Data: Papalexi et al.,
Nat Genet 2021 (GSE153056, demo); Frangieh C.J., Melms J.C., Thakore P.I., et al. *Multimodal
pooled Perturb-CITE-seq screens in patient models define mechanisms of cancer immune evasion.*
Nat Genet 53, 332–341 (2021) (Single Cell Portal SCP1064). Engineering conventions follow
weili-lab/perturbseq-pipeline.

## License

MIT (`LICENSE`, as declared in `pyproject.toml`). The
bundled demo data are derived from the public GEO accession GSE153056; cite the original
publication when using them.
