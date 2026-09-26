# Default processing parameters and their rationale

Every value below is a config key with a documented default in
`src/petrubseq_protein/config.py`; `config/scp1064.yaml` shows them all. The
priority for choosing a default was: (1) the Wei Li Perturb-seq pipeline where
transferable, (2) Scanpy / single-cell practice, (3) CITE-seq practice,
(4) validation on SCP1064. Defaults are *stable general-purpose* choices, not
values tuned for a pleasing UMAP.

## Comparison with weili-lab/perturbseq-pipeline

| stage | Wei Li Perturb-seq pipeline | petrubseq-pipeline-protein (this project) |
|---|---|---|
| inputs | 10x MTX / H5 / h5ad, guide-count matrices, guide FASTQs | dense text matrices (SCP format), 10x MTX / H5 (combined or per modality, multi-lane), h5ad (feature types or slots), guide-count matrices, provided guide assignments; no FASTQ processing |
| guide calling | dominant-guide rule (`min_umi 3`, `dominance_ratio 2`) → `targeting` / `non-targeting` / `ambiguous` / `unassigned` | same dominant rule (plus a `threshold` rule for multi-guide lists) or a provided assignment, chosen explicitly (`perturbation.assignment.source`); classes `single_targeting` / `single_control` / `multi_targeting` / `multi_control` / `mixed_control_targeting` / `ambiguous` / `unassigned`; detected guide lists kept separately |
| RNA QC | `calculate_qc_metrics` (mt, ribo, hb), MAD 3 flags, prefilter 200/3, strict filter ≥ 1000 genes / < 20 % mt, Scrublet doublets | same metrics (mt, ribo), MAD 5 flags, prefilter 200/3, strict filter ≥ 500 genes / ≤ 20 % mt with an audit table; before/after figures; no doublet detection |
| protein QC | – | total / targeting / isotype ADT, proteins detected, isotype %, MAD + extreme flags, per-antibody background flag |
| perturbation QC | guide QC report; `perturbation.min_cells_per_target: 10` | `guide_counts`, `target_counts`, condition × guide/target coverage with `low_coverage` flags at 10 cells |
| normalization | `normalize_total(target_sum=null)` + `log1p`; raw kept in `layers['counts']`, lognorm in `layers['lognorm']` | state detection; raw input → `normalize_total(1e4)` + `log1p` with `layers['counts']`; log-normalized input preserved with `layers['reconstructed_counts']` when recoverable |
| HVG | 3000, Scanpy default flavor (`seurat`) on log data | 3000, `seurat`, on log data (`auto`: only when genes > 3000) |
| scaling | scaled HVG working copy, clip 10; `X` untouched | same |
| PCA | 50 PCs, arpack, zero-centred | same (`rna.n_pcs: 50`) |
| neighbors | 15 neighbours on PCA | same (`neighbors.n_neighbors: 15`) |
| UMAP | `min_dist 0.5`, seed `run.seed = 0` | same (`umap.min_dist 0.5`, `spread 1.0`, `compute.seed 0`) |
| clustering | Leiden 1.0 (used by enrichment / lochNESS stages) | **none** (no preprocessing/QC purpose in this phase) |
| batch correction | optional Harmony | none (SCP1064 has no batch labels) |
| protein representation | – | CLR (per protein across cells) of raw counts, targeting antibodies only, scaled, 10 PCs, 15 neighbours, UMAP |
| joint representation | – | optional `concat_pcs`, off by default (see below) |
| outputs | processed h5ad, HTML + MD report, tables, figures, run manifest | processed h5ad, self-contained report.html (+ report.md), CSV tables, staged figures with a manifest, run manifest, resolved config, optional archive |

## Inputs (formats, multiplexed matrices, lanes)

Every input layout is converted by `io/adapters.py` into one canonical model
(RNA, raw ADT counts, provided normalized protein, guide counts, provided
guide assignments, metadata, embedding); nothing after the "load data" stage
knows the source format. `config/examples/*.yaml` show each layout.

| key | default | rationale |
|---|---|---|
| `inputs.<modality>.format` | `dense_csv` | `dense_csv` (SCP-style text, unchanged v0.1 reader), `mtx` (10x directory), `10x_h5` (Cell Ranger v3 feature-barcode HDF5; v2 read as gene expression), `h5ad` (`slot`: `X` / `layers` / `obsm` + `key`). Separate inputs: `rna`, `protein` (normalized), `protein_counts`, `guide_counts`. |
| `inputs.multiplexed` | unset | one combined matrix read **once** and split by `feature_types` (`rna` / `protein` / `guide` label lists; empty list = modality absent) or, for h5ad, by `slots`. A feature type present in the file but claimed by no modality is an error (nothing is dropped silently). Declaring a modality both here and as a separate input is a config error. Antibody and guide features from a combined matrix are always raw counts. |
| `inputs.multiplexed.var_names` / `inputs.<m>.var_names` | `id` | RNA feature index: `id` (Ensembl, unique) or `name` (gene symbols, needed for `rna.mito_prefix`; duplicates are made unique `-1`, `-2` ... and counted). Antibodies and guides always keep their unique feature IDs. |
| `inputs.*.lanes` | `[]` | `[{id, path}]`; lanes must share one feature set. Cell IDs become `<barcode>__<lane id>` (single-lane inputs are never suffixed); `obs['barcode_original']` and `obs['lane_id']` are written and the transformation is recorded in provenance. Lane IDs must be unique and the global IDs must be unique. |
| `inputs.lane_metadata` | unset | per-lane sample sheet joined onto cells by `lane_column` (`lane_id`); every lane needs exactly one row; its columns become metadata columns usable through `columns.*`. |
| `inputs.*.state` | `auto` | declared value state `auto` / `raw_counts` / `normalized`; feature-barcode matrices are always `raw_counts`. `normalized` maps to the RNA-side `log_normalized` hint of `rna.input_state` at the adapter boundary (the two vocabularies were kept as in v0.1). Declared raw counts with negative or non-integer values are an error; declared normalized values that are all non-negative integers give a warning. |
| `inputs.protein.required` | `true` | *a protein representation must exist*: raw ADT counts **or** normalized values satisfy it, whatever the format; neither is an error. `inputs.protein_counts.required` demands raw counts specifically. `alignment.required: [.., protein, ..]` is likewise satisfied by counts alone. Raw counts -> CLR, count-based QC; normalized-only -> values preserved, normalized-only QC (never treated as counts). |
| `protein.feature_table` | unset | csv/tsv with `feature_id` (canonical antibody name) and any of `antibody_name`, `protein_name`, `gene_symbol`, `clone`, `isotype`, `isotype_control`; merged into `uns['protein_features']` on top of the input's own feature metadata (10x `features.tsv`). Values never supplied stay empty (= not available); nothing is inferred. |

## Guide assignment and guide calling

| key | default | rationale |
|---|---|---|
| `perturbation.assignment.source` | `auto` | `provided` (guide list / metadata column drives `obs['guides']`; guide counts feed diagnostics and a comparison), `guide_counts` (calls from the count matrix drive it; a provided list is kept and compared), `auto` (the only available source; **fails** when both exist, because silently choosing one would hide a real decision). Both sources are always preserved; provenance records the configured and effective source. |
| `perturbation.assignment.method` | `dominant` | weili-lab/perturbseq-pipeline rule: top guide with `>= min_umi` UMIs and `> dominance_ratio x` the runner-up (`max_second_umi >= 0` adds a multiplet gate). Top guide below `min_umi` -> `unassigned`; reaches `min_umi` but no dominance -> `ambiguous` (the reference labels any non-zero non-dominant cell ambiguous; here ambiguity requires evidence of a real guide). `threshold`: every guide with `>= detection_min_umi` UMIs is assigned (multi-guide lists for high-MOI screens). |
| `min_umi` / `dominance_ratio` / `max_second_umi` | 3 / 2.0 / -1 | reference defaults. |
| `detection_min_umi` | 3 | a guide is *detected* at >= 3 UMIs: `obs['n_guides_detected']`, `obs['guides_detected']` (always recorded, never collapsed into the dominant call) and the `threshold` rule. |
| `ambiguous_label` | `ambiguous` | `obs['perturbation']` / `perturbation_class` value for ambiguous cells under `guide_counts` + `dominant`. |

With guide counts the per-cell diagnostics `guide_total_counts`,
`n_guides_detected`, `guide_top`, `guide_top_count`, `guide_second_count`,
`guide_dominant_call` are always written and the `guide_count_diagnostics`
figure is drawn. When a provided assignment and guide counts both exist,
`uns['petrubseq_protein']['perturbations']['assignment']['agreement']` holds
the set agreement (provided guides vs detected guides) and the single-guide
vs dominant-call agreement; disagreements are a warning, never a silent fix.

## RNA

| key | default | rationale |
|---|---|---|
| `rna.input_state` | `auto` | integer values → raw counts; constant `sum(expm1(x))` per cell → log-normalized. SCP1064 detects as ln(TPM+1) with scale 1e6. |
| `rna.normalize` | `auto` | normalize only raw counts; log-normalized input is preserved (never normalized twice). `always` is refused on non-raw input. |
| `rna.target_sum` | 10000 | Scanpy convention for raw-count input (Wei Li uses the median library size, `null`); irrelevant for SCP1064. |
| `rna.log1p` | true | standard. |
| `rna.reconstruct_counts` | `auto` | when a per-cell library size is available, integer counts are recovered as `round(expm1(x) * total / scale)` and accepted only if all entries are within `reconstruct_tolerance` of integers. Stored as `layers['reconstructed_counts']`, never as `counts`. |
| `rna.hvg.enabled` / `n_top_genes` / `flavor` | `auto` / 3000 / `seurat` | Wei Li `cluster.n_top_genes: 3000` with Scanpy's default dispersion-based flavor, which expects log-normalized data and therefore works on both raw-derived and provided ln(TPM+1) values. `seurat_v3` (needs counts) is available. |
| `rna.scale` / `scale_max_value` | true / 10 | Wei Li scales the HVG working copy (clip 10) before PCA; `X` is never modified. |
| `rna.n_pcs` | 50 | Wei Li / Scanpy default. On the full SCP1064 run 50 PCs explain 16.5 % of the scaled-HVG variance (PC1 2.6 %), the UMAP separates the three conditions, and PC1 correlates with log library size (|r| 0.85) but not with mito % (0.08) - the usual depth axis of log-normalized data, which the Wei Li pipeline also leaves uncorrected. |

## Protein (ADT)

| key | default | rationale |
|---|---|---|
| `protein.normalization` | `auto` → `provided` if a normalized matrix is shipped, else `clr` | never overwrite a published normalization; SCP1064's values are the authors' `max(0, ln((x+1)/(isotype+1)))`, reproduced exactly from the raw counts. |
| `protein.extra_representations` | `[clr]` | CLR is the standard CITE-seq transform (Stoeckius 2017; Seurat; muon). Stored as `obsm['protein_clr']` so both representations are available. |
| `protein.embedding_representation` | `auto` → `clr` when raw counts exist, else `provided` | a generic default must not depend on a dataset-specific isotype map; CLR only needs counts. On SCP1064 CLR and provided values behave similarly (protein PC1 vs log ADT depth: |r| 0.84 vs 0.76; condition silhouette 0.078 vs 0.070), so the generic choice costs nothing. |
| `protein.clr_axis` | `cells` | per-protein CLR across cells (Seurat `margin = 2`), the recommendation for small panels where per-cell compositional CLR is unstable (Hao et al. 2021). `features` (per-cell, Stoeckius 2017) is available. |
| `protein.exclude_isotypes_from_embedding` | true | isotype controls are technical background, not biology; they stay in `obsm['protein_counts']` and the QC columns (`use_isotypes_for_qc: true`). |
| `protein.scale` / `scale_max_value` | true / 10 | put antibodies with very different dynamic ranges (CD44 median 523 vs CD202b median 0) on a comparable scale before PCA. |
| `protein.n_pcs` | 10 | capped at `n_features - 1`; with ~20 targeting antibodies 10 PCs retain ~80 % of the variance (SCP1064: 78 %); 50 would be meaningless. |

Known properties, documented rather than "fixed":

* in every ADT representation the first protein PC correlates with total ADT
  counts per cell (|r| = 0.83 on the full SCP1064 run, PC1 = 36 % of variance,
  all loadings positive). This is the usual CITE-seq depth/background effect;
  denoising (e.g. dsb) is a candidate for a later version and is not applied
  by default.
* the protein UMAP of a 20-plex panel with low integer counts shows banding
  (discrete count levels of the dominant antibodies produce many tied
  distances). It is a diagnostic picture only; neighbourhood-based analyses
  should use `X_pca_protein`, not the UMAP coordinates.

## Neighbors / UMAP

| key | default | rationale |
|---|---|---|
| `neighbors.n_neighbors` | 15 | Wei Li / Scanpy default. |
| `neighbors.n_pcs` | null (= `rna.n_pcs`) | Wei Li uses all PCs. |
| `umap.min_dist` / `spread` | 0.5 / 1.0 | Wei Li `cluster.umap_min_dist: 0.5`; UMAP default spread. |
| `compute.seed` | 0 | single seed for PCA (arpack), neighbors and UMAP; Wei Li `run.seed: 0`. Same config + seed → identical PCA/UMAP/tables (tested). |
| `umap.min_features` | 5 | skip a modality UMAP when the feature space is too small to be meaningful. |

## Multimodal representation (`multimodal.enabled: false`)

Evaluated on SCP1064 (4,500-cell smoke subset, confirmed on the full run's
*Cross-modality diagnostics*: 20,000 sampled cells share 0.6 % of their 15
nearest neighbours between the RNA and protein graphs, random expectation
0.08 %): RNA and protein 15-NN graphs share only ~1-2 % of neighbours, protein
neighbourhoods are dominated by ADT depth, and block-normalized concatenation
of RNA + protein PCs (weight 1) lowers condition separation (silhouette
0.146 → 0.121) while mixing the two graphs (28 % / 21 % neighbour overlap). A
joint embedding therefore adds no diagnostic value by default. The simple
`concat_pcs` method (RNA PCs and protein PCs each divided by the square root
of their total variance, protein block times `protein_weight`, concatenated
into `obsm['X_multimodal']` → neighbors → `obsm['X_umap_multimodal']`) is
implemented and can be switched on per run. WNN-style weighting and MuData
were deliberately not introduced.

## QC: prefilter, strict filtering, flags

The run has three clearly separated QC concepts. All thresholds below are
general defaults derived from the Wei Li pipeline and common practice; they
are not universal biological truths and should be reviewed per dataset (the
before-filter figures exist so that this review is possible).

| stage | keys | default | rationale |
|---|---|---|---|
| **prefilter** (permissive) | `qc.prefilter.enabled`, `min_genes_per_cell`, `min_cells_per_gene` | true, 200, 3 | Wei Li `qc.min_genes_per_cell` / `min_cells_per_gene`; removes obvious empty droplets and never-detected genes *before* metrics so that the QC figures show real cells. Recorded as `prefilter_*` audit steps. |
| **strict filtering** | `qc.filter.enabled` | true | configured cell removal, after the before-filter figures; every step is an audit row (`rna_*`, `protein_*`, `perturbation_*`). |
| | `qc.filter.rna.min_genes` | 500 | between Wei Li's `min_genes_final: 1000` (deep ESC data) and the prefilter; conservative for most 10x data. |
| | `qc.filter.rna.min_counts` | null | off unless a dataset needs it. |
| | `qc.filter.rna.max_pct_mt` | 20 | Wei Li `qc.max_pct_mt`. |
| | `qc.filter.protein.*` | all null / false | protein thresholds are dataset-specific (panel size, depth); shipped off, available: `min_total_counts`, `min_proteins_detected`, `max_pct_isotype`, `remove_extreme_counts`. |
| | `qc.filter.perturbation.cells` | `all` | keep every cell; `assigned` (>= 1 guide) and `single_guide` are explicit choices. Coverage is never a filter. |
| **flags** (never remove) | `qc.flags.rna_n_mads`, `protein_n_mads` | 5 | MAD outlier flags on log1p depth / genes; Wei Li filters at 3 MADs, 5 marks only clear outliers. |
| | `qc.flags.extreme_fold` | 10 | `protein_extreme_counts`: total ADT > 10 x the 99th percentile (antibody aggregates; SCP1064 max 878k vs median 829). |
| | `qc.flags.background_min_fraction_above_isotype` | 0.5 | antibody is `background_dominated` when its signal exceeds the matched isotype in fewer than half of the cells. |

The same `qc.filter.*` thresholds define the `rna_qc_fail` / `protein_qc_fail`
flags, so a run with `filter.enabled: false` records exactly what *would* have
been removed. Outputs: `tables/qc_filtering_steps.csv` (step, category,
threshold, cells/genes before/after/removed), `tables/cell_qc_prefilter.csv.gz`
(every cell that entered strict filtering with all metrics, flags,
`qc_retained` and `removed_by` = first removing step), before/after figure
sets, and `uns['petrubseq_protein']['qc']`.

**SCP1064 regression configs** (`config/scp1064.yaml`, `config/scp1064_smoke.yaml`)
set `prefilter.enabled: false` and `filter.enabled: false` explicitly: the
authors already filtered the data (>= 200 genes, <= 18 % mito, genes in >= 200
cells), one retained gene is detected in only 2 cells (so the prefilter would
change the gene set) and the v0.1 object is the frozen benchmark. These are
dataset-specific regression settings, not the public defaults.

## Perturbation QC

| key | default | rationale |
|---|---|---|
| `perturbation.min_cells_per_guide` | 10 | flag threshold; no per-guide minimum exists in Wei Li, 10 mirrors its per-target minimum. |
| `perturbation.min_cells_per_target` | 10 | Wei Li `perturbation.min_cells_per_target: 10`. |
| `perturbation.preserve_multiguide` / `preserve_unassigned` | true / true | keep and label every cell; later analyses decide (the SCP1064 authors modelled all guides jointly). |
| `perturbation.control_classes` | `non_targeting` patterns | dataset configs add classes (SCP1064: `intergenic` for `ONE_NON-GENE_SITE`). |

Flags are applied per guide, per target and per (guide, condition) /
(target, condition) pair; nothing is removed.

## Clustering

Not included. Leiden clustering in the Wei Li pipeline serves downstream
stages (cluster enrichment, lochNESS). It has no preprocessing or QC role in
this pipeline, adds `leidenalg`/`igraph` dependencies and a resolution
parameter with no principled default, so it is left to the analysis phase.

## Perturbation effects (`analysis.perturbation_effects`, Stage E; off by default)

| key | default | rationale |
|---|---|---|
| `enabled` | `false` | downstream biology is opt-in; preprocessing runs are unchanged without it |
| `control_classes` | `[non_targeting]` | controls = single-guide cells of these control classes; ambiguous / multi-guide cells are never used |
| `target_gene_map` | `{}` | target label → RNA gene symbol where they differ (dataset config only) |
| `ps.top_n_genes` / `scale_factor` / `ps_threshold` | 100 / 3.0 / 0.5 | PS_python (`pertps`) and reference pipeline defaults |
| `ps.expression_cut` | `mean` | reference pipeline: the control median is degenerate (0) for most genes |
| `ps.min_cells_per_target` / `min_control_cells` / `min_pct_expressing_control` | 10 / 10 / 1 % | reference defaults |
| `lochness.n_neighbors` / `n_pcs` | 300 / 20 | pertTF / reference defaults, neighbours in `X_pca` |
| `lochness.max_k_fraction` | 0.1 | caps k at 10 % of the cells so small objects keep a local neighbourhood (added; recorded when applied) |
| `lochness.n_permutations` | 200 | label-permutation null (added, following Huang et al. 2023) |
| `modules.min_cells_per_perturbation` / `min_perturbations` / `min_genes` | 20 / 5 / 10 | reference defaults |
| `modules.program_correlation` / `module_correlation` / `linkage_method` | pearson / spearman / average | reference |
| `modules.n_programs` / `n_modules` / `cluster_distance_threshold` | 4 / 9 / 0.7 | reference |
| `modules.gene_selection` | `response` | union of each target's top 100 response genes at FDR < 0.05 (the reference uses Leiden markers; this pipeline has no clustering) |
| `modules.log2fc_pseudocount` | 1.0 | Seurat `FoldChange` convention; the reference 1e-9 gives |log2FC| > 20 for genes undetected in a small group |
| `modules.min_pct_cells_expressing` | 5 % | panel genes detected in ≥ 5 % of analysed cells (added, same reason) |
| `modules.de_lfc_threshold` / `de_fdr_alpha` | 0.5 / 0.05 | reference DE gate |
| `protein.representation` | `protein` | the primary normalized protein matrix (CLR of counts by default); never counts |
| `protein.min_cells_per_target` / `min_control_cells` / `fdr_alpha` | 10 / 10 / 0.05 | same support rule as PS |
| `concordance.min_cells` / `fdr_alpha` | 20 / 0.05 | minimum cells for a within-target Spearman correlation |

## Cell states and perturbation × cluster enrichment (`analysis.clustering`; on by default)

| key | default | rationale |
|---|---|---|
| `enabled` | `true` | part of the standard run; `false` skips clustering and enrichment (the SCP1064 regression configs do, so that the frozen v0.1 reference stays comparable) |
| `key` | `leiden` | obs column for the labels; an existing column is an error unless `overwrite: true` |
| `resolution` | 1.0 | reference (weili-lab/perturbseq-pipeline) default; results are conditional on it |
| `n_iterations` | 2 | reference; `-1` iterates until convergence |
| `neighbors_key` | `rna` | the RNA neighbour graph of the representation stage (no second preprocessing) |
| `enrichment.control` | `non_targeting` | same control rule as every perturbation analysis; `other` (all other targets) is the reference default, available explicitly |
| `enrichment.control_classes` | `[non_targeting]` | control classes counted as controls |
| `enrichment.min_cells_per_target` / `min_cells_per_cluster` / `min_control_cells` | 10 / 20 / 10 | reference |
| `enrichment.min_control_cells_in_cluster` | 10 | reference `min_reference_cells` (low-power flag) |
| `enrichment.min_cells_per_guide` | 5 | reference guide-concordance minimum |
| `enrichment.odds_pseudocount` | 0.5 | Haldane-Anscombe, display/ranking only; the Fisher odds ratio is stored unchanged |
| `enrichment.fdr_alpha` | 0.05 | BH over all tested (target, cluster) pairs of the run |
| `enrichment.stratify_by` | null | any obs column; adds a CMH test beside Fisher (own BH family) |
| `enrichment.n_permutations` | 1000 | omnibus permutation p-value (reference) |
