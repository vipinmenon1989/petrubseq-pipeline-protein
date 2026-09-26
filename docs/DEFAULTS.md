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
| RNA QC | `calculate_qc_metrics` (mt, ribo, hb), prefilter 200/3, strict filter ≥ 1000 genes / < 20 % mt / genes in ≥ 3 retained cells | same metrics (mt, ribo, hb), MAD 5 flags, prefilter 200/3, **same strict filter** (≥ 1000 genes, < 20 % mt, optional % hb, genes in ≥ 3 cells re-applied) with an audit table; before/after figures |
| protein QC | – | total / targeting / isotype ADT, proteins detected, isotype %, MAD + extreme flags, per-antibody background flag |
| perturbation QC | guide QC report; `perturbation.min_cells_per_target: 10` | `guide_counts`, `target_counts`, condition × guide/target coverage with `low_coverage` flags at 10 cells |
| normalization | `normalize_total(target_sum=null)` + `log1p`; raw kept in `layers['counts']`, lognorm in `layers['lognorm']` | state detection; raw input → `normalize_total(target_sum=null` = median library size`)` + `log1p` with `layers['counts']` (same as the reference); log-normalized input preserved with `layers['reconstructed_counts']` when recoverable |
| HVG | 3000, Scanpy default flavor (`seurat`) on log data | 3000, `seurat`, on log data (`auto`: only when genes > 3000) |
| scaling | scaled HVG working copy, clip 10 (scanpy densifies the sparse block in float64); `X` untouched | same code path (bit-identical PCs, docs/reference/PARITY_RESULTS.md) |
| PCA | 50 PCs, arpack, zero-centred | same (`rna.n_pcs: 50`) |
| neighbors | 15 neighbours on PCA | same (`neighbors.n_neighbors: 15`) |
| UMAP | `min_dist 0.5`, seed `run.seed = 0` | same (`umap.min_dist 0.5`, `spread 1.0`, `compute.seed 0`) |
| clustering | Leiden 1.0 (igraph, 2 iterations, seed) before every analysis | same, as the first analysis stage (`analysis.clustering`) |
| perturbation strength | target's own expression vs ntc / other controls, KS + MWU, BH, hit call | same (`analysis.perturbation_effects.strength`) |
| cluster enrichment | Fisher / CMH per target x cluster, both arms, primary `other` | same |
| modules / programs | Leiden-marker panel, log2FC vs ntc (pseudocount 1e-9), Pearson/Spearman clustering, score_genes | same |
| PS | pertps score, quadrants, LDA embedding | same maths re-implemented (no `pertps` dependency), same order / skipping / LDA |
| lochNESS | k = 300 in `X_pca`, no null | same; permutation null and control delta are optional extras |
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
| `rna.target_sum` | `null` | median library size, the reference `cluster.target_sum: null`; `10000` is the Scanpy convention (the SCP1064 configs set it explicitly; irrelevant there because the input is log-normalized). |
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
| | `qc.filter.rna.min_genes` | 1000 | reference `qc.min_genes_final: 1000`. |
| | `qc.filter.rna.max_pct_mt` / `max_pct_hb` | 20 / null | reference: keep `pct_counts_mt < 20` (strict), optional haemoglobin cut. |
| | `qc.filter.rna.min_cells_per_gene` | 3 | reference re-applies `filter_genes(min_cells=3)` after the cell filters, so genes that lost their last cells leave the matrix (affects HVG/PCA). |
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

## Perturbation analyses (`analysis.perturbation_effects`; on by default)

The reference Perturb-seq analyses run in the reference order once the representations exist:
Leiden clustering → **perturbation strength** → perturbation × cluster enrichment → **modules /
programs** → **PS** → **lochNESS**, then the protein extension (protein effects, RNA–protein
concordance). Every default below is the reference value unless marked *extension*
(docs/reference/REFERENCE_PIPELINE_COMPLETE_AUDIT.md, numerical parity in
docs/reference/PARITY_RESULTS.md).

| key | default | rationale |
|---|---|---|
| `enabled` | `true` | the reference runs every analysis by default; `false` gives a preprocessing-only run |
| `control_classes` | `[non_targeting]` | the `ntc` arm = single-guide cells of these control classes; ambiguous / multi-guide cells are never used |
| `target_gene_map` | `{}` | target label → RNA gene symbol where they differ (dataset config only; the reference assumes label = symbol) |
| `strength.controls` / `primary_control` | `[ntc, other]` / `ntc` | reference: both arms reported, `ntc` drives ranking and the effective-knockdown call |
| `strength.min_cells_per_target` / `min_control_cells` / `min_pct_expressing_control` | 10 / 10 / 1 % | reference |
| `strength.fdr_alpha` / `max_log2fc_for_hit` | 0.05 / 0 | hit = KS BH-FDR < α **and** log2FC < 0 (reference) |
| `strength.top_n_report` / `umap_background_fraction` | 12 / 0.1 | reference figure settings |
| `ps.top_n_genes` / `scale_factor` / `ps_threshold` | 100 / 3.0 / 0.5 | PS_python (`pertps`) and reference defaults |
| `ps.expression_cut` | `mean` | reference: the control median is degenerate (0) for most genes |
| `ps.min_cells_per_target` / `min_control_cells` / `min_pct_expressing_control` | 10 / 10 / 1 % | reference |
| `ps.score_targets_without_gene` | `false` | reference skips targets whose gene is not in `var`; `true` scores them with `not_applicable` quadrants (*extension*) |
| `ps.compute_lda_umap` / `lda_n_pcs` / `lda_max_genes` / `lda_highlight_threshold` | true / 40 / 5000 / 0.8 | reference supervised LDA embedding (PS_python `compute_lda_umap`) |
| `lochness.n_neighbors` / `n_pcs` | 300 / 20 | pertTF / reference defaults, neighbours in `X_pca` (`X_pca_harmony` when present) |
| `lochness.max_k_fraction` | 1.0 | reference: k = min(300, n − 1); a value < 1 caps k on small objects (*extension*, recorded when applied) |
| `lochness.n_permutations` | 0 | reference has no null; > 0 adds a seeded label-permutation z / p / FDR (*extension*) |
| `lochness.noise_delta` | 0 | reference option (pertTF adds 1e-4 for model training) |
| `modules.gene_selection` | `cluster_markers` | reference panel: union of the top 100 positive Wilcoxon markers of each Leiden cluster; `response` (top response genes per target) and `hvg` are options |
| `modules.cluster_key` / `n_marker_genes_per_cluster` / `marker_method` | leiden / 100 / wilcoxon | reference |
| `modules.control` | `ntc` | reference (falls back to `other`, the leave-one-target-out targeting cells, when no controls exist) |
| `modules.min_cells_per_perturbation` / `min_perturbations` / `min_genes` | 20 / 5 / 10 | reference |
| `modules.program_correlation` / `module_correlation` / `linkage_method` | pearson / spearman / average | reference |
| `modules.n_programs` / `n_modules` / `cluster_distance_threshold` | 4 / 9 / 0.7 | reference |
| `modules.log2fc_pseudocount` | 1e-9 | reference; 1.0 (Seurat `FoldChange` convention) bounds log2FC for genes undetected in one group (*option*) |
| `modules.min_pct_cells_expressing` | 0 | reference has no detection filter (*extension* when > 0) |
| `modules.de_lfc_threshold` / `de_fdr_alpha` | 0.5 / 0.05 | reference DE gate (Welch t, BH within perturbation over the panel) |
| `modules.program_scoring` / `score_programs` / `draw_networks` | score_genes / true / true | reference `sc.tl.score_genes` (ctrl_size 50, seeded) and the network graphs (needs networkx) |
| `protein.representation` | `protein` | the primary normalized protein matrix (CLR of counts by default); never counts (*extension*) |
| `protein.min_cells_per_target` / `min_control_cells` / `fdr_alpha` | 10 / 10 / 0.05 | same support rule as PS (*extension*) |
| `concordance.min_cells` / `fdr_alpha` | 20 / 0.05 | minimum cells for a within-target Spearman correlation (*extension*) |
| `top_n_report` | 12 | targets whose per-target figures are embedded in the report (all are written to disk) |

## Cell states and perturbation × cluster enrichment (`analysis.clustering`; on by default)

| key | default | rationale |
|---|---|---|
| `enabled` | `true` | part of the standard run; `false` skips clustering and enrichment (the SCP1064 regression configs do, so that the frozen v0.1 reference stays comparable) |
| `key` | `leiden` | obs column for the labels; an existing column is an error unless `overwrite: true` |
| `resolution` | 1.0 | reference default; results are conditional on it |
| `n_iterations` | 2 | reference (igraph flavour, undirected, seeded); `-1` iterates until convergence |
| `neighbors_key` | `rna` | the RNA neighbour graph of the representation stage (no second preprocessing) |
| `enrichment.controls` | `[ntc, other]` | reference: both arms computed and reported |
| `enrichment.primary_control` | `other` | reference default: the arm that drives `significant`, the ranking and the figures (non-targeting cells are few in rare clusters) |
| `enrichment.control_classes` | `[non_targeting]` | control classes forming the `ntc` arm |
| `enrichment.min_cells_per_target` / `min_cells_per_cluster` / `min_reference_cells` | 10 / 20 / 10 | reference (an arm needs ≥ 10 cells in total; pairs with < 10 reference cells in the cluster are `low_power`) |
| `enrichment.min_cells_per_guide` / `guide_concordance` | 5 / true | reference guide-concordance rule, significant pairs only |
| `enrichment.odds_pseudocount` | 0.5 | Haldane-Anscombe odds ratio (reference); the uncorrected sample odds ratio is kept as `sample_odds_ratio` |
| `enrichment.fdr_alpha` | 0.05 | BH within each control arm (reference) |
| `enrichment.stratify_by` | null | any obs column with ≥ 2 levels: CMH across its levels replaces the pooled Fisher p-value (reference); `pval_fisher` is kept |
| `enrichment.n_permutations` | 1000 | omnibus permutation p-value: every target's row resampled from a multinomial with the pooled cluster proportions (reference) |
