# Processed object contract

`petrubseq-protein run` writes one AnnData (`${RESULTS_ROOT}/<dataset>/processed/<name>.h5ad`, i.e. `../results/...` next to the repository by default).
The same description is stored in `uns['petrubseq_protein']['schema']` and
printed in the run report, so an analyst never needs the source code. The
object is identical whatever the input format (dense text, 10x MTX / H5,
h5ad); format-specific facts live only in the provenance records.

Shapes below use `n_cells` (aligned cells), `n_genes`, `n_targeting` (biological
antibodies), `n_antibodies` (targeting + isotype controls) and `n_guides`.
Cells and genes keep the order of the RNA input; cell IDs are the original
identifiers, untouched, except for multi-lane inputs where they are
`<barcode>__<lane id>` (the parts are kept in `obs`).

## `X` and `layers`

| slot | shape | dtype | content | status | provenance |
|---|---|---|---|---|---|
| `X` | n_cells × n_genes | float32 CSR | normalized, log-transformed RNA expression | **normalized** | either the provided matrix preserved as-is (`uns[...]['rna']['method']` starts with `none`) or `normalize_total(target_sum)` + `log1p` of raw counts |
| `layers['counts']` | n_cells × n_genes | int32 CSR | **observed** raw RNA UMI counts | raw | present only when the input was raw counts |
| `layers['reconstructed_counts']` | n_cells × n_genes | int32 CSR | integer counts recovered as `round(expm1(X) * library_size / scale)` | **reconstructed, not observed** | present only for log-normalized input with a per-cell library size (`columns.rna_total_counts`); acceptance stats in `uns[...]['rna']['counts_reconstruction']`. Genes removed upstream are absent, so row sums may be smaller than the original library. |

Exactly one of the two count layers exists; `uns['petrubseq_protein']['rna']['counts_layer']` names it.

## `obs` (per cell)

| column | meaning |
|---|---|
| `condition`, `sample`, `donor`, `batch`, `replicate`, `lane` | design variables mapped from the metadata (`columns.*`, which may name columns of a per-lane sample sheet); absent when `null` |
| `lane_id`, `barcode_original` | multi-lane inputs only: the lane the cell came from and its barcode before the `__<lane>` suffix |
| `moi_provided`, `rna_total_counts_provided` | metadata MOI / library size (`columns.moi`, `columns.rna_total_counts`) |
| extra columns from `columns.keep` | copied unchanged |
| `guides` | the assigned guides, `;`-joined (`''` = none): from the provided assignment or from the guide-calling rule, per `perturbation.assignment.source` |
| `n_guides` | number of assigned guides |
| `targets`, `n_targets` | unique targets among the guides |
| `guide`, `target` | the guide/target for single-guide cells, else `''` |
| `control_class` | class of a single control guide (`non_targeting`, `intergenic`, …), else `''` |
| `perturbation` | `target` for single-guide cells; `multi` / `ambiguous` / `unassigned` otherwise |
| `perturbation_class` | `single_targeting` · `single_control` · `multi_targeting` · `multi_control` · `mixed_control_targeting` · `ambiguous` (guide counts present, dominant rule undecided) · `unassigned` |
| `is_control`, `is_targeting`, `is_single_guide` | convenience booleans |
| `has_guide_counts`, `guide_total_counts`, `n_guides_detected`, `guides_detected`, `guide_top`, `guide_top_count`, `guide_second_count`, `guide_dominant_call` | guide-count diagnostics (only when a guide-count matrix exists): guide UMIs per cell, guides at `>= detection_min_umi` UMIs (`;`-joined list, never collapsed), top / runner-up guide and counts, and the dominant-rule call (guide ID, `ambiguous` or `unassigned`) regardless of which source drives `guides` |
| `total_counts`, `n_genes_by_counts`, `total_counts_mt`, `pct_counts_mt`, `total_counts_ribo`, `pct_counts_ribo` | RNA QC (Scanpy `calculate_qc_metrics` on the count layer) |
| `rna_outlier_low_counts`, `rna_outlier_high_counts`, `rna_outlier_low_genes`, `rna_outlier_high_genes` | MAD flags (`qc.flags.rna_n_mads`) |
| `rna_low_genes`, `rna_high_mt`, `rna_low_counts`, `rna_qc_fail` | threshold flags; `rna_qc_fail` is what `qc.filter` removes when enabled |
| `protein_total_counts`, `protein_total_counts_targeting`, `protein_isotype_counts`, `protein_pct_isotype`, `protein_n_detected` | protein QC from raw ADT counts (normalized-only inputs: detection-based QC on the normalized values) |
| `has_protein_counts` | false for cells absent from the raw ADT matrix (their counts are NaN) |
| `protein_outlier_low_counts`, `protein_outlier_high_counts`, `protein_extreme_counts`, `protein_qc_fail` | MAD / extreme / threshold flags |

Stage E adds, only when `analysis.perturbation_effects.enabled` (docs/PERTURBATION_EFFECTS.md):
`ps_score` (each single-guide targeting cell's PS for its own target, NaN otherwise),
`ps_quadrant` (successful_knockdown / escaper / non_responder / low_signal / not_applicable) and
`lochness_self` (each cell's lochNESS for its own label, including the control label).

Stage F adds, only when `analysis.clustering.enabled`: `obs[analysis.clustering.key]`
(default `leiden`), the categorical Leiden cluster of every cell on the RNA neighbour
graph (numbered states, not cell types; docs/CELL_STATES.md).

## `var` (per gene)

`mt`, `ribo` (gene-class flags), `n_cells_by_counts`, `mean_counts`, `total_counts`, `pct_dropout_by_counts` (Scanpy), `highly_variable`, `means`, `dispersions`, `dispersions_norm` (HVG selection). `varm['PCs']` holds the RNA loadings (zero for non-HVGs). For 10x inputs read with `var_names: name`, the gene symbol is the index (duplicates made unique with `-1`, `-2` …).

## `obsm`

| key | shape | type | content | status |
|---|---|---|---|---|
| `protein` | n_cells × n_targeting | DataFrame (named columns) | **primary** normalized protein values: provided matrix when shipped, else CLR of counts | normalized |
| `protein_counts` | n_cells × n_antibodies | DataFrame | raw ADT UMI counts incl. isotype controls; NaN rows for cells without counts | raw |
| `protein_clr` | n_cells × n_targeting | DataFrame | CLR of raw counts (`protein.clr_axis`); absent when CLR is already the primary matrix | normalized (raw-derived) |
| `protein_<name>` | n_cells × n_targeting | DataFrame | other `protein.extra_representations`, or `protein_provided` when the provided matrix is not the primary one | normalized |
| `guide_counts` | n_cells × n_guides | float32 CSR (sparse) | guide (sgRNA) UMI counts; columns = `uns['guide_features'].index`; zero rows where `obs['has_guide_counts']` is false. Only when a guide-count matrix was supplied, never fabricated | raw |
| `X_pca` | n_cells × rna.n_pcs | float32 | RNA PCA on scaled HVGs | derived |
| `X_umap_rna` | n_cells × 2 | float32 | UMAP of the RNA neighbor graph | derived |
| `X_pca_protein` | n_cells × protein.n_pcs | float32 | PCA of `obsm[embedding_key]` (targeting antibodies only, scaled); NaN rows for cells without protein values | derived |
| `X_umap_protein` | n_cells × 2 | float32 | UMAP of the protein neighbor graph | derived |
| `X_multimodal`, `X_umap_multimodal` | n_cells × (rna + protein PCs), n_cells × 2 | float32 | only when `multimodal.enabled: true` | derived |
| `ps_scores`, `ps_scores_raw` | n_cells × n_scored_targets | DataFrame (float32) | Stage E: PS per target, max-normalized to [0, 1] / clipped to [0, scale_factor]; NaN outside the target's perturbed + control cells | derived |
| `lochness` | n_cells × n_targets | DataFrame (float32) | Stage E: lochNESS per target for every cell (kNN in `X_pca`) | derived |
| `program_activity` | n_cells × n_programs | DataFrame (float32) | Stage E: mean z-scored log expression of each gene program's genes | derived |
| `X_umap_provided` (configurable) | n_cells × 2 | float32 | precomputed embedding from the input (`inputs.embedding`) | input |

`uns['petrubseq_protein']['protein']['embedding_key']` names the matrix that fed `X_pca_protein`.

## `obsp` / neighbors

`obsp['rna_connectivities']`, `obsp['rna_distances']` and `uns['rna']` (Scanpy neighbors, key `rna`); likewise `protein` and `multimodal` when every cell has protein values (otherwise the protein graph is computed on the subset and not stored).

## `uns`

| key | content |
|---|---|
| `protein_features` | DataFrame, one row per antibody (index = canonical name): `role` (`target` / `isotype_control`), `is_isotype`, `isotype_control` (matched control), `in_counts`, `in_normalized`, `name_in_normalized_file`, `in_primary_matrix`, `used_in_embedding`, `used_in_qc`, plus the structured annotation `feature_id`, `antibody_name`, `protein_name`, `gene_symbol`, `clone`, `feature_type`, `isotype`, `annotation_source` (`none`, `input_features`, `feature_table`, or both). An empty string means *not available*; nothing is inferred |
| `guide_features` | only with guide counts: DataFrame, one row per guide (index = guide ID, the `obsm['guide_counts']` columns): `target`, `control_class`, `is_control`, `total_umis`, `n_cells_detected`, `n_cells_dominant` and any 10x feature columns |
| `ps_signatures` | Stage E: long table target, rank, gene, beta, t_score (PS response signatures) |
| `perturbation_effect_matrix` | Stage E: targets × response genes log2FC vs non-targeting controls |
| `gene_programs`, `perturbation_modules` | Stage E: gene → program (`P1..`) and target → module (`M1..`) with sizes / DE counts |
| `protein_effects` | Stage E: target × protein effect, Cohen's d, Mann–Whitney p, BH-FDR, sample / guide sign support |
| `perturbation_summary` | Stage E: integrated per-target table (PS, lochNESS, module, RNA and protein effect magnitudes) |
| `perturbation_cluster_enrichment` | Stage F: target × cluster Fisher enrichment vs controls (counts, fractions, odds ratio, Haldane log2 OR, p, BH-FDR, direction, low_power, guide support, optional CMH columns) |
| `leiden` (= clustering key) | Stage F: scanpy's Leiden parameters |
| `pca`, `pca_protein` | variance ratios, parameters; `pca_protein['loadings']` is an antibody × PC DataFrame |
| `hvg`, `rna`, `protein`, `multimodal`, `umap` | Scanpy bookkeeping |
| `petrubseq_protein` | see below |

`uns['petrubseq_protein']` keys: `version`, `schema_version` (`0.2`), `run_name`; `dataset`; `inputs` (modality audit table `modalities`, input `files`, `protein_embedding_representation`, `isotype_map`, and `provenance`: per-modality format / source / state, the `multiplexed` record with its feature-type or slot split, `lanes` and the `cell_id_transformation`, `guide_counts_available`, `provided_assignments_available`); `warnings` (the run's consolidated warning list); `rna` (detected state, method, `counts_layer`, `counts_source`, reconstruction stats); `protein` (primary/extra methods, embedding key, isotype map, counts, `feature_annotation_sources`); `protein_check` (provided-vs-counts formula verification); `qc` (`rna`/`protein` settings and flag counts, `prefilter` and `filter` cell counts, `filtering_steps` audit table); `embeddings` (HVG/PCA/neighbors/UMAP parameters, multimodal, cross-modality diagnostics); `alignment` (cells per modality, duplicates, kept); `perturbations` (`source`, class counts, coverage QC, `guide_counts_available`, and `assignment`: `provided_assignments_available`, `guide_counts_available`, `configured_source`, `effective_source`, `guide_calling` parameters and per-rule counts, `agreement` between provided and count-derived assignments when both exist); `analysis` (Stage E: `perturbation_effects` with the parameters, counts and status of `ps`, `lochness`, `modules`, `protein`, `concordance`; Stage F: `cell_states` with `clustering` (method, resolution, seed, cluster sizes) and `enrichment` (control, thresholds, BH family, counts, omnibus test)); `provenance` (timestamp, host, conda env, git, SLURM, package versions, inputs); `config` (the resolved run config; list-valued keys such as `lanes` are stored as tables); `schema` (this contract, filled with actual shapes).

## Tables written next to the object

See `docs/OUTPUTS.md` (CSV files under `tables/`, cell-level tables gzipped).
