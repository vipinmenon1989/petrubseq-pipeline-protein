# Perturbation analyses: methods

Enabled with `analysis.perturbation_effects.enabled: true` (the default, as in the
reference pipeline). They run after the representations and the Leiden clustering, in
the reference order: **perturbation strength** → (perturbation × cluster enrichment,
docs/CELL_STATES.md) → **modules / gene programs** → **PS** → **lochNESS** → the
protein extension (**protein effects**, **RNA–protein concordance**). The behavioural
reference is `weili-lab/perturbseq-pipeline` (local checkout audited stage by stage in
docs/reference/REFERENCE_PIPELINE_COMPLETE_AUDIT.md); the numerical parity of every
stage on the same cells and counts is in docs/reference/PARITY_RESULTS.md
(`scripts/reference_parity.py`), and `scripts/compare_reference_perturbation.py`
cross-checks the PS / lochNESS / module mathematics against the reference functions
directly.

## Concepts

| term | meaning |
|---|---|
| perturbation target | the gene the guide acts on (`obs['target']`; the reference calls it `target_gene`) |
| downstream response genes | genes whose expression changes in cells carrying the perturbation |
| gene program (`P1..`) | a cluster of panel genes that move together across perturbations |
| perturbation module (`M1..`) | a cluster of perturbations with similar effect profiles |

Cell groups used by every analysis: **perturbed** = `single_targeting` cells of one
target; **ntc** = `single_control` cells whose `control_class` is in `control_classes`
(default `non_targeting`); **other** = `single_targeting` cells of every other target.
Ambiguous, multi-guide, mixed and unassigned cells enter no group. `target_gene_map`
maps a target label to its RNA gene symbol when they differ (dataset-specific).

## Perturbation strength (`strength`)

Reference `perturbation.py`. For every target whose gene is measured: the gene's own
log-normalized expression in perturbed vs control cells, for each arm in `controls`
(`[ntc, other]`; `primary_control: ntc` drives ranking, hit calls and the figures):
`log2fc = log2((mean(expm1 P) + 0.01) / (mean(expm1 C) + 0.01))` on de-logged means,
`pct_knockdown`, `mean_lognorm_*`, `pct_cells_expressing_*`, two-sided Kolmogorov–Smirnov
(`ks_stat`, `ks_pval`), one-sided Mann–Whitney (`mwu_pval_less`), BH-FDR across the tested
targets per arm (`ks_fdr`, `mwu_fdr`), `is_hit = ks_fdr < 0.05 & log2fc < 0`. Targets are
skipped, with the reason, when the gene is absent from `var`, when fewer than 10 perturbed
cells exist, or when < 1 % of primary-control cells express the gene; an arm with < 10
cells gives NaN statistics. Ranking: hits first, then most negative log2FC.
Tables `tables/perturbation_strength/perturbation_full.csv` (all columns), `perturbation.csv`
(reader view), `skipped.csv`; `uns['perturbation_strength']`. Figures: `perturbation_volcano`,
`perturbation_waterfall`, `perturbation_control_comparison`, `perturbation_<target>` (ECDF,
violin, UMAP; the 12 strongest embedded).

## Co-functional modules and gene programs (`modules`)

Reference `modules.py` (Zhou et al. 2023, Nature 624:154; the reference README calls it
Chen et al.). Perturbations with ≥ 20 single-guide cells (≥ 5 of them). **Gene panel**
(`gene_selection`): `cluster_markers` (reference) = union of the top 100 positive
Wilcoxon markers of each Leiden cluster (`rank_genes_groups`, ranked by log fold change);
`response` = union of each target's top 100 response genes at FDR < 0.05 (option);
`hvg`. **Effect matrix** `E[t, g] = log2((mean_P + 1e-9) / (mean_C + 1e-9))` on de-logged
means; control `ntc` (falls back to `other` = leave-one-target-out targeting cells).
Welch t-test on the log values, BH within each perturbation **over the panel**, gives the
DE mask (|log2FC| > 0.5, FDR < 0.05) used for hub sizes and TF→TF edges. **Programs**:
average-linkage clustering of genes on 1 − Pearson (`n_programs` 4, or a cut at 0.7);
**modules**: perturbations on 1 − Spearman (`n_modules` 9); labels numbered along the
dendrogram. Module × program strength = mean program-gene log2FC over the module's
perturbations; TF hubs (DE genes per perturbation), TF→TF edges (panel genes that are
themselves perturbed), module–module connectivity. Per-cell program activity =
`sc.tl.score_genes` (ctrl_size 50, seeded; `program_scoring: zscore` is the older
deterministic alternative) and its mean per cluster.
Tables: `perturbation_effect_matrix`, `perturbation_de_mask`, `gene_programs`,
`perturbation_modules`, `module_program_strength`, `perturbation_program_effects`
(extension), `program_activity_by_cluster`, `tf_hubs`, `tf_edges`, `module_connectivity`;
`obsm['program_activity']`, `uns[...]`. Figures: `regulome_heatmap`, `module_program_strength`,
`module_program_alluvial`, `module_correlation`, `program_activity_by_cluster`,
`program_<P>_umap`, `module_connectivity`, `module_network`, `tf_hub_network`, plus the
extension `perturbation_program_heatmap`.

## PS — perturbation-response score (`ps`)

Song B. et al., Nat Cell Biol 27, 493 (2025); PS_python (`pertps`) re-implemented
(identical to 4e-16). Per target with ≥ 10 cells and ≥ 10 controls: (1) rank genes
perturbed vs control with a t-test and keep the top 100 (the signature); (2) `beta_g` =
OLS slope on the perturbed indicator; (3) every cell's score = projection of its centred
expression on beta; (4) shift so the control mean is 0, clip to [0, 3], divide by the
maximum → [0, 1]. Targets whose gene is not in `var` are skipped (reference; set
`score_targets_without_gene: true` to score them anyway) and targets expressed in < 1 %
of control cells are skipped. Quadrants against the target's own expression (cut at the
control mean): `successful_knockdown`, `escaper`, `non_responder`, `low_signal`; the
same classification of the control cells gives `pct_controls_called_kd` and
`net_pct_kd`. The summary is sorted by `pct_successful_kd` (reference order).
Extensions kept as extra columns: `auc_vs_control`, `d_vs_control` (in-sample
separation of perturbed from control cells), spread statistics, per-target signatures.
`ps_vs_perturbation.csv` joins the summary to the perturbation-strength log2FC and hit.

**Supervised LDA embedding** (`compute_lda_umap`, reference): the run's HVGs (≤
`lda_max_genes`) → 2000 HVGs → scale (max 10) → PCA (`lda_n_pcs` 40) → linear
discriminant analysis (eigen solver, shrinkage auto) on the scored targets (> 5 cells)
plus the control label → UMAP (30 neighbours, min_dist 0.01, cosine, seed 42) →
`obsm['X_lda_umap']` (NaN outside the trained classes), `obs['lda_label']`.

Outputs: `obs['ps_score']`, `obs['ps_quadrant']`, `obsm['ps_scores']`,
`obsm['ps_scores_raw']`, `uns['ps_signatures']`; tables `ps_targets`, `ps_skipped`,
`ps_signatures`, `ps_vs_perturbation`. Figures: `ps_outcome_by_target`,
`ps_escaper_fraction`, `ps_vs_perturbation_strength`, `ps_quadrant_<target>`,
`ps_lda_overview`, `ps_lda_high_confidence`, `ps_lda_<target>`.

## lochNESS — local neighbourhood enrichment (`lochness`)

Huang X. et al., Nature 623, 772 (2023); pertTF port as in the reference. For every cell
and target g: `lochNESS = local_fraction(g) / overall_fraction(g) − 1`, `local_fraction`
over the cell's k nearest neighbours in `obsm['X_pca']` (first 20 PCs; `X_pca_harmony`
when present; never UMAP), `overall_fraction` over all cells. `k = min(300, n − 1)`
(reference); `max_k_fraction` < 1 caps k on small objects (extension, recorded). Summary
per target: own-cell mean (sort key), mean over all cells, maximum, `pct_cells_enriched`
(% of all cells > 0.5), mean per Leiden cluster (`by_cluster`, `top_cluster`); `self_score`
= each cell's score for its own label including the control label. Extensions (extra
columns, off by default): `n_permutations` > 0 adds a seeded label-permutation null
(z, empirical p, BH-FDR); `mean_lochness_in_control_cells` / `delta_own_vs_control`;
per-sample own-cell means.
Outputs: `obsm['lochness']`, `obs['lochness_self']`; tables `lochness_targets`,
`lochness_by_cluster`, `lochness_skipped`, `lochness_by_sample`. Figures:
`lochness_self_enrichment`, `lochness_distributions`, `lochness_by_cluster`,
`lochness_self_umap`, `lochness_<target>` (symlog map + the target's own cells with the
top cluster highlighted).

## Protein effects (extension)

For each target and each protein in `obsm['protein']` (CLR of ADT counts by default;
never counts): perturbed vs control `effect` = difference of means (on CLR a natural-log
fold change of the geometric-mean-normalized abundance), Cohen's d, two-sided
Mann–Whitney U p-value, BH-FDR over all target × protein tests, per-sample sign
consistency (`n_samples_tested`, `n_samples_same_sign`) and per-guide sign consistency
(`n_guides_tested`, `n_guides_same_sign`; guides with ≥ 5 cells). Tables `protein_effects.csv`,
`protein_effect_matrix.csv`; `uns['protein_effects']`. Normalized-only protein inputs work
unchanged.

## RNA–protein concordance (extension)

* **PS ↔ protein** (cell level, within target): Spearman rho between PS and protein value
  across a target's perturbed cells (≥ `min_cells` 20), BH-FDR over the within-target tests;
  a pooled row is descriptive only — `ps_protein_association.csv`.
* **lochNESS ↔ protein** (target level): per protein, Spearman across targets between
  own-cell mean lochNESS and the protein effect (and |effect|) — `lochness_protein_summary.csv`,
  `lochness_protein_association.csv`.
* **gene program ↔ protein**: cell level over single-guide and control cells (program
  activity vs protein value) — `program_protein_association_cells.csv`; target level
  (perturbation × program effect vs perturbation × protein effect) — `program_protein_association_targets.csv`.
* **Integrated summary** `perturbation_summary.csv` (`uns['perturbation_summary']`).

These are associations. None establishes mediation or causality.

## Report

Sections in the reference order: QC (incl. guide QC) · Clustering · Perturbation strength ·
Perturbation enrichment across clusters · Per-cell perturbation response (PS, LDA) ·
lochNESS · Co-functional modules & gene programs; then the extension: Protein QC ·
Protein effects · RNA–protein concordance · Multimodal summaries · Outputs. Per-target
figures beyond `top_n_report` (12) are written to disk and listed in the report.
