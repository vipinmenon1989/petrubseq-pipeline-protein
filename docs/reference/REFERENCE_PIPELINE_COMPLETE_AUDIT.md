# Complete audit of the reference Perturb-seq pipeline

Reference: `/local/projects-t3/lilab/vmenon/PertTF-Virtual-Challeng-Weilab/perturbseq-pipeline`
(commit `b486d65`, package `perturbseq_pipeline` 0.1.0) plus the `pertps` 0.1.0 package it
delegates the PS score to (`/autofs/projects-t3/lilab/vmenon/Validation_Pert-TF/PertTF_Combined_Python_Analysis/pertps_project`).
Audited 2026-09-25 from the code, not from figure names. Every statement below cites the
file and function it was read from. The "current" pipeline is `petrubseq_protein` in this
repository as it stood before the corrections of the same date (commit `382a9f7`).

Scope note: Energy Distance, perturbation distance and phenotype-space modules are not part
of the reference and are OUT_OF_SCOPE here.

## A. Reference repository architecture

Source (`src/perturbseq_pipeline/`, 9,900 lines):

| file | lines | role |
|---|---|---|
| `cli.py` | 484 | `run_pipeline` driver (stage order below), argparse CLI |
| `config.py` | 886 | dataclass schema, defaults, validation, YAML round trip |
| `io.py` | 1002 | 10x MTX / h5ad loading, feature splitting, guide tables, h5ad writing, archive |
| `guides.py` | 404 | target parsing, NTC detection, dominant-guide assignment, summaries |
| `guide_counting.py`, `guide_design.py` | 656, 179 | streaming guide counting from FASTQ (not exercised by the analyses audited here) |
| `qc.py` | 281 | gene classes, QC metrics, prefilter, strict filter, QC tables, guide QC warnings |
| `cluster.py` | 235 | normalize, HVG/PCA/(Harmony)/neighbors/UMAP/Leiden, cluster summary |
| `perturbation.py` | 368 | perturbation strength (target's own expression, KS/MWU, BH, hit call) |
| `enrichment.py` | 628 | perturbation x cluster Fisher/CMH enrichment, omnibus, concordance |
| `modules.py` | 522 | perturbation x gene log2FC matrix, gene programs, co-functional modules, TF network |
| `ps_score.py` | 563 | adapter around `pertps` (PS score), quadrants, LDA embedding |
| `lochness.py` | 315 | lochNESS neighbourhood enrichment |
| `plots.py` | 2146 | every figure, through `FigureRegistry` |
| `report.py` + `templates/report.html` | 342 + 492 | jinja2 report |
| `tests/` | 2,800 | pipeline, modules, guide-counting tests; synthetic generators |
| `config/{default,demo,context_a,full_4lane}.yaml`, `docs/methods.md`, `README.md` | | documentation and worked configs |

Vocabulary: `obs['target_gene']`, `obs['guide_id']`, `obs['perturbation_class']` in
{`targeting`, `non-targeting`, `ambiguous`, `unassigned`}, `obs['lane_id']`,
`layers['counts']`, `layers['lognorm']` (= `X`), `obsm['X_pca']`, optional `obsm['X_pca_harmony']`,
`obsm['X_umap']`, `obs['leiden']`.

## B. Reference execution graph (actual order in `cli.run_pipeline`)

```
load (io.load_data)                                   stage 1
  -> qc.prefilter (min_genes 200, min_cells 3)         stage 2
  -> qc.compute_qc_metrics (mt/ribo/hb; log1p=True)
  -> plot_qc("before filtering")
  -> qc.filter_cells_and_genes (min_genes_final 1000, [min_counts], pct_mt < 20, [pct_hb], min_cells 3 again)
  -> plot_qc("after filtering"); tables qc_steps, qc_summary
  -> guides.assign_guides (dominant rule)              stage 3
  -> guide QC tables/warnings; plot_guide_qc
  -> cluster.normalize (normalize_total(None=median) + log1p; layers counts/lognorm)   stage 4
  -> [cluster.assigned_only: embed all cells, write *_all_cells.h5ad, subset to singlets, reset_embedding]
  -> cluster.embed_and_cluster (HVG 3000 -> scale HVG copy (max 10) -> PCA 50 arpack -> [Harmony] -> neighbors 15 -> UMAP 0.5 -> Leiden 1.0 igraph 2 it)
  -> tables clusters; plot_clustering
  -> perturbation.test_all_targets                     stage 5   (perturbation strength)
  -> plot_perturbation_overview; plot_per_target
  -> enrichment.test_cluster_enrichment (if enabled)   stage 6
  -> plot_enrichment; plot_enrichment_per_target
  -> modules.compute_modules (if enabled)              stage 7   (needs obs['leiden'] for the gene panel)
  -> plot_modules
  -> ps_score.compute_ps_scores; attach_scores; LDA embedding; ps_vs_perturbation   stage 8
  -> plot_ps_scores; plot_ps_lda
  -> lochness.compute_lochness; attach_scores          stage 9   (needs obs['leiden'] for by_cluster)
  -> plot_lochness
  -> write tables, figure manifest, processed h5ad (+guide obsm), guide table   stage 10
  -> report.build_report; archive                      stage 11
```

Differences from the diagram in the task statement: guide assignment happens *after* QC
filtering (QC never sees guide labels); modules run *before* PS and lochNESS; clustering
runs before every analysis, and both modules (gene panel) and lochNESS (cluster summary)
consume `obs['leiden']`.

## C. Stage detail

### C.1 Input (`io.py`)

* `mtx` mode: `sc.read_10x_mtx(var_names='gene_symbols', gex_only=False)` per lane; several
  lanes concatenated with `index_unique='-'` (`<barcode>-<lane>`) or `<lane>_<barcode>` with
  `cell_id_format: prefix`; `var_names_make_unique()`. `split_features`: Gene Expression ->
  `expr`; `guide_feature_types` (`Custom`, `CRISPR Guide Capture`) -> `guides` re-indexed by
  `gene_ids`; any other feature type (e.g. Antibody Capture) is dropped. `layers['counts'] = X`.
* `h5ad` mode: guides from `var` feature types, a companion guide h5ad, a guide table, or a
  per-cell label column (`input.guide_obs_column`); `input.counts_layer` / `normalized_layer`
  honoured (`apply_layer_choices`).
* Lane id in `obs['lane_id']`; multi-lane runs require a metadata file.

### C.2 QC (`qc.py`)

| step | function | rule | notes |
|---|---|---|---|
| gene classes | `annotate_gene_classes` | `var.mt`: name starts with `MT-` (case-insensitive); `ribo`: `RPS`/`RPL`; `hb`: regex `^HB[^(P)]` | warning when no mt gene matches |
| prefilter | `prefilter` | `sc.pp.filter_cells(min_genes=200)`, `sc.pp.filter_genes(min_cells=3)` | before metrics, so figures show the tail |
| metrics | `compute_qc_metrics` | `sc.pp.calculate_qc_metrics(qc_vars=[mt,ribo,hb], log1p=True, percent_top=None)` on raw counts in `X` | `n_genes_by_counts`, `total_counts`, `pct_counts_mt/ribo/hb` |
| strict filter | `filter_cells_and_genes` | (1) `filter_cells(min_genes=min_genes_final=1000)`; (2) `filter_cells(min_counts)` if set (default None); (3) keep `pct_counts_mt < 20.0` (strict less-than); (4) keep `pct_counts_hb < max_pct_hb` if set (default None); (5) `filter_genes(min_cells=3)` **again** after the cell filters | every step recorded in `qc_steps` (step, threshold, cells/genes before/after/removed); error if 0 cells remain |
| summary | `qc_summary_table` | per lane + ALL: n_cells, median genes / UMIs / % mt / % ribo | rounded to 2 dp |
| guide QC | `guide_qc_summary`, `check_guide_qc` | class counts, median guide UMIs, guides detected (> detection_threshold 3), MOI, top:second ratio; warnings at < 50 % assigned, 0 NTC, > 30 % ambiguous, > 30 % multi-guide | |

Figures (`plots.plot_qc`, called before and after filtering): `qc_violin_<stage>` (violin per
metric, split by lane when > 1 lane; `total_counts` log axis), `qc_scatter_<stage>` (UMIs vs
genes, coloured by % mt), `qc_cells_per_lane_<stage>` (multi-lane only).

### C.3 Guide assignment (`guides.py`) — adapter stage, audited for the control definition

`top_two_guides` per cell; assigned when `top >= max(min_umi,1)` and `top > dominance_ratio *
second` (and `second <= max_second_umi` when that gate is on). `has_counts & ~assigned` ->
`ambiguous`; no counts -> `unassigned`. Target = `target_feature_column` value, else
`target_regex` group 1, else first field of the id split on `_ - .`. NTC = target label matches
`ntc_patterns` (case-insensitive: `^non[-_.]?targeting`, `^non$`, `^ntc`, `scramble`,
`^safe[-_.]?harbor`); all NTC labels collapse to `non-targeting`. `perturbation_class` follows.
Guide-QC figures (`plot_guide_qc`): `guide_umi_depth`, `guide_multiplicity`, `guide_top_vs_second`
(coloured by class, dominance line), `guide_assignment_classes`, `guide_assignment_per_lane`
(multi-lane), `target_representation` (cells per target, NTC green, `min_cells_per_target` line),
`guide_representation` (cells per guide, ranked, symlog).

### C.4 Normalization, embedding, clustering (`cluster.py`)

* `normalize`: `layers['counts']` kept; `sc.pp.normalize_total(target_sum=None)` (median library
  size) + `sc.pp.log1p`; `layers['lognorm'] = X`. Skipped when a `lognorm` layer already exists.
* `embed_and_cluster`: `sc.pp.highly_variable_genes(n_top_genes=min(3000, n_vars))` (flavor
  `seurat` default, on log data); optional `regress_out`; `n_pcs = min(50, n_obs-1, n_vars-1)`;
  when `scale_max_value` (10) is set, a **copy of the HVG block** is scaled (`sc.pp.scale(max_value=10)`)
  and `sc.tl.pca(n_comps, svd_solver='arpack', random_state=seed)` runs on it; `X_pca` and `uns['pca']`
  copied back, `X` untouched; optional Harmony (`batch_key`; `harmonypy.run_harmony` called
  directly) -> `X_pca_harmony`; `sc.pp.neighbors(n_neighbors=15, n_pcs, use_rep, random_state=seed)`;
  `sc.tl.umap(min_dist=0.5, random_state=seed)`; `sc.tl.leiden(key_added='leiden', resolution=1.0,
  flavor='igraph', n_iterations=2, directed=False, random_state=seed)`.
* Population: every QC-passing cell (default `assigned_only: false`). With `assigned_only` the
  all-cell embedding is written separately and the singlets (`targeting` + `non-targeting`) are
  re-embedded from scratch (`reset_embedding`).
* `cluster_summary`: per cluster n_cells, pct_of_total, median genes, median % mt, top lane share.
* Figures (`plot_clustering`): `pca_variance` (scree, log y), `umap_clusters`, `umap_qc_metrics`
  (genes, UMIs, % mt), `umap_lane` (multi-lane), `umap_assignment_class`, `umap_target_gene`
  (legend if <= 25 targets), `cluster_lane_composition` (multi-lane).

### C.5 Perturbation strength (`perturbation.py`) — traced exactly

* Perturbed cells of gene g: `perturbation_class == 'targeting' & target_gene == g`.
* Controls, both computed (`controls: [ntc, other]`, `primary_control: ntc`):
  `ntc` = `class == 'non-targeting'`; `other` = `class == 'targeting' & target_gene != g`
  (ambiguous/unassigned never used). An arm is dropped when it has < `min_control_cells` (10)
  cells; if the primary is unavailable the first usable arm becomes primary.
* Expression: `layers['lognorm']` column of gene g (never scaled values).
* Eligibility per target (in this order, each recorded in `skipped` with its reason): gene not in
  `var_names`; `n_perturbed < min_cells_per_target` (10); percent of **primary-control** cells with
  non-zero expression `< min_pct_expressing_control` (1.0 %). Per arm: `n_control < 10` -> NaN
  statistics for that arm.
* `compare_groups`: `mean_lognorm_*` (mean of log values); `log2fc = log2((mean(expm1 P)+0.01) /
  (mean(expm1 C)+0.01))` (pseudocount 0.01 on de-logged means); `pct_knockdown = 100(1 - ratio)`;
  `pct_cells_expressing_*`; `ks_2samp` two-sided (`ks_stat`, `ks_pval`); `mannwhitneyu(alternative='less')`
  (`mwu_pval_less`). Means, not medians.
* FDR: `benjamini_hochberg` (NaN-tolerant) separately per arm over the tested targets
  (`ks_fdr_<arm>`, `mwu_fdr_<arm>`).
* Hit: `is_hit_<arm> = ks_fdr_<arm> < 0.05 & log2fc_<arm> < 0.0` (KS drives the call; MWU is
  reported only).
* Ranking: sort by (`is_hit_<primary>` desc, `log2fc_<primary>` asc); `rank` column 1..n.
* Multiple guides / constructs: pooled per target (guide identity is not used here); lanes not used.
* Outputs: `tables/perturbation.csv` (reader view: Rank, Target, Perturbed cells, Control cells,
  log2FC, % knockdown, KS stat, KS FDR, Effective, + log2FC/KS FDR of the other arm),
  `perturbation_full.csv` (all columns), `skipped.csv`; report cards (targets tested, effective,
  control cells, not testable).
* Figures (`plot_perturbation_overview`, `plot_per_target`):
  * `perturbation_volcano`: x = `log2fc_<primary>`, y = -log10(`ks_fdr_<primary>`), hits red,
    dashed lines at FDR 0.05 and log2FC 0, the 12 most negative log2FC annotated.
  * `perturbation_waterfall`: bars of `log2fc_<primary>` sorted ascending, hits red, title with
    n effective / n tested.
  * `perturbation_control_comparison` (both arms): scatter log2FC ntc vs other, diagonal, Pearson r.
  * `perturbation_<gene>` (one per tested target, `in_report` for the top 12 by rank): ECDF of
    the gene's log-normalized expression in perturbed vs each control arm (n in legend, title
    log2FC/FDR), violin by group, UMAP with 10 % background and perturbed cells coloured by
    expression (Reds, outlined).

### C.6 Cluster enrichment (`enrichment.py`)

* Clusters tested: `obs['leiden']` levels with >= `min_cells_per_cluster` (20) cells, ordered
  numerically. Targets tested: targeting cells with >= `min_cells_per_target` (10); others in `skipped`.
* Reference arms `controls: [ntc, other]`; an arm is usable when its **total** cells >=
  `min_reference_cells` (10); `primary_control: other` (fallback to the first usable arm).
  `ntc` = non-targeting cells; `other` = targeting cells of other targets.
* Composition: `composition[gene, cluster]` = % of the target's cells in the cluster;
  `reference_composition[arm]` = % of all cells of that arm (for `other`: all targeting cells).
* Omnibus (`omnibus_test`): targeting-cell target x cluster table, zero margins pruned,
  `chi2_contingency`; % expected < 5; permutation p = (n_ge+1)/(N+1) with N = 1000 draws where
  each target's row is resampled from a multinomial with the column proportions (row margins
  fixed, columns random), seeded with `run.seed`.
* Per (target, arm, cluster): a,b,c,d counts; `pct_of_target`, `pct_of_reference`;
  `odds_ratio` = Haldane-Anscombe ((a+.5)(d+.5))/((b+.5)(c+.5)); `log2_odds_ratio`; `direction`
  enriched when `pct_t > pct_r` else depleted; `pval` from two-sided `fisher_exact`, or, when
  `stratify_by` is set and has >= 2 levels, the CMH p-value from `statsmodels StratifiedTable`
  over the strata (and the pooled OR replaces the Haldane OR when finite and > 0);
  `low_power = c < min_reference_cells`.
* FDR: BH **within each arm** (`fdr`); `significant = fdr < 0.05 & control == primary`.
* Guide concordance (`_guide_concordance`), significant pairs only: guides of the target with
  >= `min_cells_per_guide` (5) cells; a guide agrees when its in-cluster fraction lies on the
  same side of `pct_of_reference` as the direction -> `guides_concordant`, `guides_tested`.
* `effect_magnitude`: per target `composition_shift_pct` = total variation distance
  (sum |comp - ref| / 2) from the primary reference composition, `n_significant_clusters`.
* Table sorted by (`significant` desc, `log2_odds_ratio` desc). `format_enrichment_table`:
  primary arm rows with fdr < 1 sorted by |log2 OR|.
* Derived: `enrichment_matrix` (target x cluster log2 OR, primary arm), `significance_matrix`,
  `phenocopy_similarity` (Pearson correlation of composition profiles between targets).
* Figures: `enrichment_heatmap` (rows ordered by average-linkage correlation clustering of the
  log2 OR profiles, `*` at FDR < 0.05, symmetric 98th-percentile colour limit),
  `enrichment_phenocopy` (target x target Pearson r, similarity-ordered, >= 3 targets),
  `enrichment_composition` (stacked bars: reference then the top 30 targets by shift),
  `enrichment_volcano` (all primary-arm pairs, top 10 |log2 OR| significant pairs labelled),
  `enrichment_effect_magnitude` (bars of composition shift, red when a significant cluster
  exists), `enrichment_<gene>` per target (composition bars vs reference with `*`, log2 OR bars
  red when significant, UMAP of the target's cells coloured by cluster; report shows the 12
  strongest with hits first).
* Report context: n hits, targets with hits, n tests, control labels, chi2/dof/permutation p,
  % small expected, n low-power, top shift.

### C.7 Modules and gene programs (`modules.py`) — traced exactly

* Perturbations: targeting targets with >= `min_cells_per_perturbation` (20) cells, sorted; the
  stage is skipped (returns None) below `min_perturbations` (5).
* Gene panel (`select_genes`, `gene_selection: cluster_markers` default): `sc.tl.rank_genes_groups(
  groupby='leiden', method='wilcoxon', n_genes=100)` on `X` (log-normalized), keep rows with
  `logfoldchanges > 0`, per cluster the top 100 by `logfoldchanges`, union in first-seen order
  (clusters in group order), restricted to `var_names`. Falls back to HVGs when `leiden` has < 2
  levels; `gene_selection: hvg` uses the HVG flags directly. No detection filter. Skip when
  fewer than `min_genes` (10).
* Effect matrix (`build_effect_matrix`), control `ntc` (falls back to `other` when there are no
  NTC cells): perturbation indicator over cells with `target_gene == t`; `lin = expm1(lognorm)`;
  `mean_perturbed = mean(lin)` per target; `mean_control` = mean over NTC cells (or leave-one-
  target-out mean over all targeting cells for `other`); `log2fc = log2((mean_p + 1e-9) /
  (mean_c + 1e-9))` (**pseudocount 1e-9**). Welch t-test per (perturbation, gene) on the
  log values (unbiased variances, Welch df), two-sided p, BH within each perturbation; DE mask
  = `|log2fc| > 0.5 & fdr < 0.05`.
* `cluster_axis` for both axes: `1 - corr` (Pearson for genes = programs, Spearman for
  perturbations = modules), NaN -> 1, symmetrised, diagonal 0, negatives clipped, condensed,
  `scipy linkage(method='average')`; `fcluster(maxclust=k)` with k = `n_programs` 4 / `n_modules`
  9 (clamped to [2, n]) or `criterion='distance'` at `cluster_distance_threshold` 0.7 when the
  count is null; labels renumbered `P1..`/`M1..` in dendrogram-leaf order; leaf order returned.
* `module_program_strength`: mean log2FC of the program's genes over the module's perturbations.
* `tf_network`: hub size = number of DE genes per perturbation (`hubs`, sorted); TF->TF edges
  for panel genes that are themselves perturbations when the DE mask is true (`log2fc`, sign,
  modules); module x module connectivity = edges between modules / (size_i * size_j).
* Program scoring (`score_programs: true`): `sc.tl.score_genes(expr, genes, score_name='program_P<k>_score', ctrl_size=50)`
  per program -> obs columns; `program_activity` = mean score per cluster (programs x clusters).
* Tables: `effect_matrix`, `gene_programs` (gene, program, program_size), `cofunctional_modules`
  (target_gene, module, n_cells, n_de_genes), `module_program_strength`,
  `program_activity_by_cluster`, `tf_hubs`, `tf_edges`, `module_connectivity`.
* Figures (`plot_modules`): `regulome_heatmap` (genes x perturbations log2FC, both axes grouped
  by cluster in leaf order, colour bars for programs/modules, 98th-percentile limit),
  `module_program_strength` (annotated heatmap), `module_program_alluvial`,
  `module_correlation` (perturbation x perturbation Spearman, module blocks), `program_activity_by_cluster`,
  `program_<P>_umap` (per program score on the UMAP, report top 12), `module_connectivity`
  (magma heatmap), `module_network` and `tf_hub_network` (networkx spring layouts, top 40 hubs).

### C.8 PS score (`ps_score.py` + `pertps.analyzer.PerturbAnalyzer`) — traced exactly

* Input: `layers['lognorm']` copied into a working AnnData; `obs['gene']` = `target_gene` for
  targeting cells, `Non-Targeting` for NTC, `Other` for ambiguous/unassigned.
* Skipped when NTC cells < `min_control_cells` (10). Candidates: targeting targets with >=
  `min_cells_per_target` (10), sorted. Per target, skipped (with reason) when: gene not in
  `var_names`; % of NTC cells expressing the gene < `perturbation.min_pct_expressing_control`
  (1 %); `pertps` raises or returns None/empty; no perturbed cell carries a score.
* `calculate_ps_score(target, top_n=100)`: subset to `gene in {target, Non-Targeting}`; None if
  < 10 cells; `sc.tl.rank_genes_groups(groupby='gene', reference='Non-Targeting', method='t-test')`,
  top 100 names of the target group; `Y` = those genes (dense); `x` = 1 for target cells;
  `beta = (x - mean x) . (Y - mean Y) / var_x` (OLS slope = mean difference); score =
  `(Y - mean Y) . beta / (beta . beta)`; shift so the NTC mean is 0; clip to `[0, scale_factor=3]`;
  divide by the maximum (-> [0, 1]); Series indexed by the subset cells. Deterministic.
* Quadrants (`_summarize_target`): expression cut = `mean` (default; `median` / `quantile`
  options) of the NTC cells' log-normalized expression of the target gene (all cells if no
  control); `high_ps = score >= ps_threshold (0.5)`, `high_expr = expr > cut`;
  KD = high_ps & ~high_expr, escaper = high_ps & high_expr, non-responder = ~high_ps & high_expr,
  low signal otherwise. The same classification is applied to the control cells
  (`pct_controls_called_kd`) and `net_pct_kd` = pct KD - control pct KD.
* Summary row: target_gene, n_perturbed_cells, n_control_cells, mean_ps, median_ps, pct_high_ps,
  pct_successful_kd, pct_escaper, pct_non_responder, pct_low_signal, pct_controls_called_kd,
  net_pct_kd, expression_cut, expression_cut_method. **Sorted by `pct_successful_kd` desc**;
  `top_targets(n)` = head of that order.
* `attach_scores`: obs `ps_<gene>` (NaN outside the comparison), `ps_quadrant_<gene>`, `ps_score`
  (own-target score), `ps_quadrant`.
* `compare_with_perturbation_strength`: summary joined to `log2fc_<primary>` / `is_hit_<primary>`
  of the perturbation-strength table (`ps_vs_perturbation` table).
* LDA embedding (`_compute_lda_embedding`, `compute_lda_umap: true`): genes capped at
  `lda_max_genes` 5000 using the pipeline's HVG flags (i.e. the 3000 HVGs); then
  `pertps.compute_lda_umap(targets, n_pcs=40)`: copy, (normalize only if `X.max() > 20`),
  `highly_variable_genes(n_top_genes=2000)`, `scale(max_value=10)`, `pca(n_comps=40)`; labels
  `NT` for controls, target for scored targets with > 5 cells, else `Other`; subset to targets +
  NT; `LinearDiscriminantAnalysis(solver='eigen', shrinkage='auto').fit_transform(X_pca, label)`;
  `umap.UMAP(n_neighbors=30, min_dist=0.01, metric='cosine', random_state=42)`; NaN for
  unplaced cells -> `obsm['X_lda_umap']`, `obs['lda_label']`.
* Figures: `ps_outcome_by_target` (stacked quadrant % per target in summary order),
  `ps_escaper_fraction` (bars sorted by % escaper), `ps_vs_perturbation_strength` (log2FC of the
  group test vs % confirmed KD, hits red, Pearson r), `ps_quadrant_<gene>` per scored target
  (control cells grey subsampled to 2000, target cells coloured by quadrant, cuts, quadrant
  percentages; report top 12), `ps_lda_overview` (LDA-UMAP by label), `ps_lda_high_confidence`
  (cells with own score >= 0.8), `ps_lda_<gene>` (target cells shaded by score; report top 12).

### C.9 lochNESS (`lochness.py`) — traced exactly

* Representation: `use_rep` null -> `X_pca_harmony` when present else `X_pca`; `n_pcs` 20;
  `k = min(300, n_obs - 1)` (warning when reduced); `sc.pp.neighbors(n_neighbors=k, n_pcs,
  use_rep, key_added='lochness_nn', random_state=seed)`; adjacency = binarised
  `obsp['lochness_nn_distances']` (self excluded; k-1 stored neighbours); denominator = actual
  neighbour count per cell.
* Labels: `obs[genotype_key='target_gene']` over **all cells** (ambiguous/unassigned/NTC included);
  `overall = value_counts(normalize=True)`. Candidates: targeting targets with >=
  `min_cells_per_target` (10); others `skipped`.
* Score per target: `local = adj @ indicator(label == g)`; `lochness = (local / n_neighbours) /
  overall[g] - 1`. Optional Gaussian noise `noise_delta` (0 = off). No permutation null, no
  significance, no k cap other than n-1.
* Summary per target: n_cells, overall_fraction_pct, mean_lochness_all_cells,
  mean_lochness_in_own_cells, max_lochness, pct_cells_enriched (% of all cells > `enrichment_cut`
  0.5), top_cluster / top_cluster_mean (mean per `leiden` cluster); sorted by
  `mean_lochness_in_own_cells` desc. `by_cluster` table (targets x clusters mean score).
* `self_score`: each cell's score for its own label; NTC cells scored against the NTC label.
* `attach_scores`: obs `lochness_<gene>`, `lochness_self`.
* Figures: `lochness_self_enrichment` (bars of own-cell mean, red above the cut, cut line),
  `lochness_distributions` (violins over all cells for the top 30 targets),
  `lochness_by_cluster` (targets x clusters mean, similarity-ordered rows), `lochness_self_umap`
  (self score on the UMAP), `lochness_<gene>` per target (symlog-coloured score map + the
  target's own cells with cluster labels, top cluster highlighted; report top 12).

### C.10 Outputs and report

Tables written as `tables/<name>.csv` (index dropped) for every non-empty table plus
`figure_manifest.csv`; processed h5ad with `obsm['guide_counts']` + `uns['guide_names']`;
optional guide table. Report sections in order: Run summary (cards, warnings), 1 QC (filtering
steps, per-lane summary, QC figures, guide QC table + figures), 2 Clustering (cluster table,
embedding figures), 3 Perturbation strength (cards, overview figures, ranked table, skipped
table, per-target figures, extras list), 4 Enrichment (omnibus, overview figures, significant
pairs, per-target figures), 5 Per-cell response (PS figures, table, quadrants, LDA), 6 lochNESS
(figures, table, per-target maps), 7 Modules & programs (figures, modules, programs,
module x program, hubs), Outputs (deliverables, manifest, config, versions).

## D. Reference figure inventory versus the current pipeline (state before the corrections)

Statuses: MATCH, CURRENT_FIGURE_DIFFERS, MISSING_CURRENT_FIGURE, CURRENT_ANALYSIS_DIFFERS,
NOT_APPLICABLE_MULTIMODAL, OUT_OF_SCOPE.

| reference section | reference figure | plotting function | source table / data | analysis function | purpose | current equivalent (before) | status before |
|---|---|---|---|---|---|---|---|
| 1 QC | qc_violin_before/after | `plot_qc` | obs QC metrics | `qc.compute_qc_metrics` | per-cell metric distributions (per lane) | `rna_qc_distributions_*` (histograms) + `rna_*_by_sample_*` violins | CURRENT_FIGURE_DIFFERS (same information) |
| 1 QC | qc_scatter_before/after | `plot_qc` | obs | same | UMIs vs genes coloured by % mt | `rna_counts_vs_genes_*` | MATCH |
| 1 QC | qc_cells_per_lane_* | `plot_qc` | obs lane | same | lane balance | `rna_cells_per_<group>_*` | MATCH |
| 1.4 Guide QC | guide_umi_depth | `plot_guide_qc` | obs total_guide_counts | `assign_guides` | guide capture depth | `guide_count_diagnostics` panel 1 | MATCH |
| 1.4 | guide_multiplicity | `plot_guide_qc` | n_guides_detected | same | MOI | `guide_count_diagnostics` panel 2 + `guides_per_cell` | MATCH |
| 1.4 | guide_top_vs_second | `plot_guide_qc` | top/second counts | same | dominance rule | `guide_count_diagnostics` panel 3 (no class colours, no dominance line) | CURRENT_FIGURE_DIFFERS |
| 1.4 | guide_assignment_classes | `plot_guide_qc` | perturbation_class | same | assignment outcome | `perturbation_class_composition` | MATCH |
| 1.4 | guide_assignment_per_lane | `plot_guide_qc` | class x lane | same | lane failure | `perturbation_class_by_condition` (condition, not lane) | CURRENT_FIGURE_DIFFERS |
| 1.4 | target_representation | `_plot_representation` | cells per target | same | library balance, min-cells line | `cells_per_perturbation` (top 60, no threshold line) | CURRENT_FIGURE_DIFFERS |
| 1.4 | guide_representation | `_plot_representation` | `guide_representation` | same | cells per guide ranked | `coverage_histograms` (histogram, not ranked curve) | CURRENT_FIGURE_DIFFERS |
| 2 Clustering | pca_variance | `plot_clustering` | uns pca | `embed_and_cluster` | scree | `rna_pca_variance` | MATCH |
| 2 | umap_clusters | `plot_clustering` | obs leiden | same | clusters | `umap_leiden` (cell_states) | MATCH |
| 2 | umap_qc_metrics | `plot_clustering` | obs | same | technical clusters | `rna_umap_pct_counts_mt`, `rna_umap_total_counts` (no genes panel) | CURRENT_FIGURE_DIFFERS |
| 2 | umap_lane | `plot_clustering` | obs lane | same | batch | `rna_umap_<sample>` when in color_by | CURRENT_FIGURE_DIFFERS |
| 2 | umap_assignment_class | `plot_clustering` | obs class | same | unassigned islands | `rna_umap_perturbation_class` | MATCH |
| 2 | umap_target_gene | `plot_clustering` | obs target | same | perturbation islands | none | MISSING_CURRENT_FIGURE |
| 2 | cluster_lane_composition | `plot_clustering` | leiden x lane | same | batch clusters | `cluster_composition` (class + one design column) | MATCH |
| 3 Perturbation strength | perturbation_volcano | `plot_perturbation_overview` | `perturbation_full` | `test_all_targets` | knockdown significance | none | MISSING_CURRENT_FIGURE (analysis missing) |
| 3 | perturbation_waterfall | same | same | same | knockdown ranking | none | MISSING_CURRENT_FIGURE |
| 3 | perturbation_control_comparison | same | same | same | control agreement | none | MISSING_CURRENT_FIGURE |
| 3 | perturbation_<gene> (per target) | `plot_per_target` | lognorm, table | same | per-target ECDF/violin/UMAP | none | MISSING_CURRENT_FIGURE |
| 4 Enrichment | enrichment_heatmap | `plot_enrichment` | `enrichment_matrix` | `test_cluster_enrichment` | target x cluster log2 OR | `perturbation_cluster_enrichment` (rows by min FDR, single arm, default NTC reference) | CURRENT_ANALYSIS_DIFFERS |
| 4 | enrichment_phenocopy | same | `phenocopy_similarity` | same | targets that phenocopy | none | MISSING_CURRENT_FIGURE |
| 4 | enrichment_composition | same | `composition`, `effect_magnitude` | same | composition vs reference | none | MISSING_CURRENT_FIGURE |
| 4 | enrichment_volcano | same | table | same | all pairs | none | MISSING_CURRENT_FIGURE |
| 4 | enrichment_effect_magnitude | same | `effect_magnitude` | same | composition shift ranking | none | MISSING_CURRENT_FIGURE |
| 4 | enrichment_<gene> | `plot_enrichment_per_target` | table, composition, UMAP | same | per-target composition | none | MISSING_CURRENT_FIGURE |
| 5 PS | ps_outcome_by_target | `plot_ps_scores` | ps summary | `compute_ps_scores` | quadrant % per target | `ps_target_median` (median PS bars) | CURRENT_FIGURE_DIFFERS |
| 5 | ps_escaper_fraction | same | summary | same | escapers | none | MISSING_CURRENT_FIGURE |
| 5 | ps_vs_perturbation_strength | same | `ps_vs_perturbation` | same + strength | method check | none | MISSING_CURRENT_FIGURE |
| 5 | ps_quadrant_<gene> | `_plot_ps_quadrants` | scores, lognorm | same | per-cell outcome | `ps_distributions` (violins, no quadrants) | CURRENT_FIGURE_DIFFERS |
| 5.4 | ps_lda_overview | `plot_ps_lda` | X_lda_umap | `_compute_lda_embedding` | supervised embedding | none | MISSING_CURRENT_FIGURE |
| 5.4 | ps_lda_high_confidence | same | same | same | strong responders | none | MISSING_CURRENT_FIGURE |
| 5.4 | ps_lda_<gene> | same | same | same | per-target LDA map | none | MISSING_CURRENT_FIGURE |
| 6 lochNESS | lochness_self_enrichment | `plot_lochness` | summary | `compute_lochness` | self-clustering ranking | `lochness_targets` (with permutation null, k capped) | CURRENT_ANALYSIS_DIFFERS |
| 6 | lochness_distributions | same | scores | same | tails | none | MISSING_CURRENT_FIGURE |
| 6 | lochness_by_cluster | same | `by_cluster` | same | cluster counterpart | none | MISSING_CURRENT_FIGURE |
| 6 | lochness_self_umap | same | self_score | same | where identity organises | none | MISSING_CURRENT_FIGURE |
| 6 | lochness_<gene> | same | scores | same | per-target map | `lochness_umap` (6 panels, linear colour, no own-cell panel) | CURRENT_FIGURE_DIFFERS |
| 7 Modules | regulome_heatmap | `_plot_effect_heatmap` | effect matrix | `compute_modules` | modules x programs | `effect_matrix_programs` (response-gene panel, pseudocount 1) | CURRENT_ANALYSIS_DIFFERS |
| 7 | module_program_strength | `_plot_module_program` | `module_program` | same | signed strength | `module_program_strength` | CURRENT_ANALYSIS_DIFFERS (same figure, different matrix) |
| 7 | module_program_alluvial | `_plot_alluvial` | same | same | flow | none | MISSING_CURRENT_FIGURE |
| 7 | module_correlation | `_plot_module_correlation` | effect matrix | same | perturbation similarity | none | MISSING_CURRENT_FIGURE |
| 7 | program_activity_by_cluster | `_plot_program_activity` | `program_activity` | same | programs vs states | none | MISSING_CURRENT_FIGURE |
| 7 | program_<P>_umap | `_plot_program_umaps` | obs program scores | same | program on embedding | none | MISSING_CURRENT_FIGURE |
| 7 | module_connectivity | `_plot_networks` | `module_connectivity` | `tf_network` | module links | none | MISSING_CURRENT_FIGURE |
| 7 | module_network, tf_hub_network | `_plot_networks` | edges, hubs | `tf_network` | network views | none | MISSING_CURRENT_FIGURE |
| — | (current only) perturbation_program_heatmap, protein_effect_heatmap, ps_vs_protein, program_activity_vs_protein, rna_vs_protein_effects, perturbation_overview, protein QC / representation figures | | | | protein extension | present | NOT_APPLICABLE_MULTIMODAL (kept) |
| — | Energy distance, distance space | | | | | | OUT_OF_SCOPE |

## E. Analysis differences, current (before) versus reference

### E.1 QC / normalization / embedding

| item | reference | current before | classification |
|---|---|---|---|
| prefilter | 200 genes / 3 cells | same | IDENTICAL |
| strict min genes | 1000 | 500 | UNINTENTIONAL_DIVERGENCE (default) |
| % mt cut | keep `< 20` | remove `> 20` (keeps == 20) | UNINTENTIONAL_DIVERGENCE (boundary) |
| gene filter after cell filter | `filter_genes(min_cells=3)` again | none | UNINTENTIONAL_DIVERGENCE (changes the gene set fed to HVG/PCA) |
| min counts / % hb | optional, off | min counts optional; no hb metric/filter | MISSING (option) |
| QC metrics | calculate_qc_metrics on raw counts | same (+ MAD flags) | IDENTICAL (+extension) |
| normalization | median library size + log1p | target_sum 1e4 + log1p | UNINTENTIONAL_DIVERGENCE (default) |
| HVG / scale / PCA / neighbors / UMAP / Leiden | 3000 seurat; HVG copy scaled max 10; PCA 50 arpack seed; 15 NN; UMAP 0.5; Leiden 1.0 igraph 2 it | same | IDENTICAL |
| Harmony (`batch_key`) | optional, off | not available | MISSING (option; reference default off) |
| clustering population | all QC-passing cells (default) | same | IDENTICAL |
| `assigned_only` | optional double embedding | `qc.filter.perturbation.cells` filters before embedding (no all-cell object) | INTENTIONAL_EXTENSION (different mechanism; default off in both) |
| guide dominant rule | non-zero non-dominant -> ambiguous | top < min_umi -> unassigned | INTENTIONAL_EXTENSION (documented adapter difference; both classes are excluded from every analysis) |
| stage order | clustering before every analysis | clustering after PS/lochNESS/modules | UNINTENTIONAL_DIVERGENCE (modules cannot use Leiden markers; lochNESS has no cluster summary) |

### E.2 Perturbation strength

| item | reference | current before | classification |
|---|---|---|---|
| whole analysis (target's own expression, log2FC on de-logged means, KS/MWU, BH per arm, KS-based hit call, ranking, both control arms, skipped reasons, 4 figure classes) | present | absent (only coverage counts) | MISSING |

### E.3 PS

| item | reference | current before | classification |
|---|---|---|---|
| signature / projection / scaling maths | pertps | identical (`ps_for_target`, checked to 4e-16 by `scripts/compare_reference_perturbation.py`) | IDENTICAL |
| eligible cells | targeting cells of the target + NTC | `single_targeting` + `single_control` in control_classes | IDENTICAL (vocabulary) |
| expression guard | 1 % of NTC cells | same | IDENTICAL |
| targets whose gene is not in var | skipped | scored, quadrants NA | UNINTENTIONAL_DIVERGENCE (default) |
| summary order / top targets | `pct_successful_kd` desc | `median_ps` desc | UNINTENTIONAL_DIVERGENCE |
| extra statistics (auc_vs_control, d_vs_control, sd/quantiles, raw scores, signatures) | none | present | INTENTIONAL_EXTENSION |
| storage | obs `ps_<gene>`, `ps_quadrant_<gene>` | obsm `ps_scores` (+ obs ps_score / ps_quadrant) | INTENTIONAL_EXTENSION (h5ad layout) |
| ps_vs_perturbation table | present | absent | MISSING |
| LDA embedding + 3 figure classes | present | absent | MISSING |
| outcome / escaper / quadrant figures | present | absent | MISSING |

### E.4 lochNESS

| item | reference | current before | classification |
|---|---|---|---|
| formula, denominator, representation, n_pcs 20, k 300 | as C.9 | same formula | IDENTICAL (checked exactly by the comparison script at equal k) |
| k cap | `min(300, n-1)` | `min(300, max(15, 0.1 n), n-1)` -> k = 180 on 1,800 cells | UNINTENTIONAL_DIVERGENCE (changes every score on objects < 3,000 cells) |
| permutation null / z / p / FDR | none | 200 label permutations by default | INTENTIONAL_EXTENSION but on by default: flagged; must default off |
| `delta_own_vs_control`, `mean_lochness_in_control_cells`, per-sample means | none | present | INTENTIONAL_EXTENSION (extra columns; reference-compatible) |
| `noise_delta` | option (0) | absent | MISSING (option) |
| `by_cluster`, `top_cluster` | present | absent | MISSING |
| figures | 5 classes | 2 (different) | MISSING / CURRENT_FIGURE_DIFFERS |

### E.5 Modules and programs

| item | reference | current before | classification |
|---|---|---|---|
| gene panel | union of top-100 Leiden markers per cluster (Wilcoxon, positive logFC) | union of top-100 response genes per target at FDR < 0.05 | UNINTENTIONAL_DIVERGENCE (the substitution was made because the pipeline had no clustering; it now has, so the reference panel is restored as the default and `response` stays an option) |
| pseudocount | 1e-9 | 1.0 | UNINTENTIONAL_DIVERGENCE (changes every log2FC and the clustering) |
| detection filter | none | 5 % of analysed cells | UNINTENTIONAL_DIVERGENCE (default) |
| control | ntc, fallback other | ntc only | MISSING (fallback) |
| Welch t / BH / DE mask / cluster_axis / module x program | as C.7 | identical | IDENTICAL |
| program activity | `sc.tl.score_genes` + per-cluster table | mean z-scored expression | UNINTENTIONAL_DIVERGENCE |
| TF hubs / edges / connectivity | present | absent | MISSING |
| `perturbation_program_effects` table, up/down counts | absent | present | INTENTIONAL_EXTENSION |
| figures | 9 classes | 3 | MISSING |

### E.6 Clustering and enrichment

| item | reference | current before | classification |
|---|---|---|---|
| Leiden parameters and graph | as C.4 | identical | IDENTICAL |
| reference arms | both `ntc` and `other`, primary `other` | one arm, default `non_targeting` | UNINTENTIONAL_DIVERGENCE (default reference group) |
| odds ratio | Haldane in `odds_ratio` | Fisher sample OR + `log2_or_haldane` | CURRENT_FIGURE_DIFFERS (both kept after the fix) |
| FDR family | within arm | all pairs (one arm) | IDENTICAL once one arm is considered |
| CMH | replaces Fisher p when stratified | reported beside Fisher | UNINTENTIONAL_DIVERGENCE (kept as an extra column; the reference behaviour becomes the default) |
| omnibus permutation | multinomial row resampling | full label permutation | UNINTENTIONAL_DIVERGENCE (minor; reference scheme restored) |
| guide concordance | significant pairs, same rule | all pairs, same rule | IDENTICAL values |
| composition / reference composition / effect magnitude / phenocopy | present | absent | MISSING |
| lochNESS side-by-side table | absent | present | INTENTIONAL_EXTENSION |
| figures | 6 classes | 1 (heatmap, rows by min FDR) | MISSING / CURRENT_FIGURE_DIFFERS |

## F. Corrections applied (2026-09-25)

See `docs/reference/PARITY_RESULTS.md` for the numerical comparison and the list of files.
In summary: perturbation strength implemented (`analysis/perturbation_strength.py`); stage
order changed to clustering -> perturbation strength -> enrichment -> modules -> PS -> lochNESS
-> protein -> concordance; reference defaults restored (QC 1000 genes, `< 20 %` mt, post-filter
gene filter, median normalization, enrichment arms/primary, lochNESS k and null, modules
panel/pseudocount/detection/scoring, PS ordering and skipping); missing tables and every
reference figure class recreated; protein sections kept as additions.
