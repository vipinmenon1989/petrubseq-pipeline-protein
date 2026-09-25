# Perturbation-analysis reference audit (Stage E, written before implementation)

Audited 2026-09-25 from the group's local reference implementation

* the group's local checkout of `perturbseq-pipeline` (`src/perturbseq_pipeline/`)
  (`ps_score.py` 563 l., `lochness.py` 315 l., `modules.py` 522 l., `perturbation.py`, `cluster.py`, `guides.py`, `config.py`),
  its `README.md` ("What it does" sections 5–7) and `docs/methods.md` (sections 7–8), `docs/ps_python_proposal.md`;
* the `pertps` 0.1.0 package (PS_python, authors Vipin Menon and Wei Li) installed in the
  `perturbseq-pipeline` conda environment (`pertps/analyzer.py`, which holds the PS mathematics);
* public repository: https://github.com/vipinmenon1989/perturbseq-pipeline (the local copy is authoritative).

Reference vocabulary (used below): `obs['target_gene']`, `obs['perturbation_class']` ∈ {`targeting`,
`non-targeting`, `ambiguous`, `unassigned`}, `layers['lognorm']` (log-normalized), `obsm['X_pca']` /
`X_pca_harmony`, `obs['leiden']`. Our vocabulary: `obs['target']`, `obs['perturbation_class']` ∈
{`single_targeting`, `single_control`, `multi_*`, `mixed_control_targeting`, `ambiguous`, `unassigned`},
`obs['control_class']` (`non_targeting`, …), `X` (log-normalized), `obsm['X_pca']`, no clustering.

## Concepts (E1) — kept distinct throughout

| term | meaning | where it lives |
|---|---|---|
| perturbation target | the gene the guide is designed against (`obs['target']`) | assignment (Stage C) |
| downstream response genes | genes whose expression changes in cells carrying the perturbation | PS signature genes; effect-matrix columns |
| gene program | a coordinated group of response genes (cluster of the gene axis of the effect matrix) | `modules.py` → `P1..` |
| perturbation module | a group of perturbations with similar response profiles (cluster of the perturbation axis) | `modules.py` → `M1..` |

## Literature verification (E1)

| method | verified citation | note |
|---|---|---|
| PS | Song B., Liu D., Dai W., McMyn N.F., Wang Q., Yang D., et al. *Decoding heterogeneous single-cell perturbation responses.* **Nat Cell Biol** 27, 493–504 (2025). doi:10.1038/s41556-025-01626-9, PMID 40011559. Python port: PS_python (https://github.com/weili-lab/PS_python, `pertps`). | PS_python's README calls it a "translated scMAGeCK EM algorithm"; the installed `pertps` uses the closed-form OLS projection described below (comment: "mathematically equivalent to minimizing least squares"). scMAGeCK itself: Yang L. et al., Genome Biol 21, 19 (2020). |
| lochNESS | Huang X., Henck J., Qiu C., Sreenivasan V.K.A., Balachandran S., Amarie O.V., et al. *Single-cell, whole-embryo phenotyping of mammalian developmental disorders.* **Nature** 623, 772–781 (2023). doi:10.1038/s41586-023-06548-w, PMC10665194 — "local cellular heuristic neighbourhood enrichment specificity score": for each cell, the fold change of the observed vs expected number of mutant cells among its k nearest neighbours in aligned PC space. Implementation ported by the reference from pertTF (https://github.com/davidliwei/pertTF, `perttf.model.composition_change_analysis`; bioRxiv 10.64898/2026.03.12.711379). | Reference formula `local_fraction/overall_fraction − 1` = fold change − 1, k = 300, PCA space. |
| modules / gene programs | Zhou P., Shi H., Huang H., Sun X., Yuan S., Chapman N.M., et al. *Single-cell CRISPR screens in vivo map T cell fate regulomes in cancer.* **Nature** 624, 154–163 (2023). doi:10.1038/s41586-023-06733-x, PMID 37968405 (GSE216909). | **The reference README/docstring attribute this to "Chen et al. 2023"; the DOI resolves to Zhou et al. 2023 (Chi lab).** Our documentation cites Zhou et al. |

## 1. PS — perturbation-response score

**REFERENCE IMPLEMENTATION** `pertps.PerturbAnalyzer.calculate_ps_score` (analyzer.py 18–118), driven by
`perturbseq_pipeline/ps_score.py` (`_prepare_for_pertps`, `compute_ps_scores`, `_summarize_target`, `attach_scores`).

**PURPOSE** Per-cell strength of the transcriptional response to a perturbation, so that a perturbed population can be split into responders and non-responders (and, with the target's own expression, "successful knockdown / escaper / non-responder / low signal" quadrants).

**BIOLOGICAL INTERPRETATION** 0 = indistinguishable from non-targeting controls along the perturbation's downstream signature; 1 = the strongest response observed for that target. Scores are comparable *within* a target, not across targets (each target is max-normalized).

**INPUT** log-normalized expression (`layers['lognorm']`, explicitly not scaled values), `obs['gene']` with the target label per cell and a fixed control label (`Non-Targeting`); cells outside {target, controls} are excluded from the fit. Requires the target gene in `var_names` for the expressed-in-controls guard and the quadrants.

**ALGORITHM** (per target, deterministic, no random component)
1. Subset to target cells ∪ non-targeting cells; skip if `< 10` cells.
2. Guard: skip when the target gene is expressed (> 0) in `< min_pct_expressing_control` (1 %) of control cells.
3. Response-gene selection: `sc.tl.rank_genes_groups(groupby='gene', reference=control, method='t-test')`; take the **top 100** genes (`top_n`) ranked for the target group (highest t-statistic; both up and down regulated genes are not balanced — the ranking is by score, descending, so it favours up-regulated genes; this is the reference behaviour).
4. OLS signature: `x` = 1 for target cells, 0 for controls; `Y` = expression of the selected genes; `beta_g = cov(x, Y_g)/var(x)` (= mean difference target − control for a binary design).
5. Score per cell: projection of the centred expression onto beta, `s_i = (Y_i − mean(Y)) · beta / (beta · beta)`.
6. Shift so the control mean is 0; clip to `[0, scale_factor]` (3.0); divide by the maximum (→ `[0, 1]`).
7. Return a Series over the subset cells (target + control cells only).

**DEFAULT PARAMETERS** `top_n_biomarkers 100`, `scale_factor 3.0`, `ps_threshold 0.5` (quadrant cut), `expression_cut mean` (control mean of the target's expression; the reference found the median degenerate for most genes), `min_cells_per_target 10`, `min_control_cells 10`, `min_pct_expressing_control 1.0`. The optional LDA/UMAP embedding (`compute_lda_umap`) is a visualization extra, not part of the score.

**OUTPUT** `obs['ps_<TARGET>']` (score for target + control cells, NaN elsewhere), `obs['ps_quadrant_<TARGET>']`, `obs['ps_score']` (each cell's score for its *own* target; NaN for controls/others), `obs['ps_quadrant']`; per-target summary (`n_perturbed_cells`, `n_control_cells`, `mean_ps`, `median_ps`, `pct_high_ps`, `pct_successful_kd`, `pct_escaper`, `pct_non_responder`, `pct_low_signal`, `pct_controls_called_kd`, `net_pct_kd`, `expression_cut`); skipped-target table with reasons.

**ASSUMPTIONS** one target per cell; non-targeting cells are the right baseline; log-normalized values; a linear signature captures the response; the t-test ranking with `top_n` genes is a reasonable response-gene set; the target gene is measured.

**LIMITATIONS** max-normalization makes scores target-specific (documented in the reference proposal); the OLS projection uses the *same* cells to fit and to score (in-sample), so target cells are biased upward relative to a held-out estimate; no confidence interval; genes chosen by t-test favour high-expression genes; scores for cells that are neither target nor control are not defined (the reference leaves them NaN).

**DEPENDENCIES** numpy, pandas, scanpy (`rank_genes_groups`), scipy. No `pertps` import is needed for a faithful re-implementation (~60 lines); the LDA extra needs scikit-learn + umap.

**WHAT SHOULD BE PRESERVED** the definition end-to-end (subset, guard, t-test top-100 selection, OLS beta, projection, control-centering, clip at 3, max-normalization), the per-target summary columns, the quadrant logic with the control-mean cut, the skip reasons, `obs['ps_score']` as the cell's own-target score, determinism.

**WHAT SHOULD BE IMPROVED FOR THE NEW PIPELINE** (a) use `obs['control_class'] == 'non_targeting'` single-control cells as the baseline and only `single_targeting` cells as perturbed (never `ambiguous`, `multi_*`); (b) store per-target scores in one `obsm['ps_scores']` DataFrame (cells × targets, NaN outside each comparison) instead of dozens of `obs` columns, keep `obs['ps_score']`/`obs['ps_quadrant']`; (c) record `n_signature_genes`, the signature itself (`uns` table target × gene beta) and the raw (unnormalized) control-shifted scale so targets can be compared on the `[0, 3]` scale as well; (d) provenance in `uns['petrubseq_protein']['analysis']['ps']`; (e) no LDA embedding (out of scope; UMAP already exists).

## 2. lochNESS — local neighbourhood enrichment

**REFERENCE IMPLEMENTATION** `perturbseq_pipeline/lochness.py` (`_build_neighbor_graph`, `lochness_score`, `compute_lochness`, `attach_scores`), a vectorized port of pertTF's `composition_change_analysis`.

**PURPOSE** For every cell and every perturbation g: is g over- or under-represented among the cell's nearest transcriptional neighbours relative to its overall frequency? Maps *where* in the cell-state manifold a perturbation accumulates, without clusters.

**BIOLOGICAL INTERPRETATION** `lochNESS(cell, g) = local_fraction(g) / overall_fraction(g) − 1`: 0 = chance, > 0 locally enriched, < 0 depleted. A perturbation whose own cells have high lochNESS occupies a distinct cell state; a value near 0 in its own cells means it is transcriptionally indistinguishable from the population.

**INPUT** a k-nearest-neighbour graph over `obsm['X_pca_harmony']` if present else `obsm['X_pca']` (first `n_pcs` = 20 components), k = 300 (`sc.pp.neighbors`, seed = `run.seed`), the per-cell label (`genotype_key = target_gene`), the class column to restrict candidates to `targeting` cells with `>= min_cells_per_target` (10). The `overall_fraction` is computed over **all** cells in the object (every label, including controls / ambiguous).

**ALGORITHM** binary adjacency from the distances graph; `neighbor_counts` = actual neighbours per row; for each candidate target: `local = adj @ indicator`, `score = local/neighbor_counts / overall_fraction − 1`; optional Gaussian jitter (`noise_delta`, default 0). Per-target summary over all cells and over the target's own cells; `lochness_self` = each cell's score for its own label (also computed for the non-targeting label); per-cluster means when Leiden clusters exist.

**DEFAULT PARAMETERS** `n_neighbors 300`, `n_pcs 20`, `use_rep None` (→ harmony PCA or PCA), `min_cells_per_target 10`, `enrichment_cut 0.5` (for `pct_cells_enriched`), `noise_delta 0`, `recompute_neighbors True`.

**OUTPUT** `obs['lochness_<TARGET>']` per target (all cells), `obs['lochness_self']`, summary (`n_cells`, `overall_fraction_pct`, `mean_lochness_all_cells`, `mean_lochness_in_own_cells`, `max_lochness`, `pct_cells_enriched`, `top_cluster`), `by_cluster` table, skipped table.

**ASSUMPTIONS** the PCA space reflects cell state; k = 300 is small relative to n (with 1,800 cells k = 300 is 17 % of the data, so the score is heavily smoothed — see improvements); labels are exclusive.

**LIMITATIONS** no significance/permutation test (Huang et al. compare against permuted labels); no batch/lane awareness beyond using a harmony embedding; `mean_lochness_all_cells` is ≈ 0 by construction; with small n and large k the dynamic range collapses.

**DEPENDENCIES** scanpy (`pp.neighbors`), scipy.sparse.

**WHAT SHOULD BE PRESERVED** the exact score definition (fold change − 1 with actual neighbour counts), PCA-space neighbours (never UMAP), k = 300 default, the per-target summary columns (`mean_lochness_in_own_cells` as the headline), candidate/skip rules, determinism (seeded kNN).

**WHAT SHOULD BE IMPROVED** (a) cap k at a fraction of n for small objects (`k = min(300, max(15, floor(0.1·n)))`, recorded) so the demo is not degenerate, with the cap documented as an intentional deviation only when it triggers; (b) a label-permutation null per target (`n_permutations`, seeded) giving a z-score / empirical p for the own-cell mean, following Huang et al.; (c) restrict candidates to `single_targeting` cells and use `single_control` non-targeting cells as the explicit reference label; (d) store per-target scores as `obsm['lochness']` (cells × targets) plus `obs['lochness_self']`; (e) per-lane / per-sample own-cell means to expose batch-driven enrichment; (f) provenance in `uns[...]['analysis']['lochness']`.

## 3. Modules and gene programs

**REFERENCE IMPLEMENTATION** `perturbseq_pipeline/modules.py` (`select_perturbations`, `select_genes`, `build_effect_matrix`, `cluster_axis`, `module_program_strength`, `tf_network`, `compute_modules`).

**PURPOSE** Build the *regulome map*: a perturbation × response-gene effect matrix, then (i) cluster genes into co-regulated **gene programs** and (ii) cluster perturbations into co-functional **perturbation modules**, and relate them.

**BIOLOGICAL INTERPRETATION** two perturbations in one module produce similar response profiles (correlated log2FC vectors); a program is a set of genes that move together across perturbations. Module × program strength (mean log2FC of a program's genes under a module's perturbations) is the signed link between the two.

**INPUT** log-normalized expression; targets with `>= min_cells_per_perturbation` (20) `targeting` cells; a gene panel: union of the top 100 positive Leiden cluster markers per cluster (`gene_selection cluster_markers`, Wilcoxon) or HVGs (`hvg`); control = non-targeting cells (`ntc`) or leave-the-target-out other-target cells (`other`).

**ALGORITHM** effect matrix `log2((mean_pert + 1e-9)/(mean_ctrl + 1e-9))` on de-logged means (pseudobulk mean-difference log2FC), plus a Welch t-test per gene on the log values with BH-FDR within each perturbation → DE mask (|log2FC| > 0.5 & FDR < 0.05). Gene programs: hierarchical clustering (average linkage) of genes on `1 − Pearson` correlation of their effect profiles across perturbations, cut into `n_programs` (4) or at distance 0.7; perturbation modules: same on perturbations with `1 − Spearman`, `n_modules` (9). Labels `P1..`, `M1..` in dendrogram-leaf order. Per-cell program activity: `sc.tl.score_genes(program genes, ctrl_size=50)` (this is *not* deterministic unless seeded — scanpy samples control genes with `random_state=0` by default, so it is reproducible). TF-hub network from the DE mask (optional).

**DEFAULT PARAMETERS** `min_cells_per_perturbation 20`, `min_perturbations 5`, `min_genes 10`, `n_marker_genes_per_cluster 100`, `control ntc`, `program_correlation pearson`, `module_correlation spearman`, `linkage average`, `n_programs 4`, `n_modules 9`, `cluster_distance_threshold 0.7`, `hub_lfc_threshold 0.5`, `de_fdr_alpha 0.05`, `score_programs True`.

**OUTPUT** `effect_matrix` (perturbations × genes), `gene_programs` (gene, program, program_size), `modules` (target_gene, module, n_cells, n_de_genes), `module_program` strength, `program_activity` (program × cluster), `obs['program_P<k>_score']`, TF edges/hubs.

**ASSUMPTIONS** enough perturbations (≥ 5) and a response-gene panel (needs clusters or HVGs); log2FC of means is an adequate effect; cluster numbers are user choices.

**LIMITATIONS** the gene panel comes from Leiden markers, which our pipeline does not compute (no clustering by design); fixed `n_programs`/`n_modules` are arbitrary; the effect matrix is unregularized (noisy for low-cell targets); programs are unsigned clusters (up- and down-regulated genes of one response can land in different programs); the network step needs networkx.

**DEPENDENCIES** scanpy, scipy.cluster.hierarchy, (networkx optional).

**WHAT SHOULD BE PRESERVED** the two-axis construction from one perturbation × gene log2FC-vs-control matrix; Pearson for genes, Spearman for perturbations, average linkage on `1 − r`; numeric labels only (`P1..`, `M1..`), never biological names; the module × program strength; skip rules when underpowered; DE gate for "responds" counts.

**WHAT SHOULD BE IMPROVED** (a) gene panel without clusters: the union of each perturbation's significant response genes (Welch t-test vs NT, BH within perturbation, top-N by |t|), falling back to HVGs — so the panel *is* "downstream response genes" by construction; (b) an explicit `perturbation_program_effects` table (perturbation × program mean log2FC, with the per-program DE fraction) so PERTURBATION → PROGRAM is a first-class output; (c) program activity per cell as the mean z-scored expression of program genes (deterministic, no control-gene sampling), stored in `obsm['program_activity']`; (d) provenance/parameters in `uns[...]['analysis']['modules']`; (e) keep the network layer out of scope.

## 4. Protein effects, PS ↔ protein, lochNESS ↔ protein, program ↔ protein (new)

Not in the reference (RNA-only). Design: CLR values in `obsm['protein']` (never counts); per (target, protein): perturbed = `single_targeting` cells of the target, control = `single_control` non-targeting cells; effect = difference of means of CLR values (a log-ratio scale, so a mean difference is a log fold change of the geometric-mean-normalized abundance), Cohen's d, Mann–Whitney U p-value, BH-FDR across all target × protein tests; `min_cells` support gate. Associations: Spearman rho within target (PS vs CLR) with BH-FDR over the tested pairs; target-level Spearman/Pearson between lochNESS own-cell mean and protein effect across targets; program activity vs CLR Spearman over cells (all single-guide cells) and target-level program effect vs protein effect. Sample/lane structure is exposed by reporting per-sample effect signs and a per-target "driven by one sample" flag.

## 5. Computational complexity (reference)

PS: per target one `rank_genes_groups` over ≤ n_target + n_ctrl cells × all genes (t-test, O(cells × genes)) plus an O(cells × 100) projection — seconds per target. lochNESS: one kNN graph (k = 300) over n cells plus one sparse mat-vec per target — the graph dominates (minutes at 100k cells, seconds at 1.8k). Modules: dense cells × panel-genes matrices (panel ≤ a few thousand), two correlation matrices and hierarchical clustering — seconds for ≤ 60 perturbations.

## 6. Cross-check plan (E17)

Run the reference `pertps.PerturbAnalyzer.calculate_ps_score`, `lochness.lochness_score` (+ a seeded `sc.pp.neighbors` graph) and `modules.build_effect_matrix`/`cluster_axis` in the `perturbseq-pipeline` environment on the same small synthetic AnnData (translated to the reference vocabulary), and compare with our implementation: PS per-cell scores (expect identical to 1e-6 given the same gene set; the gene set itself compared as a set), lochNESS per-cell scores (identical given the same graph), effect matrix (identical), program/module labels (identical partition up to label names). Any disagreement stops the work until explained.
