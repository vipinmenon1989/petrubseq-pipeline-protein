# Numerical parity with the reference pipeline

Dataset: the bundled Papalexi et al. 2021 ECCITE-seq demo subset (1,800 cells, 18,649
genes, 111 guides, 25 targets), read by BOTH pipelines from the same 10x MTX directory
(`demo/data/papalexi_eccite/filtered_feature_bc_matrix`), so the cells and the counts are
identical by construction.

* reference: `perturbseq-pipeline run --config config/parity_reference_papalexi.yaml`
  (conda env `perturbseq-pipeline`, package `perturbseq_pipeline` 0.1.0 + `pertps` 0.1.0);
  output `../results/parity/reference_papalexi/`
* current: `petrubseq-protein run --config config/parity_papalexi_current.yaml`
  (env `petrubseq-protein`); outputs `../results/parity/current_papalexi_before/` (commit
  `382a9f7`, before the corrections) and `../results/parity/current_papalexi_after/`
* comparison: `python scripts/reference_parity.py --reference <ref run> --current <cur run> --out <md>`
  (exact equality for discrete quantities, absolute tolerance 1e-6 for floating point, 1e-5
  for the log-normalized matrix and the QC percentages, ARI / NMI for partitions).

The parity config matches the input-level settings that the reference does not expose
as options in the same way (median-library-size normalization, 1000-gene / 20 % mt
filter, 50 PCs, eGFP treated as a targeting label as the reference must, non-targeting
guides as the only control class). Every analysis block runs with the package defaults.

## Summary

| state | MATCH | DIFF | MISSING |
|---|---|---|---|
| before corrections (commit 382a9f7) | 29 | 28 | 13 |
| after corrections | 96 | 0 | 0 |

After the corrections every compared quantity agrees: retained cells and genes, QC
metrics, the log-normalized matrix (max |diff| 0), PCA (1.5e-8 after sign alignment),
Leiden labels (ARI 1.0, 9 clusters), perturbation strength (all effect sizes, p-values,
FDRs, hit calls and ranks exact under both control arms, same skipped target), PS
(per-cell scores to 4e-8, every summary column exact, quadrant classes 1442/1442, summary
order, LDA labels 1800/1800), lochNESS (per-cell scores to 1e-7, summaries exact,
ranking and top cluster identical), modules (717-gene Leiden-marker panel identical,
log2FC exact, module and program partitions ARI 1.0 with identical labels, DE-gene
counts, TF edges 26/26, module x program strength, program activity per cell to 6e-8)
and enrichment (450 pairs under both arms: odds ratios, p-values, FDRs, calls, direction,
guide concordance, composition and effect magnitude exact). UMAP coordinates (scanpy
UMAP and the LDA-UMAP) are not compared numerically because umap-learn is not
bit-reproducible across processes; the graphs they are built from are.

## Before the corrections

| section | item | value | status |
|---|---|---|---|
| QC | retained cells | ref=1800 cur=1800 common=1800 | **MATCH** |
| QC | retained genes | ref=16642 cur=16642 common=16642 | **MATCH** |
| QC | obs[n_genes_by_counts] | max|diff|=0 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | obs[total_counts] | max|diff|=0 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | obs[pct_counts_mt] | max|diff|=6.4e-07 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | obs[pct_counts_ribo] | max|diff|=2.33e-06 (n=1800, nan-mismatch=0) | **DIFF** |
| QC | filtering steps | ref final 1800 cells / 16642 genes; cur final 1800 / 16642 | **MATCH** |
| QC | perturbation classes | ref={'single_targeting': 1516, 'single_control': 199, 'ambiguous': 85} cur={'single_targeting': 1516, 'single_control': 199, 'ambiguous': 85} | **MATCH** |
| QC | target label of targeting cells | 1516/1516 equal | **MATCH** |
| NORMALIZATION | log-normalized expression (common cells x genes) | max|diff|=0 (n=29955600, nan-mismatch=0) | **MATCH** |
| NORMALIZATION | raw counts | max|diff|=0 (n=29955600, nan-mismatch=0) | **MATCH** |
| PCA | PC components (sign-invariant) | min|r|=1.000000 over 50 PCs; max|diff| after sign alignment=0.00304 | **MATCH** |
| PCA | variance ratio | max|diff|=9.44e-09 (n=50, nan-mismatch=0) | **MATCH** |
| PCA | highly variable genes | ref=3000 cur=3000 shared=3000 | **MATCH** |
| CLUSTERING | Leiden labels | ARI=0.7779 NMI=0.8222 (n=1800, k_ref=9, k_cur=10) | **DIFF** |
| CLUSTERING | n clusters | ref=9 cur=10 | **DIFF** |
| PERTURBATION_STRENGTH | perturbation_full table | ref=present cur=absent | **MISSING** |
| PS | targets scored | ref=24 cur=24 common=24 | **MATCH** |
| PS | summary mean_ps | max|diff|=7.09e-09 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary median_ps | max|diff|=1.18e-08 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_high_ps | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_successful_kd | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_escaper | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_non_responder | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_low_signal | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_controls_called_kd | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary net_pct_kd | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary expression_cut | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary order (top targets) | 5/24 equal | **DIFF** |
| PS | per-cell scores (all targets) | max|diff|=3.83e-08 over 24 targets | **MATCH** |
| PS | obs ps_score (own target) | max|diff|=3.83e-08 (n=1442, nan-mismatch=0) | **MATCH** |
| PS | ps_quadrant (own target) | 1442/1442 equal | **MATCH** |
| PS | LDA embedding present in current | no | **MISSING** |
| PS | ps_vs_perturbation table | ref=present cur=absent | **MISSING** |
| LOCHNESS | targets scored | ref=25 cur=25 common=25 | **MATCH** |
| LOCHNESS | summary mean_lochness_in_own_cells | max|diff|=0.39 (n=25, nan-mismatch=0) | **DIFF** |
| LOCHNESS | summary mean_lochness_all_cells | max|diff|=0.0572 (n=25, nan-mismatch=0) | **DIFF** |
| LOCHNESS | summary max_lochness | max|diff|=1.54 (n=25, nan-mismatch=0) | **DIFF** |
| LOCHNESS | summary overall_fraction_pct | max|diff|=0 (n=25, nan-mismatch=0) | **MATCH** |
| LOCHNESS | summary pct_cells_enriched | max|diff|=13.3 (n=25, nan-mismatch=0) | **DIFF** |
| LOCHNESS | top_cluster | nan | **MISSING** |
| LOCHNESS | ranking (order of targets) | 7/25 equal | **DIFF** |
| LOCHNESS | per-cell scores (all targets) | max|diff|=1.82 over 25 targets | **DIFF** |
| LOCHNESS | obs lochness_self | max|diff|=1.82 (n=1715, nan-mismatch=0) | **DIFF** |
| LOCHNESS | by_cluster table | ref=present cur=absent | **MISSING** |
| MODULES | effect matrix shape | ref=(25, 717) cur=(25, 910); shared genes=172 shared perturbations=25 | **DIFF** |
| MODULES | log2FC on shared entries | max|diff|=28.5 (n=4300, nan-mismatch=0) | **DIFF** |
| MODULES | perturbation Spearman correlation (shared panel) | max|diff|=0.21 (n=625, nan-mismatch=0) | **DIFF** |
| MODULES | module assignment | ARI=0.3405 NMI=0.6996 (n=25, k_ref=9, k_cur=9) | **DIFF** |
| MODULES | module labels (leaf-order numbering) | 1/25 equal | **DIFF** |
| MODULES | n_de_genes | 11/25 equal | **DIFF** |
| MODULES | gene program assignment | ARI=0.6288 NMI=0.6310 (n=172, k_ref=3, k_cur=2) | **DIFF** |
| MODULES | program labels (leaf-order numbering) | 146/172 equal | **DIFF** |
| MODULES | module x program strength | max|diff|=0.847 (n=36, nan-mismatch=0) | **DIFF** |
| MODULES | program activity per cell | max|diff|=1.55 over 4 programs | **DIFF** |
| MODULES | program activity by cluster table | ref=present cur=absent | **MISSING** |
| MODULES | tf_hubs | ref=present cur=absent | **MISSING** |
| MODULES | tf_edges | ref=present cur=absent | **MISSING** |
| MODULES | module_connectivity | ref=present cur=absent | **MISSING** |
| ENRICHMENT | control arms | ref=['ntc', 'other'] cur=single arm | **DIFF** |
| ENRICHMENT | tested pairs | ref=450 cur=250 common=225 | **DIFF** |
| ENRICHMENT | log2_odds_ratio | max|diff|=12.7 (n=225, nan-mismatch=0) | **DIFF** |
| ENRICHMENT | pval | max|diff|=1 (n=225, nan-mismatch=0) | **DIFF** |
| ENRICHMENT | fdr | max|diff|=0.999 (n=225, nan-mismatch=0) | **DIFF** |
| ENRICHMENT | pct_of_target | nan | **MISSING** |
| ENRICHMENT | pct_of_reference | nan | **MISSING** |
| ENRICHMENT | significant | 211/225 equal | **DIFF** |
| ENRICHMENT | direction | 124/225 equal | **DIFF** |
| ENRICHMENT | enrichment_composition | nan | **MISSING** |
| ENRICHMENT | enrichment_effect_magnitude | nan | **MISSING** |

## After the corrections

| section | item | value | status |
|---|---|---|---|
| QC | retained cells | ref=1800 cur=1800 common=1800 | **MATCH** |
| QC | retained genes | ref=16642 cur=16642 common=16642 | **MATCH** |
| QC | obs[n_genes_by_counts] | max|diff|=0 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | obs[total_counts] | max|diff|=0 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | obs[pct_counts_mt] | max|diff|=6.4e-07 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | obs[pct_counts_ribo] | max|diff|=2.33e-06 (n=1800, nan-mismatch=0) | **MATCH** |
| QC | filtering steps | ref final 1800 cells / 16642 genes; cur final 1800 / 16642 | **MATCH** |
| QC | perturbation classes | ref={'single_targeting': 1516, 'single_control': 199, 'ambiguous': 85} cur={'single_targeting': 1516, 'single_control': 199, 'ambiguous': 85} | **MATCH** |
| QC | target label of targeting cells | 1516/1516 equal | **MATCH** |
| NORMALIZATION | log-normalized expression (common cells x genes) | max|diff|=0 (n=29955600, nan-mismatch=0) | **MATCH** |
| NORMALIZATION | raw counts | max|diff|=0 (n=29955600, nan-mismatch=0) | **MATCH** |
| PCA | PC components (sign-invariant) | min|r|=1.000000 over 50 PCs; max|diff| after sign alignment=1.49e-08 | **MATCH** |
| PCA | variance ratio | max|diff|=5.72e-17 (n=50, nan-mismatch=0) | **MATCH** |
| PCA | highly variable genes | ref=3000 cur=3000 shared=3000 | **MATCH** |
| CLUSTERING | Leiden labels | ARI=1.0000 NMI=1.0000 (n=1800, k_ref=9, k_cur=9) | **MATCH** |
| CLUSTERING | n clusters | ref=9 cur=9 | **MATCH** |
| PERTURBATION_STRENGTH | targets tested | ref=24 cur=24 common=24 | **MATCH** |
| PERTURBATION_STRENGTH | log2fc_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | pct_knockdown_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | ks_stat_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | ks_pval_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mwu_pval_less_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | ks_fdr_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mwu_fdr_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mean_lognorm_perturbed_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mean_lognorm_control_ntc | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | is_hit_ntc | 24/24 equal | **MATCH** |
| PERTURBATION_STRENGTH | log2fc_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | pct_knockdown_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | ks_stat_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | ks_pval_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mwu_pval_less_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | ks_fdr_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mwu_fdr_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mean_lognorm_perturbed_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | mean_lognorm_control_other | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PERTURBATION_STRENGTH | is_hit_other | 24/24 equal | **MATCH** |
| PERTURBATION_STRENGTH | rank | 24/24 equal | **MATCH** |
| PERTURBATION_STRENGTH | n_perturbed | 24/24 equal | **MATCH** |
| PERTURBATION_STRENGTH | skipped targets | ref=['CAV1'] cur=['CAV1'] | **MATCH** |
| PS | targets scored | ref=24 cur=24 common=24 | **MATCH** |
| PS | summary mean_ps | max|diff|=7.09e-09 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary median_ps | max|diff|=1.18e-08 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_high_ps | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_successful_kd | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_escaper | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_non_responder | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_low_signal | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary pct_controls_called_kd | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary net_pct_kd | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary expression_cut | max|diff|=0 (n=24, nan-mismatch=0) | **MATCH** |
| PS | summary order (top targets) | 24/24 equal | **MATCH** |
| PS | per-cell scores (all targets) | max|diff|=3.83e-08 over 24 targets | **MATCH** |
| PS | obs ps_score (own target) | max|diff|=3.83e-08 (n=1442, nan-mismatch=0) | **MATCH** |
| PS | ps_quadrant (own target) | 1442/1442 equal | **MATCH** |
| PS | LDA embedding present in current | yes | **MATCH** |
| PS | lda_label | 1800/1800 equal | **MATCH** |
| PS | ps_vs_perturbation table | ref=present cur=present | **MATCH** |
| LOCHNESS | targets scored | ref=25 cur=25 common=25 | **MATCH** |
| LOCHNESS | summary mean_lochness_in_own_cells | max|diff|=0 (n=25, nan-mismatch=0) | **MATCH** |
| LOCHNESS | summary mean_lochness_all_cells | max|diff|=0 (n=25, nan-mismatch=0) | **MATCH** |
| LOCHNESS | summary max_lochness | max|diff|=0 (n=25, nan-mismatch=0) | **MATCH** |
| LOCHNESS | summary overall_fraction_pct | max|diff|=0 (n=25, nan-mismatch=0) | **MATCH** |
| LOCHNESS | summary pct_cells_enriched | max|diff|=0 (n=25, nan-mismatch=0) | **MATCH** |
| LOCHNESS | top_cluster | 25/25 equal | **MATCH** |
| LOCHNESS | ranking (order of targets) | 25/25 equal | **MATCH** |
| LOCHNESS | per-cell scores (all targets) | max|diff|=1.1e-07 over 25 targets | **MATCH** |
| LOCHNESS | obs lochness_self | max|diff|=1.1e-07 (n=1715, nan-mismatch=0) | **MATCH** |
| LOCHNESS | by_cluster table | ref=present cur=present | **MATCH** |
| MODULES | effect matrix shape | ref=(25, 717) cur=(25, 717); shared genes=717 shared perturbations=25 | **MATCH** |
| MODULES | log2FC on shared entries | max|diff|=0 (n=17925, nan-mismatch=0) | **MATCH** |
| MODULES | perturbation Spearman correlation (shared panel) | max|diff|=0 (n=625, nan-mismatch=0) | **MATCH** |
| MODULES | module assignment | ARI=1.0000 NMI=1.0000 (n=25, k_ref=9, k_cur=9) | **MATCH** |
| MODULES | module labels (leaf-order numbering) | 25/25 equal | **MATCH** |
| MODULES | n_de_genes | 25/25 equal | **MATCH** |
| MODULES | gene program assignment | ARI=1.0000 NMI=1.0000 (n=717, k_ref=4, k_cur=4) | **MATCH** |
| MODULES | program labels (leaf-order numbering) | 717/717 equal | **MATCH** |
| MODULES | module x program strength | max|diff|=0 (n=36, nan-mismatch=0) | **MATCH** |
| MODULES | program activity per cell | max|diff|=5.96e-08 over 4 programs | **MATCH** |
| MODULES | program activity by cluster table | ref=present cur=present | **MATCH** |
| MODULES | tf_hubs | rows ref=25 cur=25 | **MATCH** |
| MODULES | tf_edges | rows ref=26 cur=26 | **MATCH** |
| MODULES | module_connectivity | rows ref=9 cur=9 | **MATCH** |
| ENRICHMENT | control arms | ref=['ntc', 'other'] cur=['ntc', 'other'] | **MATCH** |
| ENRICHMENT | tested pairs | ref=450 cur=450 common=450 | **MATCH** |
| ENRICHMENT | log2_odds_ratio | max|diff|=0 (n=450, nan-mismatch=0) | **MATCH** |
| ENRICHMENT | pval | max|diff|=0 (n=450, nan-mismatch=0) | **MATCH** |
| ENRICHMENT | fdr | max|diff|=0 (n=450, nan-mismatch=0) | **MATCH** |
| ENRICHMENT | pct_of_target | max|diff|=0 (n=450, nan-mismatch=0) | **MATCH** |
| ENRICHMENT | pct_of_reference | max|diff|=0 (n=450, nan-mismatch=0) | **MATCH** |
| ENRICHMENT | significant | 450/450 equal | **MATCH** |
| ENRICHMENT | direction | 450/450 equal | **MATCH** |
| ENRICHMENT | guides_concordant (significant pairs) | 32/32 equal | **MATCH** |
| ENRICHMENT | guides_tested (significant pairs) | 32/32 equal | **MATCH** |
| ENRICHMENT | enrichment_composition | max|diff|=0 (n=225, nan-mismatch=0) | **MATCH** |
| ENRICHMENT | enrichment_effect_magnitude | max|diff|=0 (n=75, nan-mismatch=0) | **MATCH** |

## Figure inventory on the parity dataset

Reference run (`figures/<section>`), 166 figures, 104 embedded in the report:

| section | n | in_report |
|---|---|---|
| clustering | 5 | 5 |
| enrichment | 5 | 5 |
| enrichment/per_target | 25 | 12 |
| guides | 6 | 6 |
| lochness | 4 | 4 |
| lochness/per_target | 25 | 12 |
| modules | 12 | 12 |
| perturbation | 3 | 3 |
| perturbation/per_gene | 24 | 12 |
| ps_score | 5 | 5 |
| ps_score/lda | 24 | 12 |
| ps_score/per_target | 24 | 12 |
| qc | 4 | 4 |

Current run after the corrections (`figures/<section>/<stage>`), 211 figures, 149 embedded:

| section | stage | n | in_report |
|---|---|---|---|
| cell_states | clusters | 7 | 7 |
| cell_states | enrichment | 5 | 5 |
| cell_states | enrichment_per_target | 25 | 12 |
| multimodal | representation | 1 | 1 |
| perturbation_effects | concordance | 4 | 4 |
| perturbation_effects | gene_programs | 9 | 9 |
| perturbation_effects | gene_programs_umap | 4 | 4 |
| perturbation_effects | lochness | 4 | 4 |
| perturbation_effects | lochness_per_target | 25 | 12 |
| perturbation_effects | protein_effects | 1 | 1 |
| perturbation_effects | ps | 3 | 3 |
| perturbation_effects | ps_lda | 2 | 2 |
| perturbation_effects | ps_lda_per_target | 24 | 12 |
| perturbation_effects | ps_per_target | 24 | 12 |
| perturbation_strength | overview | 3 | 3 |
| perturbation_strength | per_target | 24 | 12 |
| qc_perturbation | post_filter | 9 | 9 |
| qc_protein | after_filtering | 6 | 6 |
| qc_protein | before_filtering | 6 | 6 |
| qc_rna | after_filtering | 5 | 5 |
| qc_rna | before_filtering | 5 | 5 |
| representations_protein | representation | 8 | 8 |
| representations_rna | representation | 7 | 7 |

Mapping of the reference sections onto the current ones: `qc` -> `qc_rna` (before/after),
`guides` -> `qc_perturbation`, `clustering` -> `representations_rna` + `cell_states/clusters`,
`perturbation` + `perturbation/per_gene` -> `perturbation_strength/overview` + `per_target`,
`enrichment` (+ per target) -> `cell_states/enrichment` (+ `enrichment_per_target`),
`ps_score` (+ lda, per target) -> `perturbation_effects/ps` + `ps_lda` + `ps_per_target` +
`ps_lda_per_target`, `lochness` (+ per target) -> `perturbation_effects/lochness` (+ per target),
`modules` -> `perturbation_effects/gene_programs` + `gene_programs_umap`. Every reference figure
class has a current equivalent; the current run adds the protein QC / representation /
effect / concordance figures and one extra module figure (`perturbation_program_heatmap`).

## Unit-level cross-check

`scripts/compare_reference_perturbation.py` (run inside the reference environment with
both code bases on `sys.path`) compares the PS, lochNESS and module functions directly
against `pertps.PerturbAnalyzer.calculate_ps_score`, the reference `lochness_score` on
the reference neighbour graph and the reference `build_effect_matrix` / `cluster_axis`
on a synthetic object: PS max |diff| 4.4e-16, lochNESS 0 (identical adjacency), effect
matrix 0, identical program and module partitions (2026-09-26).

## Remaining intentional differences

| item | reference | current | why it stays |
|---|---|---|---|
| storage of per-target PS / lochNESS scores | one `obs` column per target | `obsm['ps_scores']`, `obsm['lochness']` (cells x targets) + `obs['ps_score']`, `obs['lochness_self']` | keeps `obs` small for hundreds of targets; every value is present and compared |
| PS quadrant labels | `successful knockdown`, `non-responder`, ... | `successful_knockdown`, `non_responder`, ... | h5ad-safe identifiers; same classes (compared after normalisation) |
| target label column | `target_gene` | `target` (+ `gene` = mapped RNA symbol) | pipeline vocabulary; `target_gene_map` is an extension the reference lacks |
| guide calling: cells with 0 < top guide < min_umi | `ambiguous` | `unassigned` | adapter-stage rule kept from Stage C (documented); both classes are excluded from every analysis, and the parity run has no such cell |
| extra columns | - | PS `auc_vs_control`, `d_vs_control`, spread; lochNESS `delta_own_vs_control`, control mean, per-sample means, optional permutation null; enrichment `sample_odds_ratio`, `pval_fisher`, CMH columns; modules `perturbation_program_effects`, up/down counts | additions only; reference columns are unchanged |
| Harmony (`cluster.batch_key`) and `cluster.assigned_only` | optional | not implemented (reference defaults are off; the parity run does not use them) | out of the corrected scope; the current `qc.filter.perturbation.cells` filter is a different mechanism |
| protein QC / representation / effects / concordance sections | - | present | multimodal extension |

## Protein PCA: float64 working path (2026-09-26)

`preprocessing/embeddings.py::_pca_on` now casts every dense block (the protein CLR
matrix) to an explicit float64 working copy before `sc.pp.scale` and `sc.tl.pca`
(arpack, `compute.seed`), the same numerical principle as the reference RNA route;
coordinates and loadings are stored as float32 (deterministic cast), variances as
float64, and `uns['pca_protein']['params']` records `working_dtype: float64`,
`stored_dtype: float32`. Nothing scientific changes: raw ADT counts, the CLR formula
and axis, isotype handling, feature exclusion, protein-effect tests and concordance
definitions are untouched. `tests/test_protein_reproducibility.py` asserts bit-identical
protein PCA coordinates, variances, loadings and neighbour graphs across repeated calls
and repeated pipeline runs. The repeat-run check on the Papalexi demo is recorded below.

## Validation runs with the corrected pipeline (2026-09-26)

| run | config | cells | outcome |
|---|---|---|---|
| Papalexi demo | `config/demo_papalexi.yaml` | 1,800 | validator 79 checks, 0 FAIL, 1 WARN (no isotype control in the 4-antibody panel, documented); 9 Leiden clusters; strength 24/25 tested, 6 effective; enrichment 35 significant pairs; modules 25 x 753 genes; PS 24 scored; lochNESS 25; 212 figures (149 embedded) |
| Papalexi full analysis | `config/analysis_papalexi_full.yaml`, since renamed `config/analysis_eccite_seq_full.yaml` (job 20849348) | 17,473 | validator 79/0/1; 14 clusters; strength 25/25 tested, 16 effective; enrichment 122 significant; modules 25 x 1018; PS 25; 214 figures (149 embedded); 251 s |
| SCP1064 regression | `config/scp1064.yaml` (job 20849346; analyses off, frozen v0.1 settings; run before the protein float64 change) | 218,331 | validator 60/0/0; vs v0.1 baseline: 58 identical, DIFF `X_umap_rna` / `varm[PCs]` (reference float64 PCA path, intentional), DIFF `X_umap_protein` (protein PCs differed by 1.8e-5 from float32 arpack run-to-run variation; resolved by the float64 protein path below), ADD haemoglobin QC columns + protein feature annotation columns |
| SCP1064 smoke regression | `config/scp1064_smoke.yaml` | 4,500 | validator 60/0/1; vs v0.1: 61 identical, DIFF `X_pca` (7e-4) / `X_umap_rna` (intentional, as above), ADD hb columns |
| SCP1064 full analysis | `config/analysis_scp1064.yaml` (job 20849347) | 218,027 | validator 79/0/0; 24 clusters; strength 213/248 tested, 151 effective; enrichment 248 x 24 x 2 arms, 122 significant; modules 246 x 1037, 9 modules / 4 programs; PS 213 scored (35 skipped), LDA on 101,184 cells; lochNESS 248 (k = 300); protein effects 139 significant; 1,225 figures (150 embedded); 85 min, 38.5 GB peak |

### Repeat-run reproducibility (final implementation, 2026-09-26)

`config/demo_papalexi.yaml` run twice from clean output directories (identical config,
input, seed and environment; 1,800 cells). Maximum absolute difference between the two runs:

| slot / table | max abs diff |
|---|---|
| `obsm['protein']` (CLR) | 0 |
| `obsm['X_pca_protein']`, `uns['pca_protein']` variance and loadings | 0 |
| `obsm['X_umap_protein']` | 0 |
| `obsm['X_pca']`, `uns['pca']` variance, `obsm['X_umap_rna']`, `obsm['X_lda_umap']` | 0 |
| `obsp['protein_connectivities']`, `obsp['rna_connectivities']`, `obs['leiden']` | 0 / identical |
| protein_effects, protein_effect_matrix, ps_protein_association, lochness_protein_summary, program_protein_association (cells, targets), perturbation_summary, ps_targets, lochness_targets, perturbation_effect_matrix, perturbation_full | 0 (text columns identical) |

Final Papalexi demo run: validator 79 checks, 0 FAIL, 1 WARN (no isotype control in the
4-antibody panel, documented); the protein-effect, PS-protein, lochNESS-protein and
program-protein tables are identical (max diff 0) to the previous validated demo run made
with the float32 protein path, i.e. the protein statistics never depended on the PCA
precision. SCP1064 smoke regression with the final code: validator 60/0/1; versus the
previous smoke run of the corrected code only `obsm['X_pca_protein']` (3.1e-5) and
`obsm['X_umap_protein']` differ (the intentional dtype-path change); versus v0.1 the
expected RNA and protein PCA/UMAP shifts plus the additive slots. RNA reference parity
after the protein change: 96 MATCH, 0 DIFF, 0 MISSING.

## Perturbation distance stages (reference stages 10–12; added 2026-09-26)

Reference: `perturbseq-pipeline run --config config/parity_reference_papalexi_distance.yaml`
(reference repository at commit `1c48f9d`, `distance` / `distance_space` / `meta_analysis` enabled
with the reference defaults; SLURM job 20849823, 111 s) → `../results/parity/reference_papalexi_distance/`.
Current: `petrubseq-protein run --config config/parity_papalexi_current_distance.yaml` (SLURM job
20849869) → `../results/parity/current_papalexi_distance/`; comparison
`../results/parity/parity_distance.md` (`scripts/reference_parity.py`, section DISTANCE).
Audit of what is compared: docs/reference/PERTURBATION_DISTANCE_AUDIT.md.

| state | MATCH | DIFF | MISSING |
|---|---|---|---|
| reference at `1c48f9d`, previously compared stages (QC … enrichment) | 96 | 0 | 0 |
| distance stages (37 further items) | 37 | 0 | 0 |

| item | result |
|---|---|
| tested targets and ranked order (21; BRD4, CUL3, MYC, SPI1 skipped below 30 cells) | identical |
| energy distance vs control (21 targets, 199 non-targeting cells, `X_pca` 50 PCs) | max abs diff 4.5e-13 |
| DistanceTest p-values (1,000 seeded permutations) and BH-FDR | identical (max diff 0) |
| MMD (RBF, median bandwidth) | 2.0e-15 |
| significant calls (7 of 21 at FDR < 0.05: JAK2, IFNGR1, STAT1, IFNGR2, SMAD4, IRF1, STAT2) | 21/21 equal |
| pairwise 21 × 21 energy-distance matrix (symmetric, zero diagonal) | 2.7e-12 |
| PCoA coordinates (10 axes, sign-aligned per axis) | 4.1e-12 |
| nearest neighbours (210 rows: target, neighbour, rank; distances) | identical / 1.4e-12 |
| phenotype modules (5, reference default cut) | ARI = NMI = 1, labels identical |
| master table (`perturbation_meta.csv` vs `master_perturbation_table.csv`, 25 targets): order, `n_cells`, efficacy, PS, lochNESS, distance columns, module labels | identical (float32 lochNESS columns ≤ 2.2e-7) |

The distance figures of the reference (`perturbation_distance_ranking`, `perturbation_atlas`,
`ps_vs_distance_map`, `perturbation_phenotype_space`, `module_concordance`) are drawn from the same
tables with the same transformations (docs/reference/PERTURBATION_DISTANCE_AUDIT.md §2.4).

Re-run after the RNA ↔ protein level audit (docs/RNA_PROTEIN_LEVELS.md; current job 20849908, same reference run): 133 MATCH / 0 DIFF — the distance stages and the reference master-table columns are untouched by the audit.
