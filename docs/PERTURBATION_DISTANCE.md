# Perturbation distance, distance space, phenotype modules and the master tables

Reference stages 10–12 of `weili-lab/perturbseq-pipeline` (`distance.py`, `meta.py`; audit and
numerical parity: docs/reference/PERTURBATION_DISTANCE_AUDIT.md, docs/reference/PARITY_RESULTS.md).
The reference implements these statistics itself (pertpy is credited as the inspiration; it is
not executed), so this pipeline ports the reference functions and reproduces their numbers on
the same cells. Configuration lives under `analysis.perturbation_effects` (`distance`,
`distance_space`, `master_table`; every default is the reference value, docs/DEFAULTS.md). As
in the reference, `distance` and `distance_space` are **off by default** and the master table is
on; the demo and analysis configs enable them.

## Three objects

| object | question | where |
|---|---|---|
| **control distance** | how far does a perturbation move the whole transcriptional state away from the controls, and is that shift larger than chance? | `tables/perturbation_distance/perturbation_distance.csv`, `uns['perturbation_distance']` |
| **pairwise distance** | how similar are the cell-state distributions produced by two perturbations? | `perturbation_distance_matrix.csv`, `uns['perturbation_distance_matrix']`, `perturbation_neighbors.csv` |
| **phenotype space** | a low-dimensional map of the pairwise distances (PCoA) and its partition into **phenotype modules** | `perturbation_space_coordinates.csv`, `phenotype_modules.csv`, `uns['phenotype_modules']` |

All three use the cells' PCA coordinates (`distance.representation`, `X_pca`; never a UMAP).
Perturbed cells of a target = its `single_targeting` cells; controls = `single_control` cells of
`control_classes` (the `ntc` arm; the reference falls back to all targeting cells when no
non-targeting cells exist). Targets with fewer than `min_cells` (30) cells are skipped and listed.

## Control distance (`distance`)

Energy distance between the target's cells X and the control cells Y:
`E(X, Y) = 2·mean‖x − y‖ − mean‖x − x′‖ − mean‖y − y′‖` over the full pairwise Euclidean
distance matrices (the pertpy `Edistance` convention, diagonal zeros included). It is zero only
when the two distributions coincide and responds to shifts in mean, spread and shape.
Significance: a **DistanceTest** with `n_permutations` (1000) seeded label permutations of the
pooled cells, `p = (1 + #(permuted ≥ observed)) / (1 + B)`, Benjamini–Hochberg across the tested
targets, `significant = fdr < fdr_threshold` (0.05). `secondary_metric: mmd` adds the squared
maximum mean discrepancy with a Gaussian kernel (median-distance bandwidth). Groups larger than
`max_cells_per_target` (2000) / `max_control_cells` (5000) are subsampled deterministically
(`random_seed` 123; per-target seeds derived by sha256 as the reference does; `stratify_by` names
an obs column for proportional sampling, e.g. `condition` on SCP1064). Cost per target is one
`(n_t + n_c)²` distance matrix (≤ 7,000² doubles ≈ 0.4 GB) plus `B` index gathers; targets run
in parallel over `compute.n_jobs` with identical results.

Table columns: `target, n_cells, n_control, energy_distance, pvalue, mmd_distance, fdr, significant`,
sorted by energy distance. Figures (report section *Perturbation distance vs control*):
`perturbation_distance_ranking`, `perturbation_atlas` (targets × {knockdown strength, PS median,
lochNESS mean, energy distance}, column z-scores, FDR stars, module strips) and `ps_vs_distance_map`
(penetrance vs phenotype magnitude, point size lochNESS, colour phenotype module).

## Pairwise distance, phenotype space and phenotype modules (`distance_space`)

The same energy distance between every pair of eligible targets (cells sampled per target with
seed `random_seed + 43·index`) gives a symmetric target × target matrix with zero diagonal.
**PCoA** (classical multidimensional scaling: double-centred squared distances, eigen-decomposition,
strictly positive eigenvalues, `coords = V·sqrt(λ)`, `n_components` 10) gives the phenotype-space
coordinates; each target's `nearest_neighbors` (10) closest targets are listed with rank and
distance; **phenotype modules** are average-linkage clusters of the matrix (`linkage_method`),
cut into `n_modules` groups, or at `cluster_distance_threshold`, or by the reference rule
`max(2, min(9, K // 4))` when both are null (labels `PM1..`). Co-functional modules group
perturbations by shared downstream *genes*; phenotype modules group them by the similarity of
the *cell states* they produce. The two partitions are compared with ARI / NMI in the
`module_concordance` figure and `uns['petrubseq_protein']['analysis']['perturbation_effects']['module_concordance']`.
Figures: `perturbation_phenotype_space` (PCoA1/2 by module and by energy distance),
`perturbation_distance_matrix` (extension: the matrix in dendrogram order), `module_concordance`.
The pairwise stage is O(K²) pairs; K ≈ 250 targets (SCP1064) is ~31,000 pairs of ≤ 2000 × 2000
cells, feasible in well under an hour on 8 cores.

## Master tables (`master_table`)

* `tables/master_perturbation_table.csv` (`uns['master_perturbation_table']`): one row per target,
  the reference `perturbation_meta` columns first — `n_cells`, efficacy (`target_log2fc`,
  `target_pct_kd`, `target_fdr`, `is_effective_hit`; primary control arm), penetrance (`ps_mean`,
  `ps_median`, `ps_responder_fraction`, `ps_net_responder_fraction`, `ps_escaper_fraction`),
  topology (`lochness_mean`, `lochness_median`, `lochness_peak`, `lochness_pct_enriched`),
  phenotype magnitude (`energy_distance`, `mmd_distance`, `distance_pvalue`, `distance_fdr`,
  `distance_significant`), `cofunctional_module`, `phenotype_module` — then the extension columns
  `phenotype_nearest_neighbor`, `phenotype_nearest_distance`, `PCoA1`, `PCoA2`, `strongest_cluster`,
  `cluster_enrichment_log2_or`, `cluster_enrichment_fdr`, `cluster_composition_shift_pct`,
  `n_significant_clusters`, `strongest_protein`, `strongest_protein_effect`, `strongest_protein_fdr`,
  `n_proteins_significant`, `protein_effect_<protein>`, `protein_fdr_<protein>`. Sorted by energy
  distance. No composite score is ever formed (reference principle).
* `tables/master_perturbation_protein_table.csv`: one row per target × protein with the level of
  each statistic stated: the protein effect (`protein_effect`, Cohen's d, p, FDR, call, sample /
  guide sign consistency; `protein_effect_level = TARGET_LEVEL`, perturbed vs control cells of
  that target), the within-target PS ↔ protein Spearman (`ps_protein_rho`, `ps_protein_fdr`;
  `ps_protein_level = CELL_LEVEL`), and the target's master-table columns repeated on each of its
  rows (`target_columns_level = TARGET_LEVEL`).

## Distance ↔ protein (extension; `master_table.protein_associations`)

* `perturbation_distance/distance_protein_association.csv`: per protein, Spearman across targets
  between the energy distance from control and the protein effect, for the **signed** effect and
  for **|effect|** separately (two BH families; `min_targets` 5; `support` flags < 10 targets).
  Figure `distance_vs_protein`.
* `phenotype_space_protein_mantel.csv`: **Mantel test** — Spearman between the target × target
  energy distances and the target × target distances between protein-effect vectors (all proteins
  jointly = Euclidean over the effect vector; and one protein at a time = |Δeffect|), significance
  by `mantel_permutations` (999) seeded permutations of the target labels, BH across the
  single-protein rows. Chosen because it uses the whole pairwise geometry rather than a summary of
  it. `phenotype_module_protein.csv` / `phenotype_module_protein_means.csv`: the module-wise view
  — Kruskal–Wallis of each protein's effect across phenotype modules with ≥ 2 targets, and the
  module × protein means. Figures `phenotype_space_protein`, `phenotype_module_protein`.
* Levels: distance ↔ protein is `TARGET_LEVEL` (n = targets; the distance is never copied onto
  cells), the Mantel test is `TARGET_PAIR_LEVEL` (pairs share targets, hence permutation
  significance), the module test is `PHENOTYPE_MODULE_LEVEL` (modules with fewer than
  `min_targets_per_module` = 3 targets are excluded and counted). The explicit tables
  `tables/distance_protein_associations.csv`, `tables/phenotype_module_protein_associations.csv`,
  `tables/rna_protein_geometry_concordance.csv` carry an `analysis_level` column
  (docs/RNA_PROTEIN_LEVELS.md).
* With a 4-antibody panel and ~20 targets these are **underpowered, descriptive associations**;
  the tables and the report say so (`support`, `power_note`). None of them establishes mediation or causality.

## Outputs

Tables under `tables/perturbation_distance/` (`perturbation_distance`, `distance_skipped`,
`perturbation_distance_matrix`, `perturbation_space_coordinates`, `perturbation_neighbors`,
`phenotype_modules`, `distance_space_skipped`, `distance_protein_association`,
`distance_protein_targets`, `phenotype_space_protein_mantel`, `phenotype_module_protein`,
`phenotype_module_protein_means`) and the two master tables at the top level of `tables/`.
Figures under `figures/perturbation_distance/{distance, distance_space, distance_protein}/`.
Object: `uns['perturbation_distance']`, `uns['perturbation_distance_matrix']` (targets × targets),
`uns['phenotype_modules']`, `uns['master_perturbation_table']`; provenance in
`uns['petrubseq_protein']['analysis']['perturbation_effects']['distance' | 'distance_space' | 'master']`.
Report: sections *Perturbation distance vs control*, *Perturbation distance space & phenotype
modules*, *Master perturbation table* after the modules section; the distance ↔ protein results
in *RNA–protein concordance*; the protein master table in *Multimodal summaries*. Parity harness:
`scripts/reference_parity.py` (section DISTANCE) with `config/parity_reference_papalexi_distance.yaml`
/ `config/parity_papalexi_current_distance.yaml`.
