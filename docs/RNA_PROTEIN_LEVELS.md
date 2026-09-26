# Statistical level of every RNA ↔ protein analysis

Audited 2026-09-26 (`analysis/concordance.py`, `analysis/distance.py`, `analysis/master_table.py`,
`analysis/protein_effects.py`, the pipeline orchestration, report, plots and tests). The point of
this document is auditability: for every cross-modal statistic, which quantity is compared with
which, at which level each quantity exists, what the unit of observation of the test is, and
whether any value is repeated across cells. Every output table carries an `analysis_level`
column with one of these values:

| level | unit of observation | example |
|---|---|---|
| `CELL_LEVEL` | one cell | per-cell PS vs per-cell CLR protein value |
| `TARGET_LEVEL` | one perturbation target | energy distance to control vs target-level protein effect |
| `TARGET_PROTEIN_LEVEL` | one target × protein pair | the protein effect itself (rows of `protein_effects.csv`; computed from that target's cells vs controls, i.e. a target-level estimate for one protein) |
| `TARGET_PAIR_LEVEL` | one pair of targets | Mantel test of RNA pairwise geometry vs protein-profile geometry |
| `PHENOTYPE_MODULE_LEVEL` | one phenotype module (group of targets) | Kruskal–Wallis of protein effects across modules |

## The audit table

| analysis | RNA quantity (level) | protein quantity (level) | unit of the test | test | n | any value repeated across cells? | valid? | table |
|---|---|---|---|---|---|---|---|---|
| PS ↔ protein | per-cell PS of the cell's own target, `obs['ps_score']` (CELL) | per-cell CLR value `obsm['protein'][p]` (CELL) | cell, within one target's single-guide cells | Spearman ρ, two-sided p, BH over all within-target (target × protein) tests; `min_cells` 20 | that target's cells | no | **yes** | `tables/ps_protein_associations.csv`, `perturbation_effects/ps_protein_association.csv` (also holds the pooled-all-targets rows, `CELL_LEVEL_POOLED_DESCRIPTIVE`, never tested) |
| lochNESS ↔ protein, cell level | per-cell own-target lochNESS `obsm['lochness'][t]` (CELL; the kNN enrichment of target t around that cell, varies cell by cell) | per-cell CLR value (CELL) | cell, within one target's cells | Spearman, BH over within-target tests, `min_cells` 20 | that target's cells | no | **yes** | `tables/lochness_protein_associations.csv`, `perturbation_effects/lochness_protein_association_cells.csv` |
| lochNESS ↔ protein, target level | own-cell mean lochNESS per target (TARGET, a summary of the cell-level scores) | protein effect `E_t,p` (TARGET_PROTEIN) | target | Spearman across targets per protein, signed and |effect|, BH across proteins per family | targets | no (one row per target) | yes, descriptive companion | `perturbation_effects/lochness_protein_summary.csv`, `..._association_targets.csv` |
| gene program ↔ protein, cell level | per-cell program activity `obsm['program_activity']` (CELL) | per-cell CLR value (CELL) | cell, over all single-guide + control cells | Spearman, BH over program × protein | cells | no | yes (association across a mixed population; not within target) | `perturbation_effects/program_protein_association_cells.csv` |
| gene program ↔ protein, target level | target × program effect (mean log2FC of program genes; TARGET) | `E_t,p` (TARGET_PROTEIN) | target | Spearman across targets per program × protein, BH | targets (≥ 5) | no | yes | `perturbation_effects/program_protein_association_targets.csv` |
| protein effect | – | CLR difference of means, perturbed vs control cells of target t, Cohen's d, Mann–Whitney p, BH over target × protein | cell (the two-sample test) → one estimate per target × protein | Mann–Whitney U | cells of t vs controls | no | yes | `perturbation_effects/protein_effects.csv` (level `TARGET_LEVEL` in the master protein table: one estimate per target) |
| energy distance ↔ protein | energy distance to control `D_t` (TARGET: one number per target) | `E_t,p` (TARGET_PROTEIN) | **target** (n = targets, never cells) | Spearman across targets per protein; signed effect and |effect| as two BH families; `min_targets` 5, `support` flags < 10 targets | targets | **no** — `D_t` is never copied onto cells | **yes** | `tables/distance_protein_associations.csv`, `perturbation_distance/distance_protein_association.csv`, `..._targets.csv` (one row per target × protein) |
| phenotype module ↔ protein | phenotype module membership (TARGET; from the pairwise RNA distance matrix) | `E_t,p` (TARGET_PROTEIN) | target grouped by module | Kruskal–Wallis across modules with ≥ `min_targets_per_module` (3) targets (two eligible groups = Mann–Whitney), BH across proteins; singletons / small modules excluded and counted | eligible targets | no | yes when ≥ 2 eligible modules; otherwise `insufficient_support` | `tables/phenotype_module_protein_associations.csv`, `perturbation_distance/phenotype_module_protein_means.csv` |
| RNA geometry ↔ protein geometry | pairwise energy distance matrix `R_ij` (TARGET_PAIR) | protein-profile distance `P_ij` = Euclidean distance between the CLR-effect vectors of targets i and j (joint), or `|E_i,p − E_j,p|` (one protein) (TARGET_PAIR) | target pair (upper triangle, each pair once) | **Mantel**: Spearman of the upper triangles; significance by seeded permutation of the protein matrix's target labels, `p = (1 + #(ρ_perm ≥ ρ_obs)) / (1 + B)`; BH across the single-protein rows; the joint row stands alone | K(K−1)/2 pairs, which are *not* independent (they share targets), hence the permutation p rather than a naive pair-level p | no | **yes** | `tables/rna_protein_geometry_concordance.csv` |
| PCoA coordinates ↔ protein | PCoA1/2 per target (TARGET) | – | – | **not tested**: the coordinates are drawn coloured by protein effect only (`phenotype_space_protein` figure); no correlation of coordinates copied onto cells | – | never copied to `obsm` | – | – |

Checks performed on the code (2026-09-26):

* `distance_protein_association` joins the distance table (one row per target) with
  `protein.matrix` (targets × proteins) on the target key and calls Spearman on those target
  vectors; `n_targets` is reported; no cell arrays enter. `test_rna_protein_levels.py` duplicates
  every cell and shows the target-level statistics unchanged while the cell-level PS ↔ protein
  p-values shrink (the expected behaviour of a genuinely cell-level test).
* `mantel_test` requires two K × K matrices over the same targets in the same order (the caller
  builds both from one `common` list), checks symmetry and zero diagonals, uses each pair once,
  permutes target labels with `numpy.random.default_rng(seed)` (rows and columns together), and
  returns NaN rather than a naive p when no permutations are requested.
* The energy distance to control, the pairwise matrix, the PCoA coordinates and the phenotype
  modules are stored separately (`uns['perturbation_distance']`, `uns['perturbation_distance_matrix']`,
  `tables/perturbation_distance/perturbation_space_coordinates.csv`, `uns['phenotype_modules']`);
  nothing target-level is written to `obsm`.
* Master tables: `tables/master_perturbation_table.csv` is one row per target (`TARGET_LEVEL`;
  PS / lochNESS entries are target summaries of the cell-level analyses).
  `tables/master_perturbation_protein_table.csv` is one row per target × protein and every group
  of columns carries its level column (`protein_effect_level`, `ps_protein_level`,
  `lochness_protein_level`, `distance_protein_level`, `phenotype_module_protein_level`,
  `phenotype_geometry_protein_level`, `program_protein_level`, `target_columns_level`). Per-protein
  summaries (distance ↔ protein, module test, single-protein Mantel) are repeated on each of that
  protein's rows and labelled; the joint all-protein Mantel statistic is only in the geometry table.

## What was wrong before 2026-09-26 and what changed

* No pseudoreplication was found: the distance ↔ protein and geometry tests were already
  target- and target-pair-level. The demo values are therefore unchanged by design
  (PD-L1 signed ρ −0.59 / |effect| ρ 0.60, joint Mantel ρ 0.54 on the demo).
* lochNESS ↔ protein existed only at the target level; the cell-level within-target analysis was
  added (`lochness_protein_associations.csv`) and the target-level one kept and labelled.
* `analysis_level` columns, the explicit top-level tables, the module-eligibility rule
  (`min_targets_per_module`, default 3, replacing "≥ 2"), the matrix checks in the Mantel test,
  `n_guides` and the strongest gene program in the master table, and the full set of level
  columns in the master protein table were added.

## Limitations

Four antibodies (Papalexi) make every target-level and pair-level association low-dimensional
and underpowered; the `support` / `power_note` fields say so and nothing here is causal.
