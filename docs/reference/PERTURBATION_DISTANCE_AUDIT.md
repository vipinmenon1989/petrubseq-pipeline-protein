# Perturbation-distance reference audit (stages 10–12 of the reference pipeline)

Audited 2026-09-26 from the group's local reference checkout
`/local/projects-t3/lilab/vmenon/PertTF-Virtual-Challeng-Weilab/perturbseq-pipeline`, branch
`main`, commit `1c48f9d` (2026-09-24, "Merge pull request #6 … feature/unified-perturbseq-pipeline").
The distance stages entered the reference in commit `590f24f` (2026-09-23, "Integrate scalable
execution and optional perturbation modules from dev"); commit `f2f5019` added the Pertpy
attribution to the README. The earlier parity audit (`REFERENCE_PIPELINE_COMPLETE_AUDIT.md`,
`PARITY_RESULTS.md`) was made against the pre-merge working tree; re-running the parity harness
against the reference at `1c48f9d` (`../results/parity/reference_papalexi_distance/`) gives
**96 MATCH / 0 DIFF / 0 MISSING** for every previously compared quantity, so the merge changed
none of the stages already reproduced.

## 1. Where the reference implements the distance analyses

| file | what it holds |
|---|---|
| `src/perturbseq_pipeline/distance.py` (1,006 l.) | `energy_distance_from_cdist`, `compute_energy_distance`, `compute_mmd`, `_resolve_representation`, `_sample_cell_indices`, `distance_test_permutation`, `_eval_target_dist_worker`, `compute_perturbation_distance` (→ `DistanceResults`), `compute_pcoa_coordinates`, `_eval_pair_dist_worker`, `compute_distance_space` (→ `DistanceSpaceResults`) |
| `src/perturbseq_pipeline/meta.py` (259 l.) | `build_perturbation_meta` → `tables/perturbation_meta.csv` (the MASTER PERTURBATION META TABLE) |
| `src/perturbseq_pipeline/cli.py` 2094–2360 | Stage 10 `perturbation distance vs control`, Stage 11 `perturbation distance space`, Stage 12 `master perturbation meta table & plots`; tables written: `perturbation_distance.csv`, `distance_skipped.csv`, `perturbation_distance_matrix.tsv`, `perturbation_space_coordinates.csv`, `perturbation_neighbors.csv`, `phenotype_modules.csv`, `distance_space_skipped.csv`, `perturbation_meta.csv`; then `plots_mod.plot_distance_figures` and `plot_distance_space_figures` |
| `src/perturbseq_pipeline/plots.py` 8372–8870 | `_significance_stars`, `plot_perturbation_atlas`, `plot_ps_vs_distance`, `plot_perturbation_space`, `plot_module_concordance`, `plot_distance_overview`, `plot_distance_figures`, `plot_distance_space_figures`; sections `SECTION_DISTANCE = "distance"`, `SECTION_DISTANCE_SPACE = "distance_space"` |
| `src/perturbseq_pipeline/templates/report.html` 565–624, `report.py` 361–383 | report sections **8. Perturbation distance vs control** (cards: targets tested, significant, max distance, primary metric; 8.1 overview and atlas figures; 8.2 `tables.perturbation_distance`), **9. Perturbation distance space & phenotype modules** (cards: perturbations mapped, phenotype modules, PCoA dimensions, linkage; 9.1 figures; 9.2 `phenotype_modules`; 9.3 `perturbation_neighbors`), **10. Master perturbation meta table** |
| `src/perturbseq_pipeline/config.py` 833–930, 2239–2400 | `DistanceConfig`, `DistanceSpaceConfig`, `MetaAnalysisConfig`, `VisualizationConfig` and their validation |
| `src/perturbseq_pipeline/data_access.py::get_embedding` | `obsm[rep]` as **float32** (fallback `X_pca_harmony`, then `X_pca`) |
| `src/perturbseq_pipeline/compute.py::derive_seed`, `run_parallel` | sha256-derived per-target seeds; order-preserving joblib/loky execution |
| `tests/test_distance.py` (16 tests), `tests/test_meta.py` (3) | mathematical properties, permutation test, sampling, skipping, symmetry / zero diagonal, PCoA, neighbours, pipeline toggles, meta merge, lean-h5ad invariant |
| `docs/methods.md` §10–12, `README.md` "9 · Perturbation distance", "10 · Perturbation distance space", "11 · Master perturbation table" | prose; the README states the analyses were "conceptually inspired by … scverse Pertpy" |
| `config/default.yaml`, `config/examples/all_modules.yaml` | defaults (both stages **off**; meta on) and the example that turns them on |
| `results/replogle/` (K562 essential, 310,385 cells) | the reference's own real run with `distance.enabled: true` (`min_cells 10`, `secondary_metric null`, `random_seed 0`, 5,000 sampled controls, 2,049 targets tested in 5.5 h, `distance_space` off) — the source of the "8. Perturbation distance vs control" section seen in the old report |

**pertpy is not called at run time.** `grep -rn pertpy src/ tests/` finds nothing; the only
mention is the declared dependency `pertpy>=0.8.0` in `src/perturbseq_pipeline.egg-info` and the
README attribution. `pertpy` 1.0.3 is installed in the `perturbseq-pipeline` conda environment and
not in `petrubseq-protein`; nothing in either pipeline imports it. Every statistic below is the
reference's own code, so "reproduce the reference" means reproducing `distance.py`, not pertpy.

## 2. The reference statistics, line by line

### 2.1 Control distance (`compute_perturbation_distance`, stage 10)

* **representation**: `cfg.distance.representation` (`X_pca`), read as float32; `n_pcs` = all
  stored PCs (50 in the parity run; 30 in the demo config).
* **groups**: `all_targets = sorted(set(target_gene[perturbation_class == "targeting"]))`; the
  control is `perturbation.primary_control` (`ntc` = the `non-targeting` class) when it has
  ≥ `min_cells` cells, else `ntc`, else `other` (= **all** targeting cells), else any non-empty
  one; targets with < `min_cells` (30) cells are written to `distance_skipped.csv`.
* **sampling**: `_sample_cell_indices(indices, max_cells, rng, strata)` — nothing happens when
  `n ≤ max_cells`; otherwise proportional per stratum (`round(max_cells · n_s / n)`, top-up, sort),
  `rng.choice(replace=False)`. Controls: `rng = default_rng(random_seed)` (123) once; target *i*
  (index in `all_targets`, skipped targets included): `derive_seed(random_seed, f"{i}_{target}")` =
  `int(sha256("123_i_target")[:8], 16) mod (2^31 − 1)`. Strata: `cfg.distance.stratify_by`, else the
  input lane column `lane_id` when present (one lane in the parity run → one stratum → uniform).
* **energy distance**: `2·mean(d_XY) − mean(d_XX) − mean(d_YY)` over full Euclidean `cdist`
  matrices in float64, diagonal zeros included, clipped at 0 — the pertpy `Edistance.__call__`
  formula exactly (pertpy 1.0.3 `_distances.py` 667–671).
* **DistanceTest**: the pooled `N × N` matrix is built once; for `B = n_permutations` (1000)
  seeded permutations (`default_rng(target_seed)`, `rng.permutation(arange(N))`, the smaller
  group re-indexed) the statistic is recomputed from row sums; `p = (1 + #(stat ≥ obs − 1e-12)) / (1 + B)`;
  `B = 0` → NaN.
* **secondary metric**: `mmd` = squared MMD with a Gaussian RBF kernel, `gamma = 1 / median(positive squared distances)`.
* **FDR**: `benjamini_hochberg` (NaN-preserving) across tested targets; `significant = fdr < 0.05`;
  table sorted by `energy_distance` descending: `target_gene, n_cells, n_control, energy_distance, pvalue, [mmd_distance], fdr, significant`.
* `primary_metric` is recorded but not used in the computation (energy distance is always the
  tested statistic); `n_jobs` only parallelises across targets (seeds are per target, results are
  order-preserving, so the output does not depend on it).

### 2.2 Pairwise distance space (`compute_distance_space`, stage 11)

* same targets / representation; `min_cells` 30; per-target seed `(random_seed + 43·i) mod (2^31 − 1)`
  with *i* again over the full sorted target list; `max_cells_per_target` 2000; strata `lane_id`.
* `K × K` matrix of `compute_energy_distance` (or `compute_mmd`) over the sampled cells, `M[i, j] = M[j, i]`, zero diagonal; written as `perturbation_distance_matrix.tsv` (never into the h5ad: "lean H5AD" policy, `tests/test_meta.py::test_h5ad_storage_invariants`).
* **PCoA**: `D_sq = M²`, `B = −½ (D_sq − row means − column means + grand mean)`, `numpy.linalg.eigh`,
  eigenvalues sorted descending, strictly positive ones (`> 1e-10`) kept, `p = min(n_components (10), n_pos)`,
  `coords = V[:, :p] · sqrt(λ[:p])` → `perturbation_space_coordinates.csv` (`PCoA1..PCoAp`; eigenvector signs arbitrary).
* **nearest neighbours**: per target, `argsort` of its matrix row with self = ∞, `k = min(nearest_neighbors (10), K − 1)` → `target, neighbor, distance, rank`.
* **phenotype modules**: `scipy.cluster.hierarchy.linkage(squareform(M), method="average")`,
  `fcluster(maxclust, t = n_modules)` when set, else `criterion="distance"` at `cluster_distance_threshold`,
  else the default `t = max(2, min(9, K // 4 if K ≥ 8 else K))`; labels `PM1..`.

### 2.3 Master perturbation meta table (`build_perturbation_meta`, stage 12)

Union of the targets of the perturbation-strength table, PS summary, lochNESS summary, distance
table, co-functional modules and phenotype modules; columns
`n_cells` (from `n_perturbed`, else PS `n_perturbed_cells`, else distance `n_cells`),
`target_log2fc`, `target_pct_kd`, `target_fdr` (`ks_fdr_<primary>`), `is_effective_hit`,
`ps_mean`, `ps_median`, `ps_responder_fraction` (`pct_successful_kd`), `ps_net_responder_fraction`,
`ps_escaper_fraction`, `lochness_mean` (`mean_lochness_in_own_cells`), `lochness_median`,
`lochness_peak` (`max_lochness`), `lochness_pct_enriched` (`pct_own_cells_enriched`),
`energy_distance`, `mmd_distance`, `distance_pvalue`, `distance_fdr`, `distance_significant`,
`cofunctional_module`, `phenotype_module`; sorted by `energy_distance` descending (NaN last) when
present, else by target; explicitly **no composite score**.

### 2.4 Figures (reference `plots.py`)

| figure (section) | data transformation |
|---|---|
| `perturbation_distance_ranking` (distance) | top 40 targets by energy distance, bar colour = `significant` |
| `perturbation_atlas` (distance) | meta table sorted by energy distance, top `atlas_top_n` (50); columns `−target_log2fc` (KD strength), `ps_median`, `lochness_mean`, `energy_distance`; each column z-scored over its finite values (sd < 1e-8 → 1); stars from `target_fdr` and `distance_fdr` (* < 0.05, ** < 0.01, *** < 0.001); side strips for `cofunctional_module` (tab20) and `phenotype_module` (Set2); `vlag` colormap, symmetric limit `max(2.5, max|z|)` |
| `ps_vs_distance_map` (distance) | x `energy_distance`, y `ps_median` (else responder fraction); size `40 + 160·(lochness_mean − min)/(range + 1e-6)`; colour phenotype (else co-functional) module; significant targets full, others faded; 12 largest distances labelled |
| `perturbation_phenotype_space` (distance_space) | PCoA1 × PCoA2 twice: coloured by phenotype module and by energy distance (else PS median) |
| `module_concordance` (distance_space) | crosstab `cofunctional_module × phenotype_module` over targets that have both (≥ 4, ≥ 2 levels each), ARI / NMI in the title |

## 3. Comparison with pertpy 1.0.3 (installed in the reference environment)

| element | pertpy 1.0.3 | reference `distance.py` | classification |
|---|---|---|---|
| energy distance | `Edistance.__call__`: `2·δ − σ_X − σ_Y`, means of full pairwise matrices | identical formula | MATCH |
| permutation test | `DistanceTest`: unseeded `np.random.default_rng()`, `p = clip(#(perm − obs > 0), 1, ∞) / n_perms`, `multipletests(method="holm-sidak")`, unbounded groups | seeded per-target permutations, `p = (1 + #(perm ≥ obs)) / (1 + B)`, BH-FDR, bounded stratified sampling | REFERENCE_CUSTOMIZATION |
| MMD | `MMD` linear kernel by default (rbf / poly options, fixed `gamma`) | RBF with median heuristic | REFERENCE_CUSTOMIZATION |
| pairwise matrix | `Distance.pairwise` (all cells per group) | own `K × K` loop over sampled cells | REFERENCE_CUSTOMIZATION |
| phenotype space / modules / nearest perturbations | no PCoA, no `DistanceSpace.compute`, no `obsp["distances"]`, no `nearest_perturbations` / `plot_similarity` in 1.0.3 (the `PerturbationSpace` classes cover centroid / pseudobulk / clustering / discriminator spaces) | own PCoA, `fcluster`, `argsort` neighbours | REFERENCE_CUSTOMIZATION (no pertpy counterpart) |

Conclusion: the interfaces named in the task brief (`DistanceSpace.compute`, `obsp["distances"]`,
`nearest_perturbations`, `plot_similarity`) exist neither in the installed pertpy nor in the
reference; the reference is self-contained, and the objects to reproduce are its own functions.
Because pertpy is not executed by either pipeline, no `API_VERSION_CHANGE` can arise; the current
implementation ports `distance.py` / `meta.py` and reproduces their numbers on the same cells
(section 5).

## 4. Representation audit (Papalexi parity object)

* `obsm["X_pca"]`: 50 PCs of the scaled 3,000-HVG block (arpack, seed 0); stored float32 by
  both pipelines and identical to `5.8e-11` after sign alignment (`PARITY_RESULTS.md`).
* cells entering the distance: every single-guide targeting cell of a target (reference
  `targeting`, ours `single_targeting`; 1,516 cells, 25 targets) versus the 199 non-targeting
  cells; ambiguous (85) and unassigned cells enter no group.
* no sampling is triggered on 1,800 cells (largest target 102 cells < 2,000; 199 controls < 5,000),
  so the sampling code paths are exercised by unit tests only (`tests/test_perturbation_distance.py`).
* the reference's implicit stratification column `lane_id` holds a single value (`papalexi`) for
  the parity run; the current pipeline stratifies only when `stratify_by` names an obs column.

## 5. Numerical parity on the same 1,800 cells

Reference run: `perturbseq-pipeline run --config config/parity_reference_papalexi_distance.yaml`
(reference at `1c48f9d`, SLURM job 20849823, 111 s) → `../results/parity/reference_papalexi_distance/`.
Current run: `petrubseq-protein run --config config/parity_papalexi_current_distance.yaml` →
`../results/parity/current_papalexi_distance/`, compared with `scripts/reference_parity.py`
(section DISTANCE). The full parity table and the numbers are in `PARITY_RESULTS.md`
("Perturbation distance stages").

Direct check of the ported functions on the current parity object against the reference tables
(21 targets tested, 4 skipped below 30 cells; 21 × 21 matrix, 10 PCoA axes, 10 neighbours, 5
phenotype modules):

| quantity | max abs. difference | note |
|---|---|---|
| energy distance vs control | 1.5e-11 | float32 PCs → float64 cdist on both sides |
| DistanceTest p-value | 9.9e-17 | identical permutation counts (same seeds, same generator) |
| MMD | 3.2e-13 | |
| BH-FDR / significant | 1.1e-16 / 21/21 equal | 7 significant at FDR < 0.05 |
| pairwise matrix | 2.7e-11 | |
| PCoA coordinates (sign-aligned) | 2.1e-11 | eigenvalues 55.6, 13.4, 6.6, … |
| nearest neighbours | identical targets / ranks; distances 2.0e-11 | |
| phenotype modules | identical labels (PM1–PM5) | reference default cut, 21 // 4 = 5 |

Classification: MATCH (no `API_VERSION_CHANGE`, `METRIC_DEFINITION_CHANGE`,
`REPRESENTATION_DIFFERENCE`, `PERMUTATION_DIFFERENCE` or `BUG`; the reference's own deviations
from pertpy are `REFERENCE_CUSTOMIZATION` and are reproduced as such).

## 6. Vocabulary map

| reference | current |
|---|---|
| `obs["target_gene"]`, class `targeting` / `non-targeting` | `obs["target"]`, `single_targeting` / `single_control` with `control_class` in `control_classes` |
| `cfg.distance.*`, `cfg.distance_space.*`, `cfg.meta_analysis.enabled`, `cfg.visualization.*` | `analysis.perturbation_effects.distance.*`, `.distance_space.*`, `.master_table.enabled`, `.master_table.{perturbation_atlas, ps_distance_map, perturbation_space, module_concordance, atlas_top_n}` |
| `tables/perturbation_distance.csv` (`target_gene`) | `tables/perturbation_distance/perturbation_distance.csv` (`target`) |
| `tables/perturbation_distance_matrix.tsv` | `tables/perturbation_distance/perturbation_distance_matrix.csv` and `uns["perturbation_distance_matrix"]` (the task asked for the matrix in `uns`; the reference keeps it out of the h5ad) |
| `perturbation_space_coordinates.csv`, `perturbation_neighbors.csv`, `phenotype_modules.csv`, `distance_skipped.csv`, `distance_space_skipped.csv` | same names under `tables/perturbation_distance/` |
| `tables/perturbation_meta.csv` | `tables/master_perturbation_table.csv` (reference columns first, in the reference order and names, then the extension columns) |
| figures `figures/distance/*`, `figures/distance_space/*` | `figures/perturbation_distance/{distance, distance_space, distance_protein}/*` (same names) |
| report sections 8 / 9 / 10 | "Perturbation distance vs control", "Perturbation distance space & phenotype modules", "Master perturbation table" (numbered dynamically after the modules section) |
