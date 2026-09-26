# Perturbation effects (Stage E): methods

Enabled with `analysis.perturbation_effects.enabled: true` (default **off**;
`config/demo_papalexi.yaml` enables it). Runs after the representations, on the
filtered, normalized object, as stage "perturbation effects" of `petrubseq-protein run`.
The behavioural reference is the group's `weili-lab/perturbseq-pipeline`
(`ps_score.py`, `lochness.py`, `modules.py`); the audit is
`docs/reference/PERTURBATION_ANALYSIS_AUDIT.md`, and
`scripts/compare_reference_perturbation.py` reproduces the reference numbers
(PS to 4e-16, lochNESS and the effect matrix exactly, identical program / module
partitions on the same input).

Discrete cell states (Leiden) and perturbation × cluster enrichment are a separate,
independent analysis (Stage F, `analysis.clustering`): `docs/CELL_STATES.md`. lochNESS
below is the cluster-free counterpart.

## Concepts

| term | meaning |
|---|---|
| perturbation target | the gene the guide acts on (`obs['target']`) |
| downstream response genes | genes whose expression changes in cells carrying the perturbation |
| gene program (`P1..`) | a cluster of response genes that move together across perturbations |
| perturbation module (`M1..`) | a cluster of perturbations with similar response profiles |

Cell groups used by every analysis: **perturbed** = `single_targeting` cells of one
target; **controls** = `single_control` cells whose `control_class` is in
`control_classes` (default `non_targeting`). Ambiguous, multi-guide, mixed and
unassigned cells enter neither group. `target_gene_map` maps a target label to its
RNA gene symbol when they differ (dataset-specific; set in the dataset config).

## PS — perturbation-response score

**Reference.** Song B., Liu D., Dai W., et al. *Decoding heterogeneous single-cell
perturbation responses.* Nat Cell Biol 27, 493–504 (2025), doi:10.1038/s41556-025-01626-9;
Python implementation PS_python (`pertps`, https://github.com/weili-lab/PS_python).

**Definition.** Per target: (1) rank genes perturbed vs control with a t-test and
keep the top `top_n_genes` (100) — the response signature; (2) `beta_g` = OLS slope of
gene g on the perturbed indicator (= mean difference); (3) each cell's score is the
projection of its centred expression on beta, `((Y_i − Ȳ)·beta)/(beta·beta)`;
(4) shift so the control mean is 0, clip to `[0, scale_factor]` (3), divide by the
maximum → `[0, 1]`. Quadrants with the target gene's own expression (cut at the
control mean): successful knockdown / escaper / non-responder / low signal. Targets
expressed in < 1 % of control cells are skipped (knockdown unmeasurable).

**Outputs.** `obs['ps_score']` (each perturbed cell's score for its own target),
`obs['ps_quadrant']`, `obsm['ps_scores']` (cells × targets, NaN outside the target's
comparison), `obsm['ps_scores_raw']` (clipped, not max-normalized),
`uns['ps_signatures']`; `tables/perturbation_effects/ps_targets.csv`
(n cells, median / mean / IQR PS, % high PS in perturbed and control cells,
`auc_vs_control`, `d_vs_control`, knockdown quadrants, warnings), `ps_skipped.csv`,
`ps_signatures.csv`.

**Interpretation and limitations.** PS is comparable *within* a target (max-normalized).
By construction the mean unnormalized score of perturbed minus control cells is 1 for
every target, so strength across targets is reported as the separation of perturbed
from control scores (`auc_vs_control`, `d_vs_control`; additions to the reference).
The signature is fitted and applied on the same cells (in-sample), so even a
perturbation without effect separates somewhat (AUC > 0.5). Controls scored against a
target's signature give the background (`pct_high_ps_control`).

## lochNESS — local neighbourhood enrichment

**Reference.** Huang X., Henck J., Qiu C., et al. *Single-cell, whole-embryo phenotyping
of mammalian developmental disorders.* Nature 623, 772–781 (2023),
doi:10.1038/s41586-023-06548-w (lochNESS: "local cellular heuristic neighbourhood
enrichment specificity score"); implementation as ported from pertTF
(https://github.com/davidliwei/pertTF) in the reference pipeline.

**Definition.** For every cell and target g: `lochNESS = local_fraction(g) / overall_fraction(g) − 1`,
with `local_fraction` the share of the cell's k nearest neighbours (in `obsm['X_pca']`,
first `n_pcs` = 20 PCs; never UMAP) that are perturbed cells of g and `overall_fraction`
g's share of all cells. 0 = chance, > 0 enriched, < 0 depleted. k = 300 (reference),
capped at `max_k_fraction` × n cells (10 %) for small objects (the cap is recorded).

**Outputs.** `obsm['lochness']` (cells × targets), `obs['lochness_self']` (each
cell's score for its own label, including the control label);
`lochness_targets.csv`: own-cell mean / median, % own cells enriched,
`mean_lochness_in_control_cells`, `delta_own_vs_control`, label-permutation null
(`n_permutations`, seeded) mean / SD, `z_score`, `p_empirical`, BH `fdr`, samples and
largest sample share; `lochness_by_sample.csv`.

**Interpretation and limitations.** lochNESS is a *population / state* quantity: it says
whether a target's cells occupy a distinct region of transcriptional state. Because the
denominator is all cells, a perturbation without effect is still > 0 wherever control
cells live when other perturbations occupy separate regions, and the permutation null
over all cells flags that too; `delta_own_vs_control` (own-cell mean minus the target's
score in control cells, ≈ 0 for a target that sits with the controls) is the
control-referenced reading. With small objects and large k the dynamic range shrinks.

## Gene programs and perturbation modules

**Reference.** Zhou P., Shi H., Huang H., et al. *Single-cell CRISPR screens in vivo map
T cell fate regulomes in cancer.* Nature 624, 154–163 (2023), doi:10.1038/s41586-023-06733-x
(the reference pipeline's README attributes this design to "Chen et al. 2023"; the DOI it
links resolves to Zhou et al.). Implementation: reference `modules.py`.

**Definition.** Effect matrix `E[t, g] = log2((mean_P + pc) / (mean_C + pc))` on
de-logged means for targets with ≥ 20 cells, over a response-gene panel; Welch t-test on
the log values with BH within each target gives the DE gate (|log2FC| > 0.5, FDR < 0.05).
Gene programs = average-linkage clustering of genes on 1 − Pearson correlation of their
effect profiles (`n_programs` 4, or a cut at distance 0.7); perturbation modules = the
same on targets with 1 − Spearman (`n_modules`, default a cut at 0.7). Module × program
strength = mean log2FC of a program's genes under a module's perturbations.
Panel: the reference uses Leiden markers; this pipeline has no clustering, so
`gene_selection: response` takes the union of each target's top 100 genes by |t| at
FDR < 0.05 (falling back to HVGs), restricted to genes detected in ≥ `min_pct_cells_expressing`
(5 %) of the analysed cells. The pseudocount `pc` = `log2fc_pseudocount` defaults to 1 on the
normalized-count scale (Seurat `FoldChange` convention); the reference uses 1e-9, which gives
|log2FC| > 20 for genes undetected in one small group — on the 1,800-cell demo such genes formed
a separate 7-gene program. `n_modules` defaults to the reference 9 (a dendrogram cut at 0.7
when null). With `log2fc_pseudocount: 1e-9` and `min_pct_cells_expressing: 0` the effect matrix
and both partitions equal the reference exactly. Per-cell program activity = mean z-scored log
expression of the program's genes (deterministic; the reference uses `sc.tl.score_genes`).

**Outputs.** `obsm['program_activity']`, `uns['perturbation_effect_matrix']`,
`uns['gene_programs']`, `uns['perturbation_modules']`; tables `perturbation_effect_matrix.csv`,
`perturbation_de_mask.csv`, `gene_programs.csv` (gene, program, size, mean log2FC,
targets up / down), `perturbation_modules.csv` (target, module, cells, DE genes up / down),
`module_program_strength.csv`, `perturbation_program_effects.csv` (target × program mean
log2FC and fraction DE). Labels are numeric only; nothing is annotated biologically.

## Protein effects

For each target and each protein in `obsm['protein']` (CLR of ADT counts by default;
never counts): perturbed vs control `effect` = difference of means (on CLR a natural-log
fold change of the geometric-mean-normalized abundance), Cohen's d, two-sided
Mann–Whitney U p-value, BH-FDR over all target × protein tests, per-sample sign
consistency (`n_samples_tested`, `n_samples_same_sign`) and per-guide sign consistency
(`n_guides_tested`, `n_guides_same_sign`; guides with ≥ 5 cells). Tables `protein_effects.csv`,
`protein_effect_matrix.csv`; `uns['protein_effects']`. Normalized-only protein inputs
work unchanged.

## RNA–protein concordance

* **PS ↔ protein** (cell level, within target): Spearman rho between PS and protein value
  across a target's perturbed cells (≥ `min_cells` 20), BH-FDR over the within-target tests.
  A pooled row over all targets is descriptive only (it mixes between- and within-target
  variation) — `ps_protein_association.csv`.
* **lochNESS ↔ protein** (target level): lochNESS is a population quantity, so no cell-level
  correlation is computed. Per protein, Spearman across targets between own-cell mean
  lochNESS and the protein effect (and |effect|) — `lochness_protein_summary.csv`,
  `lochness_protein_association.csv`.
* **gene program ↔ protein**: cell level over single-guide and control cells (program
  activity vs protein value) — `program_protein_association_cells.csv`; target level
  (perturbation × program effect vs perturbation × protein effect across targets) —
  `program_protein_association_targets.csv`.
* **Integrated summary** `perturbation_summary.csv` (`uns['perturbation_summary']`):
  per target PS summary and separation, lochNESS and its z / FDR, module, DE genes, RNA
  effect magnitude (mean |log2FC|), the strongest programs, per-protein effect and FDR,
  largest |d| and the best within-target PS ↔ protein pair.

These are associations. None establishes mediation (perturbation → program → protein)
or causality.

## Figures (report section 6, `figures/perturbation_effects/<stage>/`)

PS: median PS per target, PS distributions vs controls. lochNESS: self-enrichment per
target with permutation null, lochNESS maps on the RNA UMAP (display only). Programs:
perturbation × program heatmap, ordered response-gene effect matrix with program borders,
module × program strength. Protein: perturbation × protein Cohen's d heatmap. Concordance:
PS vs protein (≤ 6 pairs by FDR), program activity vs protein (≤ 4 pairs), RNA / lochNESS
vs protein effect magnitude, integrated overview. Panels are chosen deterministically; every
number is in the tables.
