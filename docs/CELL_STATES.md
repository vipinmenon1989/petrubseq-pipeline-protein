# Cell states and perturbation × cluster enrichment (Stage F)

Enabled with `analysis.clustering.enabled: true` (default **off**;
`config/demo_papalexi.yaml` enables it). Runs as the stage "cell states and
perturbation enrichment" after the representations. It is independent of the
perturbation-effect analyses (`docs/PERTURBATION_EFFECTS.md`): none of them uses
the clusters. Reference audit: `docs/reference/CLUSTERING_ENRICHMENT_AUDIT.md`.

## Three complementary views of cell state

| question | analysis | needs clusters? |
|---|---|---|
| What discrete transcriptional states exist? | Leiden clustering | defines them |
| Does a perturbation preferentially occupy a defined state? | perturbation × cluster enrichment | yes |
| Does a perturbation concentrate locally in expression space? | lochNESS (Stage E) | no |

Cluster enrichment is easy to read ("60 % of JAK2 cells sit in cluster 4, versus 0 %
of controls") but depends on where the cluster boundaries fall, i.e. on the Leiden
resolution. lochNESS is cluster-free: it can see a perturbation concentrating inside
one large cluster or straddling two. A perturbation can show both, either or neither.
The two are shown side by side (`cluster_enrichment_vs_lochness.csv`) and never
combined into one score.

## Leiden clustering

Leiden (Traag V.A., Waltman L., van Eck N.J. *From Louvain to Leiden: guaranteeing
well-connected communities.* Sci Rep 9, 5233 (2019)) via `scanpy.tl.leiden` with the
igraph backend, `n_iterations` 2, undirected, `resolution` (1.0), seed `compute.seed`,
on the **existing** RNA neighbour graph `obsp['rna_connectivities']` (15 neighbours on
the RNA PCA). Nothing is re-normalized or re-embedded. Labels are written to
`obs[analysis.clustering.key]` (default `leiden`) as a categorical, numbered from the
largest cluster; an existing obs column of that name is an error unless
`overwrite: true`.

**Population.** Every QC-passing cell in the processed object, the cells the RNA graph
was built on (the reference default). Guide labels do not enter the clustering, so the
states are defined by expression alone. Ambiguous, multi-guide and unassigned cells are
clustered (they stay in the object) but are excluded from the enrichment tests; the
cluster summary reports their share per cluster so a cluster made of them is visible.

**Cluster QC** (`cluster_summary.csv`): size and fraction, cells per perturbation class,
fraction ambiguous/unassigned, median UMIs / genes / % mitochondrial, the dominant level
of each design column (sample, lane, batch, donor, condition) and its share, and the
ratio of the cluster's median library size to the overall median. Descriptive flags
(also listed as run warnings): a design level holding ≥ 80 % of a cluster while
holding < 50 % of all cells; > 50 % ambiguous/unassigned; median library size ≥ 1.5 × or
≤ 0.67 × the overall median. Nothing is removed. Composition tables:
`cluster_composition_by_class.csv`, `…_by_target.csv`, `…_by_<design column>.csv`.

## Perturbation × cluster enrichment

For every target with ≥ `min_cells_per_target` (10) single-guide cells and every cluster
with ≥ `min_cells_per_cluster` (20) cells:

```
                   in cluster k    not in k
    perturbed t          a             b          a + b = n_target_total
    control              c             d          c + d = n_control_total
```

* **perturbed** = single-guide targeting cells of t; **control** = single-guide cells
  of `control_classes` (default non-targeting). `control: other` uses the single-guide
  targeting cells of all other targets instead (explicit alternative; the reference
  pipeline's default). Ambiguous, multi-guide, mixed and unassigned cells are never counted.
* **test**: two-sided Fisher exact test (`scipy.stats.fisher_exact`); `odds_ratio` is the
  sample odds ratio a·d / (b·c) it returns (0 or ∞ when a count is 0 — kept as is);
  `log2_or_haldane` = log2((a+½)(d+½) / ((b+½)(c+½))), finite, for ranking and the heatmap.
* **multiple testing**: Benjamini–Hochberg over **all tested (target, cluster) pairs of
  the run** (one family); `significant` = FDR < `fdr_alpha` (0.05). Odds ratios alone are
  never called significant.
* **direction**: enriched when a/(a+b) > c/(c+d), depleted when lower.
* **low_power**: fewer than `min_control_cells_in_cluster` (10) control cells in the
  cluster. When the cluster holds no controls at all the test can still be decisive (a
  state occupied only by perturbed cells); the flag then says the control reference is
  thin, not that the result is weak — read `n_control_in_cluster`.
* **guide support**: for each guide of t with ≥ `min_cells_per_guide` (5) cells, the guide
  supports the pair when its own in-cluster fraction lies on the same side of the control
  fraction as the target-level direction (so depletions can be supported).
  `n_guides_observed`, `n_guides_supporting_direction`, `guide_support_fraction`. Guides
  are not tested individually and are not independent replicates; this is a consistency
  check against single-guide artefacts.
* **stratified test** (`stratify_by`, any obs column, default none): a Cochran–Mantel–
  Haenszel test across the column's levels (`statsmodels` `StratifiedTable`; strata with
  an empty margin dropped) reported *beside* the Fisher result (`cmh_odds_ratio`,
  `cmh_p_value`, `cmh_fdr`, own BH family), so a cluster that simply differs in size
  between samples or lanes cannot pass as a perturbation effect.
* **omnibus** (`uns[...]['analysis']['cell_states']['enrichment']['omnibus']`): chi-square
  of the perturbed-target × cluster table with a seeded 1,000-permutation p-value, and the
  share of expected counts below 5 (a screen, not the inferential result).

Outputs: `perturbation_cluster_enrichment.csv` (and `uns['perturbation_cluster_enrichment']`),
`cluster_enrichment_vs_lochness.csv` (per target: number of significant clusters, the
strongest enriched and depleted cluster with log2 OR and FDR, lochNESS own-cell mean,
delta and FDR when Stage E ran), `cluster_enrichment_skipped.csv`. Figures (report
section "Cell states and perturbation enrichment"): RNA UMAP by cluster; cluster sizes
with composition by perturbation class and by sample/lane/batch; target × cluster
heatmap (colour = Haldane log2 OR, * = FDR < alpha, o = low power).

## Limitations

* The number and boundaries of clusters depend on the resolution; enrichment results are
  conditional on that choice. Compare with lochNESS, which has no boundaries.
* Clusters are built from the same cells that are tested. Guide labels are not used, so
  there is no label leakage, but strong perturbations shape the clusters they are then
  enriched in — that is the phenotype being described, not an independent validation.
* Cells are treated as independent in Fisher's test; sample / lane structure is only
  addressed by the optional stratified test and the composition flags.
* Small perturbation groups and clusters with few controls limit power; untested clusters
  are listed in the provenance.
* No batch correction is applied (it could remove real perturbation biology); clusters
  dominated by one sample, lane or library depth are flagged rather than corrected.
* Clusters are numbered, never named.
