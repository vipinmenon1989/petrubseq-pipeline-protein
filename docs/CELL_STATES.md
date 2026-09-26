# Cell states and perturbation × cluster enrichment (Stage F)

Part of the standard run (`analysis.clustering.enabled`, default **on**; the
SCP1064 regression configs set it to `false` to keep the frozen v0.1 reference
comparable). Clustering is the first analysis stage after the representations, as in
the reference pipeline, because the module gene panel (Leiden markers) and the lochNESS
cluster summary consume the clusters; enrichment runs after the perturbation-strength
test. Reference audit: `docs/reference/REFERENCE_PIPELINE_COMPLETE_AUDIT.md`.

## Three complementary views of cell state

| question | analysis | needs clusters? |
|---|---|---|
| What discrete transcriptional states exist? | Leiden clustering | defines them |
| Does a perturbation preferentially occupy a defined state? | perturbation × cluster enrichment | yes |
| Does a perturbation concentrate locally in expression space? | lochNESS | no |

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

Port of the reference `enrichment.py` (audit: docs/reference/REFERENCE_PIPELINE_COMPLETE_AUDIT.md
section C.6; numerical parity: docs/reference/PARITY_RESULTS.md). For every target with
≥ `min_cells_per_target` (10) single-guide cells and every cluster with ≥ `min_cells_per_cluster`
(20) cells, and for **each reference arm** in `controls` (`ntc` = single-guide cells of
`control_classes`; `other` = single-guide targeting cells of every other target):

```
                   in cluster k    not in k
    perturbed t          a             b          a + b = n_target_cells
    reference            c             d          c + d = n_reference_cells
```

* **test**: two-sided Fisher exact test (`pval`, also kept as `pval_fisher`); with
  `stratify_by` set to an obs column with ≥ 2 levels, a Cochran–Mantel–Haenszel test over
  the strata **replaces** the pooled p-value and, when finite, the pooled odds ratio replaces
  the Haldane one (`cmh_odds_ratio`, `cmh_pval`, `cmh_n_strata` are added);
* `odds_ratio` = Haldane-Anscombe ((a+½)(d+½))/((b+½)(c+½)), finite at zero counts;
  `log2_odds_ratio`; `sample_odds_ratio` = a·d/(b·c) (extra column);
* `pct_of_target`, `pct_of_reference`; `direction` enriched when the target's share exceeds
  the reference's, else depleted; `low_power` = fewer than `min_reference_cells` (10)
  reference cells in the cluster;
* **multiple testing**: Benjamini–Hochberg **within each arm**; `significant` = FDR <
  `fdr_alpha` under the `primary_control` arm (reference default `other`: the non-targeting
  group is small and contributes few cells exactly in the rare clusters);
* **guide concordance** (significant pairs only): guides of *t* with ≥ `min_cells_per_guide`
  (5) cells whose in-cluster fraction lies on the same side of the reference fraction as the
  direction → `guides_concordant`, `guides_tested`;
* **omnibus**: chi-square of the targeting-cell target × cluster table, % expected counts
  below 5, and a seeded permutation p-value in which every target's row is resampled from a
  multinomial with the pooled cluster proportions (`n_permutations` 1000);
* per target: `composition` (% of the target's cells per cluster), each arm's
  `reference_composition`, and `effect_magnitude` (total variation distance from the primary
  reference composition, number of significant clusters); `phenocopy_similarity` = Pearson
  correlation of the composition profiles between targets (figure).

Outputs (`tables/cell_states/`): `perturbation_cluster_enrichment.csv` (both arms; also
`uns['perturbation_cluster_enrichment']`), `enrichment.csv` (reader view of the primary arm),
`enrichment_composition.csv`, `enrichment_reference_composition.csv`,
`enrichment_effect_magnitude.csv`, `cluster_enrichment_vs_lochness.csv` (descriptive
side-by-side, extension), `cluster_enrichment_skipped.csv`. Figures (report section
*Perturbation enrichment across clusters*): `enrichment_heatmap` (rows ordered by profile
similarity, `*` = FDR < α), `enrichment_phenocopy`, `enrichment_composition` (reference +
top 30 targets by shift), `enrichment_volcano`, `enrichment_effect_magnitude`, and
`enrichment_<target>` (composition vs reference, log2 OR per cluster, the target's cells on
the UMAP; the 12 strongest embedded, all on disk).

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
