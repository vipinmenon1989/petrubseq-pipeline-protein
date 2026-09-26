> **Superseded (2026-09-26)** by `REFERENCE_PIPELINE_COMPLETE_AUDIT.md` (complete stage-by-stage audit of the reference repository) and `PARITY_RESULTS.md` (numerical parity after the corrections). Kept for the history of the earlier, partial audits; where they disagree, the complete audit is authoritative.

# Cell-state clustering and perturbation × cluster enrichment: reference audit (Stage F)

Audited 2026-09-26 from the group's local checkout of weili-lab/perturbseq-pipeline
(`src/perturbseq_pipeline/cluster.py`, `enrichment.py`, `config.py` `ClusterConfig` /
`EnrichmentConfig`, `docs/methods.md` sections 4 and 6), before implementation.

## What the reference does

**Clustering** (`cluster.py`). Library-size normalization, log1p, 3,000 HVGs, PCA
(50, arpack), optional Harmony on `batch_key` (`X_pca_harmony`), 15-neighbour graph,
UMAP (`min_dist` 0.5), then `sc.tl.leiden(flavor="igraph", n_iterations=2,
directed=False, resolution=1.0, random_state=run.seed)` into `obs['leiden']`.
`cluster.assigned_only` (default **false**) re-runs PCA / graph / Leiden on
guide-assigned singlets only, to stop droplet multiplets fragmenting the embedding.

**Enrichment** (`enrichment.py`). For every (target, cluster) with ≥ 10 target cells
and ≥ 20 cluster cells, the 2 × 2 table [[target in, target out], [reference in,
reference out]] is tested with a two-sided Fisher exact test; BH-FDR over all pairs
within one control arm. Two reference arms are always computed: `ntc`
(non-targeting cells) and `other` (targeting cells of all other targets);
**`other` is the default** (`primary_control: other`) because on its demo lane
several clusters held only 3–4 NTC cells. Haldane-Anscombe odds ratio (+0.5 per
cell); pairs with < 10 reference cells in the cluster are flagged `low_power`.
Guide concordance: a guide with ≥ 5 cells agrees when its in-cluster fraction is on
the observed side of the reference fraction (direction-aware, so depletions are
supported too). Optional `stratify_by` replaces Fisher with Cochran–Mantel–Haenszel
(`statsmodels` `StratifiedTable`). An omnibus chi-square with a 1,000-permutation
p-value screens the whole target × cluster table. Figures: UMAP by cluster,
target × cluster heatmap, per-target bars.

## What is scientifically useful

Discrete states are interpretable and complement the cluster-free lochNESS; the
Fisher 2 × 2 design, the minimum sizes, the finite odds ratio for display, the
direction-aware guide agreement, the optional stratified test and the omnibus
screen are all sound and cheap.

## Ported / changed / not ported

| item | decision |
|---|---|
| Leiden (igraph, 2 iterations, resolution 1.0, seed) | ported, on the **existing** RNA graph (`obsp['rna_*']`) — no second preprocessing |
| clustering population | all QC-passing cells (reference default `assigned_only: false`); re-clustering a subset would need a new graph |
| Fisher 2 × 2, BH over all pairs, min 10 target / 20 cluster cells, low-power < 10 reference cells | ported |
| default reference | **changed to non-targeting controls** (project rule for every perturbation analysis); `control: other` kept as an explicit alternative |
| odds ratio | the sample odds ratio from Fisher is stored as is (0 / inf possible); the Haldane-corrected log2 ratio is a separate display/ranking column |
| guide agreement | ported (direction-aware, guides ≥ 5 cells), reported for every pair |
| CMH (`stratify_by`) | ported as an **additional** test beside Fisher with its own BH family (Fisher is never replaced) |
| omnibus chi-square + permutation | ported |
| Harmony / `batch_key` | **not ported** (batch correction is a separate decision and could remove perturbation biology) |
| `assigned_only` re-clustering | not ported (documented limitation) |
| per-target bar figures | not ported (one heatmap covers all targets; tables are complete) |

## Fit in the current architecture

`analysis.clustering` is a new block beside `analysis.perturbation_effects`; it
runs as its own stage after the representations (and after the perturbation
effects when both are on, only so the descriptive lochNESS comparison can be
filled). Nothing in PS, lochNESS, modules, protein effects or concordance reads
the clusters. New dependency: `igraph` (scanpy's Leiden backend).
