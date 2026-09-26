"""Stage F figures (section ``cell_states``), through the FigureRegistry.

Three figures, bounded in number whatever the size of the screen: the RNA
UMAP coloured by cluster, cluster sizes with their composition by
perturbation class and by one design column, and a target x cluster
enrichment heatmap (colour = Haldane-corrected log2 odds ratio, finite at zero
counts; * = BH-FDR < alpha, o = low power). Every value is in
``tables/cell_states/``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .figures import FigureRegistry  # noqa: E402

SECTION = "cell_states"
CLASS_COLORS = {"single_targeting": "#4C72B0", "single_control": "#C44E52", "multi_targeting": "#8172B3", "multi_control": "#DD8452", "mixed_control_targeting": "#937860", "ambiguous": "#BBBBBB", "unassigned": "#DDDDDD"}


def cell_state_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    c = res.clustering
    if c is None:
        return
    key = c.key
    labels = adata.obs[key].astype(str)
    order = list(adata.obs[key].cat.categories)
    cmap = plt.get_cmap("tab20")
    if "X_umap_rna" in adata.obsm:
        U = adata.obsm["X_umap_rna"]
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        for i, k in enumerate(order):
            m = (labels == k).to_numpy()
            ax.scatter(U[m, 0], U[m, 1], s=3, color=cmap(i % 20), rasterized=True, label=k)
            ax.text(np.median(U[m, 0]), np.median(U[m, 1]), k, fontsize=8, ha="center", va="center", weight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"Leiden clusters (resolution {cfg.analysis.clustering.resolution:g}, {len(order)} clusters)", fontsize=9)
        fig.tight_layout()
        reg.save(fig, "umap_leiden", SECTION, "clusters", "RNA UMAP by Leiden cluster", "Leiden clusters computed on the RNA neighbour graph (not on the UMAP, which is only for display). Clusters are numbered states, not cell types.")
    design = next((col for col in ("sample", "lane", "batch") if col in c.by_design), None)
    ncol = 3 if design else 2
    fig, axes = plt.subplots(1, ncol, figsize=(4.2 * ncol, 3.4))
    size = c.summary.set_index("cluster")["n_cells"].reindex(order)
    axes[0].bar(range(len(order)), size, color=[cmap(i % 20) for i in range(len(order))])
    axes[0].set_xticks(range(len(order))); axes[0].set_xticklabels(order, fontsize=7); axes[0].set_ylabel("cells"); axes[0].set_title("cluster size", fontsize=9)
    comp = c.by_class.reindex(order)
    frac = comp.div(comp.sum(axis=1), axis=0)
    bottom = np.zeros(len(order))
    for k in [k for k in CLASS_COLORS if k in frac.columns] + [k for k in frac.columns if k not in CLASS_COLORS]:
        axes[1].bar(range(len(order)), frac[k], bottom=bottom, color=CLASS_COLORS.get(k, "#999999"), label=k)
        bottom += frac[k].to_numpy()
    axes[1].set_xticks(range(len(order))); axes[1].set_xticklabels(order, fontsize=7); axes[1].set_title("perturbation class", fontsize=9); axes[1].legend(fontsize=6, loc="upper right")
    if design:
        t = c.by_design[design].reindex(order)
        f = t.div(t.sum(axis=1), axis=0)
        bottom = np.zeros(len(order))
        dmap = plt.get_cmap("Set2")
        for i, lvl in enumerate(f.columns[:20]):
            axes[2].bar(range(len(order)), f[lvl], bottom=bottom, color=dmap(i % 8), label=str(lvl))
            bottom += f[lvl].to_numpy()
        axes[2].set_xticks(range(len(order))); axes[2].set_xticklabels(order, fontsize=7); axes[2].set_title(design, fontsize=9); axes[2].legend(fontsize=6, loc="upper right")
    fig.tight_layout()
    reg.save(fig, "cluster_composition", SECTION, "clusters", "Cluster sizes and composition", f"Cells per cluster, and each cluster's composition by perturbation class{' and by ' + design if design else ''}. A cluster made mostly of one {design or 'sample'} or of ambiguous/unassigned cells is flagged in cluster_summary.csv.")
    e = res.enrichment
    if e is None or e.empty:
        return
    T = e.table
    alpha = cfg.analysis.clustering.enrichment.fdr_alpha
    L = T.pivot(index="target", columns="cluster", values="log2_or_haldane")
    S = T.pivot(index="target", columns="cluster", values="significant").fillna(False)
    W = T.pivot(index="target", columns="cluster", values="low_power").fillna(False)
    cols = [k for k in order if k in L.columns]
    rank = T.groupby("target")["fdr"].min().sort_values()
    L, S, W = L.loc[rank.index, cols], S.loc[rank.index, cols], W.loc[rank.index, cols]
    lim = float(np.nanmax(np.abs(L.to_numpy()))) if np.isfinite(L.to_numpy()).any() else 1.0
    fig, ax = plt.subplots(figsize=(1.6 + 0.45 * len(cols), 1.2 + 0.22 * len(L)))
    im = ax.imshow(L.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            if bool(S.iat[i, j]):
                ax.text(j, i, "*", ha="center", va="center", fontsize=9)
            elif bool(W.iat[i, j]):
                ax.text(j, i, "o", ha="center", va="center", fontsize=6, color="#555555")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=7); ax.set_xlabel("Leiden cluster")
    ax.set_yticks(range(len(L))); ax.set_yticklabels(L.index, fontsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=0.04); cb.set_label("log2 odds ratio (Haldane-corrected)", fontsize=8)
    fig.tight_layout()
    ctrl = "non-targeting controls" if cfg.analysis.clustering.enrichment.control == "non_targeting" else "cells of all other targets"
    reg.save(fig, "perturbation_cluster_enrichment", SECTION, "enrichment", "Perturbation x cluster enrichment", f"Enrichment (red) or depletion (blue) of each target's single-guide cells in each cluster relative to {ctrl}; * = Fisher BH-FDR < {alpha}; o = fewer than {cfg.analysis.clustering.enrichment.min_control_cells_in_cluster} control cells in the cluster (low power). Rows ordered by the target's smallest FDR. Exact odds ratios and p-values: perturbation_cluster_enrichment.csv.")
