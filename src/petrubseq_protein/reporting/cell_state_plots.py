"""Clustering figures (section ``cell_states``, stage ``clusters``), reference ``plot_clustering``
plus this pipeline's composition panel. The enrichment figures live in ``enrichment_plots``."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ._refplot import CLASS_COLORS, despine, scatter_umap, umap_coords  # noqa: E402
from .figures import FigureRegistry  # noqa: E402

SECTION = "cell_states"
ST_CLUSTERS = "clusters"


def clustering_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    c = res.clustering
    if c is None:
        return
    key = c.key
    labels = adata.obs[key].astype(str)
    order = list(adata.obs[key].cat.categories) if hasattr(adata.obs[key], "cat") else sorted(labels.unique(), key=lambda x: int(x) if str(x).isdigit() else str(x))
    cmap = plt.get_cmap("tab20")
    coords = umap_coords(adata)
    obs = adata.obs
    if coords is not None:
        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        U = coords
        for i, k in enumerate(order):
            m = (labels == k).to_numpy()
            ax.scatter(U[m, 0], U[m, 1], s=3, color=cmap(i % 20), rasterized=True, label=k, linewidths=0)
            ax.text(np.median(U[m, 0]), np.median(U[m, 1]), k, fontsize=8, ha="center", va="center", weight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
        ax.set_title(f"Leiden clusters (resolution {cfg.analysis.clustering.resolution:g})", fontsize=10)
        if len(order) <= 25:
            ax.legend(markerscale=4, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
        despine(ax, left=True, bottom=True)
        fig.tight_layout()
        reg.save(fig, "umap_clusters", SECTION, ST_CLUSTERS, "UMAP coloured by Leiden cluster", f"{len(order)} clusters at resolution {cfg.analysis.clustering.resolution:g}. Clusters are numbered states, not cell types; the graph, not the UMAP, was clustered.")
        qc_keys = [(k, l) for k, l in (("n_genes_by_counts", "Genes per cell"), ("total_counts", "Total UMIs"), ("pct_counts_mt", "% mitochondrial")) if k in obs.columns]
        if qc_keys:
            fig, axes = plt.subplots(1, len(qc_keys), figsize=(4.6 * len(qc_keys), 4.0))
            for ax, (k, l) in zip(np.atleast_1d(axes), qc_keys):
                scatter_umap(ax, coords, obs[k].to_numpy(dtype=float), False, l)
            fig.tight_layout()
            reg.save(fig, "umap_qc_metrics", SECTION, ST_CLUSTERS, "UMAP coloured by QC metrics", "A cluster driven purely by library size or mitochondrial content is a technical artefact rather than a biological state.")
        for col, title, name in (("lane", "Lane", "umap_lane"), ("sample", "Sample", "umap_sample")):
            if col in obs.columns and 1 < obs[col].astype(str).nunique() <= 25:
                fig, ax = plt.subplots(figsize=(6.0, 5.0))
                scatter_umap(ax, coords, obs[col].astype(str).to_numpy(), True, title)
                fig.tight_layout()
                reg.save(fig, name, SECTION, ST_CLUSTERS, f"UMAP coloured by {col}", f"{title}s should intermix. Separation by {col} means a batch effect (no batch correction is applied).")
        if "perturbation_class" in obs.columns:
            fig, ax = plt.subplots(figsize=(6.0, 5.0))
            scatter_umap(ax, coords, obs["perturbation_class"].astype(str).to_numpy(), True, "Guide assignment class")
            fig.tight_layout()
            reg.save(fig, "umap_assignment_class", SECTION, ST_CLUSTERS, "UMAP coloured by guide assignment", "Unassigned cells clustering together often indicates a technically distinct population (e.g. low-quality cells) rather than a biological one.")
        if "target" in obs.columns:
            tg = obs["target"].astype(str).where(obs["perturbation_class"].astype(str).isin(["single_targeting", "single_control"]), obs["perturbation_class"].astype(str))
            n_targets = tg.nunique()
            fig, ax = plt.subplots(figsize=(7.2, 5.0))
            scatter_umap(ax, coords, tg.to_numpy(), True, f"Target gene ({n_targets} levels)", legend=n_targets <= 25)
            fig.tight_layout()
            reg.save(fig, "umap_target_gene", SECTION, ST_CLUSTERS, "UMAP coloured by target gene", "Most Perturb-seq screens show perturbed cells mixed throughout the embedding; a target forming its own island has a strong phenotype." + ("" if n_targets <= 25 else " Legend omitted (too many targets)."))
    # cluster sizes + composition
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
        dmap = plt.get_cmap("tab20")
        for i, lvl in enumerate(f.columns[:20]):
            axes[2].bar(range(len(order)), f[lvl], bottom=bottom, color=dmap(i % 20), label=str(lvl))
            bottom += f[lvl].to_numpy()
        axes[2].set_xticks(range(len(order))); axes[2].set_xticklabels(order, fontsize=7); axes[2].set_title(f"{design} composition per cluster", fontsize=9); axes[2].legend(fontsize=6, loc="upper right")
    fig.tight_layout()
    reg.save(fig, "cluster_composition", SECTION, ST_CLUSTERS, "Cluster sizes and composition", f"Cells per cluster, and each cluster's composition by perturbation class{' and by ' + design if design else ''}. Clusters made up almost entirely of one {design or 'sample'} are candidates for batch effects; a cluster of ambiguous/unassigned cells is flagged in cluster_summary.csv.")


def cell_state_figures(res, adata, cfg, reg: FigureRegistry) -> None:  # backwards-compatible entry point
    from .enrichment_plots import enrichment_figures

    clustering_figures(res, adata, cfg, reg)
    enrichment_figures(res.enrichment, adata, cfg, reg)
