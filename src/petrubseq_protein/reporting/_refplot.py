"""Shared drawing helpers for the reference-style analysis figures."""

from __future__ import annotations

from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

CLASS_COLORS = {"single_targeting": "#2b6cb0", "single_control": "#38a169", "multi_targeting": "#805ad5", "multi_control": "#dd6b20", "mixed_control_targeting": "#b7791f", "ambiguous": "#dd6b20", "unassigned": "#a0aec0"}
HIT_RED, GREY, BLUE, GREEN = "#c53030", "#a0aec0", "#2b6cb0", "#38a169"


def gene_values(adata, gene: str) -> np.ndarray:
    col = adata.X[:, adata.var_names.get_loc(gene)]
    if sp.issparse(col):
        col = col.toarray()
    return np.asarray(col, dtype=float).ravel()


def umap_coords(adata) -> Optional[np.ndarray]:
    for key in ("X_umap_rna", "X_umap"):
        if key in adata.obsm:
            return np.asarray(adata.obsm[key])
    return None


def scatter_umap(ax, coords: np.ndarray, values: np.ndarray, categorical: bool, title: str, size: float = 3.0, cmap: str = "viridis", legend: bool = True, xlabel: str = "UMAP1", ylabel: str = "UMAP2") -> None:
    """Reference ``_scatter_umap``: categorical colours from tab20, numeric with a colour bar."""
    if categorical:
        cats = pd.Index(pd.unique(pd.Series(values).astype(str))).sort_values()
        palette = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(cats), 3)) if len(cats) > 20 else np.arange(max(len(cats), 3)) % 20 / 19.0)
        for i, cat in enumerate(cats):
            m = values.astype(str) == cat
            ax.scatter(coords[m, 0], coords[m, 1], s=size, color=palette[i % len(palette)], label=str(cat), linewidths=0, rasterized=True)
        if legend and len(cats) <= 25:
            ax.legend(markerscale=4, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    else:
        sca = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=size, cmap=cmap, linewidths=0, rasterized=True)
        plt.colorbar(sca, ax=ax, shrink=0.75)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)


def despine(ax, left: bool = False, bottom: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if left:
        ax.spines["left"].set_visible(False)
    if bottom:
        ax.spines["bottom"].set_visible(False)


def order_by_similarity(matrix: pd.DataFrame) -> List[str]:
    """Reference ``_order_by_similarity``: average-linkage on correlation distance."""
    if matrix.shape[0] < 3:
        return list(matrix.index)
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist

        values = np.nan_to_num(matrix.to_numpy(dtype=float))
        dist = pdist(values, metric="correlation")
        if not np.all(np.isfinite(dist)):
            return list(matrix.index)
        return [matrix.index[i] for i in leaves_list(linkage(dist, method="average"))]
    except Exception:  # pragma: no cover - ordering is cosmetic
        return list(matrix.index)


def ecdf(ax, values: np.ndarray, label: str) -> None:
    v = np.sort(np.asarray(values, dtype=float))
    y = np.arange(1, v.size + 1) / v.size
    ax.step(v, y, where="post", label=label, lw=1.2)


def violin(ax, groups: List[np.ndarray], labels: List[str], ylabel: str = "") -> None:
    data = [g if g.size else np.array([0.0]) for g in groups]
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for b in parts["bodies"]:
        b.set_alpha(0.7)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=20, fontsize=7)
    ax.set_ylabel(ylabel)


def label_clusters(ax, adata, key: str, coords: np.ndarray, highlight: str = "") -> None:
    """Reference ``_label_clusters``: cluster ids at their centroids, the highlighted one emphasised."""
    if key not in adata.obs.columns:
        return
    clusters = adata.obs[key].astype(str).to_numpy()
    for cl in np.unique(clusters):
        m = clusters == cl
        if m.sum() == 0:
            continue
        cx, cy = np.median(coords[m, 0]), np.median(coords[m, 1])
        is_top = str(cl) == str(highlight)
        ax.text(cx, cy, str(cl), fontsize=8 if is_top else 6.5, fontweight="bold" if is_top else "normal", color="#1a202c" if is_top else "#4a5568", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="#fefcbf" if is_top else "white", edgecolor="none", alpha=0.85 if is_top else 0.6))


def label_colors(labels: List[str]) -> dict:
    uniq = list(dict.fromkeys(labels))
    cmap = plt.get_cmap("tab20")
    return {lab: tuple(cmap(i % 20)) for i, lab in enumerate(uniq)}


def block_spans(seq: List[str]):
    spans = []
    if not len(seq):
        return spans
    start = 0
    for i in range(1, len(seq) + 1):
        if i == len(seq) or seq[i] != seq[start]:
            spans.append((seq[start], start, i))
            start = i
    return spans


def display_order(items, label_of, label_rank, item_rank):
    return sorted(items, key=lambda x: (label_rank[label_of[x]], item_rank[x]))
