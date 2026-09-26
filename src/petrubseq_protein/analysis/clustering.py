"""RNA cell-state clustering (Leiden) and cluster composition.

Leiden community detection (Traag, Waltman & van Eck 2019, Sci Rep 9:5233) on
the **existing** RNA neighbour graph (``obsp['rna_connectivities']``, built from
the RNA PCA in the representation stage), through ``scanpy.tl.leiden`` with the
igraph backend, ``n_iterations`` 2 and the run seed, as in
weili-lab/perturbseq-pipeline. Nothing is re-normalized or re-embedded.

Population: every QC-passing cell in the processed object, i.e. the cells the
RNA graph was built on (reference default, ``assigned_only: false``). Guide
labels never enter the clustering, so clusters are defined by expression alone;
ambiguous / multi-guide / unassigned cells are clustered but excluded later
from the enrichment tests. Clusters are numbered states, never named.

Composition tables expose clusters dominated by one sample / lane / batch
(>= 80 % of the cluster from a level holding < 50 % of all cells), by ambiguous
or unassigned cells (> 50 %), or with a median library size >= 1.5 x or
<= 0.67 x the overall median. These are descriptive flags: nothing is removed.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd

from ..config import Config

logger = logging.getLogger(__name__)

DESIGN_COLUMNS = ("sample", "lane", "batch", "donor", "condition")
DOMINANCE_SHARE = 0.8     # a level holding >= 80 % of a cluster ...
DOMINANCE_GLOBAL = 0.5    # ... while holding < 50 % of all cells is flagged
DEPTH_RATIO = 1.5         # cluster median UMIs >= 1.5 x (or <= 1/1.5 x) the overall median is flagged


@dataclass
class ClusteringResults:
    key: str
    labels: pd.Series
    summary: pd.DataFrame                   # one row per cluster
    by_class: pd.DataFrame                  # cluster x perturbation_class counts
    by_design: Dict[str, pd.DataFrame]      # column -> cluster x level counts
    by_target: pd.DataFrame                 # cluster x target counts (single-guide cells, controls as their class)
    info: Dict[str, Any] = field(default_factory=dict)
    flags: List[str] = field(default_factory=list)


def _counts(labels: pd.Series, values: pd.Series, order: List[str]) -> pd.DataFrame:
    t = pd.crosstab(labels, values.astype(str)).reindex(index=order, fill_value=0)
    t.index.name = "cluster"
    return t


def compute_clustering(adata: ad.AnnData, cfg: Config) -> ClusteringResults:
    import scanpy as sc

    cc = cfg.analysis.clustering
    key = cc.key
    nkey = cc.neighbors_key
    if f"{nkey}_connectivities" not in adata.obsp:
        raise ValueError(f"analysis.clustering needs the RNA neighbour graph obsp['{nkey}_connectivities'] (available: {sorted(adata.obsp)}); it is built in the representation stage")
    if key in adata.obs.columns and not cc.overwrite:
        raise ValueError(f"analysis.clustering.key {key!r} is already an obs column (e.g. from columns.keep); choose another key or set analysis.clustering.overwrite: true")
    try:
        import igraph  # noqa: F401
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("analysis.clustering needs the 'igraph' package (pip install igraph; listed in environment.yml)") from exc
    sc.tl.leiden(adata, resolution=cc.resolution, key_added=key, neighbors_key=nkey, flavor="igraph", n_iterations=cc.n_iterations, directed=False, random_state=cfg.compute.seed)
    # numeric, size-independent stable order: scanpy labels 0..k-1 by decreasing size
    order = sorted(adata.obs[key].astype(str).unique(), key=int)
    adata.obs[key] = pd.Categorical(adata.obs[key].astype(str), categories=order)
    labels = adata.obs[key].astype(str)
    obs = adata.obs
    n = adata.n_obs
    size = labels.value_counts().reindex(order)
    klass = obs["perturbation_class"].astype(str) if "perturbation_class" in obs.columns else pd.Series("unknown", index=obs.index)
    by_class = _counts(labels, klass, order)
    rows = []
    flags: List[str] = []
    by_design: Dict[str, pd.DataFrame] = {}
    for col in DESIGN_COLUMNS:
        if col in obs.columns and obs[col].astype(str).nunique() > 1:
            by_design[col] = _counts(labels, obs[col], order)
    for c in order:
        m = (labels == c).to_numpy()
        row: Dict[str, Any] = {"cluster": c, "n_cells": int(m.sum()), "fraction": float(m.sum() / n)}
        kc = klass[m].value_counts()
        for k in ("single_targeting", "single_control", "multi_targeting", "multi_control", "mixed_control_targeting", "ambiguous", "unassigned"):
            row[f"n_{k}"] = int(kc.get(k, 0))
        row["frac_ambiguous_or_unassigned"] = float((kc.get("ambiguous", 0) + kc.get("unassigned", 0)) / m.sum())
        for qc in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
            if qc in obs.columns:
                row[f"median_{qc}"] = float(np.nanmedian(obs.loc[m, qc].to_numpy(float)))
        warn = []
        for col, tab in by_design.items():
            share = tab.loc[c] / max(tab.loc[c].sum(), 1)
            top = share.idxmax()
            glob = float((obs[col].astype(str) == top).mean())
            row[f"top_{col}"] = top
            row[f"top_{col}_share"] = float(share.max())
            if share.max() >= DOMINANCE_SHARE and glob < DOMINANCE_GLOBAL:
                warn.append(f"{col} {top} = {100 * share.max():.0f}% of cells (vs {100 * glob:.0f}% overall)")
        if row["frac_ambiguous_or_unassigned"] > 0.5:
            warn.append(f"{100 * row['frac_ambiguous_or_unassigned']:.0f}% ambiguous/unassigned")
        row["flags"] = "; ".join(warn)
        if warn:
            flags.append(f"cluster {c}: " + "; ".join(warn))
        rows.append(row)
    summary = pd.DataFrame(rows)
    if "median_total_counts" in summary.columns and len(summary) > 1:
        med = float(np.nanmedian(obs["total_counts"].to_numpy(float)))
        summary["depth_ratio_to_overall"] = summary["median_total_counts"] / med
        for i, r in summary.iterrows():
            ratio = r["depth_ratio_to_overall"]
            if ratio >= DEPTH_RATIO or ratio <= 1 / DEPTH_RATIO:
                msg = f"median library size {ratio:.2f} x the overall median (possible depth-driven state)"
                summary.at[i, "flags"] = "; ".join([x for x in (r["flags"], msg) if x])
                flags.append(f"cluster {r['cluster']}: {msg}")
    tgt = obs["target"].astype(str).where(klass == "single_targeting", klass) if "target" in obs.columns else klass
    by_target = _counts(labels, tgt, order)
    info = {"method": "Leiden (igraph backend via scanpy.tl.leiden) on the existing RNA neighbour graph", "key": key, "neighbors_key": nkey, "resolution": cc.resolution,
            "n_iterations": cc.n_iterations, "random_state": cfg.compute.seed, "directed": False, "population": "all QC-passing cells in the processed object",
            "n_clusters": len(order), "cluster_sizes": {c: int(size[c]) for c in order}, "min_cluster_size": int(size.min()), "max_cluster_size": int(size.max()),
            "design_columns_checked": list(by_design), "n_flagged_clusters": len(flags)}
    logger.info("Leiden (resolution %g): %d clusters (sizes %d-%d)", cc.resolution, len(order), size.min(), size.max())
    return ClusteringResults(key, labels, summary, by_class, by_design, by_target, info, flags)
