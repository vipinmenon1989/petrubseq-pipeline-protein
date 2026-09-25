"""lochNESS: local neighbourhood enrichment of each perturbation.

Definition (Huang et al. 2023, Nature 623:772, "local cellular heuristic
neighbourhood enrichment specificity score"; implementation as ported from
pertTF by weili-lab/perturbseq-pipeline): for every cell and every target *g*,

    lochNESS(cell, g) = local_fraction(g) / overall_fraction(g) - 1

with ``local_fraction`` the share of the cell's k nearest neighbours (in PCA
space, never UMAP) that are perturbed cells of *g*, and ``overall_fraction``
*g*'s share of all cells in the object. 0 = chance, > 0 locally enriched,
< 0 depleted. The value in a target's own cells ("self" score) says whether
the perturbation occupies a distinct region of cell-state space.

Additions over the reference (recorded in the provenance): k is capped at
``max_k_fraction`` x n cells for small objects; a seeded label-permutation
null gives a z-score and empirical p-value for each target's own-cell mean
(Huang et al. compare against permuted labels); per-sample own-cell means
expose enrichment driven by one sample; scores live in ``obsm['lochness']``.

Control reference: because ``overall_fraction`` is taken over all cells, a
perturbation without effect is still enriched (> 0) wherever the unperturbed
cells live when other perturbations occupy separate regions, and the label
permutation (over all cells) flags it too. ``mean_lochness_in_control_cells``
(the target's score averaged over control cells) and ``delta_own_vs_control``
(own-cell mean minus that) separate "the target's cells sit where controls sit"
(delta ~ 0) from "the target's cells occupy their own state" (delta > 0).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from ..config import Config
from ._common import cell_groups

logger = logging.getLogger(__name__)


@dataclass
class LochnessResults:
    scores: pd.DataFrame           # cells x targets (all cells scored)
    self_score: pd.Series          # each cell's score for its own label (targets and the control label)
    summary: pd.DataFrame
    skipped: pd.DataFrame
    by_sample: pd.DataFrame        # targets x samples: own-cell mean per sample (empty if no sample column)
    info: Dict[str, Any] = field(default_factory=dict)


def neighbor_adjacency(adata: ad.AnnData, use_rep: str, k: int, n_pcs: int, seed: int) -> sp.csr_matrix:
    """Binary kNN adjacency (self excluded) from ``sc.pp.neighbors`` on ``obsm[use_rep]``."""
    tmp = ad.AnnData(obs=pd.DataFrame(index=adata.obs_names))
    tmp.obsm[use_rep] = np.asarray(adata.obsm[use_rep], dtype=np.float32)
    npc = min(n_pcs, tmp.obsm[use_rep].shape[1]) if n_pcs else None
    sc.pp.neighbors(tmp, n_neighbors=k, n_pcs=npc, use_rep=use_rep, key_added="lochness_nn", random_state=seed)
    adj = sp.csr_matrix(tmp.obsp["lochness_nn_distances"])
    adj.data = np.ones_like(adj.data)
    return adj


def lochness_score(adj: sp.csr_matrix, neighbor_counts: np.ndarray, indicator: np.ndarray, overall_fraction: float) -> np.ndarray:
    if overall_fraction <= 0:
        return np.full(adj.shape[0], np.nan)
    local = adj @ indicator.astype(np.float64)
    return (local / neighbor_counts) / overall_fraction - 1.0


def compute_lochness(adata: ad.AnnData, cfg: Config) -> LochnessResults:
    pe = cfg.analysis.perturbation_effects
    lc = pe.lochness
    groups = cell_groups(adata, pe.control_classes, lc.min_cells_per_target, 1)
    n = adata.n_obs
    use_rep = lc.use_rep or "X_pca"
    if use_rep not in adata.obsm:
        raise ValueError(f"lochNESS needs obsm[{use_rep!r}] (available: {sorted(adata.obsm)})")
    k_cap = max(15, int(np.floor(lc.max_k_fraction * n)))
    k = int(min(lc.n_neighbors, k_cap, max(n - 1, 2)))
    adj = neighbor_adjacency(adata, use_rep, k, lc.n_pcs, cfg.compute.seed)
    counts = np.asarray(adj.sum(axis=1)).ravel().astype(float)
    counts[counts == 0] = np.nan
    rng = np.random.default_rng(cfg.compute.seed)
    scores = pd.DataFrame(index=adata.obs_names, dtype=float)
    self_score = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    sample = adata.obs["sample"].astype(str).to_numpy() if "sample" in adata.obs.columns else None
    rows, by_sample = [], {}
    for t in groups.targets:
        ind = groups.mask(t).astype(np.float64)
        n_t = int(ind.sum())
        frac = n_t / n
        s = lochness_score(adj, counts, ind, frac)
        scores[t] = s
        own = ind.astype(bool)
        self_score[own] = s[own]
        obs_mean = float(np.nanmean(s[own]))
        row = {"target": t, "n_cells": n_t, "overall_fraction_pct": 100 * frac, "mean_lochness_in_own_cells": obs_mean,
               "median_lochness_in_own_cells": float(np.nanmedian(s[own])), "mean_lochness_all_cells": float(np.nanmean(s)), "max_lochness": float(np.nanmax(s)),
               "pct_own_cells_enriched": float(100 * np.nanmean(s[own] > lc.enrichment_cut)), "pct_all_cells_enriched": float(100 * np.nanmean(s > lc.enrichment_cut)),
               "mean_lochness_in_control_cells": float(np.nanmean(s[groups.control])) if groups.n_control else float("nan")}
        row["delta_own_vs_control"] = row["mean_lochness_in_own_cells"] - row["mean_lochness_in_control_cells"]
        if lc.n_permutations > 0:
            null = np.empty(lc.n_permutations)
            for i in range(lc.n_permutations):
                pick = rng.choice(n, size=n_t, replace=False)
                pind = np.zeros(n)
                pind[pick] = 1.0
                ps_ = lochness_score(adj, counts, pind, frac)
                null[i] = np.nanmean(ps_[pick])
            sd = float(null.std(ddof=1)) if lc.n_permutations > 1 else float("nan")
            row["null_mean"] = float(null.mean())
            row["null_sd"] = sd
            row["z_score"] = float((obs_mean - null.mean()) / sd) if sd and sd > 0 else float("nan")
            row["p_empirical"] = float((1 + np.sum(null >= obs_mean)) / (lc.n_permutations + 1))
        if sample is not None:
            per = pd.Series(s[own]).groupby(sample[own]).mean()
            by_sample[t] = per.to_dict()
            row["n_samples"] = int(per.size)
            row["max_sample_share"] = float(pd.Series(sample[own]).value_counts(normalize=True).max())
        rows.append(row)
    # control label self score
    if groups.n_control > 0:
        cind = groups.control.astype(np.float64)
        cs = lochness_score(adj, counts, cind, groups.n_control / n)
        self_score[groups.control] = cs[groups.control]
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values("mean_lochness_in_own_cells", ascending=False).reset_index(drop=True)
        if "p_empirical" in summary.columns:
            from ._common import bh_fdr
            summary["fdr"] = bh_fdr(summary["p_empirical"].to_numpy())
    info = {"method": "lochNESS = local_fraction/overall_fraction - 1 over k nearest neighbours in PCA space (Huang et al. 2023; pertTF port)",
            "use_rep": use_rep, "n_pcs": lc.n_pcs, "k_requested": lc.n_neighbors, "k_used": k, "k_capped": bool(k < lc.n_neighbors), "max_k_fraction": lc.max_k_fraction,
            "median_neighbors": float(np.nanmedian(counts)), "n_permutations": lc.n_permutations, "enrichment_cut": lc.enrichment_cut,
            "overall_fraction_denominator": "all cells in the object", "control_classes": groups.control_classes, "n_control_cells": groups.n_control,
            "n_targets_scored": int(len(summary)), "n_targets_skipped": int(len(groups.skipped)), "seed": cfg.compute.seed}
    logger.info("lochNESS: %d targets, k=%d (%s), %d permutations", len(summary), k, "capped" if info["k_capped"] else "as configured", lc.n_permutations)
    return LochnessResults(scores, self_score, summary, pd.DataFrame(groups.skipped), pd.DataFrame(by_sample).T if by_sample else pd.DataFrame(), info)
