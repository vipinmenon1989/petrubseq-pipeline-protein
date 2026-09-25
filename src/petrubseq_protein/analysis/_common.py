"""Shared helpers for the Stage E analyses: cell groups, statistics, dense access.

Group definitions used by every analysis (docs/PERTURBATION_EFFECTS.md):

* **perturbed cells of target t**: ``obs['perturbation_class'] == 'single_targeting'``
  and ``obs['target'] == t``;
* **control cells**: ``obs['perturbation_class'] == 'single_control'`` and
  ``obs['control_class']`` in ``analysis.perturbation_effects.control_classes``
  (default ``non_targeting``).

Ambiguous, multi-guide, mixed and unassigned cells belong to neither group.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

CLASS_TARGETING = "single_targeting"
CLASS_CONTROL = "single_control"


@dataclass
class Groups:
    targets: List[str]
    perturbed: Dict[str, np.ndarray]
    control: np.ndarray
    counts: pd.Series
    n_control: int
    control_classes: List[str]
    skipped: List[dict] = field(default_factory=list)

    def mask(self, target: str) -> np.ndarray:
        return self.perturbed[target]


def cell_groups(adata: ad.AnnData, control_classes: Sequence[str], min_cells: int, min_control: int = 1) -> Groups:
    obs = adata.obs
    klass = obs["perturbation_class"].astype(str).to_numpy()
    target = obs["target"].astype(str).to_numpy()
    cclass = obs["control_class"].astype(str).to_numpy() if "control_class" in obs.columns else np.full(adata.n_obs, "")
    control = (klass == CLASS_CONTROL) & np.isin(cclass, list(control_classes))
    tmask = klass == CLASS_TARGETING
    counts = pd.Series(target[tmask]).value_counts().sort_index()
    perturbed, targets, skipped = {}, [], []
    for t, n in counts.items():
        if n >= min_cells:
            perturbed[t] = tmask & (target == t)
            targets.append(t)
        else:
            skipped.append({"target": t, "n_cells": int(n), "reason": f"fewer than {min_cells} single-guide targeting cells"})
    if int(control.sum()) < min_control:
        skipped.append({"target": "*", "n_cells": int(control.sum()), "reason": f"fewer than {min_control} control cells ({', '.join(control_classes)})"})
        targets, perturbed = [], {}
    return Groups(targets, perturbed, control, counts, int(control.sum()), list(control_classes), skipped)


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    v = p[ok]
    n = v.size
    order = np.argsort(v)
    ranked = v[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


def dense(X, rows: Optional[np.ndarray] = None, cols: Optional[np.ndarray] = None, dtype=np.float64) -> np.ndarray:
    M = X
    if rows is not None:
        M = M[rows]
    if cols is not None:
        M = M[:, cols]
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    return M.astype(dtype, copy=False)


def group_mean_var(X, mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """Per-column mean and unbiased variance over the rows in ``mask`` (sparse-safe)."""
    n = int(mask.sum())
    sub = X[mask]
    if sp.issparse(sub):
        mean = np.asarray(sub.mean(axis=0)).ravel()
        sq = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
    else:
        sub = np.asarray(sub, dtype=np.float64)
        mean = sub.mean(axis=0)
        sq = (sub * sub).mean(axis=0)
    var = (sq - mean * mean) * (n / max(n - 1, 1))
    return mean, np.clip(var, 0, None), n


def welch_t(mean_a, var_a, n_a, mean_b, var_b, n_b) -> Tuple[np.ndarray, np.ndarray]:
    """Two-sided Welch t-test from summary statistics (vectorized)."""
    se2 = var_a / n_a + var_b / n_b
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (mean_a - mean_b) / np.sqrt(se2)
        df = se2 ** 2 / ((var_a / n_a) ** 2 / max(n_a - 1, 1) + (var_b / n_b) ** 2 / max(n_b - 1, 1))
        p = 2.0 * stats.t.sf(np.abs(t), np.clip(df, 1, None))
    t = np.where(np.isfinite(t), t, 0.0)
    p = np.where(np.isfinite(p) & (se2 > 0), p, 1.0)
    return t, p


def spearman(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return float("nan"), float("nan")
    r = stats.spearmanr(x[ok], y[ok])
    return float(r.statistic if hasattr(r, "statistic") else r[0]), float(r.pvalue if hasattr(r, "pvalue") else r[1])


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float("nan")
    s = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / max(na + nb - 2, 1))
    return float((a.mean() - b.mean()) / s) if s > 0 else float("nan")


def protein_matrix(adata: ad.AnnData, key: str) -> Optional[pd.DataFrame]:
    """The normalized protein DataFrame (cells x proteins) or None."""
    if key not in adata.obsm:
        return None
    P = adata.obsm[key]
    if not isinstance(P, pd.DataFrame):
        return None
    return P.astype(float)
