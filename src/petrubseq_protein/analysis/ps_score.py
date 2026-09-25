"""Perturbation-response score (PS) per cell.

Definition (Song et al. 2025, Nat Cell Biol 27:493; PS_python / ``pertps``
``PerturbAnalyzer.calculate_ps_score``, re-implemented here so the pipeline
has no optional dependency; see docs/reference/PERTURBATION_ANALYSIS_AUDIT.md):

for one target *t* with perturbed cells P and control cells C,

1. rank genes P vs C with a t-test (``scanpy.tl.rank_genes_groups``) and keep
   the top ``top_n_genes`` (the response signature);
2. ``beta_g`` = OLS slope of gene g on the binary indicator x (1 for P, 0 for C)
   = mean_P(g) - mean_C(g);
3. score_i = ((Y_i - mean(Y)) . beta) / (beta . beta) for every cell in P u C;
4. shift so the control mean is 0, clip to [0, scale_factor], divide by the
   maximum -> [0, 1] (``ps_scores``); the clipped-but-unnormalized value is
   kept as ``ps_raw`` (comparable across targets on the [0, scale_factor] scale).

Quadrants combine the score with the target gene's own expression (cut at the
control mean by default): ``successful_knockdown`` (high PS, low expression),
``escaper`` (high PS, high expression), ``non_responder``, ``low_signal``.

Departures from the reference, all recorded in the provenance: perturbed and
control cells are ``single_targeting`` / ``single_control`` cells only; a target
whose gene symbol is absent from ``var`` (or mapped through
``target_gene_map``) is still scored, only the quadrants are ``not_applicable``;
per-target scores are stored in ``obsm['ps_scores']`` (cells x targets) and
the signatures in ``uns['ps_signatures']``.

Response strength across targets: because beta is the mean difference, the
mean unnormalized score of perturbed minus control cells is exactly 1 for every
target, and the max-normalization makes [0, 1] scores target-specific. The
summary therefore reports the *separation* of perturbed from control cells on
the control-centred projection: ``auc_vs_control`` (Mann-Whitney AUC) and
``d_vs_control`` (Cohen's d). These are additions for cross-target comparison;
they are in-sample (the signature is fitted on the same cells), so a null target
still shows AUC somewhat above 0.5.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from ..config import Config
from scipy import stats

from ._common import Groups, cell_groups, cohen_d, dense

logger = logging.getLogger(__name__)

QUADRANT_KD, QUADRANT_ESCAPER, QUADRANT_NONRESP, QUADRANT_LOW, QUADRANT_NA = "successful_knockdown", "escaper", "non_responder", "low_signal", "not_applicable"


@dataclass
class PSResults:
    scores: pd.DataFrame            # cells x targets, NaN outside the target's comparison
    raw: pd.DataFrame               # same, clipped to [0, scale_factor] but not max-normalized
    own: pd.Series                  # each perturbed cell's score for its own target
    own_quadrant: pd.Series
    summary: pd.DataFrame
    skipped: pd.DataFrame
    signatures: pd.DataFrame        # long: target, rank, gene, beta, t_score
    info: Dict[str, Any] = field(default_factory=dict)


def _expression_cut(values: np.ndarray, method: str, q: float) -> float:
    if method == "median":
        return float(np.median(values))
    if method == "quantile":
        return float(np.quantile(values, q))
    return float(np.mean(values))


def _auc(a: np.ndarray, b: np.ndarray) -> float:
    """P(score of a perturbed cell > score of a control cell), Mann-Whitney AUC."""
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(u / (len(a) * len(b)))


def ps_for_target(X_sub, is_target: np.ndarray, var_names: pd.Index, top_n: int, scale_factor: float, rank_method: str = "t-test") -> Optional[Dict[str, Any]]:
    """PS for one target on an already subset matrix (rows = P u C). Returns None when no signature."""
    work = ad.AnnData(X=X_sub, obs=pd.DataFrame({"grp": np.where(is_target, "target", "control")}, index=[str(i) for i in range(X_sub.shape[0])]))
    work.var_names = var_names
    work.obs["grp"] = work.obs["grp"].astype("category")
    sc.tl.rank_genes_groups(work, groupby="grp", groups=["target"], reference="control", method=rank_method)
    names = pd.DataFrame(work.uns["rank_genes_groups"]["names"])["target"].head(top_n).tolist()
    scores_t = pd.DataFrame(work.uns["rank_genes_groups"]["scores"])["target"].head(top_n).to_numpy(float)
    cols = var_names.get_indexer(names)
    Y = dense(X_sub, cols=cols)
    x = is_target.astype(np.float64)
    xc = x - x.mean()
    var_x = float(xc @ xc)
    if var_x == 0:
        return None
    Yc = Y - Y.mean(axis=0)
    beta = (xc @ Yc) / var_x
    bb = float(beta @ beta)
    proj = Yc @ beta / bb if bb > 0 else np.zeros(Y.shape[0])
    proj = proj - proj[~is_target].mean()
    raw = np.clip(proj, 0, scale_factor)
    mx = float(raw.max())
    norm = raw / mx if mx > 0 else raw.copy()
    return {"genes": names, "beta": beta, "t_scores": scores_t, "raw": raw, "scores": norm, "max_raw": mx, "proj": proj}


def compute_ps(adata: ad.AnnData, cfg: Config) -> PSResults:
    pe = cfg.analysis.perturbation_effects
    pc = pe.ps
    groups: Groups = cell_groups(adata, pe.control_classes, pc.min_cells_per_target, pc.min_control_cells)
    n = adata.n_obs
    obs_names = adata.obs_names
    scores = pd.DataFrame(index=obs_names, dtype=float)
    raw = pd.DataFrame(index=obs_names, dtype=float)
    own = pd.Series(np.nan, index=obs_names, dtype=float)
    own_q = pd.Series(QUADRANT_NA, index=obs_names, dtype=object)
    rows, sig_rows = [], []
    skipped = list(groups.skipped)
    gene_map = dict(pe.target_gene_map)
    X = adata.X
    for t in groups.targets:
        pm = groups.mask(t)
        sub = pm | groups.control
        n_t, n_c = int(pm.sum()), int(groups.n_control)
        if n_t + n_c < 10:
            skipped.append({"target": t, "n_cells": n_t, "reason": "fewer than 10 cells in target + control"})
            continue
        gene = gene_map.get(t, t)
        gene_present = gene in adata.var_names
        expr = dense(X, cols=np.array([adata.var_names.get_loc(gene)])).ravel() if gene_present else None
        status = "scored"
        if gene_present:
            pct = 100.0 * float(np.mean(expr[groups.control] > 0))
            if pct < pc.min_pct_expressing_control:
                skipped.append({"target": t, "n_cells": n_t, "reason": f"target gene {gene} expressed in {pct:.2f}% of control cells (< {pc.min_pct_expressing_control}%)"})
                continue
        else:
            status = "scored_no_target_gene"
        idx = np.flatnonzero(sub)
        res = ps_for_target(X[idx], pm[idx], adata.var_names, pc.top_n_genes, pc.scale_factor)
        if res is None:
            skipped.append({"target": t, "n_cells": n_t, "reason": "no signature (zero variance)"})
            continue
        s = pd.Series(np.nan, index=obs_names, dtype=float)
        s.iloc[idx] = res["scores"]
        r = pd.Series(np.nan, index=obs_names, dtype=float)
        r.iloc[idx] = res["raw"]
        scores[t], raw[t] = s, r
        own[pm] = s[pm]
        high = s.iloc[idx] >= pc.ps_threshold
        is_t = pm[idx]
        row = {"target": t, "gene": gene if gene_present else "", "n_perturbed_cells": n_t, "n_control_cells": n_c, "n_signature_genes": len(res["genes"]),
               "mean_ps": float(s[pm].mean()), "median_ps": float(s[pm].median()), "sd_ps": float(s[pm].std(ddof=1)) if n_t > 1 else float("nan"),
               "q25_ps": float(s[pm].quantile(0.25)), "q75_ps": float(s[pm].quantile(0.75)), "pct_high_ps": float(100 * high[is_t].mean()),
               "mean_ps_control": float(s[groups.control].mean()), "pct_high_ps_control": float(100 * high[~is_t].mean()),
               "auc_vs_control": _auc(res["proj"][is_t], res["proj"][~is_t]), "d_vs_control": cohen_d(res["proj"][is_t], res["proj"][~is_t]), "status": status}
        if gene_present:
            cut = _expression_cut(expr[groups.control], pc.expression_cut, pc.expression_cut_quantile)
            hx = expr[idx] > cut
            quad = np.where(high & ~hx, QUADRANT_KD, np.where(high & hx, QUADRANT_ESCAPER, np.where(~high & hx, QUADRANT_NONRESP, QUADRANT_LOW)))
            tq, cq = quad[is_t], quad[~is_t]
            own_q[pm] = tq
            kd_c = float(100 * np.mean(cq == QUADRANT_KD))
            row.update({"expression_cut": cut, "expression_cut_method": pc.expression_cut, "pct_successful_kd": float(100 * np.mean(tq == QUADRANT_KD)), "pct_escaper": float(100 * np.mean(tq == QUADRANT_ESCAPER)),
                        "pct_non_responder": float(100 * np.mean(tq == QUADRANT_NONRESP)), "pct_low_signal": float(100 * np.mean(tq == QUADRANT_LOW)), "pct_controls_called_kd": kd_c, "net_pct_kd": float(100 * np.mean(tq == QUADRANT_KD) - kd_c)})
        else:
            row.update({"expression_cut": np.nan, "expression_cut_method": "", "pct_successful_kd": np.nan, "pct_escaper": np.nan, "pct_non_responder": np.nan, "pct_low_signal": np.nan, "pct_controls_called_kd": np.nan, "net_pct_kd": np.nan})
        rows.append(row)
        for k, (g, b, ts) in enumerate(zip(res["genes"], res["beta"], res["t_scores"])):
            sig_rows.append({"target": t, "rank": k + 1, "gene": g, "beta": float(b), "t_score": float(ts)})
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values("median_ps", ascending=False).reset_index(drop=True)
        summary.insert(0, "warnings", np.where(summary["n_perturbed_cells"] < 30, "low_support(<30 cells)", ""))
    info = {"method": "PS (Song et al. 2025; pertps definition): t-test top-N signature, OLS projection, control-centred, clipped [0, scale_factor], max-normalized",
            "top_n_genes": pc.top_n_genes, "scale_factor": pc.scale_factor, "ps_threshold": pc.ps_threshold, "expression_cut": pc.expression_cut,
            "min_cells_per_target": pc.min_cells_per_target, "min_control_cells": pc.min_control_cells, "min_pct_expressing_control": pc.min_pct_expressing_control,
            "control_classes": groups.control_classes, "n_control_cells": groups.n_control, "n_targets_scored": int(len(summary)), "n_targets_skipped": int(len(skipped)),
            "perturbed_cells": "single_targeting", "control_cells": "single_control", "target_gene_map": gene_map}
    logger.info("PS: %d targets scored, %d skipped (%d control cells)", len(summary), len(skipped), groups.n_control)
    return PSResults(scores, raw, own, own_q, summary, pd.DataFrame(skipped), pd.DataFrame(sig_rows), info)
