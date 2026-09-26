"""Perturbation-response score (PS) per cell.

Definition (Song et al. 2025, Nat Cell Biol 27:493; PS_python / ``pertps``
``PerturbAnalyzer.calculate_ps_score``, re-implemented here so the pipeline has no
optional dependency; verified to 4e-16 against ``pertps`` by
``scripts/compare_reference_perturbation.py``): for one target *t* with perturbed
cells P and control cells C,

1. rank genes P vs C with a t-test (``scanpy.tl.rank_genes_groups``) and keep the
   top ``top_n_genes`` (the response signature);
2. ``beta_g`` = OLS slope of gene g on the binary indicator x (1 for P, 0 for C);
3. score_i = ((Y_i - mean(Y)) . beta) / (beta . beta) for every cell in P u C;
4. shift so the control mean is 0, clip to [0, scale_factor], divide by the
   maximum -> [0, 1] (``ps_scores``); the clipped-but-unnormalized value is kept
   as ``ps_scores_raw``.

Quadrants combine the score with the target gene's own expression (cut at the
control mean by default): ``successful_knockdown`` (high PS, low expression),
``escaper`` (high PS, high expression), ``non_responder``, ``low_signal``; the
same classification of the control cells gives ``pct_controls_called_kd`` and
``net_pct_kd``. Targets expressed in < ``min_pct_expressing_control`` % of the
control cells are skipped (knockdown unmeasurable); targets whose gene symbol is
absent from ``var`` are skipped as in the reference unless
``score_targets_without_gene`` is set. The summary is sorted by
``pct_successful_kd`` (reference order; ``top_targets`` follows it).

Supervised LDA embedding (reference ``_compute_lda_embedding`` / PS_python
``compute_lda_umap``): the HVGs (at most ``lda_max_genes``) -> 2000 HVGs -> scale
(max 10) -> PCA ``lda_n_pcs`` -> LDA (eigen solver, shrinkage auto) trained on
the scored targets (> 5 cells) plus the control label -> UMAP (30 neighbours,
min_dist 0.01, cosine, seed 42) -> ``obsm['X_lda_umap']`` (NaN for cells outside
the trained classes), ``obs['lda_label']``.

Extensions kept as additional columns: ``auc_vs_control`` / ``d_vs_control``
(in-sample separation of perturbed from control cells on the control-centred
projection), spread statistics, per-target signatures.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

from ..config import Config
from ._common import Groups, cell_groups, cohen_d, dense

logger = logging.getLogger(__name__)

QUADRANT_KD, QUADRANT_ESCAPER, QUADRANT_NONRESP, QUADRANT_LOW, QUADRANT_NA = "successful_knockdown", "escaper", "non_responder", "low_signal", "not_applicable"
QUADRANT_LABELS = {QUADRANT_KD: "successful knockdown", QUADRANT_ESCAPER: "escaper", QUADRANT_NONRESP: "non-responder", QUADRANT_LOW: "low signal"}
QUADRANT_COLORS = {QUADRANT_KD: "#2f855a", QUADRANT_ESCAPER: "#c53030", QUADRANT_NONRESP: "#2b6cb0", QUADRANT_LOW: "#a0aec0"}
LDA_CONTROL_LABEL = "NT"
LDA_OTHER_LABEL = "Other"


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
    quadrants: Dict[str, pd.Series] = field(default_factory=dict)     # target -> quadrant of every scored cell
    expression_cut: Dict[str, float] = field(default_factory=dict)
    genes: Dict[str, str] = field(default_factory=dict)               # target -> gene symbol used
    lda_umap: Optional[np.ndarray] = None
    lda_label: Optional[pd.Series] = None
    lda_note: str = ""
    ps_threshold: float = 0.5

    @property
    def targets(self) -> List[str]:
        return list(self.summary["target"]) if not self.summary.empty else []

    def top_targets(self, n: int) -> List[str]:
        return list(self.summary.head(n)["target"]) if not self.summary.empty else []


def _expression_cut(values: np.ndarray, method: str, q: float) -> float:
    if values.size == 0:
        return 0.0
    if method == "median":
        return float(np.median(values))
    if method == "quantile":
        return float(np.quantile(values, q))
    return float(np.mean(values))


def _auc(a: np.ndarray, b: np.ndarray) -> float:
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
    obs_names = adata.obs_names
    scores = pd.DataFrame(index=obs_names, dtype=float)
    raw = pd.DataFrame(index=obs_names, dtype=float)
    own = pd.Series(np.nan, index=obs_names, dtype=float)
    own_q = pd.Series(QUADRANT_NA, index=obs_names, dtype=object)
    rows, sig_rows = [], []
    skipped = list(groups.skipped)
    gene_map = dict(pe.target_gene_map)
    quadrants: Dict[str, pd.Series] = {}
    cuts: Dict[str, float] = {}
    genes_used: Dict[str, str] = {}
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
        if not gene_present and not pc.score_targets_without_gene:
            skipped.append({"target": t, "n_cells": n_t, "reason": f"target gene {gene} not in the expression matrix"})
            continue
        expr = dense(X, cols=np.array([adata.var_names.get_loc(gene)])).ravel() if gene_present else None
        status = "scored"
        if gene_present:
            pct = 100.0 * float(np.mean(expr[groups.control] > 0))
            if pct < pc.min_pct_expressing_control:
                skipped.append({"target": t, "n_cells": n_t, "reason": f"not detectably expressed in control cells ({pct:.2f}% of control cells, threshold {pc.min_pct_expressing_control}%)"})
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
        genes_used[t] = gene if gene_present else ""
        high = s.iloc[idx] >= pc.ps_threshold
        is_t = pm[idx]
        row = {"target": t, "gene": gene if gene_present else "", "n_perturbed_cells": n_t, "n_control_cells": n_c, "n_signature_genes": len(res["genes"]),
               "mean_ps": float(s[pm].mean()), "median_ps": float(s[pm].median()), "pct_high_ps": float(100 * high[is_t].mean())}
        if gene_present:
            cut = _expression_cut(expr[groups.control], pc.expression_cut, pc.expression_cut_quantile)
            hx = expr[idx] > cut
            quad = np.where(high & ~hx, QUADRANT_KD, np.where(high & hx, QUADRANT_ESCAPER, np.where(~high & hx, QUADRANT_NONRESP, QUADRANT_LOW)))
            tq, cq = quad[is_t], quad[~is_t]
            own_q[pm] = tq
            quadrants[t] = pd.Series(quad, index=obs_names[idx], dtype=object)
            cuts[t] = cut
            kd_c = float(100 * np.mean(cq == QUADRANT_KD))
            row.update({"pct_successful_kd": float(100 * np.mean(tq == QUADRANT_KD)), "pct_escaper": float(100 * np.mean(tq == QUADRANT_ESCAPER)), "pct_non_responder": float(100 * np.mean(tq == QUADRANT_NONRESP)), "pct_low_signal": float(100 * np.mean(tq == QUADRANT_LOW)),
                        "pct_controls_called_kd": kd_c, "net_pct_kd": float(100 * np.mean(tq == QUADRANT_KD) - kd_c), "expression_cut": cut, "expression_cut_method": pc.expression_cut})
        else:
            row.update({"pct_successful_kd": np.nan, "pct_escaper": np.nan, "pct_non_responder": np.nan, "pct_low_signal": np.nan, "pct_controls_called_kd": np.nan, "net_pct_kd": np.nan, "expression_cut": np.nan, "expression_cut_method": ""})
        row.update({"sd_ps": float(s[pm].std(ddof=1)) if n_t > 1 else float("nan"), "q25_ps": float(s[pm].quantile(0.25)), "q75_ps": float(s[pm].quantile(0.75)),
                    "mean_ps_control": float(s[groups.control].mean()), "pct_high_ps_control": float(100 * high[~is_t].mean()),
                    "auc_vs_control": _auc(res["proj"][is_t], res["proj"][~is_t]), "d_vs_control": cohen_d(res["proj"][is_t], res["proj"][~is_t]), "status": status})
        rows.append(row)
        for k, (g, b, ts) in enumerate(zip(res["genes"], res["beta"], res["t_scores"])):
            sig_rows.append({"target": t, "rank": k + 1, "gene": g, "beta": float(b), "t_score": float(ts)})
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values("pct_successful_kd", ascending=False, na_position="last").reset_index(drop=True)
        summary["warnings"] = np.where(summary["n_perturbed_cells"] < 30, "low_support(<30 cells)", "")
    info = {"method": "PS (Song et al. 2025; pertps definition): t-test top-N signature, OLS projection, control-centred, clipped [0, scale_factor], max-normalized",
            "top_n_genes": pc.top_n_genes, "scale_factor": pc.scale_factor, "ps_threshold": pc.ps_threshold, "expression_cut": pc.expression_cut,
            "min_cells_per_target": pc.min_cells_per_target, "min_control_cells": pc.min_control_cells, "min_pct_expressing_control": pc.min_pct_expressing_control, "score_targets_without_gene": pc.score_targets_without_gene,
            "control_classes": groups.control_classes, "n_control_cells": groups.n_control, "n_targets_scored": int(len(summary)), "n_targets_skipped": int(len(skipped)),
            "perturbed_cells": "single_targeting", "control_cells": "single_control", "target_gene_map": gene_map, "summary_order": "pct_successful_kd descending (reference)",
            "median_pct_successful_kd": float(summary["pct_successful_kd"].median()) if not summary.empty else float("nan")}
    logger.info("PS: %d targets scored, %d skipped (%d control cells); median %.0f%% of perturbed cells classed as successful knockdown", len(summary), len(skipped), groups.n_control, info["median_pct_successful_kd"] if not summary.empty else float("nan"))
    res_ = PSResults(scores, raw, own, own_q, summary, pd.DataFrame(skipped), pd.DataFrame(sig_rows), info, quadrants, cuts, genes_used, None, None, "", pc.ps_threshold)
    if pc.compute_lda_umap and not summary.empty:
        res_.lda_umap, res_.lda_label, res_.lda_note = compute_lda_embedding(adata, groups, list(summary["target"]), cfg)
        info["lda"] = {"computed": res_.lda_umap is not None, "note": res_.lda_note, "n_pcs": pc.lda_n_pcs, "max_genes": pc.lda_max_genes, "n_cells_placed": int(np.isfinite(res_.lda_umap).all(axis=1).sum()) if res_.lda_umap is not None else 0}
    else:
        info["lda"] = {"computed": False, "note": "ps.compute_lda_umap is false" if not pc.compute_lda_umap else "no scored targets"}
    return res_


def compute_lda_embedding(adata: ad.AnnData, groups: Groups, targets: List[str], cfg: Config) -> Tuple[Optional[np.ndarray], Optional[pd.Series], str]:
    """PS_python's supervised LDA embedding (``PerturbAnalyzer.compute_lda_umap``), on the run's HVGs."""
    pc = cfg.analysis.perturbation_effects.ps
    try:
        import umap
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    except ImportError as exc:  # pragma: no cover - environment dependent
        return None, None, f"the LDA embedding needs umap-learn and scikit-learn ({exc})"
    logger.info("Building the supervised LDA embedding over %d target(s) (HVG 2000 -> scale -> PCA %d -> LDA -> UMAP)", len(targets), pc.lda_n_pcs)
    work = ad.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names))
    if pc.lda_max_genes and work.n_vars > pc.lda_max_genes:
        if "highly_variable" in adata.var.columns:
            keep = adata.var_names[adata.var["highly_variable"].to_numpy()]
            if len(keep) > pc.lda_max_genes:
                keep = keep[: pc.lda_max_genes]
        else:
            tmp = work.copy()
            sc.pp.highly_variable_genes(tmp, n_top_genes=pc.lda_max_genes)
            keep = tmp.var_names[tmp.var["highly_variable"].to_numpy()]
        work = work[:, list(keep)].copy()
    try:
        work = work.copy()
        if work.X.max() > 20:  # pertps: normalize raw counts (never the case for a log-normalized X)
            sc.pp.normalize_total(work, target_sum=1e4)
            sc.pp.log1p(work)
        sc.pp.highly_variable_genes(work, n_top_genes=min(2000, work.n_vars))
        sc.pp.scale(work, max_value=10)
        sc.tl.pca(work, n_comps=min(pc.lda_n_pcs, work.n_obs - 1, work.n_vars - 1))
        label = np.full(adata.n_obs, LDA_OTHER_LABEL, dtype=object)
        label[groups.control] = LDA_CONTROL_LABEL
        valid = []
        for t in targets:
            m = groups.mask(t)
            if m.sum() > 5:
                label[m] = t
                valid.append(t)
        valid_mask = np.isin(label, valid + [LDA_CONTROL_LABEL])
        lda = LinearDiscriminantAnalysis(solver="eigen", shrinkage="auto")
        X_lda = lda.fit_transform(work.obsm["X_pca"][valid_mask], label[valid_mask])
        reducer = umap.UMAP(n_neighbors=30, min_dist=0.01, metric="cosine", random_state=42)
        emb = reducer.fit_transform(X_lda)
    except Exception as exc:
        logger.warning("Skipping the LDA embedding: %s: %s", type(exc).__name__, exc)
        return None, None, f"the LDA embedding could not be built ({type(exc).__name__}: {exc})"
    full = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    full[np.flatnonzero(valid_mask)] = emb
    logger.info("LDA embedding: %d/%d cells placed (cells outside the trained classes have no coordinates)", int(valid_mask.sum()), adata.n_obs)
    return full, pd.Series(label.astype(str), index=adata.obs_names), ""


def compare_with_strength(ps: Optional[PSResults], strength) -> pd.DataFrame:
    """Join the per-cell score summary to the group-level knockdown test (reference ``compare_with_perturbation_strength``)."""
    if ps is None or ps.summary.empty or strength is None or getattr(strength, "empty", True) or strength.table.empty:
        return pd.DataFrame()
    primary = strength.primary_control
    keep = [c for c in ("target", f"log2fc_{primary}", f"ks_fdr_{primary}", f"is_hit_{primary}") if c in strength.table.columns]
    if len(keep) < 2:
        return pd.DataFrame()
    return ps.summary.merge(strength.table[keep], on="target", how="inner")
