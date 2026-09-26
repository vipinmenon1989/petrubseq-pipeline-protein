"""Perturbation x cluster enrichment (after weili-lab/perturbseq-pipeline ``enrichment.py``).

For every target t with >= ``min_cells_per_target`` perturbed cells and every
cluster k with >= ``min_cells_per_cluster`` cells, the 2 x 2 table

                   in cluster k    not in k
    perturbed t          a             b
    control              c             d

is tested with a two-sided Fisher exact test (``scipy.stats.fisher_exact``,
conditional MLE-free sample odds ratio a*d / (b*c)). **Perturbed** = single-guide
targeting cells of t; **control** = single-guide cells of the configured
control classes (default non-targeting), or with ``control: other`` the
single-guide targeting cells of every other target (explicit alternative; the
reference pipeline's default). Ambiguous, multi-guide, mixed and unassigned
cells are never counted.

* ``odds_ratio``: the sample odds ratio returned by Fisher's test (0 or inf
  when a cell is 0); ``log2_or_haldane``: log2 of the Haldane-Anscombe
  corrected ratio ((a+.5)(d+.5))/((b+.5)(c+.5)), finite, for ranking/plots;
* ``fdr``: Benjamini-Hochberg over **all tested (target, cluster) pairs of the
  run** (one family per run and control definition);
* ``direction``: enriched when a/(a+b) > c/(c+d), depleted when lower, none when equal;
* ``low_power``: fewer than ``min_control_cells_in_cluster`` control cells in the cluster;
* guide support: each guide of t with >= ``min_cells_per_guide`` cells *supports*
  the pair when its own in-cluster fraction lies on the same side of the control
  fraction as the target-level direction (no per-guide test; guides are not
  independent replicates of the whole screen, this is a consistency check);
* optional ``stratify_by`` (any obs column, e.g. sample or lane): a
  Cochran-Mantel-Haenszel test over the strata (``statsmodels`` StratifiedTable,
  strata with an empty margin dropped) is reported **alongside** the Fisher
  result, with its own BH family;
* omnibus: chi-square of the perturbed-target x cluster table plus a seeded
  label-permutation p-value (a screen; many expected counts are small).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from ..config import Config
from ._common import CLASS_CONTROL, CLASS_TARGETING, bh_fdr

logger = logging.getLogger(__name__)


@dataclass
class EnrichmentResults:
    table: pd.DataFrame
    skipped: pd.DataFrame
    lochness_comparison: pd.DataFrame
    info: Dict[str, Any] = field(default_factory=dict)
    empty: bool = False


def haldane_log2_or(a, b, c, d, pc: float = 0.5) -> float:
    return float(np.log2(((a + pc) * (d + pc)) / ((b + pc) * (c + pc))))


def cmh(tables: List[np.ndarray]):
    """Cochran-Mantel-Haenszel pooled OR and p-value over 2x2 strata (None if no usable stratum)."""
    from statsmodels.stats.contingency_tables import StratifiedTable

    usable = [t for t in tables if t.sum() > 0 and t.sum(axis=1).min() > 0 and t.sum(axis=0).min() > 0]
    if len(usable) == 0:
        return float("nan"), float("nan"), 0
    st = StratifiedTable([t.astype(float) for t in usable])
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(st.oddsratio_pooled), float(st.test_null_odds().pvalue), len(usable)


def omnibus(contingency: pd.DataFrame, n_perm: int, seed: int) -> Dict[str, Any]:
    t = contingency.loc[contingency.sum(axis=1) > 0, contingency.sum(axis=0) > 0]
    if t.shape[0] < 2 or t.shape[1] < 2:
        return {"status": "not computed (fewer than 2 targets or clusters with cells)"}
    chi2, p, dof, exp = stats.chi2_contingency(t.to_numpy())
    out = {"chi2": float(chi2), "dof": int(dof), "p_chi2": float(p), "pct_expected_below_5": float(100 * np.mean(exp < 5))}
    if n_perm > 0:
        rows = np.repeat(np.arange(t.shape[0]), t.sum(axis=1).to_numpy())
        cols = np.concatenate([np.repeat(np.arange(t.shape[1]), t.iloc[i].to_numpy()) for i in range(t.shape[0])])
        rng = np.random.default_rng(seed)
        exceed = 0
        for _ in range(n_perm):
            perm = rng.permutation(cols)
            tab = np.zeros(t.shape)
            np.add.at(tab, (rows, perm), 1)
            exceed += stats.chi2_contingency(tab)[0] >= chi2
        out["n_permutations"] = n_perm
        out["p_permutation"] = float((1 + exceed) / (n_perm + 1))
    return out


def compute_cluster_enrichment(adata: ad.AnnData, cfg: Config, key: str, lochness_summary: Optional[pd.DataFrame] = None) -> EnrichmentResults:
    ec = cfg.analysis.clustering.enrichment
    obs = adata.obs
    for col in ("perturbation_class", "target"):
        if col not in obs.columns:
            return EnrichmentResults(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {"status": f"skipped: obs['{col}'] missing (no perturbation assignment)"}, True)
    labels = obs[key].astype(str).to_numpy()
    klass = obs["perturbation_class"].astype(str).to_numpy()
    target = obs["target"].astype(str).to_numpy()
    cclass = obs["control_class"].astype(str).to_numpy() if "control_class" in obs.columns else np.full(adata.n_obs, "")
    guide = obs["guide"].astype(str).to_numpy() if "guide" in obs.columns else None
    targeting = klass == CLASS_TARGETING
    ntc = (klass == CLASS_CONTROL) & np.isin(cclass, list(ec.control_classes))
    counts = pd.Series(target[targeting]).value_counts().sort_index()
    skipped = [{"target": t, "n_cells": int(n), "reason": f"fewer than {ec.min_cells_per_target} single-guide cells"} for t, n in counts.items() if n < ec.min_cells_per_target]
    targets = [t for t, n in counts.items() if n >= ec.min_cells_per_target]
    csize = pd.Series(labels).value_counts()
    clusters = sorted([c for c in csize.index if csize[c] >= ec.min_cells_per_cluster], key=lambda x: int(x) if str(x).isdigit() else str(x))
    untested_clusters = sorted(set(csize.index) - set(clusters), key=lambda x: int(x) if str(x).isdigit() else str(x))
    if ec.control == "non_targeting" and int(ntc.sum()) < ec.min_control_cells:
        return EnrichmentResults(pd.DataFrame(), pd.DataFrame(skipped), pd.DataFrame(), {"status": f"skipped: {int(ntc.sum())} control cells (< {ec.min_control_cells})"}, True)
    strat = None
    if ec.stratify_by:
        if ec.stratify_by not in obs.columns:
            raise ValueError(f"analysis.clustering.enrichment.stratify_by {ec.stratify_by!r} is not an obs column (available: {sorted(obs.columns)[:15]} ...)")
        strat = obs[ec.stratify_by].astype(str).to_numpy()
    rows = []
    for t in targets:
        tm = targeting & (target == t)
        cm_ = ntc if ec.control == "non_targeting" else targeting & (target != t)
        nt, nc = int(tm.sum()), int(cm_.sum())
        for k in clusters:
            ink = labels == k
            a, c = int((tm & ink).sum()), int((cm_ & ink).sum())
            b, d = nt - a, nc - c
            orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            ft, fc = a / nt, c / nc if nc else float("nan")
            direction = "enriched" if ft > fc else ("depleted" if ft < fc else "none")
            row = {"target": t, "cluster": k, "n_target_total": nt, "n_target_in_cluster": a, "n_control_total": nc, "n_control_in_cluster": c,
                   "target_cluster_fraction": ft, "control_cluster_fraction": fc, "odds_ratio": float(orr), "log2_or_haldane": haldane_log2_or(a, b, c, d, ec.odds_pseudocount),
                   "p_value": float(p), "direction": direction, "low_power": bool(c < ec.min_control_cells_in_cluster)}
            if guide is not None:
                g_obs = g_sup = 0
                for gd in np.unique(guide[tm]):
                    gm = tm & (guide == gd)
                    if gm.sum() < ec.min_cells_per_guide:
                        continue
                    g_obs += 1
                    gf = float((gm & ink).sum() / gm.sum())
                    if (direction == "enriched" and gf > fc) or (direction == "depleted" and gf < fc):
                        g_sup += 1
                row.update({"n_guides_observed": g_obs, "n_guides_supporting_direction": g_sup, "guide_support_fraction": g_sup / g_obs if g_obs else float("nan")})
            if strat is not None:
                tabs = []
                for s in np.unique(strat):
                    sm = strat == s
                    tabs.append(np.array([[int((tm & ink & sm).sum()), int((tm & ~ink & sm).sum())], [int((cm_ & ink & sm).sum()), int((cm_ & ~ink & sm).sum())]]))
                por, pp, nstr = cmh(tabs)
                row.update({"cmh_odds_ratio": por, "cmh_p_value": pp, "cmh_n_strata": nstr})
            rows.append(row)
    table = pd.DataFrame(rows)
    info: Dict[str, Any] = {"method": "two-sided Fisher exact test per (target, cluster) 2x2 table; BH-FDR over all tested pairs", "control": ec.control,
                            "control_classes": list(ec.control_classes) if ec.control == "non_targeting" else [], "perturbed": "single_targeting cells of the target",
                            "min_cells_per_target": ec.min_cells_per_target, "min_cells_per_cluster": ec.min_cells_per_cluster, "min_control_cells_in_cluster": ec.min_control_cells_in_cluster,
                            "min_cells_per_guide": ec.min_cells_per_guide, "odds_pseudocount_for_log2": ec.odds_pseudocount, "fdr_alpha": ec.fdr_alpha, "stratify_by": ec.stratify_by,
                            "n_targets_tested": len(targets), "n_clusters_tested": len(clusters), "clusters_not_tested": untested_clusters, "n_control_cells": int(ntc.sum()),
                            "bh_family": "all tested (target, cluster) pairs of this run", "status": "computed"}
    if table.empty:
        info["status"] = "no testable (target, cluster) pair"
        return EnrichmentResults(table, pd.DataFrame(skipped), pd.DataFrame(), info, True)
    table["fdr"] = bh_fdr(table["p_value"].to_numpy())
    table["significant"] = table["fdr"] < ec.fdr_alpha
    if strat is not None:
        table["cmh_fdr"] = bh_fdr(table["cmh_p_value"].to_numpy())
        table["cmh_significant"] = table["cmh_fdr"] < ec.fdr_alpha
    table = table.sort_values(["fdr", "p_value", "target", "cluster"]).reset_index(drop=True)
    cont = pd.crosstab(target[targeting], labels[targeting]).reindex(index=targets, columns=clusters, fill_value=0)
    info["omnibus"] = omnibus(cont, ec.n_permutations, cfg.compute.seed)
    info["n_tests"] = int(len(table))
    info["n_significant"] = int(table["significant"].sum())
    info["n_significant_enriched"] = int((table["significant"] & (table["direction"] == "enriched")).sum())
    info["n_significant_depleted"] = int((table["significant"] & (table["direction"] == "depleted")).sum())
    if strat is not None:
        info["n_cmh_significant"] = int(table["cmh_significant"].sum())
    # descriptive comparison with lochNESS (no combined statistic) --------------------------------
    comp = []
    for t in targets:
        sub = table[table["target"] == t]
        row: Dict[str, Any] = {"target": t, "n_cells": int(counts[t]), "n_clusters_significant": int(sub["significant"].sum())}
        for dirn in ("enriched", "depleted"):
            s = sub[sub["direction"] == dirn].sort_values(["fdr", "p_value"])
            if len(s):
                r = s.iloc[0]
                row.update({f"top_{dirn}_cluster": r["cluster"], f"top_{dirn}_log2_or": r["log2_or_haldane"], f"top_{dirn}_fdr": r["fdr"]})
        if lochness_summary is not None and not lochness_summary.empty and t in set(lochness_summary["target"]):
            ls = lochness_summary.set_index("target").loc[t]
            row.update({"lochness_own_mean": float(ls["mean_lochness_in_own_cells"]), "lochness_delta_own_vs_control": float(ls.get("delta_own_vs_control", np.nan)), "lochness_fdr": float(ls.get("fdr", np.nan))})
        comp.append(row)
    comp_df = pd.DataFrame(comp).sort_values(["n_clusters_significant", "target"], ascending=[False, True]).reset_index(drop=True)
    logger.info("cluster enrichment: %d targets x %d clusters, %d significant (%d enriched, %d depleted) at FDR < %g", len(targets), len(clusters), info["n_significant"], info["n_significant_enriched"], info["n_significant_depleted"], ec.fdr_alpha)
    return EnrichmentResults(table, pd.DataFrame(skipped), comp_df, info, False)
