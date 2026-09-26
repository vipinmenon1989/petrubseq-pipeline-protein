"""Perturbation x cluster enrichment (port of the reference ``enrichment.py``).

For every target *t* with >= ``min_cells_per_target`` single-guide cells and every
cluster *k* with >= ``min_cells_per_cluster`` cells, the 2 x 2 table

                   in cluster k    not in k
    target t             a             b
    reference            c             d

is tested with a two-sided Fisher exact test, under two reference definitions
(``controls``): ``ntc`` = single-guide cells of ``control_classes``, ``other`` =
single-guide targeting cells of every other target. ``primary_control`` (reference
default ``other``) drives ``significant``, the ranking and the figures; both arms
are in the table. Ambiguous, multi-guide, mixed and unassigned cells are never
counted.

* ``odds_ratio``: Haldane-Anscombe corrected ((a+.5)(d+.5))/((b+.5)(c+.5)), finite
  at zero counts; ``log2_odds_ratio``; ``sample_odds_ratio`` (a*d)/(b*c) is kept as
  an extra column;
* ``fdr``: Benjamini-Hochberg **within each control arm**; ``significant`` =
  ``fdr < fdr_alpha`` and ``control == primary_control``;
* ``direction``: enriched when the target's in-cluster share exceeds the
  reference's, else depleted;
* ``low_power``: fewer than ``min_reference_cells`` reference cells in the cluster;
* ``stratify_by`` (any obs column with >= 2 levels): a Cochran-Mantel-Haenszel test
  over the strata replaces the pooled Fisher p-value and, when finite, the pooled
  OR replaces the Haldane OR (reference); the Fisher p-value is kept in ``pval_fisher``;
* guide concordance for the significant pairs: guides of *t* with >=
  ``min_cells_per_guide`` cells whose own in-cluster fraction lies on the same side
  of the reference fraction as the direction (``guides_concordant`` / ``guides_tested``);
* omnibus: chi-square of the targeting-cell target x cluster table, % of expected
  counts < 5 and a seeded permutation p-value in which every target's row is
  resampled from a multinomial with the pooled cluster proportions (reference);
* per target: ``composition`` (% of the target's cells per cluster), the arm's
  ``reference_composition`` and ``effect_magnitude`` (total variation distance from
  the primary reference, number of significant clusters).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

from ..config import Config
from ._common import CLASS_CONTROL, CLASS_TARGETING, bh_fdr
from .perturbation_strength import CONTROL_LABELS, CONTROL_NTC, CONTROL_OTHER

logger = logging.getLogger(__name__)


@dataclass
class EnrichmentResults:
    table: pd.DataFrame                                   # long: one row per (target, control, cluster)
    composition: pd.DataFrame                             # targets x clusters, % of the target's cells
    reference_composition: Dict[str, pd.Series]           # control -> % per cluster
    effect_magnitude: pd.DataFrame                        # per target: composition_shift_pct, n_significant_clusters
    omnibus: Dict[str, Any] = field(default_factory=dict)
    controls_used: List[str] = field(default_factory=list)
    primary_control: str = CONTROL_OTHER
    skipped: pd.DataFrame = field(default_factory=pd.DataFrame)
    cluster_key: str = "leiden"
    stratified: bool = False
    stratify_by: Optional[str] = None
    lochness_comparison: pd.DataFrame = field(default_factory=pd.DataFrame)
    info: Dict[str, Any] = field(default_factory=dict)
    empty: bool = False

    @property
    def hits(self) -> pd.DataFrame:
        if self.table.empty:
            return self.table
        return self.table[self.table["significant"].astype(bool)]

    def top_hits(self, n: int) -> pd.DataFrame:
        h = self.hits
        if h.empty:
            return h
        return h.reindex(h["log2_odds_ratio"].abs().sort_values(ascending=False).index).head(n)

    def targets_with_hits(self) -> List[str]:
        return sorted(self.hits["target"].unique()) if not self.hits.empty else []


def cluster_order(values: Sequence[str]) -> List[str]:
    uniq = list(dict.fromkeys(str(v) for v in values))
    try:
        return sorted(uniq, key=lambda x: (float(x), x))
    except ValueError:
        return sorted(uniq)


def haldane_odds_ratio(a, b, c, d, pc: float = 0.5) -> float:
    return float(((a + pc) * (d + pc)) / ((b + pc) * (c + pc)))


def haldane_log2_or(a, b, c, d, pc: float = 0.5) -> float:
    return float(np.log2(haldane_odds_ratio(a, b, c, d, pc)))


def cmh(tables: List[np.ndarray]) -> Tuple[float, float, int]:
    """Cochran-Mantel-Haenszel pooled OR and p-value over 2x2 strata (reference ``_cmh_test``)."""
    from statsmodels.stats.contingency_tables import StratifiedTable

    usable = [t for t in tables if t.sum() > 0 and t.sum(axis=1).min() > 0 and t.sum(axis=0).min() > 0]
    if not usable:
        return float("nan"), float("nan"), 0
    st = StratifiedTable([t.T.astype(float) for t in usable])
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(st.oddsratio_pooled), float(st.test_null_odds().pvalue), len(usable)


def omnibus(contingency: pd.DataFrame, n_permutations: int, seed: int) -> Dict[str, Any]:
    """Global target x cluster association (reference ``omnibus_test``)."""
    if contingency.empty or contingency.shape[0] < 2 or contingency.shape[1] < 2:
        return {}
    table = contingency.loc[contingency.sum(axis=1) > 0, contingency.sum(axis=0) > 0]
    if table.shape[0] < 2 or table.shape[1] < 2:
        return {"status": "not computed (fewer than 2 targets or clusters with cells after pruning)"}
    observed = table.to_numpy(dtype=float)
    chi2, p_chi2, dof, expected = chi2_contingency(observed)
    out: Dict[str, Any] = {"chi2": float(chi2), "dof": int(dof), "p_chi2": float(p_chi2), "pct_expected_below_5": float(100 * (expected < 5).mean()), "n_permutations": 0, "p_permutation": float("nan"),
                           "n_dropped_rows": int(contingency.shape[0] - table.shape[0]), "n_dropped_cols": int(contingency.shape[1] - table.shape[1])}
    if n_permutations and n_permutations > 0:
        rng = np.random.default_rng(seed)
        row_counts = observed.sum(axis=1).astype(int)
        col_probs = observed.sum(axis=0) / observed.sum()
        n_ge = 0
        for _ in range(n_permutations):
            sim = np.vstack([rng.multinomial(n, col_probs) for n in row_counts])
            keep = sim.sum(axis=0) > 0
            stat = chi2_contingency(sim[:, keep])[0] if keep.sum() > 1 else 0.0
            if stat >= chi2:
                n_ge += 1
        out["n_permutations"] = int(n_permutations)
        out["p_permutation"] = (n_ge + 1) / (n_permutations + 1)
    return out


def _guide_concordance(obs: pd.DataFrame, t: str, cluster: str, cluster_key: str, ref_fraction: float, min_cells: int, direction: str) -> Tuple[int, int]:
    if "guide" not in obs.columns:
        return 0, 0
    sub = obs[(obs["target"].astype(str) == t) & (obs["perturbation_class"].astype(str) == CLASS_TARGETING)]
    if sub.empty:
        return 0, 0
    tested = concordant = 0
    for _, cells in sub.groupby(sub["guide"].astype(str), observed=True):
        if len(cells) < min_cells:
            continue
        tested += 1
        frac = float((cells[cluster_key].astype(str) == cluster).mean())
        agrees = frac < ref_fraction if direction == "depleted" else frac > ref_fraction
        if agrees:
            concordant += 1
    return concordant, tested


def compute_cluster_enrichment(adata: ad.AnnData, cfg: Config, key: str) -> EnrichmentResults:
    ec = cfg.analysis.clustering.enrichment
    obs = adata.obs
    empty = EnrichmentResults(pd.DataFrame(), pd.DataFrame(), {}, pd.DataFrame(), {}, [], ec.primary_control, pd.DataFrame(), key, False, None, pd.DataFrame(), {}, True)
    for col in ("perturbation_class", "target"):
        if col not in obs.columns:
            empty.info = {"status": f"skipped: obs['{col}'] missing (no perturbation assignment)"}
            return empty
    if key not in obs.columns:
        empty.info = {"status": f"skipped: obs['{key}'] missing (clustering must run before enrichment)"}
        return empty
    clusters_all = obs[key].astype(str).to_numpy()
    targets_col = obs["target"].astype(str).to_numpy()
    klass = obs["perturbation_class"].astype(str).to_numpy()
    cclass = obs["control_class"].astype(str).to_numpy() if "control_class" in obs.columns else np.full(adata.n_obs, "")
    ntc_mask = (klass == CLASS_CONTROL) & np.isin(cclass, list(ec.control_classes))
    targeting = klass == CLASS_TARGETING

    cluster_sizes = pd.Series(clusters_all).value_counts()
    clusters = [c for c in cluster_order(clusters_all) if cluster_sizes.get(c, 0) >= ec.min_cells_per_cluster]
    dropped_clusters = [c for c in cluster_order(clusters_all) if c not in clusters]
    if not clusters:
        empty.info = {"status": f"skipped: no cluster has at least min_cells_per_cluster={ec.min_cells_per_cluster} cells"}
        return empty
    controls_used: List[str] = []
    for control in ec.controls:
        n_ref = int(ntc_mask.sum()) if control == CONTROL_NTC else int(targeting.sum())
        if n_ref < ec.min_reference_cells:
            logger.warning("enrichment: control %r has only %d cells; dropped", control, n_ref)
            continue
        controls_used.append(control)
    if not controls_used:
        empty.info = {"status": "skipped: no usable control group for the enrichment test"}
        return empty
    primary = ec.primary_control if ec.primary_control in controls_used else controls_used[0]
    if primary != ec.primary_control:
        logger.warning("enrichment: requested primary control %r unavailable; using %r", ec.primary_control, primary)
    strat_values = None
    stratified = False
    if ec.stratify_by:
        if ec.stratify_by not in obs.columns:
            raise ValueError(f"analysis.clustering.enrichment.stratify_by {ec.stratify_by!r} is not an obs column (available: {sorted(obs.columns)[:15]} ...)")
        strat_values = obs[ec.stratify_by].astype(str).to_numpy()
        if len(set(strat_values)) < 2:
            strat_values = None
        else:
            stratified = True
    target_counts = pd.Series(targets_col[targeting]).value_counts()
    testable = sorted(target_counts[target_counts >= ec.min_cells_per_target].index)
    skipped = pd.DataFrame([{"target": t, "n_cells": int(n), "reason": f"fewer than {ec.min_cells_per_target} single-guide cells"} for t, n in target_counts.items() if n < ec.min_cells_per_target])
    if not testable:
        empty.info = {"status": "skipped: no target has enough single-guide cells"}
        empty.skipped = skipped
        return empty
    in_cluster = {c: clusters_all == c for c in clusters}
    comp_rows = {}
    for t in testable:
        m = targeting & (targets_col == t)
        n = int(m.sum())
        comp_rows[t] = {c: 100 * float((m & in_cluster[c]).sum()) / n for c in clusters}
    composition = pd.DataFrame.from_dict(comp_rows, orient="index")[clusters]
    composition.index.name = "target"
    reference_composition: Dict[str, pd.Series] = {}
    for control in controls_used:
        m = ntc_mask if control == CONTROL_NTC else targeting
        n = max(int(m.sum()), 1)
        reference_composition[control] = pd.Series({c: 100 * float((m & in_cluster[c]).sum()) / n for c in clusters})
    contingency = pd.DataFrame({c: [int((targeting & (targets_col == t) & in_cluster[c]).sum()) for t in testable] for c in clusters}, index=testable)
    omni = omnibus(contingency, ec.n_permutations, cfg.compute.seed)
    rows: List[dict] = []
    for t in testable:
        tmask = targeting & (targets_col == t)
        n_target = int(tmask.sum())
        for control in controls_used:
            rmask = ntc_mask if control == CONTROL_NTC else (targeting & (targets_col != t))
            n_ref = int(rmask.sum())
            for cluster in clusters:
                cm_ = in_cluster[cluster]
                a = int((tmask & cm_).sum()); b = n_target - a
                c_ = int((rmask & cm_).sum()); d = n_ref - c_
                pct_t = 100 * a / max(n_target, 1)
                pct_r = 100 * c_ / max(n_ref, 1)
                odds = haldane_odds_ratio(a, b, c_, d, ec.odds_pseudocount)
                _, p_fisher = fisher_exact([[a, b], [c_, d]])
                pval = float(p_fisher)
                row = {"target": t, "cluster": cluster, "control": control, "n_target_cells": n_target, "n_in_cluster": a, "pct_of_target": pct_t, "n_reference_cells": n_ref, "n_reference_in_cluster": c_, "pct_of_reference": pct_r,
                       "odds_ratio": odds, "log2_odds_ratio": float(np.log2(odds)) if odds > 0 else np.nan, "sample_odds_ratio": float((a * d) / (b * c_)) if b * c_ > 0 else (np.inf if a * d > 0 else np.nan),
                       "direction": "enriched" if pct_t > pct_r else "depleted", "pval_fisher": float(p_fisher), "low_power": bool(c_ < ec.min_reference_cells)}
                if stratified:
                    tables = []
                    for s in sorted(set(strat_values)):
                        sm = strat_values == s
                        tables.append(np.array([[int((tmask & cm_ & sm).sum()), int((tmask & ~cm_ & sm).sum())], [int((rmask & cm_ & sm).sum()), int((rmask & ~cm_ & sm).sum())]], dtype=float))
                    pooled_or, p_cmh, n_strata = cmh(tables)
                    if np.isfinite(pooled_or) and pooled_or > 0:
                        row["odds_ratio"] = pooled_or
                        row["log2_odds_ratio"] = float(np.log2(pooled_or))
                    pval = p_cmh
                    row.update({"cmh_odds_ratio": pooled_or, "cmh_pval": p_cmh, "cmh_n_strata": n_strata})
                row["pval"] = float(pval)
                rows.append(row)
    table = pd.DataFrame(rows)
    table["fdr"] = np.nan
    for control in controls_used:
        m = (table["control"] == control).to_numpy()
        table.loc[m, "fdr"] = bh_fdr(table.loc[m, "pval"].to_numpy())
    table["significant"] = (table["fdr"] < ec.fdr_alpha) & (table["control"] == primary)
    table["guides_concordant"] = np.nan
    table["guides_tested"] = np.nan
    if ec.guide_concordance and "guide" in obs.columns:
        view = obs[["guide", "target", "perturbation_class", key]].copy()
        view[key] = view[key].astype(str)
        for idx in table.index[table["significant"]]:
            conc, tested = _guide_concordance(view, str(table.at[idx, "target"]), str(table.at[idx, "cluster"]), key, float(table.at[idx, "pct_of_reference"]) / 100.0, ec.min_cells_per_guide, str(table.at[idx, "direction"]))
            table.at[idx, "guides_concordant"] = conc
            table.at[idx, "guides_tested"] = tested
    ref = reference_composition[primary]
    n_sig = table[table["significant"]]["target"].value_counts()
    magnitude = pd.DataFrame([{"target": t, "n_cells": int(target_counts[t]), "composition_shift_pct": float((composition.loc[t] - ref).abs().sum() / 2.0), "n_significant_clusters": int(n_sig.get(t, 0))} for t in testable])
    magnitude = magnitude.sort_values("composition_shift_pct", ascending=False).reset_index(drop=True)
    table = table.sort_values(["significant", "log2_odds_ratio"], ascending=[False, False]).reset_index(drop=True)
    n_hits = int(table["significant"].sum())
    info: Dict[str, Any] = {"method": "two-sided Fisher exact test per (target, control, cluster) 2x2 table" + (" replaced by a Cochran-Mantel-Haenszel test across strata" if stratified else "") + "; Haldane-Anscombe odds ratios; BH-FDR within each control arm; significance under the primary control",
                            "controls_used": controls_used, "primary_control": primary, "control_labels": {c: CONTROL_LABELS[c] for c in controls_used}, "control_classes": list(ec.control_classes), "perturbed": "single_targeting cells of the target",
                            "min_cells_per_target": ec.min_cells_per_target, "min_cells_per_cluster": ec.min_cells_per_cluster, "min_reference_cells": ec.min_reference_cells, "min_cells_per_guide": ec.min_cells_per_guide, "odds_pseudocount": ec.odds_pseudocount, "fdr_alpha": ec.fdr_alpha,
                            "stratify_by": ec.stratify_by if stratified else None, "stratified": stratified, "n_targets_tested": len(testable), "n_clusters_tested": len(clusters), "clusters_not_tested": dropped_clusters, "n_control_cells": int(ntc_mask.sum()), "n_other_cells": int(targeting.sum()),
                            "n_tests_per_control": len(testable) * len(clusters), "n_significant": n_hits, "n_targets_with_hits": int(len(n_sig)), "n_low_power": int((table["low_power"] & (table["control"] == primary)).sum()), "omnibus": omni, "status": "computed"}
    logger.info("cluster enrichment: %d targets x %d clusters = %d tests per control; %d significant at FDR < %g (control: %s)", len(testable), len(clusters), len(testable) * len(clusters), n_hits, ec.fdr_alpha, CONTROL_LABELS[primary])
    return EnrichmentResults(table, composition, reference_composition, magnitude, omni, controls_used, primary, skipped, key, stratified, ec.stratify_by if stratified else None, pd.DataFrame(), info, False)


# ---------------------------------------------------------------------------
# Derived views (reference ``enrichment_matrix`` / ``significance_matrix`` / ``phenocopy_similarity``)
# ---------------------------------------------------------------------------


def enrichment_matrix(res: EnrichmentResults, control: Optional[str] = None) -> pd.DataFrame:
    if res.table.empty:
        return pd.DataFrame()
    sub = res.table[res.table["control"] == (control or res.primary_control)]
    return sub.pivot(index="target", columns="cluster", values="log2_odds_ratio")[list(res.composition.columns)]


def significance_matrix(res: EnrichmentResults, control: Optional[str] = None) -> pd.DataFrame:
    if res.table.empty:
        return pd.DataFrame()
    sub = res.table[res.table["control"] == (control or res.primary_control)]
    return sub.pivot(index="target", columns="cluster", values="fdr")[list(res.composition.columns)]


def phenocopy_similarity(res: EnrichmentResults) -> pd.DataFrame:
    comp = res.composition
    if comp.empty or comp.shape[0] < 2:
        return pd.DataFrame()
    return comp.T.corr(method="pearson")


def format_enrichment_table(res: EnrichmentResults) -> pd.DataFrame:
    """Reader-friendly view of the primary-arm pairs with FDR < 1 (reference ``format_enrichment_table``)."""
    if res.table.empty:
        return res.table
    sub = res.table[res.table["control"] == res.primary_control].copy()
    sub = sub.reindex(sub["log2_odds_ratio"].abs().sort_values(ascending=False).index)
    sub = sub[sub["fdr"] < 1]
    cols = {"target": "Target", "cluster": "Cluster", "n_in_cluster": "Cells in cluster", "pct_of_target": "% of target", "pct_of_reference": "% of reference", "odds_ratio": "Odds ratio", "fdr": "FDR", "direction": "Direction",
            "guides_concordant": "Guides agreeing", "guides_tested": "Guides tested", "low_power": "Low power", "significant": "Significant"}
    out = sub[[c for c in cols if c in sub.columns]].rename(columns=cols).reset_index(drop=True)
    for c in out.columns:
        if out[c].dtype.kind == "f":
            out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
    return out


def compare_with_lochness(res: EnrichmentResults, lochness_summary: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Descriptive side-by-side of the strongest cluster result and the lochNESS summary per target (extension; no combined statistic)."""
    if res.empty or res.table.empty:
        return pd.DataFrame()
    prim = res.table[res.table["control"] == res.primary_control]
    counts = res.effect_magnitude.set_index("target")["n_cells"] if not res.effect_magnitude.empty else pd.Series(dtype=int)
    comp = []
    for t in sorted(prim["target"].unique()):
        sub = prim[prim["target"] == t]
        row: Dict[str, Any] = {"target": t, "n_cells": int(counts.get(t, 0)), "n_clusters_significant": int(sub["significant"].sum())}
        for dirn in ("enriched", "depleted"):
            s = sub[sub["direction"] == dirn].sort_values(["fdr", "pval"])
            if len(s):
                r = s.iloc[0]
                row.update({f"top_{dirn}_cluster": r["cluster"], f"top_{dirn}_log2_or": r["log2_odds_ratio"], f"top_{dirn}_fdr": r["fdr"]})
        if lochness_summary is not None and not lochness_summary.empty and t in set(lochness_summary["target"]):
            ls = lochness_summary.set_index("target").loc[t]
            row.update({"lochness_own_mean": float(ls["mean_lochness_in_own_cells"]), "lochness_top_cluster": str(ls.get("top_cluster", "")), "lochness_delta_own_vs_control": float(ls.get("delta_own_vs_control", np.nan)), "lochness_fdr": float(ls.get("fdr", np.nan))})
        comp.append(row)
    return pd.DataFrame(comp).sort_values(["n_clusters_significant", "target"], ascending=[False, True]).reset_index(drop=True)
