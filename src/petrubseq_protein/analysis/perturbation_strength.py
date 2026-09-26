"""Perturbation strength: did each guide knock its own target down?

Port of the reference ``perturbseq_pipeline/perturbation.py`` (weili-lab
perturbseq-pipeline). For every target *g* whose gene is measured in the
expression matrix, the gene's **own** log-normalized expression is compared
between the target's perturbed cells and control cells, under two control
definitions reported side by side:

``ntc``
    single-guide control cells of the configured control classes (non-targeting).
``other``
    single-guide targeting cells assigned to a *different* target.

Per (target, control):

* ``mean_lognorm_*``: mean of the log-normalized values;
* ``log2fc = log2((mean(expm1 P) + 0.01) / (mean(expm1 C) + 0.01))`` on de-logged
  means (pseudocount 0.01 in normalized-expression units), ``pct_knockdown`` =
  the same as a percentage drop;
* ``ks_stat`` / ``ks_pval``: two-sample Kolmogorov-Smirnov (two-sided);
* ``mwu_pval_less``: one-sided Mann-Whitney U, perturbed < control;
* ``ks_fdr`` / ``mwu_fdr``: Benjamini-Hochberg across the tested targets, per arm;
* ``is_hit``: ``ks_fdr < fdr_alpha`` **and** ``log2fc < max_log2fc_for_hit`` (0).

Eligibility, each recorded in ``skipped`` with its reason: gene absent from
``var`` (after ``target_gene_map``); fewer than ``min_cells_per_target`` perturbed
cells; gene expressed in fewer than ``min_pct_expressing_control`` % of the
primary-control cells. An arm with fewer than ``min_control_cells`` cells gives
NaN statistics for that target. Targets are ranked by (hit under the primary
control, most negative log2FC). Ambiguous, multi-guide and unassigned cells are
never perturbed cells nor controls.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu

from ..config import Config
from ._common import CLASS_CONTROL, CLASS_TARGETING, bh_fdr, dense

logger = logging.getLogger(__name__)

CONTROL_NTC = "ntc"
CONTROL_OTHER = "other"
CONTROL_LABELS = {CONTROL_NTC: "non-targeting control cells", CONTROL_OTHER: "cells assigned to other target genes"}
#: pseudocount (normalized-expression units) added to both de-logged group means (reference ``_PSEUDOCOUNT``)
PSEUDOCOUNT = 0.01


@dataclass
class StrengthResults:
    table: pd.DataFrame                      # one row per tested target, reference columns, ``rank`` first
    display: pd.DataFrame                    # reader-friendly view for the report
    controls_used: List[str]
    primary_control: str
    skipped: pd.DataFrame
    n_control_cells: Dict[str, int]
    info: Dict[str, Any] = field(default_factory=dict)
    empty: bool = False

    @property
    def hits(self) -> pd.DataFrame:
        col = f"is_hit_{self.primary_control}"
        if self.table.empty or col not in self.table.columns:
            return self.table.iloc[0:0]
        return self.table[self.table[col].astype(bool)]

    def top_effects(self, n: int) -> pd.DataFrame:
        return self.table.head(n)


def control_masks(adata: ad.AnnData, control_classes) -> Dict[str, np.ndarray]:
    obs = adata.obs
    klass = obs["perturbation_class"].astype(str).to_numpy()
    cclass = obs["control_class"].astype(str).to_numpy() if "control_class" in obs.columns else np.full(adata.n_obs, "")
    return {CONTROL_NTC: (klass == CLASS_CONTROL) & np.isin(cclass, list(control_classes)), CONTROL_OTHER: klass == CLASS_TARGETING}


def control_mask_for_target(control: str, base: Dict[str, np.ndarray], targets: np.ndarray, t: str) -> np.ndarray:
    if control == CONTROL_NTC:
        return base[CONTROL_NTC]
    return base[CONTROL_OTHER] & (targets != t)


def compare_groups(perturbed: np.ndarray, control: np.ndarray) -> Dict[str, float]:
    """Effect size and significance for one target / control pair (reference ``compare_groups``)."""
    mean_p_log = float(np.mean(perturbed)) if perturbed.size else np.nan
    mean_c_log = float(np.mean(control)) if control.size else np.nan
    mean_p = float(np.mean(np.expm1(perturbed))) if perturbed.size else np.nan
    mean_c = float(np.mean(np.expm1(control))) if control.size else np.nan
    log2fc = float(np.log2((mean_p + PSEUDOCOUNT) / (mean_c + PSEUDOCOUNT)))
    pct_kd = float(100.0 * (1.0 - (mean_p + PSEUDOCOUNT) / (mean_c + PSEUDOCOUNT)))
    ks_stat, ks_p = ks_2samp(perturbed, control)
    try:
        mwu_p = float(mannwhitneyu(perturbed, control, alternative="less").pvalue)
    except ValueError:
        mwu_p = np.nan
    return {
        "mean_lognorm_perturbed": mean_p_log,
        "mean_lognorm_control": mean_c_log,
        "log2fc": log2fc,
        "pct_knockdown": pct_kd,
        "pct_cells_expressing_perturbed": float(100 * np.mean(perturbed > 0)) if perturbed.size else np.nan,
        "pct_cells_expressing_control": float(100 * np.mean(control > 0)) if control.size else np.nan,
        "ks_stat": float(ks_stat),
        "ks_pval": float(ks_p),
        "mwu_pval_less": mwu_p,
    }


def compute_perturbation_strength(adata: ad.AnnData, cfg: Config) -> StrengthResults:
    pe = cfg.analysis.perturbation_effects
    sc_ = pe.strength
    obs = adata.obs
    empty = StrengthResults(pd.DataFrame(), pd.DataFrame(), [], sc_.primary_control, pd.DataFrame(), {}, {}, True)
    for col in ("perturbation_class", "target"):
        if col not in obs.columns:
            empty.info = {"status": f"skipped: obs['{col}'] missing (no perturbation assignment)"}
            return empty
    targets_col = obs["target"].astype(str).to_numpy()
    klass = obs["perturbation_class"].astype(str).to_numpy()
    base = control_masks(adata, pe.control_classes)
    n_control_cells = {CONTROL_NTC: int(base[CONTROL_NTC].sum()), CONTROL_OTHER: int(base[CONTROL_OTHER].sum())}
    controls_used: List[str] = []
    for control in sc_.controls:
        if control == CONTROL_NTC and n_control_cells[CONTROL_NTC] < sc_.min_control_cells:
            logger.warning("Only %d non-targeting control cells (need %d); the 'ntc' arm of the perturbation-strength test is skipped", n_control_cells[CONTROL_NTC], sc_.min_control_cells)
            continue
        controls_used.append(control)
    if not controls_used:
        empty.info = {"status": f"skipped: no usable control group ({n_control_cells[CONTROL_NTC]} non-targeting, {n_control_cells[CONTROL_OTHER]} other-target cells; min_control_cells {sc_.min_control_cells})", "n_control_cells": n_control_cells}
        logger.warning("perturbation strength: %s", empty.info["status"])
        return empty
    primary = sc_.primary_control if sc_.primary_control in controls_used else controls_used[0]
    if primary != sc_.primary_control:
        logger.warning("Requested primary control %r is unavailable; using %r", sc_.primary_control, primary)
    gene_map = dict(pe.target_gene_map)
    all_targets = sorted(set(targets_col[klass == CLASS_TARGETING]))
    measured = set(adata.var_names)
    rows: List[dict] = []
    skipped: List[dict] = []
    for t in all_targets:
        gene = gene_map.get(t, t)
        pert_mask = (targets_col == t) & (klass == CLASS_TARGETING)
        n_pert = int(pert_mask.sum())
        if gene not in measured:
            skipped.append({"target": t, "gene": gene, "n_perturbed": n_pert, "reason": "target gene not present in the expression matrix"})
            continue
        if n_pert < sc_.min_cells_per_target:
            skipped.append({"target": t, "gene": gene, "n_perturbed": n_pert, "reason": f"fewer than {sc_.min_cells_per_target} perturbed cells"})
            continue
        values = dense(adata.X, cols=np.array([adata.var_names.get_loc(gene)])).ravel()
        primary_mask = control_mask_for_target(primary, base, targets_col, t)
        pct_expressing = float(100 * np.mean(values[primary_mask] > 0)) if primary_mask.any() else 0.0
        if pct_expressing < sc_.min_pct_expressing_control:
            skipped.append({"target": t, "gene": gene, "n_perturbed": n_pert, "reason": f"not detectably expressed in control cells ({pct_expressing:.2f}% of control cells, threshold {sc_.min_pct_expressing_control}%)"})
            continue
        row: Dict[str, object] = {"target": t, "gene": gene, "n_perturbed": n_pert}
        for control in controls_used:
            cmask = control_mask_for_target(control, base, targets_col, t)
            n_ctrl = int(cmask.sum())
            row[f"n_control_{control}"] = n_ctrl
            if n_ctrl < sc_.min_control_cells:
                for key in ("mean_lognorm_perturbed", "mean_lognorm_control", "log2fc", "pct_knockdown", "pct_cells_expressing_perturbed", "pct_cells_expressing_control", "ks_stat", "ks_pval", "mwu_pval_less"):
                    row[f"{key}_{control}"] = np.nan
                continue
            for key, val in compare_groups(values[pert_mask], values[cmask]).items():
                row[f"{key}_{control}"] = val
        rows.append(row)
    info: Dict[str, Any] = {"method": "target's own log-normalized expression, perturbed vs control: log2FC of de-logged means (pseudocount 0.01), two-sided KS and one-sided MWU, BH-FDR per control arm; hit = KS FDR < alpha and log2FC < max_log2fc_for_hit",
                            "controls_used": controls_used, "primary_control": primary, "control_labels": {c: CONTROL_LABELS[c] for c in controls_used}, "n_control_cells": n_control_cells,
                            "min_cells_per_target": sc_.min_cells_per_target, "min_control_cells": sc_.min_control_cells, "min_pct_expressing_control": sc_.min_pct_expressing_control,
                            "fdr_alpha": sc_.fdr_alpha, "max_log2fc_for_hit": sc_.max_log2fc_for_hit, "n_targets": len(all_targets), "n_tested": len(rows), "n_skipped": len(skipped), "target_gene_map": gene_map}
    if not rows:
        info["status"] = "no target gene was testable"
        logger.warning("perturbation strength: %s", info["status"])
        return StrengthResults(pd.DataFrame(), pd.DataFrame(), controls_used, primary, pd.DataFrame(skipped), n_control_cells, info, True)
    table = pd.DataFrame(rows)
    for control in controls_used:
        table[f"ks_fdr_{control}"] = bh_fdr(table[f"ks_pval_{control}"].to_numpy())
        table[f"mwu_fdr_{control}"] = bh_fdr(table[f"mwu_pval_less_{control}"].to_numpy())
        table[f"is_hit_{control}"] = ((table[f"ks_fdr_{control}"] < sc_.fdr_alpha) & (table[f"log2fc_{control}"] < sc_.max_log2fc_for_hit)).fillna(False).astype(bool)
    table = table.sort_values([f"is_hit_{primary}", f"log2fc_{primary}"], ascending=[False, True]).reset_index(drop=True)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    n_hits = int(table[f"is_hit_{primary}"].sum())
    info.update({"status": "computed", "n_effective": n_hits})
    logger.info("perturbation strength: %d/%d targets tested, %d effective at FDR < %.2f (control: %s)", len(table), len(all_targets), n_hits, sc_.fdr_alpha, CONTROL_LABELS[primary])
    res = StrengthResults(table, pd.DataFrame(), controls_used, primary, pd.DataFrame(skipped), n_control_cells, info, False)
    res.display = format_results_table(res)
    return res


def format_results_table(res: StrengthResults) -> pd.DataFrame:
    """Reader-friendly view (reference ``format_results_table``)."""
    if res.table.empty:
        return res.table
    primary = res.primary_control
    cols = {"rank": "Rank", "target": "Target", "gene": "Gene", "n_perturbed": "Perturbed cells", f"n_control_{primary}": "Control cells", f"log2fc_{primary}": "log2FC", f"pct_knockdown_{primary}": "% knockdown",
            f"ks_stat_{primary}": "KS stat", f"ks_fdr_{primary}": "KS FDR", f"is_hit_{primary}": "Effective"}
    other = CONTROL_OTHER if primary == CONTROL_NTC else CONTROL_NTC
    if other in res.controls_used:
        cols[f"log2fc_{other}"] = f"log2FC ({other})"
        cols[f"ks_fdr_{other}"] = f"KS FDR ({other})"
    present = {k: v for k, v in cols.items() if k in res.table.columns}
    out = res.table[list(present)].rename(columns=present).copy()
    if "Gene" in out.columns and (out["Gene"] == out["Target"]).all():
        out = out.drop(columns=["Gene"])
    for c in out.columns:
        if out[c].dtype.kind == "f":
            out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
    return out
