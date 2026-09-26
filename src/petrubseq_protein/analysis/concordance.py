"""RNA-protein concordance and the integrated multimodal summary.

Three associations, kept at the level where each quantity is defined:

* **PS <-> protein** (cell level, *within* a target): Spearman rho between a
  target's perturbed cells' PS and their normalized protein value, for every
  (target, protein) with >= ``min_cells`` cells; BH-FDR over the within-target
  tests. A separate ``pooled`` row (all perturbed cells of all targets) is
  reported only as a descriptive number and flagged, because pooling unrelated
  perturbations mixes between-target with within-target variation.
* **lochNESS <-> protein** at two levels, kept apart: (cell level, *within* a
  target) Spearman between each perturbed cell's own-target lochNESS score
  (``obsm['lochness'][t]``, the local enrichment of target t around that cell)
  and its protein value, BH-FDR over the within-target tests; and (target
  level) per protein, Spearman across targets between the own-cell mean
  lochNESS and the target's protein effect (and |effect|). The target-level
  summary is never copied onto cells.
* **gene program <-> protein**: cell level over all single-guide cells
  (Spearman of program activity vs protein value, per program x protein),
  and target level (Spearman across targets of the perturbation x program
  effect vs the perturbation x protein effect).

None of these establishes mediation or causality; they describe whether the
quantities co-vary in the analysed cells. Every table carries an
``analysis_level`` column (CELL_LEVEL / TARGET_LEVEL; docs/RNA_PROTEIN_LEVELS.md).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd

from ..config import Config
from ._common import bh_fdr, cell_groups, protein_matrix, spearman
from .lochness import LochnessResults
from .modules import ModulesResults
from .protein_effects import ProteinEffectsResults
from .ps_score import PSResults

logger = logging.getLogger(__name__)

CELL_LEVEL = "CELL_LEVEL"
TARGET_LEVEL = "TARGET_LEVEL"


@dataclass
class ConcordanceResults:
    ps_protein: pd.DataFrame                # CELL_LEVEL within target (+ descriptive pooled rows)
    lochness_protein: pd.DataFrame          # TARGET_LEVEL long table (own-cell mean lochNESS, protein effect)
    lochness_protein_summary: pd.DataFrame  # TARGET_LEVEL Spearman across targets per protein
    program_protein_cells: pd.DataFrame     # CELL_LEVEL
    program_protein_targets: pd.DataFrame   # TARGET_LEVEL
    summary: pd.DataFrame                   # integrated target-level table
    info: Dict[str, Any] = field(default_factory=dict)
    lochness_protein_cells: pd.DataFrame = field(default_factory=pd.DataFrame)   # CELL_LEVEL within target


def _status(n: int, rho: float, fdr: float, alpha: float, min_n: int) -> str:
    if n < min_n or not np.isfinite(rho):
        return "insufficient_support"
    if np.isfinite(fdr) and fdr < alpha:
        return "significant"
    return "not_significant"


def compute_concordance(adata: ad.AnnData, cfg: Config, ps: Optional[PSResults], loch: Optional[LochnessResults], mods: Optional[ModulesResults], prot: Optional[ProteinEffectsResults]) -> ConcordanceResults:
    pe = cfg.analysis.perturbation_effects
    cc = pe.concordance
    P = protein_matrix(adata, pe.protein.representation)
    proteins = list(P.columns) if P is not None else []
    has = P.notna().all(axis=1).to_numpy() if P is not None else np.zeros(adata.n_obs, bool)
    groups = cell_groups(adata, pe.control_classes, 1, 1)
    targets_all = groups.targets
    # --- PS <-> protein --------------------------------------------------------------
    rows = []
    if ps is not None and not ps.summary.empty and P is not None:
        for t in ps.summary["target"]:
            pm = groups.mask(t) & has
            x = ps.own[pm].to_numpy(float)
            for p in proteins:
                rho, pv = spearman(x, P.loc[pm, p].to_numpy(float))
                rows.append({"scope": "within_target", "target": t, "protein": p, "n_cells": int(pm.sum()), "spearman_rho": rho, "p_value": pv})
        pooled = groups.control.copy() & False
        for t in ps.summary["target"]:
            pooled |= groups.mask(t)
        pooled &= has
        x = ps.own[pooled].to_numpy(float)
        for p in proteins:
            rho, pv = spearman(x, P.loc[pooled, p].to_numpy(float))
            rows.append({"scope": "pooled_all_targets", "target": "*", "protein": p, "n_cells": int(pooled.sum()), "spearman_rho": rho, "p_value": pv})
    ps_prot = pd.DataFrame(rows)
    if not ps_prot.empty:
        w = ps_prot["scope"] == "within_target"
        ps_prot["fdr"] = np.nan
        ps_prot.loc[w, "fdr"] = bh_fdr(ps_prot.loc[w, "p_value"].to_numpy())
        ps_prot["status"] = [_status(n, r, f, cc.fdr_alpha, cc.min_cells) if s == "within_target" else "descriptive_pooled" for n, r, f, s in zip(ps_prot["n_cells"], ps_prot["spearman_rho"], ps_prot["fdr"], ps_prot["scope"])]
        ps_prot["analysis_level"] = np.where(ps_prot["scope"] == "within_target", CELL_LEVEL, CELL_LEVEL + "_POOLED_DESCRIPTIVE")
        ps_prot = ps_prot.sort_values(["scope", "fdr", "target"], na_position="last").reset_index(drop=True)
    # --- lochNESS <-> protein (cell level, within target) --------------------------------
    lc_rows = []
    if loch is not None and not loch.summary.empty and P is not None and not loch.scores.empty:
        for t in loch.summary["target"]:
            if t not in loch.scores.columns or t not in groups.perturbed:
                continue
            pm = groups.mask(t) & has
            x = loch.scores.loc[pm, t].to_numpy(float)          # the cell's own-target lochNESS score (varies cell by cell)
            for p in proteins:
                rho, pv = spearman(x, P.loc[pm, p].to_numpy(float))
                lc_rows.append({"scope": "within_target", "target": t, "protein": p, "n_cells": int(pm.sum()), "spearman_rho": rho, "p_value": pv})
    loch_cells = pd.DataFrame(lc_rows)
    if not loch_cells.empty:
        loch_cells["fdr"] = bh_fdr(loch_cells["p_value"].to_numpy())
        loch_cells["status"] = [_status(n, r, f, cc.fdr_alpha, cc.min_cells) for n, r, f in zip(loch_cells["n_cells"], loch_cells["spearman_rho"], loch_cells["fdr"])]
        loch_cells["analysis_level"] = CELL_LEVEL
        loch_cells = loch_cells.sort_values(["fdr", "target"], na_position="last").reset_index(drop=True)
    # --- lochNESS <-> protein (target level) -----------------------------------------
    lp_rows, lps_rows = [], []
    if loch is not None and not loch.summary.empty and prot is not None and not prot.empty:
        ls = loch.summary.set_index("target")
        for p in proteins:
            eff = prot.matrix[p] if p in prot.matrix.columns else pd.Series(dtype=float)
            common = [t for t in eff.index if t in ls.index]
            for t in common:
                lp_rows.append({"target": t, "protein": p, "n_cells": int(ls.loc[t, "n_cells"]), "lochness_own_mean": float(ls.loc[t, "mean_lochness_in_own_cells"]),
                                "lochness_z": float(ls.loc[t, "z_score"]) if "z_score" in ls.columns else np.nan, "protein_effect": float(eff[t]),
                                "protein_fdr": float(prot.fdr_matrix.loc[t, p]) if t in prot.fdr_matrix.index else np.nan})
            if len(common) >= 5:
                x = ls.loc[common, "mean_lochness_in_own_cells"].to_numpy(float)
                y = eff[common].to_numpy(float)
                r1, p1 = spearman(x, y)
                r2, p2 = spearman(x, np.abs(y))
                lps_rows.append({"protein": p, "n_targets": len(common), "rho_lochness_vs_effect": r1, "p_lochness_vs_effect": p1, "rho_lochness_vs_abs_effect": r2, "p_lochness_vs_abs_effect": p2})
    lp = pd.DataFrame(lp_rows)
    if not lp.empty:
        lp["analysis_level"] = TARGET_LEVEL
    lps = pd.DataFrame(lps_rows)
    if not lps.empty:
        lps["fdr_signed_effect"] = bh_fdr(lps["p_lochness_vs_effect"].to_numpy())
        lps["fdr_abs_effect"] = bh_fdr(lps["p_lochness_vs_abs_effect"].to_numpy())
        lps["support"] = np.where(lps["n_targets"] >= 10, "ok", "few_targets(<10)")
        lps["analysis_level"] = TARGET_LEVEL
    # --- program <-> protein ---------------------------------------------------------
    pp_cells, pp_targets = [], []
    if mods is not None and not mods.empty and mods.program_activity is not None and P is not None:
        cells = groups.control.copy()
        for t in targets_all:
            cells |= groups.mask(t)
        cells &= has
        A = mods.program_activity
        for prog in mods.program_labels:
            for p in proteins:
                rho, pv = spearman(A.loc[cells, prog].to_numpy(float), P.loc[cells, p].to_numpy(float))
                pp_cells.append({"gene_program": prog, "protein": p, "n_cells": int(cells.sum()), "spearman_rho": rho, "p_value": pv})
        if prot is not None and not prot.empty:
            ppe = mods.perturbation_program_effects.pivot(index="target", columns="program", values="mean_log2fc")
            common = [t for t in ppe.index if t in prot.matrix.index]
            for prog in mods.program_labels:
                for p in proteins:
                    if len(common) >= 5:
                        rho, pv = spearman(ppe.loc[common, prog].to_numpy(float), prot.matrix.loc[common, p].to_numpy(float))
                    else:
                        rho, pv = np.nan, np.nan
                    pp_targets.append({"gene_program": prog, "protein": p, "n_targets": len(common), "spearman_rho": rho, "p_value": pv})
    ppc = pd.DataFrame(pp_cells)
    if not ppc.empty:
        ppc["fdr"] = bh_fdr(ppc["p_value"].to_numpy())
        ppc["status"] = [_status(n, r, f, cc.fdr_alpha, cc.min_cells) for n, r, f in zip(ppc["n_cells"], ppc["spearman_rho"], ppc["fdr"])]
        ppc["analysis_level"] = CELL_LEVEL
        ppc = ppc.sort_values("fdr", na_position="last").reset_index(drop=True)
    ppt = pd.DataFrame(pp_targets)
    if not ppt.empty:
        ppt["fdr"] = bh_fdr(ppt["p_value"].to_numpy())
        ppt["status"] = [_status(n, r, f, cc.fdr_alpha, 5) for n, r, f in zip(ppt["n_targets"], ppt["spearman_rho"], ppt["fdr"])]
        ppt["analysis_level"] = TARGET_LEVEL
        ppt = ppt.sort_values("fdr", na_position="last").reset_index(drop=True)
    # --- integrated target-level summary ---------------------------------------------
    srows = []
    for t in targets_all:
        row: Dict[str, Any] = {"target": t, "n_cells": int(groups.counts[t])}
        if ps is not None and not ps.summary.empty and t in set(ps.summary["target"]):
            r = ps.summary.set_index("target").loc[t]
            row.update({"ps_median": float(r["median_ps"]), "ps_mean": float(r["mean_ps"]), "ps_pct_high": float(r["pct_high_ps"]), "ps_auc_vs_control": float(r["auc_vs_control"]), "ps_net_pct_kd": float(r["net_pct_kd"]) if pd.notna(r["net_pct_kd"]) else np.nan})
        if loch is not None and not loch.summary.empty and t in set(loch.summary["target"]):
            r = loch.summary.set_index("target").loc[t]
            row.update({"lochness_own_mean": float(r["mean_lochness_in_own_cells"]), "lochness_z": float(r["z_score"]) if "z_score" in r else np.nan, "lochness_fdr": float(r["fdr"]) if "fdr" in r else np.nan})
        if mods is not None and not mods.empty and t in set(mods.perturbation_modules["target"]):
            r = mods.perturbation_modules.set_index("target").loc[t]
            row.update({"module": str(r["module"]), "n_de_genes": int(r["n_de_genes"]), "rna_effect_magnitude": float(r["mean_abs_log2fc"])})
            pe_t = mods.perturbation_program_effects[mods.perturbation_program_effects["target"] == t].copy()
            pe_t["abs"] = pe_t["mean_log2fc"].abs()
            top = pe_t.sort_values("abs", ascending=False).head(2)
            row["strongest_programs"] = "; ".join(f"{a}:{b:+.2f}" for a, b in zip(top["program"], top["mean_log2fc"]))
        if prot is not None and not prot.empty and t in prot.matrix.index:
            for p in proteins:
                row[f"protein_effect_{p}"] = float(prot.matrix.loc[t, p])
                row[f"protein_fdr_{p}"] = float(prot.fdr_matrix.loc[t, p])
            d = prot.d_matrix.loc[t].abs()
            row["protein_effect_magnitude"] = float(d.max())
            row["strongest_protein"] = str(d.idxmax())
            row["n_proteins_significant"] = int((prot.fdr_matrix.loc[t] < pe.protein.fdr_alpha).sum())
        if not ps_prot.empty:
            w = ps_prot[(ps_prot["scope"] == "within_target") & (ps_prot["target"] == t)]
            if not w.empty:
                best = w.loc[w["spearman_rho"].abs().idxmax()] if w["spearman_rho"].notna().any() else None
                if best is not None:
                    row["ps_protein_best"] = f"{best['protein']} rho={best['spearman_rho']:+.2f} fdr={best['fdr']:.2g}"
        if not loch_cells.empty:
            w = loch_cells[loch_cells["target"] == t]
            if not w.empty and w["spearman_rho"].notna().any():
                best = w.loc[w["spearman_rho"].abs().idxmax()]
                row["lochness_protein_best"] = f"{best['protein']} rho={best['spearman_rho']:+.2f} fdr={best['fdr']:.2g}"
        srows.append(row)
    summary = pd.DataFrame(srows)
    if not summary.empty:
        summary = summary.sort_values("n_cells", ascending=False).reset_index(drop=True)
    info = {"ps_protein": "CELL_LEVEL: Spearman within target (perturbed cells: per-cell PS vs per-cell protein), BH-FDR; pooled rows descriptive only",
            "lochness_protein_cells": "CELL_LEVEL: Spearman within target (per-cell own-target lochNESS vs per-cell protein), BH-FDR",
            "lochness_protein": "TARGET_LEVEL: Spearman across targets per protein (own-cell mean lochNESS vs protein effect and |effect|)",
            "n_lochness_protein_cell_tests": int(len(loch_cells)), "n_lochness_protein_cell_significant": int((loch_cells["status"] == "significant").sum()) if not loch_cells.empty else 0,
            "program_protein": "cell-level Spearman over single-guide + control cells; target-level Spearman of program effect vs protein effect", "min_cells": cc.min_cells, "fdr_alpha": cc.fdr_alpha,
            "n_ps_protein_tests": int((ps_prot["scope"] == "within_target").sum()) if not ps_prot.empty else 0,
            "n_ps_protein_significant": int(((ps_prot["scope"] == "within_target") & (ps_prot["status"] == "significant")).sum()) if not ps_prot.empty else 0,
            "n_program_protein_cell_tests": int(len(ppc)), "n_program_protein_cell_significant": int((ppc["status"] == "significant").sum()) if not ppc.empty else 0,
            "n_program_protein_target_significant": int((ppt["status"] == "significant").sum()) if not ppt.empty else 0, "caveat": "associations, not mediation or causality"}
    logger.info("concordance: %d PS-protein tests (%d significant), %d program-protein cell tests (%d significant)", info["n_ps_protein_tests"], info["n_ps_protein_significant"], info["n_program_protein_cell_tests"], info["n_program_protein_cell_significant"])
    return ConcordanceResults(ps_prot, lp, lps, ppc, ppt, summary, info, loch_cells)
