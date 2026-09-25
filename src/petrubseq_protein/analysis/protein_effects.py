"""Perturbation effects on each measured protein (ADT).

For every target t with enough perturbed cells and every protein p in the
normalized protein matrix (``obsm['protein']``, CLR by default; **never** the
raw counts), compare the perturbed cells of t with the control cells:

* ``effect`` = mean_P - mean_C of the normalized values. On the CLR scale a
  difference of means is a log fold change of the geometric-mean-normalized
  abundance (natural log), so it is the natural effect size here;
* ``cohen_d`` = standardized mean difference (pooled SD);
* ``p_value`` = two-sided Mann-Whitney U (no distributional assumption);
  ``fdr`` = Benjamini-Hochberg over every (target, protein) test in the run;
* sample structure: when ``obs['sample']`` exists the effect is recomputed per
  sample (samples with >= 5 perturbed and >= 5 control cells);
  ``n_samples_tested`` / ``n_samples_same_sign`` and ``single_sample_support``
  say whether one sample carries the result; likewise per guide (guides with
  >= 5 cells): ``n_guides_tested`` / ``n_guides_same_sign`` — a real effect
  should appear across independent guides of one target.
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
from ._common import bh_fdr, cell_groups, cohen_d, protein_matrix

logger = logging.getLogger(__name__)


@dataclass
class ProteinEffectsResults:
    table: pd.DataFrame            # long: target x protein
    matrix: pd.DataFrame           # targets x proteins effect (mean difference)
    d_matrix: pd.DataFrame         # targets x proteins Cohen's d
    fdr_matrix: pd.DataFrame
    skipped: pd.DataFrame
    proteins: List[str]
    info: Dict[str, Any] = field(default_factory=dict)
    empty: bool = False


def compute_protein_effects(adata: ad.AnnData, cfg: Config) -> ProteinEffectsResults:
    pe = cfg.analysis.perturbation_effects
    pc = pe.protein
    P = protein_matrix(adata, pc.representation)
    empty = ProteinEffectsResults(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], {}, True)
    if P is None:
        empty.info = {"status": f"skipped: no normalized protein matrix obsm[{pc.representation!r}]"}
        logger.info("protein effects: %s", empty.info["status"])
        return empty
    groups = cell_groups(adata, pe.control_classes, pc.min_cells_per_target, pc.min_control_cells)
    has = P.notna().all(axis=1).to_numpy()
    sample = adata.obs["sample"].astype(str).to_numpy() if "sample" in adata.obs.columns else None
    guide = adata.obs["guide"].astype(str).to_numpy() if "guide" in adata.obs.columns else None
    proteins = list(P.columns)
    ctrl = groups.control & has
    rows: List[dict] = []
    for t in groups.targets:
        pm = groups.mask(t) & has
        n_t, n_c = int(pm.sum()), int(ctrl.sum())
        if n_t < pc.min_cells_per_target:
            groups.skipped.append({"target": t, "n_cells": n_t, "reason": f"fewer than {pc.min_cells_per_target} perturbed cells with protein values"})
            continue
        for p in proteins:
            a, b = P.loc[pm, p].to_numpy(float), P.loc[ctrl, p].to_numpy(float)
            try:
                pval = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
            except ValueError:
                pval = float("nan")
            row = {"target": t, "protein": p, "n_perturbed": n_t, "n_control": n_c, "mean_perturbed": float(a.mean()), "mean_control": float(b.mean()),
                   "median_perturbed": float(np.median(a)), "median_control": float(np.median(b)), "effect": float(a.mean() - b.mean()), "cohen_d": cohen_d(a, b), "p_value": pval}
            if sample is not None:
                signs, n_tested = [], 0
                for smp in np.unique(sample[pm | ctrl]):
                    ma, mb = pm & (sample == smp), ctrl & (sample == smp)
                    if ma.sum() >= 5 and mb.sum() >= 5:
                        n_tested += 1
                        signs.append(np.sign(P.loc[ma, p].mean() - P.loc[mb, p].mean()))
                row["n_samples_tested"] = n_tested
                row["n_samples_same_sign"] = int(sum(s == np.sign(row["effect"]) for s in signs)) if n_tested else 0
                row["single_sample_support"] = bool(n_tested <= 1)
            if guide is not None:
                gsign, g_tested = [], 0
                for gd in np.unique(guide[pm]):
                    mg = pm & (guide == gd)
                    if mg.sum() >= 5:
                        g_tested += 1
                        gsign.append(np.sign(P.loc[mg, p].mean() - b.mean()))
                row["n_guides_tested"] = g_tested
                row["n_guides_same_sign"] = int(sum(s == np.sign(row["effect"]) for s in gsign))
            rows.append(row)
    table = pd.DataFrame(rows)
    if table.empty:
        empty.info = {"status": "skipped: no target with enough perturbed cells carrying protein values", "proteins": proteins}
        empty.skipped = pd.DataFrame(groups.skipped)
        return empty
    table["fdr"] = bh_fdr(table["p_value"].to_numpy())
    table["significant"] = table["fdr"] < pc.fdr_alpha
    table["status"] = np.where(table["n_perturbed"] < 30, "low_support(<30 cells)", "ok")
    table = table.sort_values(["fdr", "target", "protein"]).reset_index(drop=True)
    mat = table.pivot(index="target", columns="protein", values="effect").reindex(columns=proteins)
    dmat = table.pivot(index="target", columns="protein", values="cohen_d").reindex(columns=proteins)
    fmat = table.pivot(index="target", columns="protein", values="fdr").reindex(columns=proteins)
    info = {"method": "mean difference of normalized protein values (CLR: log fold change) with Cohen's d and two-sided Mann-Whitney U; BH-FDR over all target x protein tests",
            "representation": pc.representation, "proteins": proteins, "n_targets": int(mat.shape[0]), "n_tests": int(len(table)), "n_significant": int(table["significant"].sum()),
            "fdr_alpha": pc.fdr_alpha, "min_cells_per_target": pc.min_cells_per_target, "min_control_cells": pc.min_control_cells, "n_control_cells": int(ctrl.sum()),
            "control_classes": groups.control_classes, "status": "computed"}
    logger.info("protein effects: %d targets x %d proteins, %d significant at FDR < %g", mat.shape[0], len(proteins), info["n_significant"], pc.fdr_alpha)
    return ProteinEffectsResults(table, mat, dmat, fmat, pd.DataFrame(groups.skipped), proteins, info, False)
