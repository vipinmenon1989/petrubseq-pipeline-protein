"""Perturbation x response-gene effect matrix, gene programs and perturbation modules.

After Zhou et al. 2023 (Nature 624:154, "regulome" map) as implemented in
weili-lab/perturbseq-pipeline ``modules.py``:

* **effect matrix** E[t, g] = log2((mean_P(g) + pc) / (mean_C(g) + pc)) on
  de-logged (expm1) means, perturbed cells P of target t vs control cells C,
  over a *response-gene panel*; a Welch t-test on the log values (BH within
  each target) gives the DE gate |log2FC| > ``de_lfc_threshold`` and FDR < alpha;
* **gene programs** (``P1..``): hierarchical clustering (``linkage_method``) of
  the genes on 1 - Pearson correlation of their effect profiles across
  perturbations;
* **perturbation modules** (``M1..``): the same on the perturbations with
  1 - Spearman correlation of their effect profiles across genes;
* module x program strength = mean log2FC of the program's genes under the
  module's perturbations; perturbation x program effect = the same per target.

Gene panel: the reference takes Leiden cluster markers; this pipeline has no
clustering, so ``gene_selection: response`` takes the union over targets of the
top ``n_response_genes_per_perturbation`` genes by |t| among FDR < ``response_fdr``
(the genes that respond to at least one perturbation), falling back to HVGs
when the union is too small (``hvg`` selects HVGs directly).

Deliberate differences from the reference (documented in docs/PERTURBATION_EFFECTS.md):
``log2fc_pseudocount`` defaults to 1 (Seurat ``FoldChange`` convention) instead of
1e-9, and panel genes must be detected in ``min_pct_cells_expressing`` (5 %) of the
analysed cells; with the reference values genes undetected in one small group get
|log2FC| > 20 and form their own programs. Both are parameters, and
``scripts/compare_reference_perturbation.py`` sets them to the reference values.

Per-cell program activity = mean z-scored log expression of the program's
genes (deterministic; ``obsm['program_activity']``). Labels are numeric only;
nothing is annotated biologically.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd

from ..config import Config
from ._common import bh_fdr, cell_groups, dense, group_mean_var, welch_t

logger = logging.getLogger(__name__)


@dataclass
class ModulesResults:
    effect: pd.DataFrame                    # targets x genes log2FC
    de_mask: pd.DataFrame                   # targets x genes bool
    gene_programs: pd.DataFrame             # gene, program, program_size, mean_abs_log2fc, n_targets_de
    perturbation_modules: pd.DataFrame      # target, module, n_cells, n_de_genes, n_up, n_down
    module_program: pd.DataFrame            # modules x programs strength
    perturbation_program_effects: pd.DataFrame  # long: target, program, mean_log2fc, frac_de, n_genes
    program_activity: Optional[pd.DataFrame]    # cells x programs
    program_labels: List[str]
    module_labels: List[str]
    gene_order: List[str]
    target_order: List[str]
    skipped: pd.DataFrame
    info: Dict[str, Any] = field(default_factory=dict)
    empty: bool = False


def cluster_axis(items: pd.DataFrame, correlation: str, linkage_method: str, k: Optional[int], threshold: float, prefix: str) -> Tuple[pd.Series, List[str], List[str]]:
    """Hierarchical clustering of the rows of ``items`` on 1 - correlation; labels in leaf order."""
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
    from scipy.spatial.distance import squareform

    names = list(items.index)
    n = len(names)
    if n < 2:
        lab = pd.Series([f"{prefix}1"] * n, index=names)
        return lab, names, [f"{prefix}1"] if n else []
    corr = items.T.corr(method=correlation).to_numpy()
    dist = 1.0 - corr
    dist = np.nan_to_num(dist, nan=1.0, posinf=2.0, neginf=0.0)
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    dist[dist < 0] = 0.0
    Z = linkage(squareform(dist, checks=False), method=linkage_method)
    raw = fcluster(Z, max(2, min(k, n)), criterion="maxclust") if k is not None else fcluster(Z, threshold, criterion="distance")
    order = leaves_list(Z)
    remap: Dict[int, str] = {}
    for i in order:
        c = int(raw[i])
        if c not in remap:
            remap[c] = f"{prefix}{len(remap) + 1}"
    labels = pd.Series([remap[int(c)] for c in raw], index=names)
    ordered_labels = [remap[c] for c in dict.fromkeys(int(raw[i]) for i in order)]
    return labels, [names[i] for i in order], ordered_labels


def compute_modules(adata: ad.AnnData, cfg: Config) -> ModulesResults:
    pe = cfg.analysis.perturbation_effects
    mc = pe.modules
    groups = cell_groups(adata, pe.control_classes, mc.min_cells_per_perturbation, 1)
    skipped = list(groups.skipped)
    empty = ModulesResults(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None, [], [], [], [], pd.DataFrame(skipped), {}, True)
    if len(groups.targets) < mc.min_perturbations or groups.n_control == 0:
        empty.info = {"status": f"skipped: {len(groups.targets)} perturbations with >= {mc.min_cells_per_perturbation} cells (need {mc.min_perturbations}); {groups.n_control} control cells"}
        logger.info("modules: %s", empty.info["status"])
        return empty
    X = adata.X
    genes_all = adata.var_names
    mean_c, var_c, n_c = group_mean_var(X, groups.control)
    # per-target Welch t-test on the log values, all genes ------------------------------------
    t_rows, p_rows, fdr_rows, mean_rows = [], [], [], []
    for t in groups.targets:
        m, v, npert = group_mean_var(X, groups.mask(t))
        tt, pp = welch_t(m, v, npert, mean_c, var_c, n_c)
        t_rows.append(tt); p_rows.append(pp); fdr_rows.append(bh_fdr(pp)); mean_rows.append(m)
    T = np.vstack(t_rows); FDR = np.vstack(fdr_rows); M = np.vstack(mean_rows)
    # gene panel ---------------------------------------------------------------------------------
    panel: List[str] = []
    selection = mc.gene_selection
    if selection == "response":
        chosen: Dict[str, None] = {}
        for i, t in enumerate(groups.targets):
            sig = np.flatnonzero(FDR[i] < mc.response_fdr)
            top = sig[np.argsort(-np.abs(T[i][sig]))][: mc.n_response_genes_per_perturbation]
            for j in top:
                chosen.setdefault(str(genes_all[j]), None)
        panel = list(chosen)
        if len(panel) < mc.min_genes:
            logger.warning("modules: only %d response genes at FDR < %g; falling back to highly variable genes", len(panel), mc.response_fdr)
            selection = "hvg (fallback)"
            panel = []
    if not panel:
        panel = genes_all[adata.var["highly_variable"].to_numpy()].tolist() if "highly_variable" in adata.var else genes_all.tolist()
    if len(panel) < mc.min_genes:
        empty.info = {"status": f"skipped: only {len(panel)} panel genes (need {mc.min_genes})"}
        return empty
    if mc.min_pct_cells_expressing > 0:
        used = groups.control.copy()
        for t in groups.targets:
            used |= groups.mask(t)
        Xu = X[used]
        det = np.asarray((Xu > 0).mean(axis=0)).ravel() * 100 if hasattr(Xu, "tocsr") else (np.asarray(Xu) > 0).mean(axis=0) * 100
        keep = pd.Index(genes_all)[det >= mc.min_pct_cells_expressing]
        n_before = len(panel)
        panel = [g for g in panel if g in set(keep)]
        n_dropped_expr = n_before - len(panel)
    else:
        n_dropped_expr = 0
    if len(panel) < mc.min_genes:
        empty.info = {"status": f"skipped: only {len(panel)} panel genes after the expression filter (need {mc.min_genes})"}
        return empty
    gidx = genes_all.get_indexer(panel)
    eps = mc.log2fc_pseudocount
    # effect matrix on de-logged means -----------------------------------------------------------
    lin_c = np.asarray(np.expm1(dense(X, rows=groups.control, cols=gidx)).mean(axis=0)).ravel()
    eff = np.zeros((len(groups.targets), len(panel)))
    for i, t in enumerate(groups.targets):
        lin_p = np.expm1(dense(X, rows=groups.mask(t), cols=gidx)).mean(axis=0)
        eff[i] = np.log2((lin_p + eps) / (lin_c + eps))
    effect = pd.DataFrame(eff, index=groups.targets, columns=panel)
    de = pd.DataFrame((np.abs(eff) > mc.de_lfc_threshold) & (FDR[:, gidx] < mc.de_fdr_alpha), index=groups.targets, columns=panel)
    sign = np.sign(eff)
    # clustering -------------------------------------------------------------------------------
    gene_program, gene_order, program_labels = cluster_axis(effect.T, mc.program_correlation, mc.linkage_method, mc.n_programs, mc.cluster_distance_threshold, "P")
    pert_module, target_order, module_labels = cluster_axis(effect, mc.module_correlation, mc.linkage_method, mc.n_modules, mc.cluster_distance_threshold, "M")
    sizes = gene_program.value_counts()
    gene_programs = pd.DataFrame({"gene": gene_program.index, "program": gene_program.to_numpy(), "program_size": [int(sizes[p]) for p in gene_program.to_numpy()],
                                  "mean_abs_log2fc": np.abs(eff).mean(axis=0), "mean_log2fc": eff.mean(axis=0), "n_targets_de": de.sum(axis=0).to_numpy(),
                                  "n_targets_up": ((sign > 0) & de.to_numpy()).sum(axis=0), "n_targets_down": ((sign < 0) & de.to_numpy()).sum(axis=0)})
    gene_programs = gene_programs.sort_values(["program", "mean_abs_log2fc"], ascending=[True, False], ignore_index=True)
    modules_tbl = pd.DataFrame({"target": pert_module.index, "module": pert_module.to_numpy(), "n_cells": [int(groups.counts[t]) for t in pert_module.index],
                                "n_de_genes": de.sum(axis=1).to_numpy(), "n_up": ((sign > 0) & de.to_numpy()).sum(axis=1), "n_down": ((sign < 0) & de.to_numpy()).sum(axis=1),
                                "mean_abs_log2fc": np.abs(eff).mean(axis=1)}).sort_values(["module", "n_de_genes"], ascending=[True, False], ignore_index=True)
    mp = pd.DataFrame(index=module_labels, columns=program_labels, dtype=float)
    ppe_rows = []
    for p in program_labels:
        pg = gene_program.index[gene_program == p]
        cols = effect.columns.get_indexer(pg)
        for t in groups.targets:
            i = groups.targets.index(t)
            ppe_rows.append({"target": t, "program": p, "n_genes": int(len(pg)), "mean_log2fc": float(eff[i, cols].mean()), "frac_de": float(de.to_numpy()[i, cols].mean()), "n_de": int(de.to_numpy()[i, cols].sum())})
        for m in module_labels:
            perts = pert_module.index[pert_module == m]
            mp.loc[m, p] = float(effect.loc[perts, pg].to_numpy().mean())
    ppe = pd.DataFrame(ppe_rows)
    # per-cell program activity -------------------------------------------------------------------
    activity = None
    if mc.score_programs:
        Z = dense(X, cols=gidx, dtype=np.float32)
        mu, sd = Z.mean(axis=0), Z.std(axis=0)
        sd[sd == 0] = 1.0
        Z = (Z - mu) / sd
        act = {p: Z[:, effect.columns.get_indexer(gene_program.index[gene_program == p])].mean(axis=1) for p in program_labels}
        activity = pd.DataFrame(act, index=adata.obs_names).astype(np.float32)
    info = {"method": "perturbation x gene log2FC-of-means vs control; programs = genes clustered on 1-Pearson, modules = perturbations on 1-Spearman (Zhou et al. 2023 design as in weili-lab/perturbseq-pipeline)",
            "gene_selection": selection, "n_panel_genes": len(panel), "n_panel_genes_dropped_low_expression": int(n_dropped_expr), "min_pct_cells_expressing": mc.min_pct_cells_expressing, "log2fc_pseudocount": eps, "n_response_genes_per_perturbation": mc.n_response_genes_per_perturbation, "response_fdr": mc.response_fdr,
            "n_perturbations": len(groups.targets), "n_programs": len(program_labels), "n_modules": len(module_labels), "n_programs_requested": mc.n_programs, "n_modules_requested": mc.n_modules,
            "cluster_distance_threshold": mc.cluster_distance_threshold, "linkage_method": mc.linkage_method, "program_correlation": mc.program_correlation, "module_correlation": mc.module_correlation,
            "de_lfc_threshold": mc.de_lfc_threshold, "de_fdr_alpha": mc.de_fdr_alpha, "control_classes": groups.control_classes, "n_control_cells": groups.n_control,
            "program_activity": "mean z-scored log expression of program genes" if activity is not None else "not computed", "status": "computed"}
    logger.info("modules: %d perturbations x %d genes -> %d programs, %d modules", len(groups.targets), len(panel), len(program_labels), len(module_labels))
    return ModulesResults(effect, de, gene_programs, modules_tbl, mp, ppe, activity, program_labels, module_labels, gene_order, target_order, pd.DataFrame(skipped), info, False)
