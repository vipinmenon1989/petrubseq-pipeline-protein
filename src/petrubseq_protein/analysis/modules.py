"""Perturbation x gene effect matrix, co-functional modules and gene programs.

Port of the reference weili-lab/perturbseq-pipeline ``modules.py`` (the network
analysis of Zhou et al. 2023, Nature 624:154; the reference README calls it
"Chen et al. 2023"). Perturbations with >= ``min_cells_per_perturbation`` cells,
a gene panel, a **perturbation x gene log2FC-vs-control matrix**, then

* **gene programs** (``P1..``): hierarchical clustering (``linkage_method``) of the
  genes on 1 - Pearson correlation of their effect profiles across perturbations;
* **co-functional modules** (``M1..``): the same on the perturbations with
  1 - Spearman correlation across genes;
* a signed module x program strength matrix (mean program-gene log2FC per
  perturbation, averaged over the module) and a TF-hub / module-module network.

Gene panel (``gene_selection``): ``cluster_markers`` (reference default) = union of
the top ``n_marker_genes_per_cluster`` positive markers of each cell-state
cluster (``sc.tl.rank_genes_groups`` on ``obs[cluster_key]``, ``marker_method``,
ranked by log fold change); ``response`` = union of each perturbation's top
response genes (Welch t-test vs control, BH within perturbation; an option of
this pipeline); ``hvg`` = highly variable genes. With ``min_pct_cells_expressing``
> 0 (default 0 = reference) panel genes must be detected in that share of the
analysed cells.

Effect value: ``log2((mean_P(g) + pc) / (mean_C(g) + pc))`` on de-logged
(``expm1``) means over the log-normalized values, pseudocount ``pc`` =
``log2fc_pseudocount`` (reference 1e-9). Control ``ntc`` = single-guide cells of
the configured control classes (falls back to ``other``, the leave-one-target-out
mean over all single-guide targeting cells, when there are none). A per-gene
Welch t-test on the log values, BH-corrected within each perturbation, gives
the DE mask ``|log2FC| > de_lfc_threshold & fdr < de_fdr_alpha`` used for hub
sizes and network edges. Per-cell program activity: ``sc.tl.score_genes``
(``ctrl_size`` 50, reference) or the mean z-scored expression (``zscore``).
Labels are numeric; nothing is annotated biologically.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd

from ..config import Config
from ._common import CLASS_TARGETING, bh_fdr, cell_groups, dense, group_mean_var, welch_t
from .perturbation_strength import CONTROL_LABELS, CONTROL_NTC, CONTROL_OTHER

logger = logging.getLogger(__name__)


@dataclass
class ModulesResults:
    effect: pd.DataFrame                    # targets x genes log2FC
    de_mask: pd.DataFrame                   # targets x genes bool
    gene_programs: pd.DataFrame             # gene, program, program_size, ...
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
    control: str = CONTROL_NTC
    program_activity_by_cluster: pd.DataFrame = field(default_factory=pd.DataFrame)   # programs x clusters
    tf_edges: pd.DataFrame = field(default_factory=pd.DataFrame)
    hubs: pd.DataFrame = field(default_factory=pd.DataFrame)
    module_connectivity: pd.DataFrame = field(default_factory=pd.DataFrame)
    program_genes: Dict[str, List[str]] = field(default_factory=dict)
    module_members: Dict[str, List[str]] = field(default_factory=dict)

    @property
    def n_programs(self) -> int:
        return len(self.program_labels)

    @property
    def n_modules(self) -> int:
        return len(self.module_labels)


def cluster_axis(items: pd.DataFrame, correlation: str, linkage_method: str, k: Optional[int], threshold: float, prefix: str) -> Tuple[pd.Series, List[str], List[str]]:
    """Hierarchical clustering of the rows of ``items`` on 1 - correlation; labels in leaf order (reference ``cluster_axis``)."""
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


def select_marker_genes(adata: ad.AnnData, key: str, n_per_cluster: int, method: str) -> List[str]:
    """Union of the top positive markers of each cluster (reference ``select_genes`` / ``cluster_markers``)."""
    import scanpy as sc

    work = ad.AnnData(X=adata.X, obs=adata.obs[[key]].copy(), var=pd.DataFrame(index=adata.var_names))
    work.obs[key] = work.obs[key].astype(str).astype("category")
    sc.tl.rank_genes_groups(work, key, method=method, n_genes=n_per_cluster)
    df = sc.get.rank_genes_groups_df(work, group=None)
    df = df[df["logfoldchanges"] > 0]
    genes: List[str] = []
    for _, sub in df.groupby("group", observed=True):
        genes.extend(sub.sort_values("logfoldchanges", ascending=False).head(n_per_cluster)["names"].tolist())
    seen: Dict[str, None] = {}
    var = set(adata.var_names)
    for g in genes:
        if g in var:
            seen.setdefault(g, None)
    return list(seen)


def tf_network(effect: pd.DataFrame, de_mask: pd.DataFrame, perturbation_module: pd.Series, n_cells: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """TF->TF edges, hub sizes and module-module connectivity (reference ``tf_network``)."""
    hub_counts = de_mask.sum(axis=1)
    hubs = pd.DataFrame({"target": hub_counts.index, "module": perturbation_module.reindex(hub_counts.index).to_numpy(), "n_cells": [int(n_cells.get(t, 0)) for t in hub_counts.index], "n_de_genes": hub_counts.to_numpy().astype(int)}).sort_values("n_de_genes", ascending=False, ignore_index=True)
    tf_genes = [g for g in effect.columns if g in effect.index]
    edges = []
    for src in effect.index:
        for tgt in tf_genes:
            if src == tgt:
                continue
            if bool(de_mask.at[src, tgt]):
                val = effect.at[src, tgt]
                edges.append({"source": src, "target": tgt, "log2fc": float(val), "sign": "positive" if val > 0 else "negative", "source_module": perturbation_module.get(src, ""), "target_module": perturbation_module.get(tgt, "")})
    tf_edges = pd.DataFrame(edges, columns=["source", "target", "log2fc", "sign", "source_module", "target_module"])
    labels = sorted(perturbation_module.unique(), key=lambda s: (len(s), s))
    sizes = perturbation_module.value_counts()
    conn = pd.DataFrame(0.0, index=labels, columns=labels)
    if not tf_edges.empty:
        for _, e in tf_edges.iterrows():
            a, b = e["source_module"], e["target_module"]
            if a in labels and b in labels:
                conn.loc[a, b] += 1.0
    for a in labels:
        for b in labels:
            denom = float(sizes.get(a, 1) * sizes.get(b, 1))
            conn.loc[a, b] = conn.loc[a, b] / denom if denom else 0.0
    conn.index.name = "module"
    return tf_edges, hubs, conn


def compute_modules(adata: ad.AnnData, cfg: Config) -> ModulesResults:
    pe = cfg.analysis.perturbation_effects
    mc = pe.modules
    groups = cell_groups(adata, pe.control_classes, mc.min_cells_per_perturbation, 0)
    skipped = list(groups.skipped)
    empty = ModulesResults(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None, [], [], [], [], pd.DataFrame(skipped), {}, True)
    if len(groups.targets) < mc.min_perturbations:
        empty.info = {"status": f"skipped: {len(groups.targets)} perturbations with >= {mc.min_cells_per_perturbation} cells (need {mc.min_perturbations})"}
        logger.info("modules: %s", empty.info["status"])
        return empty
    control = mc.control
    if control == CONTROL_NTC and groups.n_control == 0:
        logger.warning("modules: no control cells of %s; using 'other' (other-target cells) as the control instead of 'ntc'", groups.control_classes)
        control = CONTROL_OTHER
    X = adata.X
    genes_all = adata.var_names
    klass = adata.obs["perturbation_class"].astype(str).to_numpy()
    targeting = klass == CLASS_TARGETING
    # per-target Welch t-test on the log values, all genes ------------------------------------
    if control == CONTROL_NTC:
        mean_c, var_c, n_c = group_mean_var(X, groups.control)
        ctrl_stats = {t: (mean_c, var_c, n_c) for t in groups.targets}
    else:
        tot_mean, tot_var, n_tot = group_mean_var(X, targeting)
        tot_sum = tot_mean * n_tot
        tot_sq = tot_var * max(n_tot - 1, 1) + n_tot * tot_mean ** 2
        ctrl_stats = {}
        for t in groups.targets:
            m, v, npert = group_mean_var(X, groups.mask(t))
            n_o = max(n_tot - npert, 1)
            s_o = tot_sum - m * npert
            sq_o = tot_sq - (v * max(npert - 1, 1) + npert * m ** 2)
            mean_o = s_o / n_o
            var_o = np.clip((sq_o - s_o * mean_o) / max(n_o - 1, 1), 0, None)
            ctrl_stats[t] = (mean_o, var_o, n_o)
    t_rows, p_rows = [], []
    for t in groups.targets:
        m, v, npert = group_mean_var(X, groups.mask(t))
        mc_, vc_, nc_ = ctrl_stats[t]
        tt, pp = welch_t(m, v, npert, mc_, vc_, nc_)
        t_rows.append(tt); p_rows.append(pp)
    T = np.vstack(t_rows); P = np.vstack(p_rows)
    # all-gene BH (within perturbation) is only needed to choose the 'response' panel
    FDR = np.vstack([bh_fdr(P[i]) for i in range(P.shape[0])]) if mc.gene_selection == "response" else None
    # gene panel ---------------------------------------------------------------------------------
    panel: List[str] = []
    selection = mc.gene_selection
    if selection == "cluster_markers":
        key = mc.cluster_key
        if key in adata.obs.columns and adata.obs[key].astype(str).nunique() >= 2:
            panel = select_marker_genes(adata, key, mc.n_marker_genes_per_cluster, mc.marker_method)
        else:
            logger.warning("modules: obs[%r] has < 2 groups or is absent; falling back to highly variable genes for the gene panel", key)
            selection = "hvg (fallback: no clusters)"
    elif selection == "response":
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
    n_dropped_expr = 0
    if mc.min_pct_cells_expressing > 0:
        used = groups.control.copy()
        for t in groups.targets:
            used |= groups.mask(t)
        Xu = X[used]
        det = np.asarray((Xu > 0).mean(axis=0)).ravel() * 100 if hasattr(Xu, "tocsr") else (np.asarray(Xu) > 0).mean(axis=0) * 100
        keep = set(pd.Index(genes_all)[det >= mc.min_pct_cells_expressing])
        n_before = len(panel)
        panel = [g for g in panel if g in keep]
        n_dropped_expr = n_before - len(panel)
    if len(panel) < mc.min_genes:
        empty.info = {"status": f"skipped: only {len(panel)} panel genes after the expression filter (need {mc.min_genes})"}
        return empty
    gidx = genes_all.get_indexer(panel)
    eps = mc.log2fc_pseudocount
    # effect matrix on de-logged means -----------------------------------------------------------
    lin_all = None
    if control == CONTROL_NTC:
        lin_c = np.asarray(np.expm1(dense(X, rows=groups.control, cols=gidx)).mean(axis=0)).ravel()
    else:
        lin_all = np.expm1(dense(X, rows=targeting, cols=gidx))
        tot_lin = lin_all.sum(axis=0)
        n_tot_t = int(targeting.sum())
    eff = np.zeros((len(groups.targets), len(panel)))
    for i, t in enumerate(groups.targets):
        lin_p_mat = np.expm1(dense(X, rows=groups.mask(t), cols=gidx))
        lin_p = lin_p_mat.mean(axis=0)
        if control == CONTROL_NTC:
            ctrl_mean = lin_c
        else:
            n_p = lin_p_mat.shape[0]
            ctrl_mean = (tot_lin - lin_p_mat.sum(axis=0)) / max(n_tot_t - n_p, 1)
        eff[i] = np.log2((lin_p + eps) / (ctrl_mean + eps))
    effect = pd.DataFrame(eff, index=groups.targets, columns=panel)
    effect.index.name = "target"
    # reference: the Welch t-test p-values are BH-corrected within each perturbation over the PANEL genes
    FDR_panel = np.vstack([bh_fdr(P[i, gidx]) for i in range(P.shape[0])])
    de = pd.DataFrame((np.abs(eff) > mc.de_lfc_threshold) & (FDR_panel < mc.de_fdr_alpha), index=groups.targets, columns=panel)
    sign = np.sign(eff)
    # clustering -------------------------------------------------------------------------------
    gene_program, gene_order, program_labels = cluster_axis(effect.T, mc.program_correlation, mc.linkage_method, mc.n_programs, mc.cluster_distance_threshold, "P")
    pert_module, target_order, module_labels = cluster_axis(effect, mc.module_correlation, mc.linkage_method, mc.n_modules, mc.cluster_distance_threshold, "M")
    sizes = gene_program.value_counts()
    gene_programs = pd.DataFrame({"gene": gene_program.index, "program": gene_program.to_numpy(), "program_size": [int(sizes[p]) for p in gene_program.to_numpy()],
                                  "mean_abs_log2fc": np.abs(eff).mean(axis=0), "mean_log2fc": eff.mean(axis=0), "n_targets_de": de.sum(axis=0).to_numpy(),
                                  "n_targets_up": ((sign > 0) & de.to_numpy()).sum(axis=0), "n_targets_down": ((sign < 0) & de.to_numpy()).sum(axis=0)})
    gene_programs = gene_programs.sort_values(["program", "mean_abs_log2fc"], ascending=[True, False], ignore_index=True)
    n_cells_all = pd.Series(adata.obs["target"].astype(str).to_numpy()).value_counts()
    modules_tbl = pd.DataFrame({"target": pert_module.index, "module": pert_module.to_numpy(), "n_cells": [int(groups.counts[t]) for t in pert_module.index],
                                "n_de_genes": de.sum(axis=1).to_numpy(), "n_up": ((sign > 0) & de.to_numpy()).sum(axis=1), "n_down": ((sign < 0) & de.to_numpy()).sum(axis=1),
                                "mean_abs_log2fc": np.abs(eff).mean(axis=1)}).sort_values(["module", "n_de_genes"], ascending=[True, False], ignore_index=True)
    mp = pd.DataFrame(index=module_labels, columns=program_labels, dtype=float)
    mp.index.name = "module"
    ppe_rows = []
    program_genes = {p: gene_program.index[gene_program == p].tolist() for p in program_labels}
    for p in program_labels:
        pg = program_genes[p]
        cols = effect.columns.get_indexer(pg)
        for i, t in enumerate(groups.targets):
            ppe_rows.append({"target": t, "program": p, "n_genes": int(len(pg)), "mean_log2fc": float(eff[i, cols].mean()), "frac_de": float(de.to_numpy()[i, cols].mean()), "n_de": int(de.to_numpy()[i, cols].sum())})
        for m in module_labels:
            perts = pert_module.index[pert_module == m]
            block = effect.loc[perts, pg].to_numpy()
            mp.loc[m, p] = float(np.nanmean(block)) if block.size else np.nan
    ppe = pd.DataFrame(ppe_rows)
    module_members = {m: pert_module.index[pert_module == m].tolist() for m in module_labels}
    tf_edges, hubs, conn = tf_network(effect, de, pert_module, n_cells_all)
    # per-cell program activity -------------------------------------------------------------------
    activity = None
    pa_cluster = pd.DataFrame()
    if mc.score_programs:
        if mc.program_scoring == "score_genes":
            import scanpy as sc

            work = ad.AnnData(X=X, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names))
            cols_ = {}
            for p in program_labels:
                sc.tl.score_genes(work, program_genes[p], score_name=f"program_{p}_score", ctrl_size=50, random_state=cfg.compute.seed)
                cols_[p] = work.obs[f"program_{p}_score"].to_numpy(dtype=np.float32)
            activity = pd.DataFrame(cols_, index=adata.obs_names)
        else:
            Z = dense(X, cols=gidx, dtype=np.float32)
            mu, sd = Z.mean(axis=0), Z.std(axis=0)
            sd[sd == 0] = 1.0
            Z = (Z - mu) / sd
            activity = pd.DataFrame({p: Z[:, effect.columns.get_indexer(program_genes[p])].mean(axis=1) for p in program_labels}, index=adata.obs_names).astype(np.float32)
        key = mc.cluster_key
        if key in adata.obs.columns:
            pa_cluster = activity.groupby(adata.obs[key].astype(str).to_numpy()).mean().T
            try:
                pa_cluster = pa_cluster[sorted(pa_cluster.columns, key=lambda c: (float(c), c))]
            except ValueError:
                pa_cluster = pa_cluster[sorted(pa_cluster.columns)]
            pa_cluster.index.name = "program"
    info = {"method": "perturbation x gene log2FC of de-logged means vs control; programs = genes clustered on 1-Pearson, modules = perturbations on 1-Spearman (Zhou et al. 2023 design as in weili-lab/perturbseq-pipeline)",
            "gene_selection": selection, "cluster_key": mc.cluster_key, "n_marker_genes_per_cluster": mc.n_marker_genes_per_cluster, "marker_method": mc.marker_method, "n_panel_genes": len(panel), "n_panel_genes_dropped_low_expression": int(n_dropped_expr),
            "min_pct_cells_expressing": mc.min_pct_cells_expressing, "log2fc_pseudocount": eps, "n_response_genes_per_perturbation": mc.n_response_genes_per_perturbation, "response_fdr": mc.response_fdr,
            "control": control, "control_label": CONTROL_LABELS[control], "n_perturbations": len(groups.targets), "n_programs": len(program_labels), "n_modules": len(module_labels), "n_programs_requested": mc.n_programs, "n_modules_requested": mc.n_modules,
            "cluster_distance_threshold": mc.cluster_distance_threshold, "linkage_method": mc.linkage_method, "program_correlation": mc.program_correlation, "module_correlation": mc.module_correlation,
            "de_lfc_threshold": mc.de_lfc_threshold, "de_fdr_alpha": mc.de_fdr_alpha, "control_classes": groups.control_classes, "n_control_cells": groups.n_control, "median_de_genes_per_perturbation": float(de.sum(axis=1).median()),
            "n_tf_edges": int(len(tf_edges)), "top_hub": str(hubs.iloc[0]["target"]) if len(hubs) else "", "top_hub_n_de": int(hubs.iloc[0]["n_de_genes"]) if len(hubs) else 0,
            "program_activity": ("sc.tl.score_genes (ctrl_size 50)" if mc.program_scoring == "score_genes" else "mean z-scored log expression of program genes") if activity is not None else "not computed", "status": "computed"}
    logger.info("modules: %d perturbations x %d genes (log2FC vs %s; panel: %s) -> %d programs, %d modules; median %d DE genes/perturbation", len(groups.targets), len(panel), control, selection, len(program_labels), len(module_labels), int(de.sum(axis=1).median()))
    return ModulesResults(effect, de, gene_programs, modules_tbl, mp, ppe, activity, program_labels, module_labels, gene_order, target_order, pd.DataFrame(skipped), info, False, control, pa_cluster, tf_edges, hubs, conn, program_genes, module_members)
