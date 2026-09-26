"""Perturbation-effect figures (section ``perturbation_effects``), through the FigureRegistry.

Reference figure classes (weili-lab/perturbseq-pipeline ``plots.py``) for the
modules / programs, PS and lochNESS stages, followed by the protein-extension
figures (protein effects, RNA-protein concordance). Per-target figures beyond
``top_n_report`` are written to disk but not embedded in the report.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..analysis._common import CLASS_CONTROL, CLASS_TARGETING  # noqa: E402
from ..analysis.perturbation_strength import CONTROL_LABELS  # noqa: E402
from ..analysis.ps_score import QUADRANT_COLORS, QUADRANT_ESCAPER, QUADRANT_KD, QUADRANT_LABELS, QUADRANT_LOW, QUADRANT_NONRESP  # noqa: E402
from ._refplot import GREY, HIT_RED, block_spans, despine, display_order, gene_values, label_clusters, label_colors, order_by_similarity, scatter_umap, umap_coords  # noqa: E402
from .figures import FigureRegistry  # noqa: E402

logger = logging.getLogger("petrubseq_protein")
SECTION = "perturbation_effects"
ST_PS, ST_PS_TARGET, ST_LDA, ST_LDA_TARGET = "ps", "ps_per_target", "ps_lda", "ps_lda_per_target"
ST_LOCH, ST_LOCH_TARGET = "lochness", "lochness_per_target"
ST_PROG, ST_PROG_UMAP = "gene_programs", "gene_programs_umap"
ST_PROT, ST_CONC = "protein_effects", "concordance"
BLUE, RED = "#4C72B0", "#C44E52"


def _trend(ax, x: np.ndarray, y: np.ndarray, color: str = RED) -> None:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 10:
        return
    edges = np.unique(np.quantile(x, np.linspace(0, 1, 9)))
    if len(edges) < 3:
        return
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    mx = [np.median(x[idx == i]) for i in range(len(edges) - 1) if (idx == i).any()]
    my = [np.median(y[idx == i]) for i in range(len(edges) - 1) if (idx == i).any()]
    ax.plot(mx, my, color=color, lw=1.5)


def _heatmap(ax, M: pd.DataFrame, cmap: str = "RdBu_r", center: bool = True, label: str = "", marks: Optional[pd.DataFrame] = None):
    v = M.to_numpy(float)
    lim = np.nanmax(np.abs(v)) if np.isfinite(v).any() else 1.0
    im = ax.imshow(v, aspect="auto", cmap=cmap, vmin=-lim if center else None, vmax=lim if center else None, interpolation="nearest")
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=6 if M.shape[0] > 30 else 7)
    if M.shape[1] <= 40:
        ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, fontsize=7, rotation=90)
    else:
        ax.set_xticks([])
    if marks is not None:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if bool(marks.iat[i, j]):
                    ax.text(j, i, "*", ha="center", va="center", fontsize=8)
    cb = plt.colorbar(im, ax=ax, fraction=0.04)
    cb.set_label(label, fontsize=8)


# ---------------------------------------------------------------------------------------- modules
def _plot_effect_heatmap(mo, reg, cfg) -> None:
    prog = mo.gene_programs.set_index("gene")["program"]
    mod = mo.perturbation_modules.set_index("target")["module"]
    prog_rank = {l: i for i, l in enumerate(mo.program_labels)}
    mod_rank = {l: i for i, l in enumerate(mo.module_labels)}
    g_rank = {g: i for i, g in enumerate(mo.gene_order)}
    p_rank = {p: i for i, p in enumerate(mo.target_order)}
    genes = display_order(mo.gene_order, prog, prog_rank, g_rank)
    perts = display_order(mo.target_order, mod, mod_rank, p_rank)
    D = mo.effect.loc[perts, genes].T
    row_lab = [prog[g] for g in genes]; col_lab = [mod[t] for t in perts]
    prog_colors = label_colors(mo.program_labels); mod_colors = label_colors(mo.module_labels)
    lim = float(np.nanpercentile(np.abs(D.to_numpy()), 98)) or 1.0
    w = max(7.0, 0.06 * len(perts) + 3); h = max(6.0, 0.02 * len(genes) + 3)
    fig = plt.figure(figsize=(w, h))
    gs = fig.add_gridspec(2, 3, width_ratios=[0.02, 1, 0.04], height_ratios=[0.03, 1], wspace=0.02, hspace=0.02)
    ax_top = fig.add_subplot(gs[0, 1]); ax_left = fig.add_subplot(gs[1, 0]); ax = fig.add_subplot(gs[1, 1]); cax = fig.add_subplot(gs[1, 2])
    im = ax.imshow(D.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_yticks([])
    if len(perts) <= 60:
        ax.set_xticks(range(len(perts))); ax.set_xticklabels(perts, rotation=90, fontsize=5)
    else:
        ax.set_xticks([])
    ax.set_xlabel(f"{len(perts)} perturbations (grouped into {mo.n_modules} modules)")
    ax.set_ylabel(f"{len(genes)} genes (grouped into {mo.n_programs} programs)")
    ax_top.imshow(np.array([[mod_colors[l][:3] for l in col_lab]]), aspect="auto"); ax_top.set_xticks([]); ax_top.set_yticks([])
    for lab, s, e in block_spans(col_lab):
        ax_top.text((s + e - 1) / 2, 0, lab, ha="center", va="center", fontsize=7, fontweight="bold")
    ax_left.imshow(np.array([[prog_colors[l][:3]] for l in row_lab]), aspect="auto"); ax_left.set_xticks([]); ax_left.set_yticks([])
    for lab, s, e in block_spans(row_lab):
        ax_left.text(0, (s + e - 1) / 2, lab, ha="center", va="center", rotation=90, fontsize=7, fontweight="bold")
    for _, _, e in block_spans(col_lab)[:-1]:
        ax.axvline(e - 0.5, color="white", lw=0.6)
    for _, _, e in block_spans(row_lab)[:-1]:
        ax.axhline(e - 0.5, color="white", lw=0.6)
    plt.colorbar(im, cax=cax, label="log2FC vs control")
    ax_top.set_title(f"Regulome map: {mo.n_modules} co-functional modules x {mo.n_programs} gene programs", fontsize=11, pad=14)
    reg.save(fig, "regulome_heatmap", SECTION, ST_PROG, "Co-functional modules and gene programs",
             f"Perturbation x gene log2FC vs {CONTROL_LABELS.get(mo.control, mo.control)}. Columns are {len(perts)} perturbations grouped into {mo.n_modules} co-functional modules ({mo.info.get('module_correlation')}-correlation clustering); rows are {len(genes)} downstream genes grouped into {mo.n_programs} programs ({mo.info.get('program_correlation')}-correlation clustering). Red = up, blue = down after perturbation. Module (M) and program (P) labels are arbitrary cluster ids.")


def _plot_module_program(mo, reg) -> None:
    mp = mo.module_program.astype(float)
    if mp.empty:
        return
    lim = float(np.nanmax(np.abs(mp.to_numpy()))) or 1.0
    fig, ax = plt.subplots(figsize=(max(4, 0.7 * mp.shape[1] + 2), max(3, 0.5 * mp.shape[0] + 1.5)))
    im = ax.imshow(mp.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(mp.shape[1])); ax.set_xticklabels(mp.columns); ax.set_yticks(range(mp.shape[0])); ax.set_yticklabels(mp.index)
    ax.set_xlabel("Gene program"); ax.set_ylabel("Co-functional module")
    for i in range(mp.shape[0]):
        for j in range(mp.shape[1]):
            v = mp.to_numpy()[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7, color="white" if abs(v) > 0.6 * lim else "black")
    plt.colorbar(im, ax=ax, shrink=0.7, label="mean log2FC")
    ax.set_title("Module -> program regulatory strength", fontsize=11)
    fig.tight_layout()
    reg.save(fig, "module_program_strength", SECTION, ST_PROG, "Module x program strength", "Signed strength of each module's regulation of each program: the mean program-gene log2FC per perturbation, averaged over the module's perturbations. Red = net activation, blue = net repression.")


def _plot_alluvial(mo, reg) -> None:
    mp = mo.module_program.astype(float)
    if mp.empty:
        return
    mod_colors = label_colors(mo.module_labels)
    mag = mp.abs().fillna(0.0)
    mod_tot = mag.sum(axis=1); prog_tot = mag.sum(axis=0)
    total = float(mag.to_numpy().sum())
    if total <= 0:
        return
    gap = 0.02

    def _stack(totals):
        pos, y = {}, 1.0
        n = len(totals)
        usable = 1.0 - gap * (n - 1)
        for name, t in totals.items():
            hgt = usable * (t / total) if total else 0.0
            pos[name] = (y - hgt, y)
            y -= hgt + gap
        return pos

    left = _stack(mod_tot); right = _stack(prog_tot)
    fig, ax = plt.subplots(figsize=(7, max(4, 0.5 * mp.shape[0] + 2)))
    left_cursor = {m: left[m][1] for m in mp.index}; right_cursor = {p: right[p][1] for p in mp.columns}
    xs = np.linspace(0, 1, 40); sm = xs * xs * (3 - 2 * xs)
    for m in mp.index:
        for p in mp.columns:
            v = mp.loc[m, p]
            if not np.isfinite(v) or v == 0:
                continue
            thick = (1.0 - gap * (max(len(mp.index), len(mp.columns)) - 1)) * (abs(v) / total)
            l_hi = left_cursor[m]; l_lo = l_hi - thick; left_cursor[m] = l_lo
            r_hi = right_cursor[p]; r_lo = r_hi - thick; right_cursor[p] = r_lo
            X = 0.08 + 0.84 * xs
            ax.fill_between(X, l_lo + (r_lo - l_lo) * sm, l_hi + (r_hi - l_hi) * sm, color=("#c0392b" if v > 0 else "#2c6fbb"), alpha=0.55, lw=0)
    for m, (lo, hi) in left.items():
        ax.add_patch(plt.Rectangle((0.04, lo), 0.04, hi - lo, color=mod_colors[m]))
        ax.text(0.02, (lo + hi) / 2, m, ha="right", va="center", fontsize=8, fontweight="bold")
    for p, (lo, hi) in right.items():
        ax.add_patch(plt.Rectangle((0.92, lo), 0.04, hi - lo, color="#666666"))
        ax.text(0.98, (lo + hi) / 2, p, ha="left", va="center", fontsize=8, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.axis("off")
    ax.set_title("Module -> program regulation (alluvial)\nred = activation, blue = repression", fontsize=11)
    reg.save(fig, "module_program_alluvial", SECTION, ST_PROG, "Module -> program alluvial", "Regulatory flow from co-functional modules (left) to gene programs (right). Ribbon width is the magnitude of the mean effect; red = activation, blue = repression. Bar heights show each module/program's total regulatory strength.")


def _plot_module_correlation(mo, reg) -> None:
    perts = mo.target_order
    if len(perts) < 3:
        return
    mod = mo.perturbation_modules.set_index("target")["module"]
    mod_rank = {l: i for i, l in enumerate(mo.module_labels)}
    p_rank = {p: i for i, p in enumerate(perts)}
    order = display_order(perts, mod, mod_rank, p_rank)
    corr = mo.effect.loc[order].T.corr(method=mo.info.get("module_correlation", "spearman"))
    fig, ax = plt.subplots(figsize=(max(5, 0.14 * len(order) + 2),) * 2)
    im = ax.imshow(corr.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1)
    show = len(order) <= 70
    ax.set_xticks(range(len(order)) if show else [])
    if show:
        ax.set_xticklabels(order, rotation=90, fontsize=4); ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=4)
    else:
        ax.set_yticks([])
    for _, _, e in block_spans([mod[t] for t in order])[:-1]:
        ax.axvline(e - 0.5, color="black", lw=0.5); ax.axhline(e - 0.5, color="black", lw=0.5)
    plt.colorbar(im, ax=ax, shrink=0.6, label=f"{mo.info.get('module_correlation', 'spearman')} r")
    ax.set_title("Perturbation similarity (co-functional modules)", fontsize=11)
    fig.tight_layout()
    reg.save(fig, "module_correlation", SECTION, ST_PROG, "Perturbation correlation (modules)", f"{str(mo.info.get('module_correlation', 'spearman')).title()} correlation between perturbations of their effect profiles over the gene panel. Black lines separate the co-functional modules; red blocks are perturbations acting together.")


def _plot_program_activity(mo, reg, cfg) -> None:
    act = mo.program_activity_by_cluster
    if act is None or act.empty:
        return
    lim = float(np.nanmax(np.abs(act.to_numpy(dtype=float)))) or 1.0
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * act.shape[1] + 2), max(2.5, 0.5 * act.shape[0] + 1.5)))
    im = ax.imshow(act.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(act.shape[1])); ax.set_xticklabels(act.columns, fontsize=7); ax.set_yticks(range(act.shape[0])); ax.set_yticklabels(act.index)
    ax.set_xlabel(f"Cluster ({cfg.analysis.perturbation_effects.modules.cluster_key})"); ax.set_ylabel("Gene program")
    plt.colorbar(im, ax=ax, shrink=0.7, label="mean program score")
    ax.set_title("Gene-program activity by cluster", fontsize=11)
    fig.tight_layout()
    reg.save(fig, "program_activity_by_cluster", SECTION, ST_PROG, "Program activity by cluster", f"Mean per-cell program score ({mo.info.get('program_activity')}) in each cell-state cluster - which transcriptional states each program marks.")


def _plot_program_umaps(mo, adata, reg, cfg) -> None:
    coords = umap_coords(adata)
    if coords is None or mo.program_activity is None:
        return
    top_n = cfg.analysis.perturbation_effects.modules.top_n_report
    for rank, label in enumerate(mo.program_labels):
        fig, ax = plt.subplots(figsize=(5, 4.2))
        scatter_umap(ax, coords, mo.program_activity[label].to_numpy(dtype=float), False, f"Program {label} activity", size=4, cmap="RdBu_r")
        fig.tight_layout()
        reg.save(fig, f"program_{label}_umap", SECTION, ST_PROG_UMAP, f"Program {label} activity (UMAP)", f"Per-cell score for program {label} ({len(mo.program_genes.get(label, []))} genes) on the embedding.", in_report=rank < top_n)


def _plot_networks(mo, reg, cfg) -> None:
    conn = mo.module_connectivity
    if conn.shape[0] >= 2:
        fig, ax = plt.subplots(figsize=(max(4, 0.5 * conn.shape[0] + 2),) * 2)
        im = ax.imshow(conn.to_numpy(dtype=float), cmap="magma", aspect="auto")
        ax.set_xticks(range(conn.shape[1])); ax.set_xticklabels(conn.columns, fontsize=7); ax.set_yticks(range(conn.shape[0])); ax.set_yticklabels(conn.index, fontsize=7)
        plt.colorbar(im, ax=ax, shrink=0.7, label="TF-TF edges / (size_i x size_j)")
        ax.set_title("Module-module connectivity", fontsize=11)
        fig.tight_layout()
        reg.save(fig, "module_connectivity", SECTION, ST_PROG, "Module-module connectivity", "Regulatory connectivity between modules: TF->TF edges spanning two modules, normalised by the product of their sizes.")
    if not cfg.analysis.perturbation_effects.modules.draw_networks:
        return
    try:
        import networkx as nx
    except ImportError:  # pragma: no cover - optional
        logger.warning("modules.draw_networks is on but networkx is not installed; skipping the network graphs")
        return
    mod_colors = label_colors(mo.module_labels)
    if conn.shape[0] >= 2:
        G = nx.Graph()
        sizes = mo.perturbation_modules["module"].value_counts()
        for m in conn.index:
            G.add_node(m, size=int(sizes.get(m, 1)))
        for i, a in enumerate(conn.index):
            for b in conn.columns[i + 1:]:
                w = conn.loc[a, b] + conn.loc[b, a]
                if w > 0:
                    G.add_edge(a, b, weight=w)
        pos = nx.spring_layout(G, seed=0, weight="weight")
        fig, ax = plt.subplots(figsize=(6, 5))
        nx.draw_networkx_edges(G, pos, ax=ax, width=[2 + 8 * G[u][v]["weight"] for u, v in G.edges], edge_color="#b0b0b0")
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[80 + 40 * G.nodes[n]["size"] for n in G.nodes], node_color=[mod_colors[n] for n in G.nodes])
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=9, font_weight="bold")
        ax.axis("off"); ax.set_title("Module interaction network", fontsize=11)
        reg.save(fig, "module_network", SECTION, ST_PROG, "Module interaction network", "Modules (nodes, sized by number of member perturbations) linked by their TF-TF regulatory connectivity (edge width).")
    edges, hubs = mo.tf_edges, mo.hubs
    if not edges.empty and not hubs.empty:
        top_hubs = set(hubs.head(40)["target"])
        H = nx.DiGraph()
        mod_of = mo.perturbation_modules.set_index("target")["module"].to_dict()
        hub_size = hubs.set_index("target")["n_de_genes"].to_dict()
        for _, e in edges.iterrows():
            if e["source"] in top_hubs and e["target"] in top_hubs:
                H.add_edge(e["source"], e["target"], sign=e["sign"])
        if H.number_of_nodes() >= 2:
            pos = nx.spring_layout(H, seed=0)
            fig, ax = plt.subplots(figsize=(7, 6))
            nx.draw_networkx_edges(H, pos, ax=ax, edge_color=["#c0392b" if H[u][v]["sign"] == "positive" else "#2c6fbb" for u, v in H.edges], alpha=0.5, arrowsize=8, width=0.8)
            nx.draw_networkx_nodes(H, pos, ax=ax, node_size=[40 + 12 * hub_size.get(n, 1) for n in H.nodes], node_color=[mod_colors.get(mod_of.get(n, ""), "#999999") for n in H.nodes])
            nx.draw_networkx_labels(H, pos, ax=ax, font_size=6)
            ax.axis("off"); ax.set_title("TF regulatory network (hubs)\nred = activation, blue = repression", fontsize=11)
            reg.save(fig, "tf_hub_network", SECTION, ST_PROG, "TF hub network", "Regulatory edges between the top hub perturbations (node size = number of DE genes it perturbs; colour = its module). Red = activating, blue = repressing edge.")


def modules_figures(mo, adata, cfg, reg: FigureRegistry) -> None:
    if mo is None or mo.empty:
        return
    _plot_effect_heatmap(mo, reg, cfg)
    _plot_module_program(mo, reg)
    _plot_alluvial(mo, reg)
    _plot_module_correlation(mo, reg)
    _plot_program_activity(mo, reg, cfg)
    _plot_program_umaps(mo, adata, reg, cfg)
    _plot_networks(mo, reg, cfg)
    # extension: perturbation x program heatmap with module labels
    ppe = mo.perturbation_program_effects.pivot(index="target", columns="program", values="mean_log2fc").reindex(index=mo.target_order, columns=mo.program_labels)
    mod = mo.perturbation_modules.set_index("target")["module"]
    ppe.index = [f"{t} [{mod[t]}]" for t in ppe.index]
    fig, ax = plt.subplots(figsize=(1.2 + 0.6 * ppe.shape[1], 1.2 + 0.2 * ppe.shape[0]))
    _heatmap(ax, ppe, label="mean log2FC of program genes")
    ax.set_title("perturbation x gene program", fontsize=9)
    fig.tight_layout()
    reg.save(fig, "perturbation_program_heatmap", SECTION, ST_PROG, "Perturbation x gene-program effects", "Mean log2FC of each gene program's genes under each perturbation; rows ordered by the perturbation-module dendrogram, module in brackets (extension of the reference set).")
    logger.info("Wrote module/program figures (%d modules, %d programs)", mo.n_modules, mo.n_programs)


# --------------------------------------------------------------------------------------------- PS
def ps_figures(ps, ps_vs_strength: Optional[pd.DataFrame], strength, adata, cfg, reg: FigureRegistry) -> None:
    if ps is None or ps.summary.empty:
        return
    summary = ps.summary
    top_n = cfg.analysis.perturbation_effects.top_n_report
    # --- 1. outcome by target ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(max(6, 0.20 * len(summary)), 4.0))
    bottom = np.zeros(len(summary))
    for key, q in (("pct_successful_kd", QUADRANT_KD), ("pct_escaper", QUADRANT_ESCAPER), ("pct_non_responder", QUADRANT_NONRESP), ("pct_low_signal", QUADRANT_LOW)):
        vals = summary[key].fillna(0).to_numpy(dtype=float)
        ax.bar(range(len(summary)), vals, bottom=bottom, color=QUADRANT_COLORS[q], label=QUADRANT_LABELS[q])
        bottom += vals
    ax.set_xticks(range(len(summary))); ax.set_xticklabels(summary["target"], rotation=90, fontsize=6)
    ax.set_ylabel("% of perturbed cells"); ax.set_title("Per-cell perturbation outcome by target", fontsize=11)
    ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "ps_outcome_by_target", SECTION, ST_PS, "Per-cell perturbation outcome", "Each perturbed cell is classified by its perturbation score and the target's own expression. Green is a confirmed knockdown; red are escapers, which carry the guide and show the signature yet still express the gene.")
    # --- 2. escaper fraction -----------------------------------------------------------------
    esc = summary.sort_values("pct_escaper", ascending=False)
    fig, ax = plt.subplots(figsize=(max(6, 0.20 * len(esc)), 3.6))
    ax.bar(range(len(esc)), esc["pct_escaper"].fillna(0), color=HIT_RED)
    ax.set_xticks(range(len(esc))); ax.set_xticklabels(esc["target"], rotation=90, fontsize=6)
    ax.set_ylabel("% escapers"); ax.set_title("Escaper fraction per target", fontsize=11)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "ps_escaper_fraction", SECTION, ST_PS, "Escaper fraction per target", "Cells carrying a guide whose target is nonetheless still expressed. A high fraction means the population-level effect understates how well the guide works in the cells where it does work.")
    # --- 3. agreement with the group-level test ---------------------------------------------
    if ps_vs_strength is not None and not ps_vs_strength.empty and strength is not None:
        primary = strength.primary_control
        lfc_col = f"log2fc_{primary}"
        merged = ps_vs_strength
        if lfc_col in merged.columns:
            ok = merged[lfc_col].notna() & merged["pct_successful_kd"].notna()
            if ok.sum() > 2:
                x = merged.loc[ok, lfc_col].to_numpy(dtype=float); y = merged.loc[ok, "pct_successful_kd"].to_numpy(dtype=float)
                r = float(np.corrcoef(x, y)[0, 1])
                fig, ax = plt.subplots(figsize=(5.4, 4.6))
                hit = merged.loc[ok].get(f"is_hit_{primary}")
                colors = [HIT_RED if h else GREY for h in hit] if hit is not None else BLUE
                ax.scatter(x, y, s=28, c=colors)
                for xi, yi, name in zip(x, y, merged.loc[ok, "target"]):
                    if yi > np.percentile(y, 85) or xi < np.percentile(x, 15):
                        ax.annotate(name, (xi, yi), fontsize=6, xytext=(3, 3), textcoords="offset points")
                ax.set_xlabel("log2FC of the target's own expression (group-level test)"); ax.set_ylabel("% cells with confirmed knockdown (per-cell score)")
                ax.set_title(f"Per-cell score vs group-level knockdown (Pearson r = {r:.2f})", fontsize=10)
                despine(ax)
                fig.tight_layout()
                reg.save(fig, "ps_vs_perturbation_strength", SECTION, ST_PS, "Per-cell scores vs the group-level test", "The group-level test uses the target's own expression, while the per-cell score projects cells onto the perturbation's whole downstream signature; a gene can be strongly knocked down yet change little downstream, or the reverse, so these need not track each other closely.")
    # --- 4. quadrant scatter per target ------------------------------------------------------
    obs = adata.obs
    targets = obs["target"].astype(str); klass = obs["perturbation_class"].astype(str)
    rng = np.random.default_rng(cfg.compute.seed)
    thr = ps.ps_threshold
    for rank, t in enumerate(summary["target"]):
        gene = ps.genes.get(t, "")
        if t not in ps.scores.columns or not gene or gene not in adata.var_names:
            continue
        series = ps.scores[t].dropna()
        cells = series.index
        expression = pd.Series(gene_values(adata, gene), index=adata.obs_names).loc[cells]
        is_target = ((targets.loc[cells] == t) & (klass.loc[cells] == CLASS_TARGETING)).to_numpy()
        cut = ps.expression_cut.get(t, float(np.median(expression)))
        fig, ax = plt.subplots(figsize=(6.4, 5.0))
        ctrl_idx = cells[~is_target]
        if len(ctrl_idx) > 2000:
            ctrl_idx = pd.Index(rng.choice(ctrl_idx, size=2000, replace=False))
        ax.scatter(series.loc[ctrl_idx], expression.loc[ctrl_idx], s=14, c="#cbd5e0", alpha=0.45, linewidths=0, label="control cells", rasterized=True)
        tgt = cells[is_target]
        quad = ps.quadrants.get(t)
        colors = [QUADRANT_COLORS.get(str(quad.get(c, QUADRANT_LOW)), GREY) for c in tgt] if quad is not None else "#e53e3e"
        ax.scatter(series.loc[tgt], expression.loc[tgt], s=26, c=colors, alpha=0.85, edgecolors="white", linewidths=0.4, label=f"{t} cells", rasterized=True)
        ax.axvline(thr, color="black", ls="--", lw=1, alpha=0.6); ax.axhline(cut, color="black", ls="--", lw=1, alpha=0.6)
        row = summary[summary["target"] == t].iloc[0]
        xmax = float(max(series.max(), thr * 2)); ymax = float(max(expression.max(), cut * 2)) or 1.0
        ax.text(thr + (xmax - thr) * 0.5, ymax * 0.95, f"ESCAPERS\n{row['pct_escaper']:.0f}%", fontsize=8, ha="center", color=QUADRANT_COLORS[QUADRANT_ESCAPER], fontweight="bold")
        ax.text(thr + (xmax - thr) * 0.5, ymax * 0.05, f"KNOCKED DOWN\n{row['pct_successful_kd']:.0f}%", fontsize=8, ha="center", color=QUADRANT_COLORS[QUADRANT_KD], fontweight="bold")
        ax.text(thr * 0.5, ymax * 0.95, f"NON-RESPONDER\n{row['pct_non_responder']:.0f}%", fontsize=8, ha="center", color=QUADRANT_COLORS[QUADRANT_NONRESP], fontweight="bold")
        ax.text(thr * 0.5, ymax * 0.05, f"LOW SIGNAL\n{row['pct_low_signal']:.0f}%", fontsize=8, ha="center", color="#718096", fontweight="bold")
        ax.set_xlabel("Perturbation score (per cell)"); ax.set_ylabel(f"{gene} expression (log-normalized)")
        ax.set_title(f"{t}: per-cell perturbation outcome (n={int(row['n_perturbed_cells']):,})", fontsize=10)
        ax.legend(fontsize=7, frameon=False, loc="upper right")
        despine(ax)
        fig.tight_layout()
        reg.save(fig, f"ps_quadrant_{t}", SECTION, ST_PS_TARGET, f"{t} perturbation score vs expression",
                 f"Cells carrying {t} guides, split by perturbation score (vertical cut at {thr}) and by {gene} expression relative to the control {row.get('expression_cut_method', 'mean')} (horizontal cut). {row['pct_successful_kd']:.0f}% show a confirmed knockdown and {row['pct_escaper']:.0f}% escape it.",
                 in_report=rank < top_n)
    logger.info("Wrote %d perturbation-score quadrant figures (%d shown in the report)", len(summary), min(top_n, len(summary)))
    # --- 5. LDA embedding ---------------------------------------------------------------------
    if ps.lda_umap is None:
        return
    coords = np.asarray(ps.lda_umap, dtype=float)
    placed = np.isfinite(coords).all(axis=1)
    if placed.sum() < 10:
        return
    labels = ps.lda_label.astype(str).to_numpy() if ps.lda_label is not None else np.full(adata.n_obs, "?", dtype=object)
    n_targets = len(summary)
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    scatter_umap(ax, coords[placed], labels[placed], True, f"Supervised LDA embedding ({n_targets} targets + control)", size=4, legend=n_targets <= 24, xlabel="LDA-UMAP1", ylabel="LDA-UMAP2")
    fig.tight_layout()
    reg.save(fig, "ps_lda_overview", SECTION, ST_LDA, "Supervised LDA embedding", "Linear discriminant analysis trained on the perturbation labels, then embedded with UMAP. Unlike the unsupervised embedding, the axes here are chosen to separate perturbations, so groups that overlap there can resolve." + ("" if n_targets <= 24 else " Legend omitted (too many targets)."))
    own = ps.own.to_numpy(dtype=float)
    thr_l = cfg.analysis.perturbation_effects.ps.lda_highlight_threshold
    strong = placed & np.isfinite(own) & (own >= thr_l)
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    ax.scatter(coords[placed, 0], coords[placed, 1], s=4, color="#e2e8f0", linewidths=0, rasterized=True, label="all cells")
    sca = ax.scatter(coords[strong, 0], coords[strong, 1], s=12, c=own[strong], cmap="viridis", vmin=thr_l, vmax=1.0, linewidths=0.2, edgecolors="#2d3748", rasterized=True)
    plt.colorbar(sca, ax=ax, shrink=0.75, label="perturbation score")
    ax.set_title(f"High-confidence responders (score >= {thr_l}): {int(strong.sum()):,} of {int(placed.sum()):,} cells", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.legend(fontsize=7, frameon=False, loc="best", markerscale=3)
    despine(ax, left=True, bottom=True)
    fig.tight_layout()
    reg.save(fig, "ps_lda_high_confidence", SECTION, ST_LDA, "High-confidence responders on the LDA map", f"Cells whose perturbation score reaches {thr_l}, coloured by score. Where these concentrate is where the screen produced its clearest phenotypes.")
    tg = obs["target"].astype(str).to_numpy()
    for rank, t in enumerate(summary["target"]):
        if t not in ps.scores.columns:
            continue
        score = ps.scores[t].to_numpy(dtype=float)
        is_target = (tg == t) & placed
        if is_target.sum() == 0:
            continue
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        bg = placed & ~is_target
        ax.scatter(coords[bg, 0], coords[bg, 1], s=4, color="#e2e8f0", alpha=0.6, linewidths=0, rasterized=True, label="other cells")
        order = np.argsort(np.nan_to_num(score[is_target]))
        idx = np.where(is_target)[0][order]
        sca = ax.scatter(coords[idx, 0], coords[idx, 1], s=18, c=np.nan_to_num(score[idx]), cmap="Blues", vmin=0, vmax=1, linewidths=0.3, edgecolors="#2d3748", rasterized=True)
        plt.colorbar(sca, ax=ax, shrink=0.75, label="perturbation score")
        row = summary[summary["target"] == t].iloc[0]
        ax.set_title(f"{t} on the LDA map (n={int(row['n_perturbed_cells']):,}, {row['pct_successful_kd']:.0f}% knocked down)", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.legend(fontsize=7, frameon=False, loc="best", markerscale=3)
        despine(ax, left=True, bottom=True)
        fig.tight_layout()
        reg.save(fig, f"ps_lda_{t}", SECTION, ST_LDA_TARGET, f"{t} on the LDA embedding", f"Cells carrying {t} guides, shaded by perturbation score, against all other cells in grey. Darker cells respond more strongly; a tight darker cluster means the perturbation drives a consistent state.", in_report=rank < top_n)
    logger.info("Wrote %d per-target LDA figures (%d shown in the report)", n_targets, min(top_n, n_targets))


# --------------------------------------------------------------------------------------- lochNESS
def _lochness_norm(vmax: float):
    import matplotlib.colors as mcolors

    vmax = max(float(vmax), 1.0)
    return mcolors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=-vmax, vmax=vmax, base=10)


def lochness_figures(lo, adata, cfg, reg: FigureRegistry) -> None:
    if lo is None or lo.summary.empty:
        return
    lcfg = cfg.analysis.perturbation_effects.lochness
    summary = lo.summary
    coords = umap_coords(adata)
    key = cfg.analysis.clustering.key
    top_n = cfg.analysis.perturbation_effects.top_n_report
    # --- 1. self-enrichment ranking -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(max(6, 0.20 * len(summary)), 4.0))
    vals = summary["mean_lochness_in_own_cells"].to_numpy(dtype=float)
    ax.bar(range(len(summary)), vals, color=[HIT_RED if v > lcfg.enrichment_cut else GREY for v in vals])
    if "null_mean" in summary.columns:
        ax.errorbar(range(len(summary)), summary["null_mean"], yerr=2 * summary["null_sd"], fmt="_", color="black", lw=0.8, label="permutation null (mean +/- 2 SD)")
        ax.legend(fontsize=7, frameon=False)
    ax.axhline(0, color="black", lw=0.8); ax.axhline(lcfg.enrichment_cut, color="#718096", ls="--", lw=1)
    ax.set_xticks(range(len(summary))); ax.set_xticklabels(summary["target"], rotation=90, fontsize=6)
    ax.set_ylabel("mean lochNESS in its own cells")
    ax.set_title(f"How strongly each perturbation clusters with itself (k = {lo.n_neighbors} neighbours)", fontsize=10)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "lochness_self_enrichment", SECTION, ST_LOCH, "Self-enrichment per perturbation", "Average lochNESS of each perturbation's own cells. 0 means its cells sit among neighbours at the background rate; a positive value means cells sharing the perturbation are neighbours far more often than chance, i.e. the perturbation drives a distinct state.")
    # --- 2. distributions ----------------------------------------------------------------------
    top = list(summary["target"].head(30))
    fig, ax = plt.subplots(figsize=(max(6, 0.34 * len(top)), 4.2))
    data = [lo.scores[t].to_numpy(dtype=float) for t in top]
    data = [d[np.isfinite(d)] if np.isfinite(d).any() else np.array([0.0]) for d in data]
    parts = ax.violinplot(data, showextrema=False)
    for b in parts["bodies"]:
        b.set_alpha(0.7)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(1, len(top) + 1)); ax.set_xticklabels(top, rotation=90, fontsize=6)
    ax.set_ylabel("lochNESS"); ax.set_title("lochNESS across all cells, per perturbation (top 30)", fontsize=10)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "lochness_distributions", SECTION, ST_LOCH, "lochNESS distribution per perturbation", "Each violin is one perturbation's score across every cell. A long upper tail means a subset of the manifold is strongly enriched for it, even when most cells sit at background.")
    # --- 3. target x cluster heatmap ------------------------------------------------------------
    if lo.by_cluster is not None and not lo.by_cluster.empty:
        mat = lo.by_cluster
        order = order_by_similarity(mat)
        mat = mat.loc[order]
        lim = float(np.nanpercentile(np.abs(mat.to_numpy(dtype=float)), 98)) or 1.0
        fig, ax = plt.subplots(figsize=(max(6, 0.55 * mat.shape[1] + 4), max(4, 0.20 * len(mat) + 1.5)))
        im = ax.imshow(mat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, fontsize=8); ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index, fontsize=6)
        ax.set_xlabel(f"Cluster ({key})")
        plt.colorbar(im, ax=ax, shrink=0.6, label="mean lochNESS")
        ax.set_title("Mean lochNESS per cluster", fontsize=11)
        fig.tight_layout()
        reg.save(fig, "lochness_by_cluster", SECTION, ST_LOCH, "lochNESS by cluster", "Average score of each perturbation within each cluster, rows ordered by similarity. This is the continuous counterpart of the cluster-enrichment test - agreement between the two is a good sign, and structure here that the cluster test missed is worth a look.")
    # --- 4. self score on the embedding -------------------------------------------------------
    if coords is not None:
        vals = lo.self_score.to_numpy(dtype=float)
        ok = np.isfinite(vals)
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        ax.scatter(coords[~ok, 0], coords[~ok, 1], s=3, color="#edf2f7", linewidths=0, rasterized=True, label="unassigned / ambiguous")
        lim = float(np.nanpercentile(np.abs(vals[ok]), 98)) if ok.any() else 1.0
        lim = lim or 1.0
        sca = ax.scatter(coords[ok, 0], coords[ok, 1], s=5, c=vals[ok], cmap="RdBu_r", vmin=-lim, vmax=lim, linewidths=0, rasterized=True)
        plt.colorbar(sca, ax=ax, shrink=0.75, label="lochNESS (own perturbation)")
        ax.set_title("Each cell scored for its own perturbation", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.legend(fontsize=7, frameon=False, loc="best", markerscale=3)
        despine(ax, left=True, bottom=True)
        fig.tight_layout()
        reg.save(fig, "lochness_self_umap", SECTION, ST_LOCH, "Self lochNESS on the embedding", "Red regions are where cells sit among others sharing their own perturbation more often than chance - the parts of the manifold that perturbation identity actually organises.")
    # --- 5. one map per perturbation ------------------------------------------------------------
    if coords is None:
        return
    own_all = adata.obs["target"].astype(str).to_numpy()
    klass = adata.obs["perturbation_class"].astype(str).to_numpy()
    for rank, t in enumerate(summary["target"]):
        score = lo.scores[t].to_numpy(dtype=float)
        ok = np.isfinite(score)
        own = (own_all == t) & (klass == CLASS_TARGETING)
        row = summary[summary["target"] == t].iloc[0]
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
        ax = axes[0]
        srt = np.argsort(np.nan_to_num(score))
        vmax = float(np.nanmax(score[ok])) if ok.any() else 1.0
        sca = ax.scatter(coords[srt, 0], coords[srt, 1], s=4, c=np.nan_to_num(score[srt]), cmap="RdBu_r", norm=_lochness_norm(vmax), linewidths=0, rasterized=True)
        cbar = plt.colorbar(sca, ax=ax, shrink=0.78, label="lochNESS (symlog)"); cbar.ax.tick_params(labelsize=7)
        ax.set_title(f"{t}: neighbourhood enrichment", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        despine(ax, left=True, bottom=True)
        ax = axes[1]
        ax.scatter(coords[~own, 0], coords[~own, 1], s=3, color="#e2e8f0", linewidths=0, rasterized=True, label="other cells")
        ax.scatter(coords[own, 0], coords[own, 1], s=14, color=HIT_RED, linewidths=0.3, edgecolors="#2d3748", rasterized=True, label=f"{t} cells")
        label_clusters(ax, adata, key, coords, highlight=str(row.get("top_cluster", "")))
        ax.set_title(f"where the {int(own.sum()):,} {t} cells actually are (top cluster {row.get('top_cluster', '?')} highlighted)", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.legend(fontsize=7, frameon=False, loc="best", markerscale=2)
        despine(ax, left=True, bottom=True)
        fig.tight_layout()
        reg.save(fig, f"lochness_{t}", SECTION, ST_LOCH_TARGET, f"{t} lochNESS map",
                 f"Left: every cell scored for how enriched {t} is among its neighbours (red = enriched, blue = depleted). Right: the cells actually carrying {t}, for comparison. Mean score in its own cells {row['mean_lochness_in_own_cells']:.2f}; {row['pct_cells_enriched']:.1f}% of all cells score above {lcfg.enrichment_cut}.",
                 in_report=rank < top_n)
    logger.info("Wrote %d per-perturbation lochNESS maps (%d shown in the report)", len(summary), min(top_n, len(summary)))


# ---------------------------------------------------------------------------- protein extension
def protein_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    pr = res.protein
    alpha = cfg.analysis.perturbation_effects.concordance.fdr_alpha
    ps, mo = res.ps, res.modules
    if pr is not None and not pr.empty:
        D = pr.d_matrix.copy()
        order = D.abs().max(axis=1).sort_values(ascending=False).index
        D = D.loc[order]
        F = pr.fdr_matrix.loc[order, D.columns] < cfg.analysis.perturbation_effects.protein.fdr_alpha
        fig, ax = plt.subplots(figsize=(1.8 + 0.7 * D.shape[1], 1.2 + 0.2 * D.shape[0]))
        _heatmap(ax, D, label="Cohen's d (perturbed - control)", marks=F)
        ax.set_title("perturbation x protein", fontsize=9)
        fig.tight_layout()
        reg.save(fig, "protein_effect_heatmap", SECTION, ST_PROT, "Protein effects", f"Standardized effect (Cohen's d) of each perturbation on each measured protein ({pr.info.get('representation')} values vs non-targeting controls); * = BH-FDR < {cfg.analysis.perturbation_effects.protein.fdr_alpha} (Mann-Whitney U).")
    co = res.concordance
    if co is None:
        return
    pp = co.ps_protein
    if pp is not None and not pp.empty and ps is not None:
        w = pp[(pp["scope"] == "within_target") & (pp["status"] != "insufficient_support")].copy()
        w["abs"] = w["spearman_rho"].abs()
        w = w.sort_values(["fdr", "abs", "target"], ascending=[True, False, True]).head(6)
        if not w.empty:
            P = adata.obsm[cfg.analysis.perturbation_effects.protein.representation]
            n = len(w)
            fig, axes = plt.subplots(1, n, figsize=(2.8 * n, 2.8), squeeze=False)
            for ax, (_, r) in zip(axes[0], w.iterrows()):
                m = (adata.obs["target"].astype(str) == r["target"]).to_numpy() & ps.own.notna().to_numpy()
                x, y = ps.own[m].to_numpy(float), P.loc[m, r["protein"]].to_numpy(float)
                ax.scatter(x, y, s=5, alpha=0.6, color=BLUE); _trend(ax, x, y)
                ax.set_title(f"{r['target']}: {r['protein']}\nrho={r['spearman_rho']:+.2f}, FDR={r['fdr']:.2g}, n={r['n_cells']}", fontsize=7)
                ax.set_xlabel("PS", fontsize=8); ax.set_ylabel(f"{r['protein']} (CLR)", fontsize=8)
            fig.tight_layout()
            reg.save(fig, "ps_vs_protein", SECTION, ST_CONC, "PS vs protein (within target)", "Per-cell PS vs protein value inside one target's perturbed cells, for the pairs with the lowest FDR; red = binned-median trend (visual aid). The statistic is the within-target Spearman rho in ps_protein_association.csv.")
    ppc = co.program_protein_cells
    if ppc is not None and not ppc.empty and mo is not None and not mo.empty and mo.program_activity is not None:
        w = ppc.copy(); w["abs"] = w["spearman_rho"].abs()
        w = w.sort_values(["fdr", "abs"], ascending=[True, False]).head(4)
        P = adata.obsm[cfg.analysis.perturbation_effects.protein.representation]
        cells = adata.obs["perturbation_class"].astype(str).isin([CLASS_TARGETING, CLASS_CONTROL]).to_numpy() & P.notna().all(axis=1).to_numpy()
        n = len(w)
        fig, axes = plt.subplots(1, n, figsize=(2.8 * n, 2.8), squeeze=False)
        for ax, (_, r) in zip(axes[0], w.iterrows()):
            x = mo.program_activity.loc[cells, r["gene_program"]].to_numpy(float); y = P.loc[cells, r["protein"]].to_numpy(float)
            ax.hexbin(x, y, gridsize=30, cmap="Blues", mincnt=1); _trend(ax, x, y)
            ax.set_title(f"{r['gene_program']} vs {r['protein']}\nrho={r['spearman_rho']:+.2f}, FDR={r['fdr']:.2g}", fontsize=7)
            ax.set_xlabel(f"{r['gene_program']} activity", fontsize=8); ax.set_ylabel(f"{r['protein']} (CLR)", fontsize=8)
        fig.tight_layout()
        reg.save(fig, "program_activity_vs_protein", SECTION, ST_CONC, "Gene-program activity vs protein", "Cell-level program activity vs protein value over single-guide and control cells, strongest pairs by FDR; association only, not mediation.")
    sm = co.summary
    if sm is not None and not sm.empty and "protein_effect_magnitude" in sm.columns:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
        for ax, (col, lab) in zip(axes, (("rna_effect_magnitude", "RNA effect (mean |log2FC| of panel genes)"), ("lochness_own_mean", "lochNESS in own cells"))):
            if col not in sm.columns:
                ax.axis("off"); continue
            ok = sm[col].notna() & sm["protein_effect_magnitude"].notna()
            ax.scatter(sm.loc[ok, col], sm.loc[ok, "protein_effect_magnitude"], s=np.clip(sm.loc[ok, "n_cells"], 10, 200), color=BLUE, alpha=0.7)
            for _, r in sm[ok].iterrows():
                ax.annotate(r["target"], (r[col], r["protein_effect_magnitude"]), fontsize=6)
            ax.set_xlabel(lab, fontsize=8); ax.set_ylabel("max |Cohen's d| over proteins", fontsize=8)
        fig.tight_layout()
        reg.save(fig, "rna_vs_protein_effects", SECTION, ST_CONC, "RNA vs protein effect magnitude", "Target-level comparison: transcriptional effect size (left) and lochNESS state shift (right) vs the largest protein effect; point size = cells. Across-target association statistics are in lochness_protein_summary.csv.")
        cols = [c for c in ("ps_auc_vs_control", "lochness_own_mean", "rna_effect_magnitude", "protein_effect_magnitude") if c in sm.columns]
        O = sm.set_index("target")[cols].astype(float)
        O = O.loc[O.abs().sum(axis=1).sort_values(ascending=False).index]
        Z = (O - O.mean()) / O.std(ddof=0).replace(0, 1)
        fig, ax = plt.subplots(figsize=(1.5 + 0.8 * Z.shape[1], 1.2 + 0.2 * Z.shape[0]))
        _heatmap(ax, Z, label="z-score across targets")
        fig.tight_layout()
        reg.save(fig, "perturbation_overview", SECTION, ST_CONC, "Integrated perturbation overview", "Per target: PS separation from controls (AUC), lochNESS self-enrichment, RNA effect magnitude and protein effect magnitude, each z-scored across targets (values in perturbation_summary.csv).")


def perturbation_effect_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    """Reference order: modules / programs, PS, lochNESS; then the protein extension."""
    modules_figures(res.modules, adata, cfg, reg)
    ps_figures(res.ps, res.ps_vs_strength, res.strength, adata, cfg, reg)
    lochness_figures(res.lochness, adata, cfg, reg)
    protein_figures(res, adata, cfg, reg)
