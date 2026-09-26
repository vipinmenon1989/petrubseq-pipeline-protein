"""Perturbation x cluster enrichment figures (reference ``plot_enrichment`` / ``plot_enrichment_per_target``)."""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..analysis._common import CLASS_TARGETING  # noqa: E402
from ..analysis.cluster_enrichment import enrichment_matrix, phenocopy_similarity, significance_matrix  # noqa: E402
from ..analysis.perturbation_strength import CONTROL_LABELS  # noqa: E402
from ._refplot import GREY, HIT_RED, despine, order_by_similarity, umap_coords  # noqa: E402
from .figures import FigureRegistry  # noqa: E402

logger = logging.getLogger("petrubseq_protein")
SECTION = "cell_states"
ST_ENRICH, ST_PER_TARGET = "enrichment", "enrichment_per_target"


def _has_hit(rows: pd.DataFrame, cluster) -> bool:
    try:
        return bool(rows.at[cluster, "significant"])
    except (KeyError, ValueError):
        return False


def enrichment_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    if res is None or res.empty or res.table.empty:
        return
    ecfg = cfg.analysis.clustering.enrichment
    key = res.cluster_key
    lor = enrichment_matrix(res)
    fdr = significance_matrix(res)
    if lor.empty:
        return
    order = order_by_similarity(lor)
    lor = lor.loc[order]; fdr = fdr.loc[order]
    # --- 1. heatmap ---------------------------------------------------------------------------
    lim = float(np.nanpercentile(np.abs(lor.to_numpy()), 98)) or 1.0
    fig, ax = plt.subplots(figsize=(max(6.0, 0.55 * lor.shape[1] + 4), max(4.0, 0.20 * len(lor) + 1.5)))
    im = ax.imshow(lor.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(lor.shape[1])); ax.set_xticklabels(lor.columns, fontsize=8)
    ax.set_yticks(range(lor.shape[0])); ax.set_yticklabels(lor.index, fontsize=6)
    ax.set_xlabel(f"Cluster ({key})")
    n_marked = 0
    for i in range(lor.shape[0]):
        for j in range(lor.shape[1]):
            v = fdr.iat[i, j]
            if pd.notna(v) and v < ecfg.fdr_alpha:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9, color="black"); n_marked += 1
    plt.colorbar(im, ax=ax, shrink=0.6, label="log2 odds ratio")
    ax.set_title(f"Perturbation enrichment across clusters\n* FDR < {ecfg.fdr_alpha} vs {CONTROL_LABELS[res.primary_control]}", fontsize=10)
    fig.tight_layout()
    reg.save(fig, "enrichment_heatmap", SECTION, ST_ENRICH, "Perturbation enrichment across clusters", f"Red means cells carrying that perturbation are over-represented in the cluster, blue under-represented; asterisks mark significant pairs. Rows are ordered by profile similarity, so perturbations with the same phenotype sit together ({n_marked} significant pair(s)).")
    # --- 2. phenocopy -------------------------------------------------------------------------
    sim = phenocopy_similarity(res)
    if not sim.empty and sim.shape[0] >= 3:
        so = order_by_similarity(sim)
        sim = sim.loc[so, so]
        size = max(5.0, 0.16 * len(sim) + 2)
        fig, ax = plt.subplots(figsize=(size, size))
        im = ax.imshow(sim.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(sim))); ax.set_xticklabels(sim.index, rotation=90, fontsize=5)
        ax.set_yticks(range(len(sim))); ax.set_yticklabels(sim.index, fontsize=5)
        plt.colorbar(im, ax=ax, shrink=0.6, label="Pearson r")
        ax.set_title("Do perturbations phenocopy each other?", fontsize=11)
        fig.tight_layout()
        reg.save(fig, "enrichment_phenocopy", SECTION, ST_ENRICH, "Perturbation similarity (phenocopy map)", "Correlation between targets of their cluster-composition profiles. Red blocks are groups of perturbations producing the same cell-state shift - subunits of one complex are expected to land together, a check that needs no prior knowledge of the complexes.")
    # --- 3. composition bars ------------------------------------------------------------------
    comp = res.composition
    mag = res.effect_magnitude
    top = list(mag["target"].head(30))
    ref = res.reference_composition[res.primary_control]
    plot_df = pd.concat([ref.to_frame("REFERENCE").T, comp.loc[top]])
    fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(plot_df) + 2), 4.4))
    bottom = np.zeros(len(plot_df))
    cmap = plt.get_cmap("tab20")
    for i, cl in enumerate(plot_df.columns):
        ax.bar(range(len(plot_df)), plot_df[cl].to_numpy(dtype=float), bottom=bottom, color=cmap(i % 20), label=str(cl))
        bottom += plot_df[cl].to_numpy(dtype=float)
    ax.set_xticks(range(len(plot_df))); ax.set_xticklabels(plot_df.index, rotation=90, fontsize=6)
    ax.set_ylabel("% of cells"); ax.set_xlim(-0.6, len(plot_df) - 0.4)
    ax.axvline(0.5, color="black", lw=1.2)
    ax.legend(title="cluster", fontsize=6, title_fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", ncol=1)
    ax.set_title("Cluster composition per perturbation (top 30 by shift)", fontsize=10)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "enrichment_composition", SECTION, ST_ENRICH, "Cluster composition per perturbation", "Each bar is one perturbation's distribution across clusters; the leftmost bar (left of the black line) is the reference. Perturbations are sorted by how far their composition sits from it.")
    # --- 4. volcano ---------------------------------------------------------------------------
    sub = res.table[res.table["control"] == res.primary_control]
    x = sub["log2_odds_ratio"].to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        y = -np.log10(np.clip(sub["fdr"].to_numpy(dtype=float), 1e-300, 1))
    sig = sub["significant"].to_numpy(dtype=bool)
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.scatter(x[~sig], y[~sig], s=12, color=GREY, label="not significant")
    ax.scatter(x[sig], y[sig], s=20, color=HIT_RED, label=f"FDR < {ecfg.fdr_alpha}")
    ax.axhline(-np.log10(ecfg.fdr_alpha), color="#718096", ls="--", lw=1); ax.axvline(0, color="#718096", ls="--", lw=1)
    labelled = sub[sig].reindex(sub[sig]["log2_odds_ratio"].abs().sort_values(ascending=False).index).head(10)
    for _, row in labelled.iterrows():
        ax.annotate(f"{row['target']}:{row['cluster']}", (row["log2_odds_ratio"], -np.log10(max(row["fdr"], 1e-300))), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("log2 odds ratio (enriched >0, depleted <0)"); ax.set_ylabel("-log10 FDR")
    ax.set_title("Enrichment across all target x cluster pairs", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "enrichment_volcano", SECTION, ST_ENRICH, "Enrichment volcano", "Every target/cluster pair. Points to the right are perturbations that accumulate in a cluster; to the left, ones depleted from it.")
    # --- 5. effect magnitude ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(max(6, 0.16 * len(mag)), 3.8))
    ax.bar(range(len(mag)), mag["composition_shift_pct"], color=[HIT_RED if n > 0 else GREY for n in mag["n_significant_clusters"]])
    ax.set_xticks(range(len(mag))); ax.set_xticklabels(mag["target"], rotation=90, fontsize=6)
    ax.set_ylabel("Composition shift (%)")
    ax.set_title("How far each perturbation moves cells between clusters (red = has a significant cluster)", fontsize=10)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "enrichment_effect_magnitude", SECTION, ST_ENRICH, "Composition shift per perturbation", "Total variation distance between each perturbation's cluster composition and the reference: 0% means indistinguishable, 100% means the cells sit in entirely different clusters.")
    # --- 6. per target ------------------------------------------------------------------------
    obs = adata.obs
    targets_col = obs["target"].astype(str).to_numpy()
    klass = obs["perturbation_class"].astype(str).to_numpy()
    clusters_col = obs[key].astype(str).to_numpy()
    coords = umap_coords(adata)
    clusters = list(res.composition.columns)
    ranked = list(mag["target"])
    hit_targets = res.targets_with_hits()
    ordered = [t for t in ranked if t in hit_targets] + [t for t in ranked if t not in hit_targets]
    top_n = ecfg.top_n_report
    sub_tbl = res.table[res.table["control"] == res.primary_control]
    palette = {c: plt.get_cmap("tab20")(i % 20) for i, c in enumerate(clusters)}
    for rank, t in enumerate(ordered):
        if t not in res.composition.index:
            continue
        comp_t = res.composition.loc[t]
        rows = sub_tbl[sub_tbl["target"] == t].set_index("cluster")
        n_panels = 2 if coords is None else 3
        fig, axes = plt.subplots(1, n_panels, figsize=(4.7 * n_panels, 4.0))
        ax = axes[0]
        idx = np.arange(len(clusters))
        ax.bar(idx - 0.2, ref[clusters].to_numpy(dtype=float), width=0.4, color=GREY, label="reference")
        ax.bar(idx + 0.2, comp_t[clusters].to_numpy(dtype=float), width=0.4, color="#2b6cb0", label=t)
        for i, cl in enumerate(clusters):
            if cl in rows.index and bool(rows.at[cl, "significant"]):
                ax.text(i, max(comp_t[cl], ref[cl]) + 1, "*", ha="center", fontsize=11)
        ax.set_xticks(idx); ax.set_xticklabels(clusters, fontsize=7)
        ax.set_xlabel(f"Cluster ({key})"); ax.set_ylabel("% of cells")
        ax.legend(fontsize=7, frameon=False); ax.set_title(f"{t} cluster composition", fontsize=10)
        despine(ax)
        ax = axes[1]
        lor_t = rows["log2_odds_ratio"].reindex(clusters)
        sig_t = rows["significant"].reindex(clusters).fillna(False).to_numpy(dtype=bool)
        ax.bar(idx, np.nan_to_num(lor_t.to_numpy(dtype=float)), color=[HIT_RED if s else GREY for s in sig_t])
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(idx); ax.set_xticklabels(clusters, fontsize=7)
        ax.set_xlabel(f"Cluster ({key})"); ax.set_ylabel("log2 odds ratio")
        ax.set_title(f"{t} enrichment (red = FDR < {ecfg.fdr_alpha})", fontsize=10)
        despine(ax)
        if coords is not None:
            ax = axes[2]
            pert = (targets_col == t) & (klass == CLASS_TARGETING)
            ax.scatter(coords[~pert, 0], coords[~pert, 1], s=2, color="#e2e8f0", linewidths=0, rasterized=True)
            ax.scatter(coords[pert, 0], coords[pert, 1], s=10, c=[palette.get(c, (0.3, 0.3, 0.3, 1)) for c in clusters_col[pert]], linewidths=0.2, edgecolors="#2d3748", rasterized=True)
            ax.set_title(f"{t} cells by cluster (n={int(pert.sum()):,})", fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
            despine(ax, left=True, bottom=True)
        best = rows["log2_odds_ratio"].abs().idxmax() if len(rows) and rows["log2_odds_ratio"].notna().any() else None
        caption = f"Cluster distribution of cells perturbed for {t}."
        if best is not None and _has_hit(rows, best):
            caption += f" Strongest association: cluster {best} ({rows.at[best, 'pct_of_target']:.1f}% of {t} cells vs {rows.at[best, 'pct_of_reference']:.1f}% of reference, FDR = {rows.at[best, 'fdr']:.2g})."
        fig.tight_layout()
        reg.save(fig, f"enrichment_{t}", SECTION, ST_PER_TARGET, f"{t} cluster enrichment", caption, in_report=rank < top_n)
    logger.info("Wrote %d per-target enrichment figures (%d shown in the report)", len(ordered), min(top_n, len(ordered)))
