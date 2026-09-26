"""Perturbation-strength figures (reference ``plot_perturbation_overview`` / ``plot_per_target``)."""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..analysis._common import CLASS_CONTROL, CLASS_TARGETING  # noqa: E402
from ..analysis.perturbation_strength import CONTROL_LABELS, CONTROL_NTC, CONTROL_OTHER, control_masks  # noqa: E402
from ._refplot import GREY, HIT_RED, despine, ecdf, gene_values, umap_coords, violin  # noqa: E402
from .figures import FigureRegistry  # noqa: E402

logger = logging.getLogger("petrubseq_protein")
SECTION = "perturbation_strength"
ST_OVERVIEW, ST_PER_TARGET = "overview", "per_target"


def perturbation_strength_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    if res is None or res.empty or res.table.empty:
        return
    tbl = res.table
    primary = res.primary_control
    scfg = cfg.analysis.perturbation_effects.strength
    lfc = tbl[f"log2fc_{primary}"].to_numpy(dtype=float)
    fdr = tbl[f"ks_fdr_{primary}"].to_numpy(dtype=float)
    hit = tbl[f"is_hit_{primary}"].to_numpy(dtype=bool)
    names = tbl["target"].to_numpy()
    # --- volcano ---------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    with np.errstate(divide="ignore"):
        y = -np.log10(np.clip(fdr, 1e-300, 1))
    ax.scatter(lfc[~hit], y[~hit], s=18, color=GREY, label="not significant")
    ax.scatter(lfc[hit], y[hit], s=22, color=HIT_RED, label="effective knockdown")
    ax.axhline(-np.log10(scfg.fdr_alpha), color="#718096", ls="--", lw=1)
    ax.axvline(scfg.max_log2fc_for_hit, color="#718096", ls="--", lw=1)
    ok = np.isfinite(lfc)
    for i in np.argsort(np.where(ok, lfc, np.inf))[: min(12, int(ok.sum()))]:
        ax.annotate(names[i], (lfc[i], y[i]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("log2 fold change (perturbed / control)")
    ax.set_ylabel("-log10 FDR (KS test)")
    ax.set_title(f"Perturbation strength - control: {CONTROL_LABELS[primary]}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "perturbation_volcano", SECTION, ST_OVERVIEW, "Volcano of perturbation strength", "Each point is a target gene, comparing its own expression in perturbed vs control cells. Effective perturbations fall in the upper-left: significantly reduced expression.")
    # --- waterfall -------------------------------------------------------------------------
    order = np.argsort(np.where(ok, lfc, np.inf))
    fig, ax = plt.subplots(figsize=(max(6, 0.18 * len(order)), 4.0))
    ax.bar(range(len(order)), np.nan_to_num(lfc[order]), color=[HIT_RED if hit[i] else GREY for i in order])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(names[order], rotation=90, fontsize=6)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("log2 fold change")
    ax.set_title(f"Knockdown per target ({int(hit.sum())}/{len(hit)} effective at FDR < {scfg.fdr_alpha})", fontsize=10)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "perturbation_waterfall", SECTION, ST_OVERVIEW, "Knockdown strength per target", "Targets sorted by fold change. Red bars are significant reductions; bars near zero mean the guide did not measurably reduce its target.")
    # --- control comparison ------------------------------------------------------------------
    if len(res.controls_used) > 1:
        other = CONTROL_OTHER if primary == CONTROL_NTC else CONTROL_NTC
        lfc2 = tbl[f"log2fc_{other}"].to_numpy(dtype=float)
        ok2 = ~(np.isnan(lfc) | np.isnan(lfc2))
        if ok2.sum() > 2:
            fig, ax = plt.subplots(figsize=(4.8, 4.6))
            ax.scatter(lfc[ok2], lfc2[ok2], s=20, color="#2b6cb0", alpha=0.8)
            lim = [min(lfc[ok2].min(), lfc2[ok2].min()) - 0.1, max(lfc[ok2].max(), lfc2[ok2].max()) + 0.1]
            ax.plot(lim, lim, ls="--", color="#718096", lw=1)
            r = float(np.corrcoef(lfc[ok2], lfc2[ok2])[0, 1])
            ax.set_xlabel(f"log2FC vs {CONTROL_LABELS[primary]}")
            ax.set_ylabel(f"log2FC vs {CONTROL_LABELS[other]}")
            ax.set_title(f"Control comparison (Pearson r = {r:.2f})", fontsize=10)
            despine(ax)
            fig.tight_layout()
            reg.save(fig, "perturbation_control_comparison", SECTION, ST_OVERVIEW, "Effect size under both control definitions", "Agreement between the two control groups. Points off the diagonal are targets whose apparent effect depends on the control used.")
    # --- per target -------------------------------------------------------------------------
    rng = np.random.default_rng(cfg.compute.seed)
    obs = adata.obs
    targets_col = obs["target"].astype(str).to_numpy()
    klass = obs["perturbation_class"].astype(str).to_numpy()
    base = control_masks(adata, cfg.analysis.perturbation_effects.control_classes)
    coords = umap_coords(adata)
    top_n = scfg.top_n_report
    for i, row in tbl.iterrows():
        t = str(row["target"]); gene = str(row.get("gene", t)) or t
        if gene not in adata.var_names:
            continue
        values = gene_values(adata, gene)
        pert = (targets_col == t) & (klass == CLASS_TARGETING)
        groups = {"perturbed": values[pert]}
        for control in res.controls_used:
            m = base[CONTROL_NTC] if control == CONTROL_NTC else (base[CONTROL_OTHER] & (targets_col != t))
            groups[CONTROL_LABELS[control]] = values[m]
        n_panels = 3 if coords is not None else 2
        fig, axes = plt.subplots(1, n_panels, figsize=(4.6 * n_panels, 4.0))
        ax = axes[0]
        for label, vals in groups.items():
            if vals.size:
                ecdf(ax, vals, f"{label} (n={vals.size:,})")
        ax.set_xlabel(f"{gene} expression (log-normalized)")
        ax.set_ylabel("Cumulative fraction of cells")
        ax.legend(fontsize=7, frameon=False, loc="lower right")
        fdr_i = float(row.get(f"ks_fdr_{primary}", np.nan)); lfc_i = float(row.get(f"log2fc_{primary}", np.nan))
        ax.set_title(f"{t}: log2FC = {lfc_i:.2f}, FDR = {fdr_i:.2g}", fontsize=10)
        despine(ax)
        ax = axes[1]
        violin(ax, list(groups.values()), list(groups.keys()), f"{gene} expression")
        ax.set_title(f"{gene} expression by group", fontsize=10)
        despine(ax)
        if coords is not None:
            ax = axes[2]
            frac = scfg.umap_background_fraction
            bg = rng.random(adata.n_obs) < frac
            keep = bg | pert
            ax.scatter(coords[keep & ~pert, 0], coords[keep & ~pert, 1], s=3, color="#e2e8f0", linewidths=0, rasterized=True, label=f"other cells ({int(frac * 100)}% shown)")
            sca = ax.scatter(coords[pert, 0], coords[pert, 1], s=10, c=values[pert], cmap="Reds", vmin=0, edgecolors="#2d3748", linewidths=0.25, rasterized=True)
            plt.colorbar(sca, ax=ax, shrink=0.75, label=f"{gene} expression")
            ax.set_title(f"Cells perturbed for {t} (n={int(pert.sum()):,})", fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
            ax.legend(fontsize=7, frameon=False, loc="best", markerscale=3)
            despine(ax, left=True, bottom=True)
        fig.tight_layout()
        is_hit = bool(row.get(f"is_hit_{primary}", False))
        reg.save(fig, f"perturbation_{t}", SECTION, ST_PER_TARGET, f"{t} perturbation effect",
                 f"Expression of {gene} in cells carrying {t} guides versus control cells. log2FC = {lfc_i:.2f}, KS FDR = {fdr_i:.2g}" + (" - effective knockdown." if is_hit else " - no significant reduction."),
                 in_report=i < top_n)
    logger.info("Wrote %d per-target perturbation figures (%d shown in the report)", len(tbl), min(top_n, len(tbl)))
