"""All figures of a run. Each ``*_figures`` orchestrator draws with matplotlib
and registers through the :class:`FigureRegistry`; nothing here knows the
final report layout. Before/after QC figures share the same drawing code and
differ only in ``stage``."""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..config import Config  # noqa: E402
from .figures import (  # noqa: E402
    SECTION_MULTIMODAL,
    SECTION_QC_PERTURBATION,
    SECTION_QC_PROTEIN,
    SECTION_QC_RNA,
    SECTION_REP_PROTEIN,
    SECTION_REP_RNA,
    STAGE_POST,
    STAGE_REP,
    FigureRegistry,
)

logger = logging.getLogger("petrubseq_protein")

MAX_CATEGORIES = 30
MAX_POINTS = 100_000
BLUE, ORANGE, RED, GREY = "#4C72B0", "#DD8452", "#C44E52", "#8C8C8C"


def _stage_label(stage: str) -> str:
    return stage.replace("_", " ")


def _sample_idx(n: int, seed: int = 0) -> np.ndarray:
    if n <= MAX_POINTS:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, MAX_POINTS, replace=False))


def choose_group(obs: pd.DataFrame, candidates: Sequence[str] = ("condition", "sample", "lane", "batch", "donor")) -> Optional[str]:
    for c in candidates:
        if c in obs.columns and 2 <= obs[c].nunique() <= MAX_CATEGORIES:
            return c
    return None


def _hist(ax, v: np.ndarray, label: str, log10: bool = False, thresholds: Optional[List] = None, color: str = BLUE) -> None:
    v = v[np.isfinite(v)]
    if log10:
        v = np.log10(v[v > 0])
        label = f"log10({label})"
    ax.hist(v, bins=60, color=color)
    ax.set_xlabel(label, fontsize=9)
    ax.set_ylabel("cells", fontsize=9)
    for thr, txt in thresholds or []:
        if thr is None:
            continue
        x = np.log10(thr) if log10 and thr > 0 else thr
        ax.axvline(x, color=RED, ls="--", lw=1)
        ax.text(x, ax.get_ylim()[1] * 0.95, txt, color=RED, fontsize=7, rotation=90, va="top", ha="right")


# ---------------------------------------------------------------------------
# RNA QC (before / after)
# ---------------------------------------------------------------------------


def rna_qc_figures(obs: pd.DataFrame, cfg: Config, reg: FigureRegistry, stage: str) -> None:
    fr = cfg.qc.filter.rna
    on = cfg.qc.filter.enabled
    sl = _stage_label(stage)
    thr_note = "" if on else " (strict filter disabled: thresholds shown for reference)"
    metrics = [
        ("n_genes_by_counts", "genes detected", False, [(fr.min_genes, f"min_genes={fr.min_genes}")]),
        ("total_counts", "total RNA counts", True, [(fr.min_counts, f"min_counts={fr.min_counts}")]),
        ("pct_counts_mt", "% mitochondrial", False, [(fr.max_pct_mt, f"max_pct_mt={fr.max_pct_mt}")]),
        ("pct_counts_ribo", "% ribosomal", False, []),
    ]
    metrics = [m for m in metrics if m[0] in obs.columns]
    if metrics:
        fig, axes = plt.subplots(1, len(metrics), figsize=(3.6 * len(metrics), 3.2))
        for ax, (col, label, lg, thr) in zip(np.atleast_1d(axes), metrics):
            _hist(ax, obs[col].to_numpy(float), label, lg, thr)
        fig.suptitle(f"RNA QC distributions ({sl}, n = {len(obs):,} cells)", fontsize=11)
        fig.tight_layout()
        reg.save(fig, f"rna_qc_distributions_{stage}", SECTION_QC_RNA, stage, f"RNA QC distributions ({sl})",
                 f"Per-cell genes detected, total counts, mitochondrial and ribosomal fractions; dashed lines are the configured qc.filter.rna thresholds{thr_note}.")
    if {"total_counts", "n_genes_by_counts"} <= set(obs.columns):
        idx = _sample_idx(len(obs))
        x, y = obs["total_counts"].to_numpy(float)[idx], obs["n_genes_by_counts"].to_numpy(float)[idx]
        c = obs["pct_counts_mt"].to_numpy(float)[idx] if "pct_counts_mt" in obs else None
        fig, ax = plt.subplots(figsize=(5, 4.2))
        sca = ax.scatter(x, y, s=2, c=c, cmap="viridis", alpha=0.5, rasterized=True)
        if c is not None:
            fig.colorbar(sca, ax=ax, label="% mitochondrial")
        ax.set_xscale("log")
        ax.set_xlabel("total RNA counts")
        ax.set_ylabel("genes detected")
        if fr.min_genes:
            ax.axhline(fr.min_genes, color=RED, ls="--", lw=1)
        if fr.min_counts:
            ax.axvline(fr.min_counts, color=RED, ls="--", lw=1)
        ax.set_title(f"Counts vs genes ({sl})", fontsize=11)
        fig.tight_layout()
        reg.save(fig, f"rna_counts_vs_genes_{stage}", SECTION_QC_RNA, stage, f"Total counts vs genes detected ({sl})",
                 "Each point is a cell (subsampled to 100k), coloured by mitochondrial fraction; dashed lines are the configured thresholds.")
    g = choose_group(obs)
    if g:
        counts = obs[g].astype(str).value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(max(4, 0.6 * len(counts)) + 1, 3.2))
        ax.bar(range(len(counts)), counts.to_numpy(), color=GREY)
        ax.set_xticks(range(len(counts)))
        ax.set_xticklabels(counts.index, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("cells")
        ax.set_title(f"Cells per {g} ({sl})", fontsize=11)
        fig.tight_layout()
        reg.save(fig, f"rna_cells_per_{g}_{stage}", SECTION_QC_RNA, stage, f"Cells per {g} ({sl})", f"Cell counts per {g} at this stage.")
        for col, label in (("n_genes_by_counts", "genes detected"), ("pct_counts_mt", "% mitochondrial")):
            if col in obs.columns:
                fig = _violin_by_group(obs, col, g, label)
                reg.save(fig, f"rna_{col}_by_{g}_{stage}", SECTION_QC_RNA, stage, f"{label} by {g} ({sl})", f"Distribution of {label} within each {g}.")


def _violin_by_group(obs: pd.DataFrame, col: str, group: str, label: str, log: bool = False):
    cats = sorted(obs[group].astype(str).unique())
    data = []
    for c in cats:
        v = obs.loc[obs[group].astype(str) == c, col].to_numpy(float)
        v = v[np.isfinite(v)]
        if log:
            v = np.log10(v[v > 0])
        data.append(v if v.size else np.array([0.0]))
    fig, ax = plt.subplots(figsize=(max(4, 0.6 * len(cats)) + 1, 3.2))
    ax.violinplot(data, showmedians=True)
    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"log10({label})" if log else label, fontsize=9)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Protein QC (before / after)
# ---------------------------------------------------------------------------


def protein_qc_figures(obs: pd.DataFrame, counts: Optional[pd.DataFrame], prot: Optional[pd.DataFrame], table: Optional[pd.DataFrame], cfg: Config, reg: FigureRegistry, stage: str, primary_desc: str = "") -> None:
    sl = _stage_label(stage)
    fp = cfg.qc.filter.protein
    if counts is not None and "protein_total_counts" in obs.columns:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.2))
        _hist(axes[0], obs["protein_total_counts"].to_numpy(float), "total ADT counts", True, [(fp.min_total_counts, f"min_total_counts={fp.min_total_counts}")], ORANGE)
        _hist(axes[1], obs["protein_n_detected"].to_numpy(float), "antibodies detected", False, [(fp.min_proteins_detected, f"min_proteins_detected={fp.min_proteins_detected}")], ORANGE)
        _hist(axes[2], obs["protein_pct_isotype"].to_numpy(float), "% isotype counts", False, [(fp.max_pct_isotype, f"max_pct_isotype={fp.max_pct_isotype}")], ORANGE)
        fig.suptitle(f"Protein QC distributions ({sl}, n = {int(obs['protein_total_counts'].notna().sum()):,} cells with ADT counts)", fontsize=11)
        fig.tight_layout()
        reg.save(fig, f"protein_qc_distributions_{stage}", SECTION_QC_PROTEIN, stage, f"Protein QC distributions ({sl})",
                 "Total ADT UMIs per cell, antibodies with >= 1 count, and the share of ADT counts on isotype controls; dashed lines are the configured qc.filter.protein thresholds.")
        # targeting vs isotype
        t = obs.get("protein_total_counts_targeting"); i = obs.get("protein_isotype_counts")
        if t is not None and i is not None:
            ok = np.isfinite(t.to_numpy(float)) & np.isfinite(i.to_numpy(float))
            idx = np.where(ok)[0][_sample_idx(int(ok.sum()))]
            ext = obs["protein_extreme_counts"].to_numpy(bool)[idx] if "protein_extreme_counts" in obs else np.zeros(len(idx), bool)
            fig, ax = plt.subplots(figsize=(5, 4.2))
            ax.scatter(np.log10(t.to_numpy(float)[idx] + 1), np.log10(i.to_numpy(float)[idx] + 1), s=2, color=ORANGE, alpha=0.4, rasterized=True, label="cells")
            if ext.any():
                ax.scatter(np.log10(t.to_numpy(float)[idx][ext] + 1), np.log10(i.to_numpy(float)[idx][ext] + 1), s=8, color=RED, label=f"extreme ADT ({int(ext.sum())})")
            ax.set_xlabel("log10(targeting ADT counts + 1)")
            ax.set_ylabel("log10(isotype ADT counts + 1)")
            ax.legend(fontsize=7, frameon=False)
            ax.set_title(f"Targeting vs isotype ADT ({sl})", fontsize=11)
            fig.tight_layout()
            reg.save(fig, f"protein_targeting_vs_isotype_{stage}", SECTION_QC_PROTEIN, stage, f"Targeting vs isotype ADT counts ({sl})",
                     "Isotype-control counts track non-specific background; cells far above the cloud in both axes are antibody aggregates (extreme ADT flag).")
        # RNA depth vs ADT depth
        if "total_counts" in obs.columns:
            r = obs["total_counts"].to_numpy(float); p = obs["protein_total_counts"].to_numpy(float)
            ok = np.isfinite(r) & np.isfinite(p) & (r > 0) & (p > 0)
            idx = np.where(ok)[0][_sample_idx(int(ok.sum()))]
            fig, ax = plt.subplots(figsize=(5, 4.2))
            c = obs["pct_counts_mt"].to_numpy(float)[idx] if "pct_counts_mt" in obs else None
            sca = ax.scatter(np.log10(r[idx]), np.log10(p[idx]), s=2, c=c, cmap="viridis", alpha=0.4, rasterized=True)
            if c is not None:
                fig.colorbar(sca, ax=ax, label="% mitochondrial")
            ax.set_xlabel("log10(total RNA counts)")
            ax.set_ylabel("log10(total ADT counts)")
            rr = np.corrcoef(np.log10(r[idx]), np.log10(p[idx]))[0, 1]
            ax.set_title(f"RNA depth vs ADT depth ({sl}), r = {rr:.2f}", fontsize=11)
            fig.tight_layout()
            reg.save(fig, f"protein_depth_vs_rna_depth_{stage}", SECTION_QC_PROTEIN, stage, f"RNA depth vs ADT depth ({sl})",
                     "Per-cell sequencing depth of the two modalities; a weak relationship is normal for CITE-seq.")
        # per-antibody distributions
        fig = _antibody_grid(counts, table, log1p=True, title=f"Raw ADT counts per antibody ({sl})")
        reg.save(fig, f"protein_antibody_distributions_{stage}", SECTION_QC_PROTEIN, stage, f"Per-antibody raw count distributions ({sl})",
                 "log1p raw counts per antibody; isotype controls in red, background-dominated antibodies marked (*).")
        fig = _antibody_summary(table, title=f"Antibody summary ({sl})")
        reg.save(fig, f"protein_antibody_summary_{stage}", SECTION_QC_PROTEIN, stage, f"Antibody summary ({sl})",
                 "Median raw counts per antibody (isotype controls red, hatched = background-dominated) and the fraction of cells above the matched isotype control.")
    if prot is not None:
        fig = _antibody_grid(prot, table, log1p=False, title=f"Normalized protein values ({sl}){' - ' + primary_desc if primary_desc else ''}", color=BLUE)
        reg.save(fig, f"protein_normalized_distributions_{stage}", SECTION_QC_PROTEIN, stage, f"Normalized protein value distributions ({sl})",
                 f"Distribution of obsm['protein'] per antibody ({primary_desc or 'primary normalized representation'}).")


def _antibody_grid(mat: pd.DataFrame, table: Optional[pd.DataFrame], log1p: bool, title: str, color: str = ORANGE):
    cols = list(mat.columns)
    iso = set(table.index[table["is_isotype"].astype(bool)]) if table is not None and "is_isotype" in table else set()
    bg = set(table.index[table["background_dominated"].fillna(False).astype(bool)]) if table is not None and "background_dominated" in table else set()
    ncol = 6
    nrow = int(np.ceil(len(cols) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, c in zip(axes, cols):
        v = mat[c].to_numpy(float)
        v = v[np.isfinite(v)]
        ax.hist(np.log1p(v) if log1p else v, bins=50, color=RED if c in iso else color)
        ax.set_title(f"{c}{' *' if c in bg else ''}", fontsize=8, color=RED if c in iso else "black")
        ax.tick_params(labelsize=6)
    for ax in axes[len(cols):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig


def _antibody_summary(table: Optional[pd.DataFrame], title: str):
    if table is None or "counts_median" not in table.columns:
        return None
    t = table.sort_values("counts_median", ascending=False)
    iso = t["is_isotype"].astype(bool).to_numpy() if "is_isotype" in t else np.zeros(len(t), bool)
    bg = t["background_dominated"].fillna(False).astype(bool).to_numpy() if "background_dominated" in t else np.zeros(len(t), bool)
    fig, axes = plt.subplots(2, 1, figsize=(max(5, 0.35 * len(t)), 5.4), sharex=True)
    axes[0].bar(range(len(t)), t["counts_median"], color=[RED if i else BLUE for i in iso], hatch=["//" if b else "" for b in bg])
    axes[0].set_ylabel("median raw counts")
    axes[0].set_yscale("symlog")
    if "frac_cells_above_isotype" in t:
        axes[1].bar(range(len(t)), t["frac_cells_above_isotype"].fillna(0), color=[RED if i else GREY for i in iso])
        axes[1].axhline(0.5, color=RED, ls="--", lw=1)
        axes[1].set_ylabel("frac. cells > isotype")
    axes[1].set_xticks(range(len(t)))
    axes[1].set_xticklabels(t.index, rotation=90, fontsize=7)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Perturbation QC (post-filter)
# ---------------------------------------------------------------------------


def perturbation_qc_figures(obs: pd.DataFrame, pc: pd.DataFrame, gc: pd.DataFrame, tc: pd.DataFrame, cond_target: Optional[pd.DataFrame], cfg: Config, reg: FigureRegistry) -> None:
    S, st = SECTION_QC_PERTURBATION, STAGE_POST
    # class composition
    vc = obs["perturbation_class"].astype(str).value_counts()
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    ax.bar(range(len(vc)), vc.to_numpy(), color=BLUE)
    ax.set_xticks(range(len(vc)))
    ax.set_xticklabels(vc.index, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("cells")
    for i, v in enumerate(vc.to_numpy()):
        ax.text(i, v, f"{100 * v / len(obs):.1f}%", ha="center", va="bottom", fontsize=7)
    ax.set_title("Perturbation class composition", fontsize=11)
    fig.tight_layout()
    reg.save(fig, "perturbation_class_composition", S, st, "Perturbation class composition", "Cells per perturbation class after filtering (single/multi-guide, control, unassigned).")
    if "condition" in obs.columns:
        reg.save(_class_by_condition(obs), "perturbation_class_by_condition", S, st, "Perturbation class by condition", "Fraction of cells in each perturbation class within each condition.")
    # guides per cell
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ng = obs["n_guides"].to_numpy(int)
    ax.hist(ng, bins=np.arange(-0.5, ng.max() + 1.5, 1), color=BLUE)
    ax.set_xlabel("guides per cell")
    ax.set_ylabel("cells")
    ax.set_title(f"Guides per cell (mean {ng.mean():.2f}; {100 * (ng == 0).mean():.1f}% unassigned, {100 * (ng > 1).mean():.1f}% multi-guide)", fontsize=9)
    fig.tight_layout()
    reg.save(fig, "guides_per_cell", S, st, "Guides per cell", "Number of assigned guides per cell.")
    reg.save(_perturbation_bar(pc), "cells_per_perturbation", S, st, "Cells per perturbation", "Top perturbation labels by cell count (control classes red, multi/unassigned grey).")
    reg.save(_coverage_hist(gc, tc), "coverage_histograms", S, st, "Guide and target coverage", "Single-guide cells per targeting guide and per targeting target; dashed lines are the low-coverage flag thresholds (nothing is removed).")
    if cond_target is not None:
        reg.save(_condition_heatmap(cond_target), "condition_target_coverage", S, st, "Condition x target coverage", "Single-guide cells per target and condition (log10); x marks pairs below the coverage threshold.")
    # guide-count diagnostics only when counts exist
    if "guide_total_counts" in obs.columns:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.2))
        _hist(axes[0], obs["guide_total_counts"].to_numpy(float), "guide UMIs per cell", True, [(cfg.perturbation.assignment.min_umi if hasattr(cfg.perturbation, "assignment") else None, "min_umi")])
        _hist(axes[1], obs["n_guides_detected"].to_numpy(float), "guides detected per cell")
        if "guide_top_count" in obs and "guide_second_count" in obs:
            axes[2].scatter(obs["guide_top_count"], obs["guide_second_count"], s=2, alpha=0.4, rasterized=True)
            axes[2].set_xscale("symlog"); axes[2].set_yscale("symlog")
            axes[2].set_xlabel("top guide UMIs"); axes[2].set_ylabel("second guide UMIs")
        fig.suptitle("Guide count diagnostics", fontsize=11)
        fig.tight_layout()
        reg.save(fig, "guide_count_diagnostics", S, st, "Guide count diagnostics", "Guide UMI depth, guides detected per cell and top-vs-second guide counts.")


def _class_by_condition(obs: pd.DataFrame):
    ct = pd.crosstab(obs["condition"], obs["perturbation_class"])
    frac = ct.div(ct.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(max(4, 0.9 * len(ct)) + 2.5, 3.4))
    bottom = np.zeros(len(frac))
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(frac.columns):
        ax.bar(range(len(frac)), frac[c], bottom=bottom, color=cmap(i % 10), label=c)
        bottom += frac[c].to_numpy()
    ax.set_xticks(range(len(frac)))
    ax.set_xticklabels(frac.index, rotation=30, ha="right")
    ax.set_ylabel("fraction of cells")
    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    fig.tight_layout()
    return fig


def _perturbation_bar(table: pd.DataFrame, top: int = 60):
    if table is None or table.empty:
        return None
    t = table.sort_values("n_cells", ascending=False).head(top)
    colors = {"single_targeting": BLUE, "multi_targeting": GREY, "multi_control": GREY, "mixed_control_targeting": GREY, "ambiguous": "#BBBBBB", "unassigned": "#BBBBBB"}
    cols = [colors.get(c, RED) for c in t["perturbation_class"].astype(str)]
    fig, ax = plt.subplots(figsize=(max(6, 0.22 * len(t)), 3.6))
    ax.bar(range(len(t)), t["n_cells"], color=cols)
    ax.set_xticks(range(len(t)))
    ax.set_xticklabels(t.index, rotation=90, fontsize=6)
    ax.set_ylabel("cells")
    ax.set_title(f"Cells per perturbation (top {len(t)})", fontsize=11)
    fig.tight_layout()
    return fig


def _coverage_hist(guides: pd.DataFrame, targets: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.4))
    for ax, (df, name) in zip(axes, ((guides, "guide"), (targets, "target"))):
        d = df.loc[df["class"] == "targeting", "n_cells_single_guide"].to_numpy(float)
        if d.size == 0:
            ax.axis("off")
            continue
        ax.hist(d, bins=50, color=BLUE)
        thr = float(df["low_coverage_threshold"].iloc[0])
        ax.axvline(thr, color=RED, ls="--", label=f"min_cells_per_{name}={int(thr)}")
        ax.set_xlabel(f"single-guide cells per {name}")
        ax.set_ylabel(f"{name}s")
        ax.set_title(f"{len(d)} targeting {name}s; {int((d < thr).sum())} flagged low coverage", fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


def _condition_heatmap(cond_target: pd.DataFrame, max_rows: int = 300):
    key, cond = cond_target.columns[0], cond_target.columns[1]
    t = cond_target.loc[cond_target["class"] == "targeting"]
    if t.empty:
        return None
    mat = t.pivot(index=key, columns=cond, values="n_cells").fillna(0)
    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index].iloc[:max_rows]
    thr = float(t["low_coverage_threshold"].iloc[0])
    fig, ax = plt.subplots(figsize=(1.1 * mat.shape[1] + 2.5, max(3, 0.09 * len(mat) + 1)))
    im = ax.imshow(np.log10(mat.to_numpy() + 1), aspect="auto", cmap="viridis")
    ys, xs = np.where(mat.to_numpy() < thr)
    ax.scatter(xs, ys, marker="x", color=RED, s=6, linewidths=0.6)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_ylabel(f"{len(mat)} targeting {key}s (sorted by total cells)")
    fig.colorbar(im, ax=ax, label="log10(cells + 1)")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Representations
# ---------------------------------------------------------------------------


def _scree(vr: np.ndarray, title: str):
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    ax.bar(range(1, len(vr) + 1), 100 * vr, color=BLUE)
    ax.plot(range(1, len(vr) + 1), 100 * np.cumsum(vr), color=RED, marker=".", ms=3, label="cumulative")
    ax.set_xlabel("PC")
    ax.set_ylabel("% variance")
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    return fig


def _pc_depth(Z: np.ndarray, covs: Dict[str, np.ndarray], title: str, n_pcs: int = 10):
    n = min(n_pcs, Z.shape[1])
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    w = 0.8 / max(1, len(covs))
    for j, (name, y) in enumerate(covs.items()):
        ok = np.isfinite(y) & np.isfinite(Z).all(axis=1)
        r = [abs(np.corrcoef(Z[ok, i], y[ok])[0, 1]) if ok.sum() > 2 else np.nan for i in range(n)]
        ax.bar(np.arange(1, n + 1) + (j - (len(covs) - 1) / 2) * w, r, width=w, label=name)
    ax.set_xlabel("PC")
    ax.set_ylabel("|Pearson r|")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title, fontsize=11)
    fig.tight_layout()
    return fig


def embedding_figure(coords: np.ndarray, obs: pd.DataFrame, color: str, title: str, seed: int = 0):
    if color not in obs.columns:
        return None
    ok = np.isfinite(coords).all(axis=1)
    idx = np.where(ok)[0]
    rng = np.random.default_rng(seed)
    if idx.size > MAX_POINTS:
        idx = np.sort(rng.choice(idx, MAX_POINTS, replace=False))
    xy = coords[idx]
    vals = obs[color].iloc[idx]
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    num = pd.to_numeric(vals, errors="coerce")
    if vals.dtype.name in ("category", "object", "bool") or num.isna().mean() > 0.5:
        cats = vals.astype(str)
        uniq = list(pd.unique(cats))
        if len(uniq) > MAX_CATEGORIES:
            top = cats.value_counts().index[: MAX_CATEGORIES - 1].tolist()
            cats = cats.where(cats.isin(top), "other")
            uniq = top + ["other"]
        cmap = plt.get_cmap("tab20" if len(uniq) > 10 else "tab10")
        order = rng.permutation(len(idx))
        for k, u in enumerate(uniq):
            m = (cats.to_numpy() == u)[order]
            ax.scatter(xy[order][m, 0], xy[order][m, 1], s=1.5, color=cmap(k % cmap.N), label=f"{u} ({m.sum()})", alpha=0.6, rasterized=True)
        ax.legend(markerscale=6, fontsize=6, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    else:
        v = num.to_numpy(dtype=float)
        sca = ax.scatter(xy[:, 0], xy[:, 1], s=1.5, c=v, cmap="viridis", alpha=0.7, rasterized=True)
        fig.colorbar(sca, ax=ax, label=color)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{title}: {color}", fontsize=11)
    return fig


def protein_loadings_figure(loadings: pd.DataFrame, variance_ratio, n_pcs: int = 4):
    if loadings is None or loadings.empty:
        return None
    n = min(n_pcs, loadings.shape[1])
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 0.22 * len(loadings) + 1.5), sharey=True)
    order = loadings.iloc[:, 0].sort_values().index
    for i, ax in enumerate(np.atleast_1d(axes)):
        v = loadings.loc[order].iloc[:, i]
        ax.barh(range(len(v)), v, color=np.where(v > 0, BLUE, RED))
        ax.set_yticks(range(len(v)))
        ax.set_yticklabels(order, fontsize=7)
        ax.set_title(f"PC{i+1} ({100*float(variance_ratio[i]):.1f}%)", fontsize=9)
        ax.axvline(0, color="k", lw=0.5)
    fig.suptitle("Protein PCA loadings")
    fig.tight_layout()
    return fig


def representation_figures(adata, cfg: Config, reg: FigureRegistry, color_by: Sequence[str], provided_key: Optional[str], diagnostics: Dict[str, Any]) -> None:
    obs = adata.obs
    seed = cfg.compute.seed
    # ---- RNA ----
    if "X_pca" in adata.obsm and "pca" in adata.uns:
        vr = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
        reg.save(_scree(vr, f"RNA PCA variance ({100*vr.sum():.1f}% in {len(vr)} PCs)"), "rna_pca_variance", SECTION_REP_RNA, STAGE_REP, "RNA PCA variance ratio", "Variance explained per principal component of the scaled HVG matrix; red = cumulative.")
        covs = {}
        if "total_counts" in obs: covs["log total RNA counts"] = np.log1p(obs["total_counts"].to_numpy(float))
        if "pct_counts_mt" in obs: covs["% mitochondrial"] = obs["pct_counts_mt"].to_numpy(float)
        if "n_genes_by_counts" in obs: covs["genes detected"] = obs["n_genes_by_counts"].to_numpy(float)
        reg.save(_pc_depth(np.asarray(adata.obsm["X_pca"]), covs, "RNA PCs vs depth / QC covariates"), "rna_pc_vs_depth", SECTION_REP_RNA, STAGE_REP, "RNA PCs vs depth", "Absolute correlation of the leading RNA PCs with library size, mitochondrial fraction and genes detected. A depth-associated PC1 is expected for log-normalized data and is reported, not corrected.")
    if "X_umap_rna" in adata.obsm:
        for c in color_by:
            reg.save(embedding_figure(np.asarray(adata.obsm["X_umap_rna"]), obs, c, "RNA UMAP", seed), f"rna_umap_{c}", SECTION_REP_RNA, STAGE_REP, f"RNA UMAP coloured by {c}", "UMAP of the RNA neighbour graph (HVG -> scaled PCA -> 15-NN).")
    if provided_key and provided_key in adata.obsm:
        for c in color_by[:2]:
            reg.save(embedding_figure(np.asarray(adata.obsm[provided_key]), obs, c, "Provided embedding", seed), f"provided_embedding_{c}", SECTION_REP_RNA, STAGE_REP, f"Provided embedding coloured by {c}", f"Precomputed embedding from the input (obsm['{provided_key}']), for comparison only.")
    # ---- protein ----
    if "X_pca_protein" in adata.obsm and "pca_protein" in adata.uns:
        vr = np.asarray(adata.uns["pca_protein"]["variance_ratio"], dtype=float)
        reg.save(_scree(vr, f"Protein PCA variance ({100*vr.sum():.1f}% in {len(vr)} PCs)"), "protein_pca_variance", SECTION_REP_PROTEIN, STAGE_REP, "Protein PCA variance ratio", "Variance explained per protein PC (targeting antibodies only, scaled).")
        lo = adata.uns["pca_protein"].get("loadings")
        if isinstance(lo, pd.DataFrame):
            reg.save(protein_loadings_figure(lo, vr), "protein_pca_loadings", SECTION_REP_PROTEIN, STAGE_REP, "Protein PCA loadings", "Antibody loadings of the first protein PCs; an all-positive PC1 is the ADT depth axis.")
        covs = {}
        if "protein_total_counts" in obs: covs["log total ADT counts"] = np.log1p(obs["protein_total_counts"].to_numpy(float))
        if "protein_pct_isotype" in obs: covs["% isotype"] = obs["protein_pct_isotype"].to_numpy(float)
        if "total_counts" in obs: covs["log total RNA counts"] = np.log1p(obs["total_counts"].to_numpy(float))
        if covs:
            reg.save(_pc_depth(np.asarray(adata.obsm["X_pca_protein"]), covs, "Protein PCs vs ADT depth / isotype"), "protein_pc_vs_depth", SECTION_REP_PROTEIN, STAGE_REP, "Protein PCs vs depth", "Absolute correlation of protein PCs with ADT depth, isotype fraction and RNA depth (documented CITE-seq depth effect).")
    if "X_umap_protein" in adata.obsm:
        for c in list(color_by) + [c for c in ("protein_total_counts",) if c in obs and c not in color_by]:
            reg.save(embedding_figure(np.asarray(adata.obsm["X_umap_protein"]), obs, c, "Protein UMAP", seed), f"protein_umap_{c}", SECTION_REP_PROTEIN, STAGE_REP, f"Protein UMAP coloured by {c}", "UMAP of the protein neighbour graph; banding reflects discrete low counts and is diagnostic only.")
    # ---- multimodal ----
    d = diagnostics or {}
    if "rna_protein_knn_overlap_mean" in d:
        fig, ax = plt.subplots(figsize=(4.6, 3.2))
        names, vals = ["RNA vs protein"], [d["rna_protein_knn_overlap_mean"]]
        if "multimodal_knn_overlap_with_rna" in d:
            names += ["multimodal vs RNA", "multimodal vs protein"]; vals += [d["multimodal_knn_overlap_with_rna"], d["multimodal_knn_overlap_with_protein"]]
        ax.bar(range(len(vals)), vals, color=BLUE)
        ax.axhline(d.get("rna_protein_knn_overlap_expected_random", 0), color=RED, ls="--", lw=1, label="random expectation")
        ax.set_xticks(range(len(vals))); ax.set_xticklabels(names, fontsize=8)
        ax.set_ylabel(f"shared {d.get('k', 15)}-NN fraction")
        ax.legend(fontsize=7, frameon=False)
        ax.set_title(f"Neighbourhood agreement ({d.get('sampled_cells', 0):,} sampled cells)", fontsize=10)
        fig.tight_layout()
        reg.save(fig, "knn_overlap", SECTION_MULTIMODAL, STAGE_REP, "RNA / protein neighbourhood agreement", "Mean fraction of shared nearest neighbours between modality graphs; near the random line means the modalities organise cells independently.")
    if "X_umap_multimodal" in adata.obsm:
        for c in color_by:
            reg.save(embedding_figure(np.asarray(adata.obsm["X_umap_multimodal"]), obs, c, "Multimodal UMAP", seed), f"multimodal_umap_{c}", SECTION_MULTIMODAL, STAGE_REP, f"Multimodal UMAP coloured by {c}", "UMAP of the concatenated RNA + protein PC space (multimodal.enabled).")
