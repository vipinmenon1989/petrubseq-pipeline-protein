"""Stage E figures (section ``perturbation_effects``), all through the FigureRegistry.

Panel selection is deterministic and bounded: per-target panels show the
``top_n_report`` targets by the relevant statistic (ties broken by name);
scatter panels show at most 6 pairs chosen by FDR then |rho|. Every value is
also in the tables under ``tables/perturbation_effects/``.
"""

from __future__ import annotations

from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .figures import FigureRegistry  # noqa: E402

SECTION = "perturbation_effects"
ST_PS, ST_LOCH, ST_PROG, ST_PROT, ST_CONC = "ps", "lochness", "gene_programs", "protein_effects", "concordance"
BLUE, RED, GREY = "#4C72B0", "#C44E52", "#8C8C8C"


def _trend(ax, x: np.ndarray, y: np.ndarray, color: str = RED) -> None:
    """Binned-median trend line (visual aid only; the statistic is the reported Spearman rho)."""
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
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(M.index, fontsize=6 if M.shape[0] > 30 else 7)
    if M.shape[1] <= 40:
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels(M.columns, fontsize=7, rotation=90)
    else:
        ax.set_xticks([])
    if marks is not None:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if bool(marks.iat[i, j]):
                    ax.text(j, i, "*", ha="center", va="center", fontsize=8)
    cb = plt.colorbar(im, ax=ax, fraction=0.04)
    cb.set_label(label, fontsize=8)


def perturbation_effect_figures(res, adata, cfg, reg: FigureRegistry) -> None:
    top_n = cfg.analysis.perturbation_effects.top_n_report
    alpha = cfg.analysis.perturbation_effects.concordance.fdr_alpha
    # ---------------------------------------------------------------- PS
    ps = res.ps
    if ps is not None and not ps.summary.empty:
        s = ps.summary.sort_values(["median_ps", "target"], ascending=[False, True])
        fig, ax = plt.subplots(figsize=(max(5, 0.28 * len(s)), 3.4))
        ax.bar(range(len(s)), s["median_ps"], color=BLUE)
        ax.errorbar(range(len(s)), s["median_ps"], yerr=[s["median_ps"] - s["q25_ps"], s["q75_ps"] - s["median_ps"]], fmt="none", ecolor=GREY, lw=0.8)
        ax.set_xticks(range(len(s))); ax.set_xticklabels([f"{t} ({n})" for t, n in zip(s["target"], s["n_perturbed_cells"])], rotation=90, fontsize=7)
        ax.set_ylabel("median PS (IQR)"); ax.axhline(cfg.analysis.perturbation_effects.ps.ps_threshold, color=RED, ls="--", lw=0.8)
        fig.tight_layout()
        reg.save(fig, "ps_target_median", SECTION, ST_PS, "Median PS per target", "Median perturbation-response score of each target's single-guide cells (IQR bars; cell count in brackets); dashed line = quadrant threshold. PS is max-normalized within each target; for cross-target comparison of response strength use auc_vs_control / d_vs_control in ps_targets.csv (in-sample separation of perturbed from control cells).")
        show = s.head(top_n)["target"].tolist()
        data = [ps.own[(adata.obs["target"].astype(str) == t).to_numpy() & ps.own.notna().to_numpy()].to_numpy() for t in show]
        ctrl = [ps.scores[t][adata.obs["perturbation_class"].astype(str).eq("single_control").to_numpy()].dropna().to_numpy() for t in show]
        fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(show)), 3.4))
        pos = np.arange(len(show))
        ax.violinplot([d if len(d) else [0] for d in data], positions=pos - 0.18, widths=0.34, showmedians=True)
        vp = ax.violinplot([c if len(c) else [0] for c in ctrl], positions=pos + 0.18, widths=0.34, showmedians=True)
        for b in vp["bodies"]:
            b.set_facecolor(GREY)
        ax.set_xticks(pos); ax.set_xticklabels(show, rotation=90, fontsize=7); ax.set_ylabel("PS")
        fig.tight_layout()
        reg.save(fig, "ps_distributions", SECTION, ST_PS, "PS distributions (perturbed vs controls)", f"Per-cell PS of the {len(show)} targets with the highest median PS (blue: perturbed cells; grey: control cells scored against the same signature).")
    # ---------------------------------------------------------------- lochNESS
    lo = res.lochness
    if lo is not None and not lo.summary.empty:
        s = lo.summary.sort_values(["mean_lochness_in_own_cells", "target"], ascending=[False, True])
        fig, ax = plt.subplots(figsize=(max(5, 0.28 * len(s)), 3.4))
        sig = s["fdr"] < alpha if "fdr" in s.columns else np.zeros(len(s), bool)
        ax.bar(range(len(s)), s["mean_lochness_in_own_cells"], color=np.where(sig, RED, BLUE))
        if "null_mean" in s.columns:
            ax.errorbar(range(len(s)), s["null_mean"], yerr=2 * s["null_sd"], fmt="_", color="black", lw=0.8, label="permutation null (mean ± 2 SD)")
            ax.legend(fontsize=7)
        ax.axhline(0, color=GREY, lw=0.8)
        ax.set_xticks(range(len(s))); ax.set_xticklabels([f"{t} ({n})" for t, n in zip(s["target"], s["n_cells"])], rotation=90, fontsize=7)
        ax.set_ylabel("mean lochNESS in own cells")
        fig.tight_layout()
        reg.save(fig, "lochness_targets", SECTION, ST_LOCH, "lochNESS self-enrichment per target", f"Mean lochNESS of each target's own cells (0 = chance; > 0 = the target's cells cluster together in PCA space); red = permutation FDR < {alpha}. k = {lo.info.get('k_used')} neighbours in {lo.info.get('use_rep')}.")
        if "X_umap_rna" in adata.obsm:
            show = s.head(min(6, top_n))["target"].tolist()
            n = len(show)
            fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 2.6), squeeze=False)
            U = adata.obsm["X_umap_rna"]
            for ax, t in zip(axes[0], show):
                v = lo.scores[t].to_numpy()
                lim = np.nanmax(np.abs(v))
                ax.scatter(U[:, 0], U[:, 1], c=v, cmap="RdBu_r", vmin=-lim, vmax=lim, s=2, rasterized=True)
                own = (adata.obs["target"].astype(str) == t).to_numpy() & adata.obs["perturbation_class"].astype(str).eq("single_targeting").to_numpy()
                ax.scatter(U[own, 0], U[own, 1], facecolors="none", edgecolors="black", s=5, lw=0.3)
                ax.set_title(t, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
            fig.tight_layout()
            reg.save(fig, "lochness_umap", SECTION, ST_LOCH, "lochNESS maps", "lochNESS of the strongest targets on the RNA UMAP (colour; computed in PCA space, UMAP only for display); black circles = the target's own cells.")
    # ---------------------------------------------------------------- programs / modules
    mo = res.modules
    if mo is not None and not mo.empty:
        ppe = mo.perturbation_program_effects.pivot(index="target", columns="program", values="mean_log2fc").reindex(index=mo.target_order, columns=mo.program_labels)
        mod = mo.perturbation_modules.set_index("target")["module"]
        ppe.index = [f"{t} [{mod[t]}]" for t in ppe.index]
        fig, ax = plt.subplots(figsize=(1.2 + 0.6 * ppe.shape[1], 1.2 + 0.2 * ppe.shape[0]))
        _heatmap(ax, ppe, label="mean log2FC of program genes")
        ax.set_title("perturbation x gene program", fontsize=9)
        fig.tight_layout()
        reg.save(fig, "perturbation_program_heatmap", SECTION, ST_PROG, "Perturbation x gene-program effects", "Mean log2FC (vs non-targeting controls) of each gene program's genes under each perturbation; rows ordered by the perturbation-module dendrogram, module in brackets. Programs (P) are gene clusters, modules (M) are perturbation clusters; neither carries a biological label.")
        E = mo.effect.reindex(index=mo.target_order, columns=mo.gene_order)
        prog = mo.gene_programs.set_index("gene")["program"].reindex(mo.gene_order)
        fig, ax = plt.subplots(figsize=(9, 1.2 + 0.2 * E.shape[0]))
        _heatmap(ax, E.clip(-3, 3), label="log2FC (clipped ±3)")
        b = np.flatnonzero(prog.to_numpy()[1:] != prog.to_numpy()[:-1]) + 0.5
        for x in b:
            ax.axvline(x, color="black", lw=0.6)
        for p in mo.program_labels:
            pos = np.flatnonzero(prog.to_numpy() == p)
            ax.text(pos.mean(), -0.8, f"{p} ({len(pos)})", ha="center", fontsize=7)
        fig.tight_layout()
        reg.save(fig, "effect_matrix_programs", SECTION, ST_PROG, "Response-gene effect matrix", "Perturbation x response-gene log2FC ordered by both dendrograms; vertical lines separate gene programs (size in brackets). Gene membership and effects: gene_programs.csv.")
        mp = mo.module_program.astype(float)
        fig, ax = plt.subplots(figsize=(1.5 + 0.6 * mp.shape[1], 1.2 + 0.35 * mp.shape[0]))
        _heatmap(ax, mp, label="mean log2FC")
        fig.tight_layout()
        reg.save(fig, "module_program_strength", SECTION, ST_PROG, "Module x program strength", "Mean log2FC of a program's genes under the perturbations of a module (signed link between perturbation modules and gene programs).")
    # ---------------------------------------------------------------- protein
    pr = res.protein
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
    # ---------------------------------------------------------------- concordance
    co = res.concordance
    if co is not None:
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
                    ax.scatter(x, y, s=5, alpha=0.6, color=BLUE)
                    _trend(ax, x, y)
                    ax.set_title(f"{r['target']}: {r['protein']}\nrho={r['spearman_rho']:+.2f}, FDR={r['fdr']:.2g}, n={r['n_cells']}", fontsize=7)
                    ax.set_xlabel("PS", fontsize=8); ax.set_ylabel(f"{r['protein']} (CLR)", fontsize=8)
                fig.tight_layout()
                reg.save(fig, "ps_vs_protein", SECTION, ST_CONC, "PS vs protein (within target)", "Per-cell PS vs protein value inside one target's perturbed cells, for the pairs with the lowest FDR; red = binned-median trend (visual aid). The statistic is the within-target Spearman rho in ps_protein_association.csv.")
        ppc = co.program_protein_cells
        if ppc is not None and not ppc.empty and mo is not None and mo.program_activity is not None:
            w = ppc.copy()
            w["abs"] = w["spearman_rho"].abs()
            w = w.sort_values(["fdr", "abs"], ascending=[True, False]).head(4)
            P = adata.obsm[cfg.analysis.perturbation_effects.protein.representation]
            cells = adata.obs["perturbation_class"].astype(str).isin(["single_targeting", "single_control"]).to_numpy() & P.notna().all(axis=1).to_numpy()
            n = len(w)
            fig, axes = plt.subplots(1, n, figsize=(2.8 * n, 2.8), squeeze=False)
            for ax, (_, r) in zip(axes[0], w.iterrows()):
                x = mo.program_activity.loc[cells, r["gene_program"]].to_numpy(float)
                y = P.loc[cells, r["protein"]].to_numpy(float)
                ax.hexbin(x, y, gridsize=30, cmap="Blues", mincnt=1)
                _trend(ax, x, y)
                ax.set_title(f"{r['gene_program']} vs {r['protein']}\nrho={r['spearman_rho']:+.2f}, FDR={r['fdr']:.2g}", fontsize=7)
                ax.set_xlabel(f"{r['gene_program']} activity", fontsize=8); ax.set_ylabel(f"{r['protein']} (CLR)", fontsize=8)
            fig.tight_layout()
            reg.save(fig, "program_activity_vs_protein", SECTION, ST_CONC, "Gene-program activity vs protein", "Cell-level program activity (mean z-scored expression of the program's genes) vs protein value over single-guide and control cells, strongest pairs by FDR; association only, not mediation.")
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
