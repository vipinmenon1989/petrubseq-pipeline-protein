"""Perturbation-distance figures: the reference ``plots.py`` distance / distance-space / atlas
figures (``perturbation_distance_ranking``, ``perturbation_atlas``, ``ps_vs_distance_map``,
``perturbation_phenotype_space``, ``module_concordance``) drawn with the same data
transformations (column z-scores, KD strength = -log2FC, FDR stars, top-N by energy
distance, PCoA1/2 coloured by module and by distance, crosstab with ARI / NMI), plus the
extension figures: the pairwise distance matrix ordered by the module dendrogram, and the
distance <-> protein and phenotype-space <-> protein views."""

from __future__ import annotations

import logging
from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ._refplot import GREY, despine  # noqa: E402
from .figures import FigureRegistry  # noqa: E402

logger = logging.getLogger("petrubseq_protein")
SECTION = "perturbation_distance"
ST_DIST, ST_SPACE, ST_PROT = "distance", "distance_space", "distance_protein"
SIG, NSIG = "#dd6b20", "#a0aec0"


def _stars(fdr: float) -> str:
    if pd.isna(fdr):
        return ""
    return "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else ""


def _palette(cats: List[str], name: str = "tab10") -> dict:
    cmap = plt.get_cmap(name)
    return {c: cmap(i % cmap.N) for i, c in enumerate(cats)}


# ---------------------------------------------------------------------------------------- reference figures


def plot_distance_ranking(dist, reg: FigureRegistry, cfg) -> None:
    if dist is None or dist.empty:
        return
    tbl = dist.table.sort_values("energy_distance", ascending=False).reset_index(drop=True)
    sub = tbl.iloc[: min(40, len(tbl))]
    alpha = cfg.analysis.perturbation_effects.distance.fdr_threshold
    fig, ax = plt.subplots(figsize=(max(6.0, 0.25 * len(sub) + 1.5), 4.5))
    ax.bar(range(len(sub)), sub["energy_distance"], color=[SIG if s else NSIG for s in sub["significant"].astype(bool)], edgecolor="none", width=0.8)
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub["target"], rotation=90, fontsize=8)
    ax.set_ylabel("Energy distance vs control", fontsize=10, fontweight="bold")
    ax.set_title("Perturbation distance vs control ranking", fontsize=11, fontweight="bold")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=SIG, label=f"FDR < {alpha}"), plt.Rectangle((0, 0), 1, 1, color=NSIG, label="not significant")], fontsize=8, loc="upper right", frameon=False)
    despine(ax)
    fig.tight_layout()
    reg.save(fig, "perturbation_distance_ranking", SECTION, ST_DIST, "Perturbation distance vs control ranking",
             f"Energy distance between each target's cells and the {dist.control_used} control cells in {dist.representation}, ranked; orange = DistanceTest permutation FDR < {alpha} ({dist.info.get('n_permutations')} permutations).")


def plot_perturbation_atlas(master: pd.DataFrame, reg: FigureRegistry, cfg) -> None:
    """Reference ``plot_perturbation_atlas``: column-z-scored heatmap of efficacy (-log2FC), PS median, lochNESS mean, energy distance; FDR stars; module side strips."""
    if master is None or master.empty or "target" not in master.columns:
        return
    df = master.copy()
    sort_col = "energy_distance" if "energy_distance" in df.columns and df["energy_distance"].notna().any() else ("ps_median" if "ps_median" in df.columns else "target")
    df = df.sort_values(sort_col, ascending=(sort_col == "target"), na_position="last").reset_index(drop=True)
    max_n = cfg.analysis.perturbation_effects.master_table.atlas_top_n
    if len(df) > max_n:
        df = df.iloc[:max_n].copy()
    specs = [("target_log2fc", "Efficacy\n(KD strength)", True, "target_fdr"), ("ps_median", "Penetrance\n(PS median)", False, None), ("lochness_mean", "Topology\n(lochNESS)", False, None), ("energy_distance", "Phenotype\n(energy dist.)", False, "distance_fdr")]
    cols, labels, fdr_cols, data = [], [], [], []
    for col, label, invert, fdr_col in specs:
        if col in df.columns and df[col].notna().any():
            vals = df[col].to_numpy(dtype=float, na_value=np.nan)
            data.append(-vals if invert else vals)
            cols.append(col)
            labels.append(label)
            fdr_cols.append(fdr_col if fdr_col in df.columns else None)
    if len(cols) < 2:
        return
    M = np.column_stack(data)
    Z = np.zeros_like(M)
    for j in range(M.shape[1]):
        v = M[:, j]
        ok = ~np.isnan(v)
        if ok.sum() > 1:
            sd = np.std(v[ok])
            Z[ok, j] = (v[ok] - np.mean(v[ok])) / (sd if sd > 1e-8 else 1.0)
    targets = df["target"].tolist()
    n = len(targets)
    has_co = "cofunctional_module" in df.columns and df["cofunctional_module"].notna().any()
    has_ph = "phenotype_module" in df.columns and df["phenotype_module"].notna().any()
    n_side = int(has_co) + int(has_ph)
    fig = plt.figure(figsize=(8.0 + 1.2 * n_side, max(6.0, 0.28 * n + 2.0)))
    gs = fig.add_gridspec(1, 1 + n_side + 1, width_ratios=[0.4] * n_side + [4.0, 0.2], wspace=0.15)
    k = 0
    for present, col, label, cmap_name in ((has_co, "cofunctional_module", "Co-func\nmodule", "tab20"), (has_ph, "phenotype_module", "Pheno\nmodule", "Set2")):
        if not present:
            continue
        ax = fig.add_subplot(gs[0, k])
        k += 1
        mods = df[col].fillna("None").astype(str).tolist()
        pal = _palette(sorted(set(mods)), cmap_name)
        ax.imshow(np.array([pal[m][:3] for m in mods])[:, None, :], aspect="auto", interpolation="nearest")
        ax.set_xticks([0])
        ax.set_xticklabels([label], rotation=90, fontsize=8)
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    ax = fig.add_subplot(gs[0, k])
    cax = fig.add_subplot(gs[0, k + 1])
    vmax = max(2.5, float(np.nanmax(np.abs(Z))))
    im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    fig.colorbar(im, cax=cax, label="column z-score")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9, fontweight="bold")
    ax.set_yticks(range(n))
    ax.set_yticklabels(targets, fontsize=8)
    for i in range(n):
        for j, fc in enumerate(fdr_cols):
            if fc:
                s = _stars(df.iloc[i][fc])
                if s:
                    ax.text(j, i, s, ha="center", va="center", color="black" if abs(Z[i, j]) < 1.2 else "white", fontsize=9, fontweight="bold")
    ax.set_title("Perturbation atlas", fontsize=12, fontweight="bold", pad=12)
    reg.save(fig, "perturbation_atlas", SECTION, ST_DIST, "Perturbation atlas",
             "Multi-dimensional overview per target (rows sorted by energy distance): efficacy (knockdown strength = -log2FC of the target's own gene), penetrance (PS median), manifold topology (own-cell mean lochNESS) and phenotype magnitude (energy distance from control). "
             "Columns are z-scored for display; stars mark FDR (* < 0.05, ** < 0.01, *** < 0.001) of the knockdown test and of the DistanceTest. Side strips: co-functional module (gene-effect similarity) and phenotype module (cell-state distance similarity).")


def plot_ps_vs_distance(master: pd.DataFrame, reg: FigureRegistry, cfg) -> None:
    """Reference ``plot_ps_vs_distance``: PS median vs energy distance, point size = lochNESS mean, colour = phenotype (or co-functional) module."""
    if master is None or master.empty or "energy_distance" not in master.columns:
        return
    y_col = "ps_median" if "ps_median" in master.columns else ("ps_responder_fraction" if "ps_responder_fraction" in master.columns else None)
    if y_col is None:
        return
    sub = master[master["energy_distance"].notna() & master[y_col].notna()].copy()
    if len(sub) < 3:
        return
    alpha = cfg.analysis.perturbation_effects.distance.fdr_threshold
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    if "lochness_mean" in sub.columns and sub["lochness_mean"].notna().any():
        l = sub["lochness_mean"].fillna(0).to_numpy(float)
        sizes = 40.0 + 160.0 * (l - l.min()) / (l.max() - l.min() + 1e-6)
    else:
        sizes = np.full(len(sub), 60.0)
    color_col = next((c for c in ("phenotype_module", "cofunctional_module") if c in sub.columns and sub[c].notna().any()), None)
    if color_col:
        cats = sub[color_col].fillna("None").astype(str)
        pal = _palette(sorted(set(cats)))
        colors = np.array([pal[c] for c in cats])
    else:
        colors = np.array(["#2b6cb0"] * len(sub))
    if "distance_significant" in sub.columns:
        is_sig = sub["distance_significant"].fillna(False).astype(bool).to_numpy()
    elif "distance_fdr" in sub.columns:
        is_sig = (sub["distance_fdr"].fillna(1.0) < alpha).to_numpy()
    else:
        is_sig = np.ones(len(sub), bool)
    x, y = sub["energy_distance"].to_numpy(float), sub[y_col].to_numpy(float)
    ax.scatter(x[is_sig], y[is_sig], s=sizes[is_sig], c=colors[is_sig], alpha=0.85, edgecolors="#1a202c", linewidths=1.2, label=f"significant (FDR < {alpha})")
    if (~is_sig).any():
        ax.scatter(x[~is_sig], y[~is_sig], s=sizes[~is_sig] * 0.7, c=colors[~is_sig], alpha=0.35, edgecolors="gray", linewidths=0.8, label="not significant")
    for _, row in sub.nlargest(min(12, len(sub)), "energy_distance").iterrows():
        ax.annotate(str(row["target"]), (row["energy_distance"], row[y_col]), xytext=(5, 5), textcoords="offset points", fontsize=8, fontweight="bold", alpha=0.9)
    ax.set_xlabel("Phenotype magnitude (energy distance from control)", fontsize=10, fontweight="bold")
    ax.set_ylabel(f"Penetrance ({y_col.replace('_', ' ')})", fontsize=10, fontweight="bold")
    ax.set_title("Perturbation penetrance vs phenotype distance map", fontsize=11, fontweight="bold")
    if color_col and len(pal) <= 10:
        ax.legend(handles=[plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=pal[c], markersize=8, label=c) for c in pal], title=color_col.replace("_", " "), loc="best", fontsize=8)
    despine(ax)
    reg.save(fig, "ps_vs_distance_map", SECTION, ST_DIST, "PS vs perturbation distance map",
             "Per-cell penetrance (PS median) against global phenotype magnitude (energy distance from control). Point size = own-cell mean lochNESS; colour = phenotype module; faded points are not significant in the DistanceTest; the 12 most distant targets are labelled.")


def plot_phenotype_space(space, master: Optional[pd.DataFrame], reg: FigureRegistry, cfg) -> None:
    """Reference ``plot_perturbation_space``: PCoA1/2 coloured by phenotype module and by energy distance."""
    if space is None or space.empty or space.coordinates.empty or "PCoA2" not in space.coordinates.columns:
        return
    co = space.coordinates.copy()
    if not space.phenotype_modules.empty:
        co = co.merge(space.phenotype_modules, on="target", how="left")
    if master is not None and not master.empty:
        co = co.merge(master[[c for c in ("target", "energy_distance", "ps_median") if c in master.columns]], on="target", how="left")
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.5))
    ax1 = axes[0]
    if "phenotype_module" in co.columns and co["phenotype_module"].notna().any():
        cats = co["phenotype_module"].fillna("None").astype(str)
        pal = _palette(sorted(set(cats), key=lambda s: (len(s), s)))
        for c in pal:
            s = co[cats == c]
            ax1.scatter(s["PCoA1"], s["PCoA2"], c=[pal[c]], label=c, s=50, alpha=0.85, edgecolors="none")
        if len(pal) <= 12:
            ax1.legend(title="phenotype module", fontsize=8, loc="best")
    else:
        ax1.scatter(co["PCoA1"], co["PCoA2"], c="#2b6cb0", s=50, alpha=0.85)
    for _, r in co.iterrows():
        ax1.annotate(str(r["target"]), (r["PCoA1"], r["PCoA2"]), fontsize=6.5, xytext=(3, 3), textcoords="offset points", alpha=0.8)
    ax1.set_xlabel("PCoA 1", fontsize=10, fontweight="bold")
    ax1.set_ylabel("PCoA 2", fontsize=10, fontweight="bold")
    ax1.set_title("Perturbation phenotype space (modules)", fontsize=11, fontweight="bold")
    despine(ax1)
    ax2 = axes[1]
    metric = next((c for c in ("energy_distance", "ps_median") if c in co.columns and co[c].notna().any()), None)
    if metric:
        sc = ax2.scatter(co["PCoA1"], co["PCoA2"], c=co[metric].to_numpy(float), cmap="viridis", s=50, alpha=0.85, edgecolors="none")
        fig.colorbar(sc, ax=ax2, label=metric.replace("_", " "))
    else:
        ax2.scatter(co["PCoA1"], co["PCoA2"], c="#4a5568", s=50, alpha=0.85)
    ax2.set_xlabel("PCoA 1", fontsize=10, fontweight="bold")
    ax2.set_ylabel("PCoA 2", fontsize=10, fontweight="bold")
    ax2.set_title(f"Phenotype space ({metric.replace('_', ' ') if metric else 'PCoA'})", fontsize=11, fontweight="bold")
    despine(ax2)
    ev = space.eigenvalues
    frac = (ev[:2] / ev.sum() * 100) if ev.size >= 2 and ev.sum() > 0 else None
    reg.save(fig, "perturbation_phenotype_space", SECTION, ST_SPACE, "Perturbation phenotype space (PCoA)",
             "Classical multidimensional scaling (PCoA) of the pairwise target x target energy distances" + (f"; PCoA 1 / 2 carry {frac[0]:.0f} % / {frac[1]:.0f} % of the positive eigenvalue mass" if frac is not None else "") +
             ". Left: phenotype modules (average-linkage clusters of the same matrix); right: energy distance from control. Targets close together produce similar cell-state distributions.")


def plot_module_concordance(master: pd.DataFrame, reg: FigureRegistry) -> Optional[dict]:
    """Reference ``plot_module_concordance``: co-functional module x phenotype module crosstab with ARI / NMI."""
    if master is None or master.empty or "cofunctional_module" not in master.columns or "phenotype_module" not in master.columns:
        return None
    ok = master["cofunctional_module"].notna() & master["phenotype_module"].notna()
    if ok.sum() < 4:
        return None
    co, ph = master.loc[ok, "cofunctional_module"].astype(str), master.loc[ok, "phenotype_module"].astype(str)
    if co.nunique() < 2 or ph.nunique() < 2:
        return None
    ct = pd.crosstab(co, ph)
    metric = {}
    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

        metric = {"ari": float(adjusted_rand_score(co, ph)), "nmi": float(normalized_mutual_info_score(co, ph)), "n_targets": int(ok.sum())}
        mstr = f" (ARI = {metric['ari']:.3f}, NMI = {metric['nmi']:.3f})"
    except Exception:  # pragma: no cover
        mstr = ""
    fig, ax = plt.subplots(figsize=(max(5.0, 0.6 * ct.shape[1] + 2.0), max(4.5, 0.5 * ct.shape[0] + 1.5)))
    im = ax.imshow(ct.to_numpy(), cmap="Blues", aspect="auto")
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = int(ct.iat[i, j])
            ax.text(j, i, str(v), ha="center", va="center", fontsize=9, color="white" if v > ct.to_numpy().max() / 2 else "black")
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels(ct.columns, fontsize=9)
    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index, fontsize=9)
    ax.set_xlabel("Phenotype modules (distance space)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Co-functional modules (gene programs)", fontsize=10, fontweight="bold")
    ax.set_title(f"Module concordance{mstr}", fontsize=11, fontweight="bold", pad=12)
    fig.colorbar(im, ax=ax, fraction=0.04, label="targets")
    reg.save(fig, "module_concordance", SECTION, ST_SPACE, "Module concordance heatmap",
             f"Cross-tabulation of the gene-effect co-functional modules against the cell-state phenotype modules{mstr}. ARI / NMI measure the agreement of the two partitions over the targets present in both.")
    return metric


# ---------------------------------------------------------------------------------------- extension figures


def plot_distance_matrix(space, reg: FigureRegistry) -> None:
    if space is None or space.empty:
        return
    M = space.distance_matrix
    order = list(M.index)
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        order = [M.index[i] for i in leaves_list(linkage(squareform(M.to_numpy(float), checks=False), method=space.linkage_method))]
    except Exception:  # pragma: no cover - ordering is cosmetic
        pass
    Mo = M.loc[order, order]
    mods = space.phenotype_modules.set_index("target")["phenotype_module"] if not space.phenotype_modules.empty else None
    n = len(order)
    fig = plt.figure(figsize=(max(6.5, 0.32 * n + 2.5), max(5.5, 0.32 * n + 2.0)))
    gs = fig.add_gridspec(1, 3 if mods is not None else 2, width_ratios=([0.25] if mods is not None else []) + [4.0, 0.2], wspace=0.08)
    k = 0
    if mods is not None:
        axm = fig.add_subplot(gs[0, 0])
        labs = [str(mods.get(t, "None")) for t in order]
        pal = _palette(sorted(set(labs), key=lambda s: (len(s), s)), "Set2")
        axm.imshow(np.array([pal[m][:3] for m in labs])[:, None, :], aspect="auto", interpolation="nearest")
        axm.set_xticks([0])
        axm.set_xticklabels(["module"], rotation=90, fontsize=8)
        axm.set_yticks(range(n))
        axm.set_yticklabels(order, fontsize=7)
        for sp in axm.spines.values():
            sp.set_visible(False)
        k = 1
    ax = fig.add_subplot(gs[0, k])
    cax = fig.add_subplot(gs[0, k + 1])
    im = ax.imshow(Mo.to_numpy(float), cmap="magma_r", aspect="auto", interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_yticks(range(n) if mods is None else [])
    if mods is None:
        ax.set_yticklabels(order, fontsize=7)
    fig.colorbar(im, cax=cax, label=f"{space.metric} between targets")
    ax.set_title("Pairwise perturbation distance matrix", fontsize=11, fontweight="bold")
    reg.save(fig, "perturbation_distance_matrix", SECTION, ST_SPACE, "Pairwise perturbation distance matrix",
             f"Target x target {space.metric} in {space.info.get('representation', 'X_pca')} (zero diagonal), rows and columns in the order of the {space.linkage_method}-linkage dendrogram that defines the phenotype modules (side strip). Extension figure; the matrix itself is the reference object.")


def plot_distance_vs_protein(mt, reg: FigureRegistry, cfg) -> None:
    if mt is None or mt.distance_protein_targets.empty:
        return
    long, summ = mt.distance_protein_targets, mt.distance_protein.set_index("protein") if not mt.distance_protein.empty else None
    proteins = list(dict.fromkeys(long["protein"]))
    ncol = min(4, len(proteins))
    nrow = int(np.ceil(len(proteins) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.3 * nrow), squeeze=False)
    for ax, p in zip(axes.ravel(), proteins):
        d = long[long["protein"] == p]
        sig_d = d["distance_significant"].astype(bool).to_numpy()
        sig_p = (d["protein_fdr"] < cfg.analysis.perturbation_effects.protein.fdr_alpha).to_numpy()
        ax.axhline(0, color="#718096", lw=0.8, ls="--")
        ax.scatter(d["energy_distance"][~sig_d], d["protein_effect"][~sig_d], s=22, color=GREY, alpha=0.7)
        ax.scatter(d["energy_distance"][sig_d], d["protein_effect"][sig_d], s=28, color=SIG, alpha=0.9, edgecolors=np.where(sig_p[sig_d], "#1a202c", "none"), linewidths=0.9)
        for _, r in d[sig_d | sig_p].iterrows():
            ax.annotate(r["target"], (r["energy_distance"], r["protein_effect"]), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
        title = p
        if summ is not None and p in summ.index:
            s = summ.loc[p]
            title += f"\nrho signed {s['rho_signed']:+.2f} (FDR {s['fdr_signed']:.2g}); |effect| {s['rho_abs']:+.2f} (FDR {s['fdr_abs']:.2g})"
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("energy distance from control", fontsize=8)
        ax.set_ylabel("protein effect (CLR difference)", fontsize=8)
        despine(ax)
    for ax in axes.ravel()[len(proteins):]:
        ax.axis("off")
    fig.suptitle("Phenotype distance vs protein effect (target level)", fontsize=10, fontweight="bold")
    fig.tight_layout()
    reg.save(fig, "distance_vs_protein", SECTION, ST_PROT, "Perturbation distance vs protein effect",
             "Per protein, one point per target: energy distance from control (x) against the target's protein effect (y, CLR mean difference vs control). Orange = DistanceTest significant; black edge = protein effect significant. Titles give the target-level Spearman for the signed effect and for |effect| with their BH-FDR (two families). Descriptive association at the target level.")


def plot_phenotype_space_protein(space, mt, reg: FigureRegistry) -> None:
    if space is None or space.empty or mt is None or mt.distance_protein_targets.empty or "PCoA2" not in space.coordinates.columns:
        return
    co = space.coordinates.set_index("target")
    long = mt.distance_protein_targets
    proteins = list(dict.fromkeys(long["protein"]))
    ncol = min(4, len(proteins))
    nrow = int(np.ceil(len(proteins) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.8 * ncol, 3.4 * nrow), squeeze=False)
    mantel = mt.phenotype_protein_mantel.set_index("protein") if not mt.phenotype_protein_mantel.empty else None
    for ax, p in zip(axes.ravel(), proteins):
        d = long[long["protein"] == p].set_index("target")
        common = [t for t in co.index if t in d.index]
        v = d.loc[common, "protein_effect"].to_numpy(float)
        lim = max(np.nanmax(np.abs(v)), 1e-9)
        sc = ax.scatter(co.loc[common, "PCoA1"], co.loc[common, "PCoA2"], c=v, cmap="RdBu_r", vmin=-lim, vmax=lim, s=48, edgecolors="#4a5568", linewidths=0.5)
        for t in common:
            ax.annotate(t, (co.loc[t, "PCoA1"], co.loc[t, "PCoA2"]), fontsize=6, xytext=(3, 3), textcoords="offset points", alpha=0.8)
        fig.colorbar(sc, ax=ax, fraction=0.05, label="protein effect")
        title = p
        if mantel is not None and p in mantel.index:
            title += f"\nMantel rho {mantel.loc[p, 'mantel_rho']:+.2f}, p {mantel.loc[p, 'p_permutation']:.3f}"
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("PCoA 1", fontsize=8)
        ax.set_ylabel("PCoA 2", fontsize=8)
        despine(ax)
    for ax in axes.ravel()[len(proteins):]:
        ax.axis("off")
    joint = ""
    if mantel is not None and "all_proteins" in mantel.index:
        joint = f" - joint Mantel rho {mantel.loc['all_proteins', 'mantel_rho']:+.2f}, p {mantel.loc['all_proteins', 'p_permutation']:.3f}"
    fig.suptitle(f"Phenotype space coloured by protein effect{joint}", fontsize=10, fontweight="bold")
    fig.tight_layout()
    reg.save(fig, "phenotype_space_protein", SECTION, ST_PROT, "Phenotype space vs protein effects",
             "The PCoA phenotype space with each target coloured by its effect on one protein. Mantel statistics: Spearman between the target x target energy distances and the target x target differences in protein effect (one protein, or the Euclidean distance over all proteins), seeded label permutations. Underpowered with a 4-antibody panel; descriptive.")


def plot_phenotype_module_protein(space, mt, reg: FigureRegistry) -> None:
    if space is None or space.empty or space.phenotype_modules.empty or mt is None or mt.phenotype_module_means.empty or mt.distance_protein_targets.empty:
        return
    means = mt.phenotype_module_means
    long = mt.distance_protein_targets
    modmap = space.phenotype_modules.set_index("target")["phenotype_module"]
    mods = sorted(means["phenotype_module"].unique(), key=lambda s: (len(str(s)), str(s)))
    proteins = list(dict.fromkeys(means["protein"]))
    kw = mt.phenotype_module_protein.set_index("protein") if not mt.phenotype_module_protein.empty else None
    pal = _palette(mods, "Set2")
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, len(proteins), figsize=(3.2 * len(proteins), 3.6), squeeze=False)
    for ax, p in zip(axes[0], proteins):
        d = long[long["protein"] == p].set_index("target")
        for i, m in enumerate(mods):
            members = [t for t in d.index if modmap.get(t) == m]
            if members:
                x = i + rng.uniform(-0.18, 0.18, len(members))
                ax.scatter(x, d.loc[members, "protein_effect"], s=22, color=pal[m], edgecolors="#4a5568", linewidths=0.4, zorder=3)
            row = means[(means["phenotype_module"] == m) & (means["protein"] == p)]
            if len(row):
                ax.hlines(float(row["mean_effect"].iloc[0]), i - 0.3, i + 0.3, color="black", lw=1.5, zorder=4)
        ax.axhline(0, color="#718096", lw=0.8, ls="--")
        ax.set_xticks(range(len(mods)))
        ax.set_xticklabels([f"{m}\n(n={int(means[(means['phenotype_module'] == m) & (means['protein'] == p)]['n_targets'].iloc[0])})" for m in mods], fontsize=7)
        title = p
        if kw is not None and p in kw.index and np.isfinite(kw.loc[p, "fdr"]):
            title += f"\nKruskal-Wallis FDR {kw.loc[p, 'fdr']:.2g}"
        ax.set_title(title, fontsize=8)
        ax.set_ylabel("protein effect (CLR difference)", fontsize=8)
        despine(ax)
    fig.suptitle("Protein effect by phenotype module", fontsize=10, fontweight="bold")
    fig.tight_layout()
    reg.save(fig, "phenotype_module_protein", SECTION, ST_PROT, "Protein effect by phenotype module",
             "Each target's protein effect grouped by its phenotype module (points; black bar = module mean). The Kruskal-Wallis test compares the effect distributions across modules with >= 2 targets (BH-FDR across proteins). Small modules make this descriptive.")


def distance_figures(res, adata, cfg, reg: FigureRegistry) -> Optional[dict]:
    """Reference order: distance ranking, atlas, PS-vs-distance map; phenotype space, module concordance; then the extension figures."""
    mtc = cfg.analysis.perturbation_effects.master_table
    master = res.master.perturbation if getattr(res, "master", None) is not None else None
    plot_distance_ranking(res.distance, reg, cfg)
    concordance_metric = None
    if master is not None and not master.empty:
        if mtc.perturbation_atlas:
            plot_perturbation_atlas(master, reg, cfg)
        if mtc.ps_distance_map:
            plot_ps_vs_distance(master, reg, cfg)
    if res.distance_space is not None and not res.distance_space.empty:
        if mtc.perturbation_space:
            plot_phenotype_space(res.distance_space, master, reg, cfg)
        plot_distance_matrix(res.distance_space, reg)
    if master is not None and not master.empty and mtc.module_concordance:
        concordance_metric = plot_module_concordance(master, reg)
    if getattr(res, "master", None) is not None:
        plot_distance_vs_protein(res.master, reg, cfg)
        plot_phenotype_space_protein(res.distance_space, res.master, reg)
        plot_phenotype_module_protein(res.distance_space, res.master, reg)
    return concordance_metric
