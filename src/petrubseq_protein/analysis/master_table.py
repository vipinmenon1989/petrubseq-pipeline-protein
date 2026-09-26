"""Master perturbation tables and the distance <-> protein associations.

* ``build_master_perturbation_table``: port of the reference ``meta.build_perturbation_meta``
  (one row per target; efficacy, PS penetrance, lochNESS topology, energy distance /
  DistanceTest, co-functional and phenotype modules; sorted by energy distance; no composite
  score), followed by extension columns (phenotype nearest neighbour, PCoA1/2, strongest
  cluster association, strongest protein effect, per-protein effects). The reference column
  names are kept so that ``tables/master_perturbation_table.csv`` compares one-to-one with the
  reference ``tables/perturbation_meta.csv`` (``target_gene`` -> ``target``).
* ``build_master_protein_table``: target x protein rows. Every statistic is labelled with the
  level it was computed at: the protein effect (perturbed vs control cells of one target) and
  every joined target statistic are ``TARGET_LEVEL``; the PS <-> protein Spearman is a
  ``CELL_LEVEL`` within-target correlation.
* ``distance_protein_association``: per protein, Spearman across targets between the energy
  distance from control and the protein effect, signed and absolute separately (two BH
  families); underpowered with few proteins / targets and reported as such.
* ``phenotype_space_protein_association``: does the pairwise phenotype geometry relate to the
  protein effects? (a) Mantel test: Spearman between the target x target energy-distance
  matrix and the target x target distance between protein-effect vectors (all proteins jointly
  and one protein at a time), seeded label permutations of the protein matrix; (b) module-wise:
  Kruskal-Wallis of each protein's effect across phenotype modules with >= 2 members, plus the
  module x protein mean effects. Both are descriptive associations; with 4 proteins and ~20
  targets they are underpowered, which the ``support`` columns state.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

from ._common import bh_fdr, spearman

logger = logging.getLogger(__name__)

TARGET_LEVEL = "TARGET_LEVEL"
CELL_LEVEL = "CELL_LEVEL"
TARGET_PAIR_LEVEL = "TARGET_PAIR_LEVEL"
PHENOTYPE_MODULE_LEVEL = "PHENOTYPE_MODULE_LEVEL"
REFERENCE_META_COLUMNS = ["target", "n_cells", "target_log2fc", "target_pct_kd", "target_fdr", "is_effective_hit", "ps_mean", "ps_median", "ps_responder_fraction", "ps_net_responder_fraction", "ps_escaper_fraction",
                          "lochness_mean", "lochness_median", "lochness_peak", "lochness_pct_enriched", "energy_distance", "mmd_distance", "distance_pvalue", "distance_fdr", "distance_significant", "cofunctional_module", "phenotype_module"]


@dataclass
class MasterTables:
    perturbation: pd.DataFrame                                  # one row per target
    protein: pd.DataFrame = field(default_factory=pd.DataFrame)  # target x protein
    distance_protein: pd.DataFrame = field(default_factory=pd.DataFrame)          # per protein: signed / absolute Spearman
    distance_protein_targets: pd.DataFrame = field(default_factory=pd.DataFrame)  # long target x protein
    phenotype_protein_mantel: pd.DataFrame = field(default_factory=pd.DataFrame)
    phenotype_module_protein: pd.DataFrame = field(default_factory=pd.DataFrame)  # per protein Kruskal-Wallis across modules
    phenotype_module_means: pd.DataFrame = field(default_factory=pd.DataFrame)    # module x protein
    info: Dict[str, Any] = field(default_factory=dict)

    @property
    def empty(self) -> bool:
        return self.perturbation.empty


def _tbl(res, attr: str = "table") -> Optional[pd.DataFrame]:
    if res is None or getattr(res, "empty", False):
        return None
    df = getattr(res, attr, None)
    return df if isinstance(df, pd.DataFrame) and not df.empty and "target" in df.columns else None


def _first(df: pd.DataFrame, *names: str) -> Optional[str]:
    return next((c for c in names if c in df.columns), None)


def build_master_perturbation_table(strength=None, ps=None, lochness=None, distance=None, modules=None, space=None, protein=None, enrichment=None, primary_control: str = "ntc", n_guides: Optional[pd.Series] = None) -> pd.DataFrame:
    """Reference ``build_perturbation_meta`` semantics, then the extension columns (one row per target, TARGET_LEVEL)."""
    st = _tbl(strength)
    pss = _tbl(ps, "summary")
    ls = _tbl(lochness, "summary")
    dt = _tbl(distance)
    cm = _tbl(modules, "perturbation_modules")
    pm = _tbl(space, "phenotype_modules")
    sets = [set(df["target"].dropna().astype(str)) for df in (st, pss, ls, dt, cm, pm) if df is not None]
    if not sets:
        return pd.DataFrame(columns=["target"])
    meta = pd.DataFrame({"target": sorted(set().union(*sets))})
    t = meta["target"]
    # 1. cell counts and efficacy (perturbation strength, primary arm)
    if st is not None:
        p = st.set_index("target")
        c = _first(p, "n_perturbed", "n_cells", "n_perturbed_cells")
        if c:
            meta["n_cells"] = t.map(p[c]).astype("Int64")
        for out, cands in (("target_log2fc", (f"log2fc_{primary_control}", "log2fc")), ("target_pct_kd", (f"pct_knockdown_{primary_control}", "pct_knockdown")), ("target_fdr", (f"ks_fdr_{primary_control}", f"mwu_fdr_{primary_control}", "fdr"))):
            c = _first(p, *cands)
            if c:
                meta[out] = t.map(p[c]).astype(float)
        c = _first(p, f"is_hit_{primary_control}", "is_hit")
        if c:
            meta["is_effective_hit"] = t.map(p[c]).astype("boolean")
    # 2. PS penetrance
    if pss is not None:
        p = pss.set_index("target")
        if "n_cells" not in meta.columns and "n_perturbed_cells" in p.columns:
            meta["n_cells"] = t.map(p["n_perturbed_cells"]).astype("Int64")
        for out, c in (("ps_mean", "mean_ps"), ("ps_median", "median_ps"), ("ps_responder_fraction", "pct_successful_kd"), ("ps_net_responder_fraction", "net_pct_kd"), ("ps_escaper_fraction", "pct_escaper")):
            if c in p.columns:
                meta[out] = t.map(p[c]).astype(float)
    # 3. lochNESS topology
    if ls is not None:
        p = ls.set_index("target")
        c = _first(p, "mean_lochness_in_own_cells", "mean_lochness")
        if c:
            meta["lochness_mean"] = t.map(p[c]).astype(float)
        for out, c in (("lochness_median", "median_lochness_in_own_cells"), ("lochness_peak", "max_lochness"), ("lochness_pct_enriched", "pct_own_cells_enriched")):
            if c in p.columns:
                meta[out] = t.map(p[c]).astype(float)
    # 4. energy distance / DistanceTest
    if dt is not None:
        p = dt.set_index("target")
        if "n_cells" not in meta.columns and "n_cells" in p.columns:
            meta["n_cells"] = t.map(p["n_cells"]).astype("Int64")
        for out, c in (("energy_distance", "energy_distance"), ("mmd_distance", "mmd_distance"), ("distance_pvalue", "pvalue"), ("distance_fdr", "fdr")):
            if c in p.columns:
                meta[out] = t.map(p[c]).astype(float)
        if "significant" in p.columns:
            meta["distance_significant"] = t.map(p["significant"]).astype("boolean")
    # 5. co-functional and phenotype modules
    if cm is not None:
        c = _first(cm, "module", "cofunctional_module")
        if c:
            meta["cofunctional_module"] = t.map(cm.set_index("target")[c]).astype(object).where(lambda s: s.notna(), None)
    if pm is not None:
        c = _first(pm, "phenotype_module", "module")
        if c:
            meta["phenotype_module"] = t.map(pm.set_index("target")[c]).astype(object).where(lambda s: s.notna(), None)
    # ---- extension columns (after the reference columns) ---------------------------------
    if n_guides is not None:
        meta["n_guides"] = t.map(n_guides).astype("Int64")
    if modules is not None and not getattr(modules, "empty", True):
        ppe = modules.perturbation_program_effects.copy()
        if not ppe.empty:
            ppe["abs"] = ppe["mean_log2fc"].abs()
            top = ppe.sort_values("abs", ascending=False).groupby("target").head(1).set_index("target")
            meta["strongest_gene_program"] = t.map(top["program"]).astype(object).where(lambda s: s.notna(), None)
            meta["gene_program_effect"] = t.map(top["mean_log2fc"]).astype(float)
    if space is not None and not getattr(space, "empty", True):
        nb = space.neighbors
        if isinstance(nb, pd.DataFrame) and not nb.empty:
            first = nb[nb["rank"] == 1].set_index("target")
            meta["phenotype_nearest_neighbor"] = t.map(first["neighbor"]).astype(object).where(lambda s: s.notna(), None)
            meta["phenotype_nearest_distance"] = t.map(first["distance"]).astype(float)
        co = space.coordinates
        if isinstance(co, pd.DataFrame) and not co.empty:
            ci = co.set_index("target")
            for ax in ("PCoA1", "PCoA2"):
                if ax in ci.columns:
                    meta[ax] = t.map(ci[ax]).astype(float)
    if enrichment is not None and not getattr(enrichment, "empty", True) and not enrichment.table.empty:
        prim = enrichment.table[enrichment.table["control"] == enrichment.primary_control].sort_values(["fdr", "pval"])
        best = prim.groupby("target").head(1).set_index("target")
        meta["strongest_cluster"] = t.map(best["cluster"].astype(str)).astype(object).where(lambda s: s.notna(), None)
        meta["cluster_enrichment_log2_or"] = t.map(best["log2_odds_ratio"]).astype(float)
        meta["cluster_enrichment_fdr"] = t.map(best["fdr"]).astype(float)
        mag = enrichment.effect_magnitude.set_index("target")
        meta["cluster_composition_shift_pct"] = t.map(mag["composition_shift_pct"]).astype(float)
        meta["n_significant_clusters"] = t.map(mag["n_significant_clusters"]).astype("Int64")
    if protein is not None and not getattr(protein, "empty", True) and not protein.matrix.empty:
        E, F, Dm = protein.matrix, protein.fdr_matrix, protein.d_matrix
        common = t[t.isin(E.index)]
        best = Dm.loc[common].abs().idxmax(axis=1)
        meta["strongest_protein"] = t.map(best).astype(object).where(lambda s: s.notna(), None)
        meta["strongest_protein_effect"] = [float(E.loc[x, best[x]]) if x in best.index else np.nan for x in t]
        meta["strongest_protein_fdr"] = [float(F.loc[x, best[x]]) if x in best.index else np.nan for x in t]
        alpha = float(protein.info.get("fdr_alpha", 0.05))
        meta["n_proteins_significant"] = t.map((F.loc[common] < alpha).sum(axis=1)).astype("Int64")
        for p_ in E.columns:
            meta[f"protein_effect_{p_}"] = t.map(E[p_]).astype(float)
            meta[f"protein_fdr_{p_}"] = t.map(F[p_]).astype(float)
    if "energy_distance" in meta.columns and meta["energy_distance"].notna().any():
        meta = meta.sort_values("energy_distance", ascending=False, na_position="last").reset_index(drop=True)
    else:
        meta = meta.sort_values("target").reset_index(drop=True)
    return meta


def build_master_protein_table(master: pd.DataFrame, protein=None, concordance=None, distance_protein: Optional[pd.DataFrame] = None, module_protein: Optional[pd.DataFrame] = None, mantel: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """One row per target x protein; every statistic carries the level it was computed at.

    ``protein_effect_level = TARGET_LEVEL`` (perturbed vs control cells of one target),
    ``ps_protein_level = CELL_LEVEL`` and ``lochness_protein_level = CELL_LEVEL`` (within-target
    correlations over that target's cells), ``distance_protein_level = TARGET_LEVEL`` (one Spearman
    across targets per protein, repeated on each of the protein's rows),
    ``phenotype_module_protein_level = PHENOTYPE_MODULE_LEVEL`` and
    ``phenotype_geometry_protein_level = TARGET_PAIR_LEVEL`` (per-protein summaries, repeated per
    protein; the joint geometry statistic stays in the separate geometry table).
    """
    if protein is None or getattr(protein, "empty", True) or protein.table.empty:
        return pd.DataFrame()
    pt = protein.table
    keep = {"target": "target", "protein": "protein", "n_perturbed": "n_perturbed", "n_control": "n_control", "effect": "protein_effect", "cohen_d": "protein_cohen_d", "p_value": "protein_p_value", "fdr": "protein_fdr",
            "significant": "protein_significant", "n_samples_tested": "n_samples_tested", "n_samples_same_sign": "n_samples_same_sign", "n_guides_tested": "n_guides_tested", "n_guides_same_sign": "n_guides_same_sign"}
    out = pt[[c for c in keep if c in pt.columns]].rename(columns=keep).copy()
    out["protein_effect_level"] = TARGET_LEVEL
    if concordance is not None:
        for attr, prefix in (("ps_protein", "ps_protein"), ("lochness_protein_cells", "lochness_protein")):
            df = getattr(concordance, attr, None)
            if isinstance(df, pd.DataFrame) and not df.empty:
                w = df[df["scope"] == "within_target"][["target", "protein", "n_cells", "spearman_rho", "p_value", "fdr", "status"]] if "scope" in df.columns else df[["target", "protein", "n_cells", "spearman_rho", "p_value", "fdr", "status"]]
                w = w.rename(columns={"n_cells": f"{prefix}_n_cells", "spearman_rho": f"{prefix}_rho", "p_value": f"{prefix}_p", "fdr": f"{prefix}_fdr", "status": f"{prefix}_status"})
                out = out.merge(w, on=["target", "protein"], how="left")
                out[f"{prefix}_level"] = CELL_LEVEL
    if master is not None and not master.empty:
        cols = [c for c in ("n_cells", "n_guides", "target_log2fc", "target_pct_kd", "target_fdr", "is_effective_hit", "ps_mean", "ps_median", "ps_responder_fraction", "lochness_mean", "lochness_peak", "energy_distance", "distance_pvalue", "distance_fdr", "distance_significant",
                            "cofunctional_module", "phenotype_module", "phenotype_nearest_neighbor", "strongest_gene_program", "gene_program_effect") if c in master.columns]
        out = out.merge(master[["target"] + cols], on="target", how="left")
        out["target_columns_level"] = TARGET_LEVEL
        # program <-> protein at the target level: the target's strongest program against this protein
        ppt = getattr(concordance, "program_protein_targets", None) if concordance is not None else None
        if isinstance(ppt, pd.DataFrame) and not ppt.empty and "strongest_gene_program" in out.columns:
            key = ppt.set_index(["gene_program", "protein"])
            idx = list(zip(out["strongest_gene_program"], out["protein"]))
            out["program_protein_rho"] = [float(key.loc[k, "spearman_rho"]) if k in key.index else np.nan for k in idx]
            out["program_protein_fdr"] = [float(key.loc[k, "fdr"]) if k in key.index else np.nan for k in idx]
            out["program_protein_level"] = TARGET_LEVEL
    if isinstance(distance_protein, pd.DataFrame) and not distance_protein.empty:
        dp = distance_protein.set_index("protein")
        for src, dst in (("n_targets", "distance_protein_n_targets"), ("rho_signed", "distance_protein_rho_signed"), ("p_signed", "distance_protein_p_signed"), ("fdr_signed", "distance_protein_fdr_signed"), ("rho_abs", "distance_protein_rho_abs"), ("p_abs", "distance_protein_p_abs"), ("fdr_abs", "distance_protein_fdr_abs")):
            out[dst] = out["protein"].map(dp[src]).astype(float)
        out["distance_protein_level"] = TARGET_LEVEL
    if isinstance(module_protein, pd.DataFrame) and not module_protein.empty:
        mp = module_protein.set_index("protein")
        for src, dst in (("n_modules_eligible", "phenotype_module_protein_n_modules"), ("statistic", "phenotype_module_protein_statistic"), ("p_value", "phenotype_module_protein_p"), ("fdr", "phenotype_module_protein_fdr")):
            if src in mp.columns:
                out[dst] = out["protein"].map(mp[src]).astype(float)
        out["phenotype_module_protein_level"] = PHENOTYPE_MODULE_LEVEL
    if isinstance(mantel, pd.DataFrame) and not mantel.empty:
        mm = mantel[mantel["protein"] != "all_proteins"].set_index("protein")
        for src, dst in (("mantel_rho", "phenotype_geometry_protein_rho"), ("p_permutation", "phenotype_geometry_protein_p"), ("fdr", "phenotype_geometry_protein_fdr")):
            if src in mm.columns:
                out[dst] = out["protein"].map(mm[src]).astype(float)
        out["phenotype_geometry_protein_level"] = TARGET_PAIR_LEVEL
    if master is not None and not master.empty:
        order = {t: i for i, t in enumerate(master["target"])}
        out["_o"] = out["target"].map(order).fillna(len(order))
        out = out.sort_values(["_o", "protein"]).drop(columns="_o").reset_index(drop=True)
    return out


def _status(n: int, rho: float, fdr: float, alpha: float, min_n: int) -> str:
    if n < min_n or not np.isfinite(rho):
        return "insufficient_support"
    return "significant" if np.isfinite(fdr) and fdr < alpha else "not_significant"


def distance_protein_association(distance, protein, min_targets: int = 5, alpha: float = 0.05):
    """Per protein: Spearman across targets of energy distance vs protein effect (signed) and vs |effect| (two BH families)."""
    dt = _tbl(distance)
    if dt is None or protein is None or getattr(protein, "empty", True) or protein.matrix.empty:
        return pd.DataFrame(), pd.DataFrame()
    d = dt.set_index("target")
    E, F = protein.matrix, protein.fdr_matrix
    common = [t for t in d.index if t in E.index]
    long_rows, rows = [], []
    for p in E.columns:
        for t in common:
            long_rows.append({"target": t, "protein": p, "energy_distance": float(d.loc[t, "energy_distance"]), "distance_fdr": float(d.loc[t, "fdr"]), "distance_significant": bool(d.loc[t, "significant"]),
                              "protein_effect": float(E.loc[t, p]), "protein_fdr": float(F.loc[t, p])})
        x = d.loc[common, "energy_distance"].to_numpy(float)
        y = E.loc[common, p].to_numpy(float)
        r_s, p_s = spearman(x, y) if len(common) >= min_targets else (np.nan, np.nan)
        r_a, p_a = spearman(x, np.abs(y)) if len(common) >= min_targets else (np.nan, np.nan)
        rows.append({"protein": p, "n_targets": len(common), "rho_signed": r_s, "p_signed": p_s, "rho_abs": r_a, "p_abs": p_a})
    summ = pd.DataFrame(rows)
    if not summ.empty:
        summ["fdr_signed"] = bh_fdr(summ["p_signed"].to_numpy())
        summ["fdr_abs"] = bh_fdr(summ["p_abs"].to_numpy())
        summ["status_signed"] = [_status(n, r, f, alpha, min_targets) for n, r, f in zip(summ["n_targets"], summ["rho_signed"], summ["fdr_signed"])]
        summ["status_abs"] = [_status(n, r, f, alpha, min_targets) for n, r, f in zip(summ["n_targets"], summ["rho_abs"], summ["fdr_abs"])]
        summ["support"] = np.where(summ["n_targets"] >= 10, "ok", "few_targets(<10)")
        summ["min_targets"] = min_targets
        summ["test"] = "spearman_across_targets"
        summ["analysis_level"] = TARGET_LEVEL
    long = pd.DataFrame(long_rows)
    if not long.empty:
        long["analysis_level"] = TARGET_LEVEL
    return summ, long


def mantel_test(D1: np.ndarray, D2: np.ndarray, n_perm: int, seed: int):
    """Mantel test between two K x K distance matrices over the SAME K targets in the SAME order.

    Statistic: Spearman rho between the upper triangles (K(K-1)/2 pairs, each pair once). The pairs
    share targets and are not independent, so no parametric p-value is used: the second matrix's
    target labels are permuted ``n_perm`` times (rows and columns together, ``default_rng(seed)``)
    and ``p = (1 + #(rho_perm >= rho_obs - 1e-12)) / (1 + n_perm)`` (one-sided, plus-one corrected).
    Both matrices must be symmetric with a zero diagonal (checked).
    """
    D1, D2 = np.asarray(D1, float), np.asarray(D2, float)
    K = D1.shape[0]
    if D1.shape != (K, K) or D2.shape != (K, K):
        raise ValueError("Mantel: both matrices must be K x K over the same targets")
    for name, D in (("RNA", D1), ("protein", D2)):
        if not np.allclose(D, D.T, atol=1e-8) or not np.allclose(np.diag(D), 0.0, atol=1e-8):
            raise ValueError(f"Mantel: the {name} distance matrix must be symmetric with a zero diagonal")
    iu = np.triu_indices(K, 1)
    obs, _ = spearman(D1[iu], D2[iu])
    if not np.isfinite(obs) or n_perm <= 0:
        return obs, np.nan, int(len(iu[0]))
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(K)
        r, _ = spearman(D1[iu], D2[np.ix_(perm, perm)][iu])
        if np.isfinite(r) and r >= obs - 1e-12:
            count += 1
    return obs, (1.0 + count) / (1.0 + n_perm), int(len(iu[0]))


def protein_profile_distances(V: np.ndarray) -> np.ndarray:
    """Target x target Euclidean distance between protein-effect vectors (targets x proteins).

    All effects are CLR mean differences on the same scale, so no per-protein rescaling is applied;
    a single protein gives |E_i - E_j|.
    """
    V = np.asarray(V, float)
    if V.ndim == 1:
        V = V[:, None]
    return np.sqrt(((V[:, None, :] - V[None, :, :]) ** 2).sum(-1))


def phenotype_space_protein_association(space, protein, min_targets: int = 5, n_perm: int = 999, seed: int = 0, alpha: float = 0.05, min_targets_per_module: int = 3):
    """(a) Mantel: RNA pairwise energy distances vs protein-effect-profile distances (TARGET_PAIR_LEVEL);
    (b) module-wise: Kruskal-Wallis of each protein's target-level effect across the phenotype modules
    with >= ``min_targets_per_module`` targets (PHENOTYPE_MODULE_LEVEL), plus the module x protein means."""
    empty = pd.DataFrame()
    if space is None or getattr(space, "empty", True) or protein is None or getattr(protein, "empty", True) or protein.matrix.empty:
        return empty, empty, empty
    M = space.distance_matrix
    E = protein.matrix
    common = [t for t in M.index if t in E.index]          # one order for both matrices
    mantel_rows: List[dict] = []
    if len(common) >= min_targets:
        D1 = M.loc[common, common].to_numpy(float)
        V = E.loc[common].to_numpy(float)
        r, p, n_pairs = mantel_test(D1, protein_profile_distances(V), n_perm, seed)
        mantel_rows.append({"protein": "all_proteins", "n_targets": len(common), "n_pairs": n_pairs, "n_proteins": int(V.shape[1]), "protein_distance_metric": "euclidean over the protein-effect vector (CLR differences)", "mantel_rho": r, "p_permutation": p})
        for k, prot in enumerate(E.columns):
            r, p, n_pairs = mantel_test(D1, protein_profile_distances(V[:, k]), n_perm, seed + 1 + k)
            mantel_rows.append({"protein": prot, "n_targets": len(common), "n_pairs": n_pairs, "n_proteins": 1, "protein_distance_metric": "|E_i,p - E_j,p|", "mantel_rho": r, "p_permutation": p})
    mantel = pd.DataFrame(mantel_rows)
    if not mantel.empty:
        single = mantel["protein"] != "all_proteins"
        mantel["fdr"] = np.nan
        mantel.loc[single, "fdr"] = bh_fdr(mantel.loc[single, "p_permutation"].to_numpy())
        mantel["status"] = [("significant" if (np.isfinite(p) and p < alpha) else "not_significant") if pr == "all_proteins" else _status(n, r, f, alpha, min_targets)
                            for pr, n, r, f, p in zip(mantel["protein"], mantel["n_targets"], mantel["mantel_rho"], mantel["fdr"], mantel["p_permutation"])]
        mantel["test"] = "mantel_spearman_label_permutation"
        mantel["n_permutations"] = n_perm
        mantel["seed"] = seed
        mantel["support"] = np.where((mantel["n_targets"] >= 10) & ((mantel["n_proteins"] >= 5) | (mantel["protein"] != "all_proteins")), "ok", "underpowered(<10 targets or <5 proteins)")
        mantel["analysis_level"] = TARGET_PAIR_LEVEL
    # module-wise
    pm = space.phenotype_modules
    kw_rows, mean_rows = [], []
    if isinstance(pm, pd.DataFrame) and not pm.empty:
        mod = pm.set_index("target")["phenotype_module"]
        tt = [t for t in common if t in mod.index]
        if tt:
            lab = mod.loc[tt]
            sizes = lab.value_counts()
            usable = sorted([m for m in sizes.index if sizes[m] >= min_targets_per_module], key=lambda s: (len(str(s)), str(s)))
            alpha_p = float(protein.info.get("fdr_alpha", 0.05))
            for prot in E.columns:
                groups = [E.loc[lab.index[lab == m], prot].to_numpy(float) for m in usable]
                if len(groups) >= 2:
                    try:
                        h, p = stats.kruskal(*groups)
                    except ValueError:
                        h, p = np.nan, np.nan
                else:
                    h, p = np.nan, np.nan
                kw_rows.append({"protein": prot, "n_targets": int(sum(len(g) for g in groups)), "n_modules_eligible": len(usable), "modules_eligible": ";".join(usable), "n_modules_excluded": int(len(sizes) - len(usable)),
                                "min_targets_per_module": min_targets_per_module, "test": "kruskal_wallis" if len(groups) != 2 else "kruskal_wallis(2 groups = Mann-Whitney U)", "statistic": float(h), "p_value": float(p)})
                for m in sizes.index:
                    v = E.loc[lab.index[lab == m], prot].to_numpy(float)
                    mean_rows.append({"phenotype_module": m, "protein": prot, "n_targets": int(len(v)), "eligible": bool(sizes[m] >= min_targets_per_module), "mean_effect": float(np.mean(v)), "median_effect": float(np.median(v)),
                                      "n_significant": int((protein.fdr_matrix.loc[lab.index[lab == m], prot] < alpha_p).sum())})
    kw = pd.DataFrame(kw_rows)
    if not kw.empty:
        kw["fdr"] = bh_fdr(kw["p_value"].to_numpy())
        kw["status"] = [("insufficient_support" if k < 2 else ("significant" if np.isfinite(f) and f < alpha else "not_significant")) for k, f in zip(kw["n_modules_eligible"], kw["fdr"])]
        kw["support"] = np.where(kw["n_targets"] >= 10, "ok", "few_targets(<10)")
        kw["analysis_level"] = PHENOTYPE_MODULE_LEVEL
    means = pd.DataFrame(mean_rows)
    if not means.empty:
        means["analysis_level"] = PHENOTYPE_MODULE_LEVEL
        means = means.sort_values(["phenotype_module", "protein"]).reset_index(drop=True)
    return mantel, kw, means


def build_master_tables(cfg, strength=None, ps=None, lochness=None, distance=None, modules=None, space=None, protein=None, enrichment=None, concordance=None, n_guides: Optional[pd.Series] = None) -> MasterTables:
    pe = cfg.analysis.perturbation_effects
    mt = pe.master_table
    master = build_master_perturbation_table(strength, ps, lochness, distance, modules, space, protein, enrichment, primary_control=getattr(strength, "primary_control", pe.strength.primary_control), n_guides=n_guides)
    res = MasterTables(master)
    res.info = {"status": "ok" if not master.empty else "no target-level results", "analysis_level": TARGET_LEVEL, "n_targets": int(len(master)), "n_columns": int(master.shape[1]), "reference_columns": [c for c in REFERENCE_META_COLUMNS if c in master.columns],
                "extension_columns": [c for c in master.columns if c not in REFERENCE_META_COLUMNS], "sort": "energy_distance descending (reference) or target", "composite_score": "none (reference principle)"}
    if master.empty:
        return res
    if mt.protein_associations and protein is not None and not getattr(protein, "empty", True):
        res.distance_protein, res.distance_protein_targets = distance_protein_association(distance, protein, mt.min_targets, pe.protein.fdr_alpha)
        res.phenotype_protein_mantel, res.phenotype_module_protein, res.phenotype_module_means = phenotype_space_protein_association(space, protein, mt.min_targets, mt.mantel_permutations, cfg.compute.seed, pe.protein.fdr_alpha, mt.min_targets_per_module)
        n_prot = int(protein.matrix.shape[1])
        res.info.update({"distance_protein": "TARGET_LEVEL: per protein Spearman across targets, energy distance vs effect (signed) and vs |effect| (BH within each family)",
                         "phenotype_protein": "TARGET_PAIR_LEVEL Mantel (Spearman of upper triangles, seeded target-label permutations) of the RNA distance matrix vs protein-effect-profile distances; PHENOTYPE_MODULE_LEVEL Kruskal-Wallis of effects across eligible modules",
                         "n_proteins": n_prot, "n_distance_protein_tests": int(len(res.distance_protein)) * 2,
                         "n_distance_protein_significant_signed": int((res.distance_protein.get("status_signed", pd.Series(dtype=str)) == "significant").sum()) if not res.distance_protein.empty else 0,
                         "n_distance_protein_significant_abs": int((res.distance_protein.get("status_abs", pd.Series(dtype=str)) == "significant").sum()) if not res.distance_protein.empty else 0,
                         "mantel_all_proteins_rho": float(res.phenotype_protein_mantel.iloc[0]["mantel_rho"]) if not res.phenotype_protein_mantel.empty else np.nan,
                         "mantel_all_proteins_p": float(res.phenotype_protein_mantel.iloc[0]["p_permutation"]) if not res.phenotype_protein_mantel.empty else np.nan,
                         "min_targets_per_module": mt.min_targets_per_module,
                         "power_note": f"{n_prot} proteins: target-level and geometry associations are underpowered and descriptive" if n_prot < 5 else "", "caveat": "associations, not mediation or causality"})
    res.protein = build_master_protein_table(master, protein, concordance, res.distance_protein, res.phenotype_module_protein, res.phenotype_protein_mantel)
    res.info["n_protein_rows"] = int(len(res.protein))
    res.info["levels"] = {"protein_effect": TARGET_LEVEL, "ps_protein": CELL_LEVEL, "lochness_protein": CELL_LEVEL, "distance_protein": TARGET_LEVEL, "phenotype_module_protein": PHENOTYPE_MODULE_LEVEL, "phenotype_geometry_protein": TARGET_PAIR_LEVEL, "program_protein": TARGET_LEVEL, "target_columns": TARGET_LEVEL}
    return res
