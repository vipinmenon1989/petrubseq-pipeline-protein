"""Perturbation / guide / target coverage QC and condition/sample summaries.

Nothing here removes cells or perturbations; low coverage is *flagged*.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from ..config import Config

CLASS_ORDER = ["single_targeting", "single_control", "multi_targeting", "multi_control", "mixed_control_targeting", "ambiguous", "unassigned"]


def _crosstab(obs: pd.DataFrame, key: str, col: Optional[str], prefix: bool = True) -> Optional[pd.DataFrame]:
    if col is None or col not in obs.columns:
        return None
    ct = pd.crosstab(obs[key], obs[col])
    if prefix:
        ct.columns = [f"{col}={v}" for v in ct.columns]
    return ct


def guide_counts(obs: pd.DataFrame, cfg: Config, condition: Optional[str] = "condition") -> pd.DataFrame:
    """One row per guide: target, class, cells (single-guide and any-guide), per condition, low-coverage flag."""
    single = obs.loc[obs["is_single_guide"]]
    df = single.groupby("guide", observed=True).agg(target=("target", "first"), control_class=("control_class", "first")).astype(str)
    df["class"] = np.where(df["control_class"].isin(["", "nan"]), "targeting", df["control_class"])
    df["n_cells_single_guide"] = single.groupby("guide", observed=True).size().reindex(df.index).fillna(0).astype(int)
    exploded = obs.loc[obs["n_guides"] > 0, ["guides"]].assign(g=lambda d: d["guides"].astype(str).str.split(";")).explode("g")
    any_counts = exploded["g"].value_counts()
    missing = [g for g in any_counts.index if g not in df.index]
    if missing:  # guides that only ever appear in multi-guide cells
        extra = pd.DataFrame({"target": [_first_target(obs, g) for g in missing], "control_class": "", "n_cells_single_guide": 0}, index=pd.Index(missing, name="guide"))
        extra["class"] = "targeting"
        df = pd.concat([df, extra])
    df["n_cells_any"] = any_counts.reindex(df.index).fillna(0).astype(int)
    ct = _crosstab(single, "guide", condition)
    if ct is not None:
        df = df.join(ct.reindex(df.index).fillna(0).astype(int))
        df["min_cells_per_condition"] = ct.reindex(df.index).fillna(0).min(axis=1).astype(int)
    thr = cfg.perturbation.min_cells_per_guide
    df["low_coverage"] = df["n_cells_single_guide"] < thr
    df["low_coverage_threshold"] = thr
    df = df.drop(columns=["control_class"])
    return df.sort_values(["class", "n_cells_single_guide"], ascending=[True, False])


def _first_target(obs: pd.DataFrame, g: str) -> str:
    m = obs["guides"].astype(str).str.split(";").map(lambda l: g in l)
    if not m.any():
        return ""
    row = obs.loc[m].iloc[0]
    gl, tl = str(row["guides"]).split(";"), str(row["targets"]).split(";")
    return tl[min(gl.index(g), len(tl) - 1)] if g in gl else ""


def target_counts(obs: pd.DataFrame, cfg: Config, condition: Optional[str] = "condition") -> pd.DataFrame:
    """One row per target (genes and control classes): guides, cells, per condition, low-coverage flag."""
    single = obs.loc[obs["is_single_guide"]]
    df = single.groupby("target", observed=True).agg(control_class=("control_class", "first")).astype(str)
    df["class"] = np.where(df["control_class"].isin(["", "nan"]), "targeting", df["control_class"])
    df["n_guides"] = single.groupby("target", observed=True)["guide"].nunique().reindex(df.index).fillna(0).astype(int)
    df["n_cells_single_guide"] = single.groupby("target", observed=True).size().reindex(df.index).fillna(0).astype(int)
    exploded = obs.loc[obs["n_guides"] > 0, ["targets"]].assign(t=lambda d: d["targets"].astype(str).str.split(";")).explode("t")
    df["n_cells_any"] = exploded["t"].value_counts().reindex(df.index).fillna(0).astype(int)
    ct = _crosstab(single, "target", condition)
    if ct is not None:
        df = df.join(ct.reindex(df.index).fillna(0).astype(int))
        df["min_cells_per_condition"] = ct.reindex(df.index).fillna(0).min(axis=1).astype(int)
    thr = cfg.perturbation.min_cells_per_target
    df["low_coverage"] = df["n_cells_single_guide"] < thr
    df["low_coverage_threshold"] = thr
    df = df.drop(columns=["control_class"])
    return df.sort_values(["class", "n_cells_single_guide"], ascending=[True, False])


def perturbation_counts(obs: pd.DataFrame, condition: Optional[str] = "condition") -> pd.DataFrame:
    """One row per ``obs['perturbation']`` label (targets + multi + unassigned)."""
    df = obs.groupby("perturbation", observed=True).size().rename("n_cells").to_frame()
    cls = obs.groupby("perturbation", observed=True)["perturbation_class"].first().astype(str)
    df.insert(0, "perturbation_class", cls.reindex(df.index))
    ct = _crosstab(obs, "perturbation", condition)
    if ct is not None:
        df = df.join(ct.reindex(df.index).fillna(0).astype(int))
    return df.sort_values(["perturbation_class", "n_cells"], ascending=[True, False])


def condition_coverage(obs: pd.DataFrame, key: str, cfg: Config, condition: Optional[str] = "condition") -> Optional[pd.DataFrame]:
    """Long table: one row per (guide|target, condition) with cells and a low-coverage flag."""
    if condition is None or condition not in obs.columns:
        return None
    single = obs.loc[obs["is_single_guide"]]
    thr = cfg.perturbation.min_cells_per_guide if key == "guide" else cfg.perturbation.min_cells_per_target
    ct = pd.crosstab(single[key], single[condition])
    long = ct.stack().rename("n_cells").reset_index()
    long.columns = [key, condition, "n_cells"]
    cls = single.groupby(key, observed=True)["control_class"].first().astype(str)
    long["class"] = long[key].map(lambda k: "targeting" if cls.get(k, "") in ("", "nan") else cls.get(k))
    long["low_coverage"] = long["n_cells"] < thr
    long["low_coverage_threshold"] = thr
    return long.sort_values([key, condition]).reset_index(drop=True)


def moi_summary(obs: pd.DataFrame, condition: Optional[str] = "condition") -> pd.DataFrame:
    if condition and condition in obs.columns:
        return pd.crosstab(obs["n_guides"], obs[condition], margins=True)
    return obs["n_guides"].value_counts().sort_index().rename("n_cells").to_frame()


def perturbation_qc(obs: pd.DataFrame, cfg: Config, guides: pd.DataFrame, targets: pd.DataFrame, cond_guide: Optional[pd.DataFrame], cond_target: Optional[pd.DataFrame]) -> Dict[str, Any]:
    """Scalar perturbation-QC summary for the report / uns."""
    n = len(obs)
    cls_counts = obs["perturbation_class"].value_counts()
    tg = targets.loc[targets["class"] == "targeting"]
    gg = guides.loc[guides["class"] == "targeting"]
    out: Dict[str, Any] = {
        "n_cells": int(n),
        "class_counts": {c: int(cls_counts.get(c, 0)) for c in CLASS_ORDER},
        "frac_unassigned": round(float((obs["n_guides"] == 0).mean()), 4),
        "frac_multi_guide": round(float((obs["n_guides"] > 1).mean()), 4),
        "frac_single_guide": round(float((obs["n_guides"] == 1).mean()), 4),
        "n_guides_total": int(len(guides)),
        "n_targets_total": int(len(targets)),
        "n_targeting_guides": int(len(gg)),
        "n_targeting_targets": int(len(tg)),
        "control_classes": {c: {"n_guides": int((guides["class"] == c).sum()), "n_cells_single_guide": int(guides.loc[guides["class"] == c, "n_cells_single_guide"].sum())} for c in sorted(set(guides["class"]) - {"targeting"})},
        "guides_per_target": {"median": float(tg["n_guides"].median()) if len(tg) else np.nan, "min": int(tg["n_guides"].min()) if len(tg) else 0, "max": int(tg["n_guides"].max()) if len(tg) else 0},
        "cells_per_targeting_guide": _dist(gg["n_cells_single_guide"]),
        "cells_per_targeting_target": _dist(tg["n_cells_single_guide"]),
        "min_cells_per_guide": cfg.perturbation.min_cells_per_guide,
        "min_cells_per_target": cfg.perturbation.min_cells_per_target,
        "n_low_coverage_guides": int(gg["low_coverage"].sum()),
        "n_low_coverage_targets": int(tg["low_coverage"].sum()),
    }
    if cond_guide is not None:
        cg = cond_guide.loc[cond_guide["class"] == "targeting"]
        out["n_low_coverage_guide_condition_pairs"] = int(cg["low_coverage"].sum())
        out["n_guide_condition_pairs"] = int(len(cg))
    if cond_target is not None:
        ctg = cond_target.loc[cond_target["class"] == "targeting"]
        out["n_low_coverage_target_condition_pairs"] = int(ctg["low_coverage"].sum())
        out["n_target_condition_pairs"] = int(len(ctg))
        out["targets_covered_in_all_conditions"] = int((tg["min_cells_per_condition"] >= cfg.perturbation.min_cells_per_target).sum()) if "min_cells_per_condition" in tg else None
    return out


def _dist(s: pd.Series) -> Dict[str, Any]:
    if len(s) == 0:
        return {}
    d = {"median": float(s.median()), "min": int(s.min()), "max": int(s.max())}
    for k in (10, 20, 30, 50):
        d[f"n_lt_{k}"] = int((s < k).sum())
    return d


def group_summary(obs: pd.DataFrame, key: str, qc_cols: Optional[List[str]] = None) -> Optional[pd.DataFrame]:
    """Per-group cell counts, perturbation-class breakdown and median QC metrics."""
    if key not in obs.columns:
        return None
    df = obs.groupby(key, observed=True).size().rename("n_cells").to_frame()
    ct = _crosstab(obs, key, "perturbation_class")
    if ct is not None:
        df = df.join(ct)
    qc_cols = qc_cols or ["total_counts", "n_genes_by_counts", "pct_counts_mt", "protein_total_counts", "protein_n_detected"]
    for c in qc_cols:
        if c in obs.columns:
            df[f"median_{c}"] = obs.groupby(key, observed=True)[c].median()
    return df
