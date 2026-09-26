"""QC metric calculation (always) and optional, conservative filtering.

Calculation and filtering are deliberately separate functions. v0.1 defaults
to ``qc.filter: false``: every metric and flag is written to ``obs`` and to
``tables/cell_qc.tsv`` so thresholds can be chosen from the data.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import Config

logger = logging.getLogger("petrubseq_protein")


def _mad_flags(x: np.ndarray, n_mads: float, log: bool = True) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    v = np.log1p(x) if log else x.astype(float)
    ok = np.isfinite(v)
    med = np.nanmedian(v[ok])
    mad = np.nanmedian(np.abs(v[ok] - med)) * 1.4826
    low = med - n_mads * mad
    high = med + n_mads * mad
    return (v < low) & ok, (v > high) & ok, {"median": float(med), "mad": float(mad), "low": float(low), "high": float(high), "log1p": log}


def rna_qc(adata: ad.AnnData, cfg: Config) -> Dict[str, Any]:
    """Add RNA QC columns to ``adata.obs`` / ``adata.var``; return a summary dict.

    Uses ``layers['counts']`` when present. Otherwise total counts come from
    ``obs['rna_total_counts_provided']`` (if any) and the gene fractions are
    computed on ``expm1(X)``, which is exact for log-normalized data.
    """
    import scanpy as sc

    rcfg = cfg.rna
    var = adata.var
    names = var.index.str.upper()
    var["mt"] = names.str.startswith(rcfg.mito_prefix.upper())
    var["ribo"] = names.str.startswith(tuple(p.upper() for p in rcfg.ribo_prefixes))
    var["hb"] = names.str.contains(rcfg.hb_pattern, regex=True)
    info: Dict[str, Any] = {"n_mito_genes": int(var["mt"].sum()), "n_ribo_genes": int(var["ribo"].sum()), "n_hb_genes": int(var["hb"].sum())}
    layer = "counts" if "counts" in adata.layers else ("reconstructed_counts" if "reconstructed_counts" in adata.layers else None)
    if layer is not None:
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], layer=layer, percent_top=None, log1p=False, inplace=True)
        info["source"] = f"layers['{layer}']"
    else:
        X = sp.csr_matrix(adata.X)
        E = X.copy()
        E.data = np.expm1(E.data)
        tot = np.asarray(E.sum(axis=1)).ravel()
        adata.obs["n_genes_by_counts"] = np.diff(X.indptr)
        if "rna_total_counts_provided" in adata.obs:
            adata.obs["total_counts"] = adata.obs["rna_total_counts_provided"].to_numpy()
            info["source"] = "obs['rna_total_counts_provided'] + expm1(X) fractions"
        else:
            adata.obs["total_counts"] = tot
            info["source"] = "expm1(X) sums (relative units)"
        with np.errstate(divide="ignore", invalid="ignore"):
            adata.obs["pct_counts_mt"] = 100 * np.asarray(E[:, var["mt"].to_numpy()].sum(axis=1)).ravel() / tot
            adata.obs["pct_counts_ribo"] = 100 * np.asarray(E[:, var["ribo"].to_numpy()].sum(axis=1)).ravel() / tot
            adata.obs["pct_counts_hb"] = 100 * np.asarray(E[:, var["hb"].to_numpy()].sum(axis=1)).ravel() / tot
        adata.var["n_cells_by_counts"] = np.diff(X.tocsc().indptr)
    obs = adata.obs
    flags = pd.DataFrame(index=obs.index)
    thr: Dict[str, Any] = {}
    lo, hi, t = _mad_flags(obs["total_counts"].to_numpy(dtype=float), cfg.qc.flags.rna_n_mads)
    flags["rna_outlier_low_counts"], flags["rna_outlier_high_counts"], thr["total_counts"] = lo, hi, t
    lo, hi, t = _mad_flags(obs["n_genes_by_counts"].to_numpy(dtype=float), cfg.qc.flags.rna_n_mads)
    flags["rna_outlier_low_genes"], flags["rna_outlier_high_genes"], thr["n_genes_by_counts"] = lo, hi, t
    fr = cfg.qc.filter.rna
    flags["rna_low_genes"] = obs["n_genes_by_counts"].to_numpy() < fr.min_genes if fr.min_genes is not None else False
    flags["rna_high_mt"] = obs["pct_counts_mt"].to_numpy() >= fr.max_pct_mt if fr.max_pct_mt is not None else False
    if fr.min_counts is not None:
        flags["rna_low_counts"] = obs["total_counts"].to_numpy() < fr.min_counts
    flags["rna_qc_fail"] = flags[["rna_low_genes", "rna_high_mt"] + (["rna_low_counts"] if "rna_low_counts" in flags else [])].any(axis=1)
    for c in flags.columns:
        adata.obs[c] = flags[c].to_numpy()
    info["thresholds"] = {"min_genes": fr.min_genes, "max_pct_mt": fr.max_pct_mt, "max_pct_hb": fr.max_pct_hb, "min_counts": fr.min_counts, "min_cells_per_gene": fr.min_cells_per_gene, "n_mads": cfg.qc.flags.rna_n_mads, "mad": thr}
    info["summary"] = {k: {"median": float(np.nanmedian(obs[k])), "min": float(np.nanmin(obs[k])), "max": float(np.nanmax(obs[k]))} for k in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb") if k in obs.columns}
    info["flag_counts"] = {c: int(flags[c].sum()) for c in flags.columns}
    return info


def protein_qc(adata: ad.AnnData, cfg: Config) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """Per-cell protein QC into ``obs`` and a per-feature summary table."""
    prot: Optional[pd.DataFrame] = adata.obsm.get("protein")
    counts: Optional[pd.DataFrame] = adata.obsm.get("protein_counts")
    feats: pd.DataFrame = adata.uns["protein_features"].copy() if "protein_features" in adata.uns else pd.DataFrame()
    info: Dict[str, Any] = {}
    obs = adata.obs
    if counts is not None:
        iso_cols = [c for c in counts.columns if feats.get("is_isotype", pd.Series(dtype=bool)).get(c, False)] if cfg.protein.use_isotypes_for_qc else []
        tgt_cols = [c for c in counts.columns if not feats.get("is_isotype", pd.Series(dtype=bool)).get(c, False)]
        obs["protein_total_counts"] = counts.sum(axis=1, min_count=1).to_numpy()
        obs["protein_total_counts_targeting"] = counts[tgt_cols].sum(axis=1, min_count=1).to_numpy()
        obs["protein_isotype_counts"] = counts[iso_cols].sum(axis=1, min_count=1).to_numpy() if iso_cols else np.nan
        obs["protein_n_detected"] = (counts[tgt_cols] > 0).sum(axis=1).where(counts[tgt_cols].notna().any(axis=1)).to_numpy()
        obs["has_protein_counts"] = counts.notna().all(axis=1).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            obs["protein_pct_isotype"] = 100 * obs["protein_isotype_counts"] / obs["protein_total_counts"]
        x = obs["protein_total_counts"].to_numpy(dtype=float)
        lo, hi, t = _mad_flags(np.where(np.isfinite(x), x, np.nan), cfg.qc.flags.protein_n_mads)
        obs["protein_outlier_low_counts"], obs["protein_outlier_high_counts"] = lo, hi
        q99 = float(np.nanquantile(x, 0.99)) if np.isfinite(x).any() else np.nan
        fl, fp = cfg.qc.flags, cfg.qc.filter.protein
        obs["protein_extreme_counts"] = np.isfinite(x) & (x > fl.extreme_fold * q99)
        info["extreme_counts_threshold"] = {"q99": q99, "fold": fl.extreme_fold, "n_cells": int(obs["protein_extreme_counts"].sum())}
        fail = np.zeros(adata.n_obs, dtype=bool)
        if fp.min_total_counts is not None:
            fail |= x < fp.min_total_counts
        if fp.min_proteins_detected is not None:
            fail |= obs["protein_n_detected"].to_numpy(dtype=float) < fp.min_proteins_detected
        if fp.max_pct_isotype is not None:
            pi = obs["protein_pct_isotype"].to_numpy(dtype=float)
            obs["protein_high_isotype"] = np.isfinite(pi) & (pi > fp.max_pct_isotype)
            fail |= obs["protein_high_isotype"].to_numpy(bool)
        if fp.remove_extreme_counts:
            fail |= obs["protein_extreme_counts"].to_numpy(bool)
        obs["protein_qc_fail"] = fail
        info["source"] = "obsm['protein_counts']"
        info["thresholds"] = {"min_total_counts": fp.min_total_counts, "min_proteins_detected": fp.min_proteins_detected, "max_pct_isotype": fp.max_pct_isotype, "remove_extreme_counts": fp.remove_extreme_counts, "n_mads": fl.protein_n_mads, "mad_total_counts": t}
        info["cells_without_counts"] = int((~obs["has_protein_counts"]).sum())
    elif prot is not None:
        obs["protein_n_detected"] = (prot > 0).sum(axis=1).to_numpy()
        obs["protein_qc_fail"] = False
        info["source"] = "obsm['protein'] (normalized; detection = value > 0)"
    # Per-feature table
    rows: List[Dict[str, Any]] = []
    names = list(feats.index) if len(feats) else list((counts if counts is not None else prot).columns)
    for n in names:
        r: Dict[str, Any] = {"protein": n}
        if len(feats):
            for c in feats.columns:
                r[c] = feats.loc[n, c]
        if counts is not None and n in counts.columns:
            v = counts[n].dropna().to_numpy(dtype=float)
            r.update(counts_mean=float(v.mean()), counts_median=float(np.median(v)), counts_q90=float(np.quantile(v, 0.9)), counts_max=float(v.max()), counts_zero_frac=float((v == 0).mean()), n_cells_with_counts=int(v.size))
        if prot is not None and n in prot.columns:
            v = prot[n].dropna().to_numpy(dtype=float)
            r.update(norm_mean=float(v.mean()), norm_median=float(np.median(v)), norm_q90=float(np.quantile(v, 0.9)), norm_max=float(v.max()), norm_zero_frac=float((v == 0).mean()))
        rows.append(r)
    table = pd.DataFrame(rows).set_index("protein")
    if counts is not None and "is_isotype" in table.columns and table["is_isotype"].any():
        # Fraction of cells in which the antibody exceeds its matched isotype
        # control (or the max over all isotypes when no mapping is known).
        iso_cols = table.index[table["is_isotype"].astype(bool)].tolist()
        iso_map = table["isotype_control"].to_dict() if "isotype_control" in table.columns else {}
        frac = {}
        for n in table.index:
            if n in iso_cols or n not in counts.columns:
                continue
            ref = counts[iso_map[n]] if iso_map.get(n) in counts.columns else counts[iso_cols].max(axis=1)
            ok = counts[n].notna() & ref.notna()
            frac[n] = float((counts.loc[ok, n] > ref[ok]).mean()) if ok.any() else np.nan
        table["frac_cells_above_isotype"] = pd.Series(frac).reindex(table.index)
        thr = cfg.qc.flags.background_min_fraction_above_isotype
        table["background_dominated"] = (~table["is_isotype"].astype(bool)) & (table["frac_cells_above_isotype"] < thr)
        info["background_dominated"] = table.index[table["background_dominated"].fillna(False)].tolist()
        info["background_criterion"] = f"fraction of cells with count > matched isotype count < {thr}"
    for k in ("protein_total_counts", "protein_n_detected", "protein_pct_isotype"):
        if k in obs:
            v = obs[k].to_numpy(dtype=float)
            info.setdefault("summary", {})[k] = {"median": float(np.nanmedian(v)), "min": float(np.nanmin(v)), "max": float(np.nanmax(v))}
    return info, table


def qc_summary_table(obs: pd.DataFrame, group_col: Optional[str] = None) -> pd.DataFrame:
    """Median QC metrics per group (condition/lane) plus an ALL row."""
    metrics = [("n_genes_by_counts", "median_genes"), ("total_counts", "median_counts"), ("pct_counts_mt", "median_pct_mt"), ("pct_counts_ribo", "median_pct_ribo"), ("protein_total_counts", "median_adt_counts"), ("protein_n_detected", "median_proteins_detected"), ("protein_pct_isotype", "median_pct_isotype")]
    metrics = [(m, l) for m, l in metrics if m in obs.columns]
    rows = []
    groups = [(g, sub) for g, sub in obs.groupby(obs[group_col].astype(str), observed=True)] if group_col and group_col in obs.columns and obs[group_col].nunique() <= 50 else []
    for g, sub in groups + [("ALL", obs)]:
        row = {"group": g, "n_cells": int(len(sub))}
        for m, l in metrics:
            row[l] = float(np.nanmedian(sub[m].to_numpy(dtype=float)))
        rows.append(row)
    return pd.DataFrame(rows).round(2)
