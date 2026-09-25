"""Standard representations (diagnostic, not analytical, in this phase).

RNA      : HVG (Scanpy ``seurat`` flavor on log data) -> scaled HVG copy ->
           PCA -> neighbors -> UMAP.  ``X`` is never scaled in place.
Protein  : selected normalized matrix (targeting antibodies only) -> scale ->
           PCA -> neighbors -> UMAP.
Joint    : optional ``concat_pcs`` (block-normalized RNA PCs + weighted
           protein PCs) -> neighbors -> UMAP.
Diagnostic: RNA-vs-protein neighbourhood agreement and PC-vs-depth
           correlations on a cell sample (always computed).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import Config

logger = logging.getLogger("petrubseq_protein")


def _pca_on(X: np.ndarray | sp.spmatrix, n_comps: int, scale: bool, max_value: Optional[float], seed: int) -> Dict[str, Any]:
    import scanpy as sc

    tmp = ad.AnnData(X=X.copy() if sp.issparse(X) else np.ascontiguousarray(X, dtype=np.float32))
    if scale:
        if sp.issparse(tmp.X):
            tmp.X = tmp.X.toarray()
        sc.pp.scale(tmp, zero_center=True, max_value=max_value)
    sc.pp.pca(tmp, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    return {"X_pca": tmp.obsm["X_pca"].astype(np.float32), "loadings": tmp.varm["PCs"].astype(np.float32), "variance_ratio": np.asarray(tmp.uns["pca"]["variance_ratio"], dtype=float), "variance": np.asarray(tmp.uns["pca"]["variance"], dtype=float)}


def _neighbors_umap(adata: ad.AnnData, rep: str, key: str, n_neighbors: int, n_pcs: Optional[int], cfg: Config, umap_key: str, do_umap: bool) -> Dict[str, Any]:
    import scanpy as sc

    info: Dict[str, Any] = {"neighbors": {"use_rep": rep, "n_neighbors": n_neighbors, "n_pcs": n_pcs, "metric": cfg.neighbors.metric, "key": key}}
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep, key_added=key, metric=cfg.neighbors.metric, random_state=cfg.compute.seed)
    if do_umap:
        sc.tl.umap(adata, neighbors_key=key, min_dist=cfg.umap.min_dist, spread=cfg.umap.spread, random_state=cfg.compute.seed)
        adata.obsm[umap_key] = adata.obsm.pop("X_umap").astype(np.float32)
        adata.uns.pop("umap", None)
        info["umap"] = {"key": umap_key, "min_dist": cfg.umap.min_dist, "spread": cfg.umap.spread, "random_state": cfg.compute.seed}
    else:
        info["umap"] = None
    return info


# ---------------------------------------------------------------------------
# RNA
# ---------------------------------------------------------------------------


def rna_embedding(adata: ad.AnnData, cfg: Config) -> Dict[str, Any]:
    import scanpy as sc

    rcfg = cfg.rna
    n_genes = adata.n_vars
    use_hvg = {"auto": n_genes > rcfg.hvg.n_top_genes, "true": True, "false": False}[rcfg.hvg.enabled]
    info: Dict[str, Any] = {}
    if use_hvg:
        kw = {"n_top_genes": min(rcfg.hvg.n_top_genes, n_genes), "flavor": rcfg.hvg.flavor, "batch_key": rcfg.hvg.batch_key}
        if rcfg.hvg.flavor == "seurat_v3":
            layer = "counts" if "counts" in adata.layers else ("reconstructed_counts" if "reconstructed_counts" in adata.layers else None)
            if layer is None:
                raise ValueError("rna.hvg.flavor=seurat_v3 needs a counts layer")
            kw["layer"] = layer
        sc.pp.highly_variable_genes(adata, **kw)
        mask = adata.var["highly_variable"].to_numpy()
        info["hvg"] = {"flavor": rcfg.hvg.flavor, "n_top_genes": int(mask.sum()), "batch_key": rcfg.hvg.batch_key, "input": "X (log-normalized)"}
    else:
        mask = np.ones(n_genes, dtype=bool)
        adata.var["highly_variable"] = True
        info["hvg"] = {"flavor": None, "n_top_genes": int(n_genes), "note": "all genes used"}
    n_comps = int(min(rcfg.n_pcs, adata.n_obs - 1, int(mask.sum()) - 1))
    res = _pca_on(adata.X[:, mask], n_comps, rcfg.scale, rcfg.scale_max_value, cfg.compute.seed)
    adata.obsm["X_pca"] = res["X_pca"]
    load = np.zeros((n_genes, n_comps), dtype=np.float32)
    load[mask] = res["loadings"]
    adata.varm["PCs"] = load
    adata.uns["pca"] = {"variance_ratio": res["variance_ratio"], "variance": res["variance"], "params": {"n_comps": n_comps, "scaled": rcfg.scale, "scale_max_value": rcfg.scale_max_value, "use_highly_variable": bool(use_hvg), "svd_solver": "arpack", "random_state": cfg.compute.seed}}
    info["pca"] = {"key": "X_pca", "n_comps": n_comps, "scaled": rcfg.scale, "scale_max_value": rcfg.scale_max_value, "variance_ratio_top10": [round(float(v), 4) for v in res["variance_ratio"][:10]], "variance_ratio_total": round(float(res["variance_ratio"].sum()), 4)}
    n_pcs = int(min(cfg.neighbors.n_pcs or n_comps, n_comps))
    info.update(_neighbors_umap(adata, "X_pca", "rna", cfg.neighbors.n_neighbors, n_pcs, cfg, "X_umap_rna", cfg.umap.enabled))
    logger.info("RNA representation: %d HVGs -> PCA(%d, scaled=%s) -> neighbors(%d) %s", int(mask.sum()), n_comps, rcfg.scale, cfg.neighbors.n_neighbors, "-> UMAP" if cfg.umap.enabled else "")
    return info


# ---------------------------------------------------------------------------
# Protein
# ---------------------------------------------------------------------------


def protein_embedding(adata: ad.AnnData, cfg: Config, emb_key: Optional[str]) -> Dict[str, Any]:
    pcfg = cfg.protein
    info: Dict[str, Any] = {}
    if emb_key is None or emb_key not in adata.obsm:
        info["skipped"] = "no protein matrix"
        return info
    P: pd.DataFrame = adata.obsm[emb_key]
    feats = adata.uns.get("protein_features")
    cols = list(P.columns)
    if pcfg.exclude_isotypes_from_embedding and feats is not None and "is_isotype" in feats:
        iso = set(feats.index[feats["is_isotype"].astype(bool)])
        cols = [c for c in cols if c not in iso]
    n_feat = len(cols)
    ok = P[cols].notna().all(axis=1).to_numpy()
    n_comps = int(min(pcfg.n_pcs, n_feat - 1, int(ok.sum()) - 1))
    res = _pca_on(P.loc[ok, cols].to_numpy(dtype=np.float32), n_comps, pcfg.scale, pcfg.scale_max_value, cfg.compute.seed)
    pca = np.full((adata.n_obs, n_comps), np.nan, dtype=np.float32)
    pca[ok] = res["X_pca"]
    adata.obsm["X_pca_protein"] = pca
    adata.uns["pca_protein"] = {"variance_ratio": res["variance_ratio"], "variance": res["variance"], "loadings": pd.DataFrame(res["loadings"], index=cols, columns=[f"PC{i+1}" for i in range(n_comps)]), "params": {"representation": emb_key, "features": cols, "n_comps": n_comps, "scaled": pcfg.scale, "scale_max_value": pcfg.scale_max_value, "random_state": cfg.compute.seed}}
    info["pca"] = {"key": "X_pca_protein", "representation": emb_key, "n_features": n_feat, "isotypes_excluded": pcfg.exclude_isotypes_from_embedding, "n_comps": n_comps, "scaled": pcfg.scale, "cells_used": int(ok.sum()), "variance_ratio": [round(float(v), 4) for v in res["variance_ratio"]], "variance_ratio_total": round(float(res["variance_ratio"].sum()), 4)}
    do_umap = cfg.umap.enabled and n_feat >= cfg.umap.min_features
    if not do_umap and cfg.umap.enabled:
        info["skipped_umap"] = f"only {n_feat} protein features (< umap.min_features={cfg.umap.min_features})"
    if ok.all():
        info.update(_neighbors_umap(adata, "X_pca_protein", "protein", cfg.neighbors.n_neighbors, n_comps, cfg, "X_umap_protein", do_umap))
    else:
        sub = ad.AnnData(obs=adata.obs.loc[ok, []].copy())
        sub.obsm["X_pca_protein"] = pca[ok]
        info.update(_neighbors_umap(sub, "X_pca_protein", "protein", cfg.neighbors.n_neighbors, n_comps, cfg, "X_umap_protein", do_umap))
        if do_umap:
            um = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
            um[ok] = sub.obsm["X_umap_protein"]
            adata.obsm["X_umap_protein"] = um
        info["neighbors"]["note"] = f"graph computed on the {int(ok.sum())} cells with protein values; not stored in obsp"
    logger.info("protein representation: %s (%d features, isotypes excluded=%s) -> PCA(%d) %s", emb_key, n_feat, pcfg.exclude_isotypes_from_embedding, n_comps, "-> UMAP" if do_umap else "")
    return info


# ---------------------------------------------------------------------------
# Multimodal (optional) + cross-modality diagnostic (always)
# ---------------------------------------------------------------------------


def _block_normalize(Z: np.ndarray) -> np.ndarray:
    v = np.nansum(np.nanvar(Z, axis=0))
    return Z / np.sqrt(v) if v > 0 else Z


def multimodal_embedding(adata: ad.AnnData, cfg: Config) -> Dict[str, Any]:
    mcfg = cfg.multimodal
    info: Dict[str, Any] = {"enabled": mcfg.enabled, "method": mcfg.method}
    if not mcfg.enabled:
        return info
    if "X_pca" not in adata.obsm or "X_pca_protein" not in adata.obsm:
        info["skipped"] = "needs both X_pca and X_pca_protein"
        return info
    R = _block_normalize(np.asarray(adata.obsm["X_pca"], dtype=np.float32))
    P = np.asarray(adata.obsm["X_pca_protein"], dtype=np.float32)
    ok = np.isfinite(P).all(axis=1)
    P = _block_normalize(np.where(ok[:, None], P, 0.0)) * mcfg.protein_weight
    M = np.hstack([R, P]).astype(np.float32)
    M[~ok] = np.nan
    adata.obsm["X_multimodal"] = M
    info["params"] = {"protein_weight": mcfg.protein_weight, "rna_dims": R.shape[1], "protein_dims": P.shape[1], "block_scaling": "each block divided by sqrt(total variance)"}
    if ok.all():
        info.update(_neighbors_umap(adata, "X_multimodal", "multimodal", cfg.neighbors.n_neighbors, M.shape[1], cfg, "X_umap_multimodal", mcfg.umap and cfg.umap.enabled))
    else:
        sub = ad.AnnData(obs=adata.obs.loc[ok, []].copy())
        sub.obsm["X_multimodal"] = M[ok]
        info.update(_neighbors_umap(sub, "X_multimodal", "multimodal", cfg.neighbors.n_neighbors, M.shape[1], cfg, "X_umap_multimodal", mcfg.umap and cfg.umap.enabled))
        if "X_umap_multimodal" in sub.obsm:
            um = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
            um[ok] = sub.obsm["X_umap_multimodal"]
            adata.obsm["X_umap_multimodal"] = um
    logger.info("multimodal representation: concat_pcs (protein_weight=%.2f) -> X_multimodal", mcfg.protein_weight)
    return info


def modality_diagnostics(adata: ad.AnnData, cfg: Config) -> Dict[str, Any]:
    """Cheap, sampled checks: how much RNA and protein neighbourhoods agree, and how
    strongly the first PCs of each modality track sequencing/ADT depth."""
    from sklearn.neighbors import NearestNeighbors

    rng = np.random.default_rng(cfg.compute.seed)
    n = adata.n_obs
    k = cfg.neighbors.n_neighbors
    out: Dict[str, Any] = {"sampled_cells": None, "k": k}

    def _corr(Z: np.ndarray, y: np.ndarray, n_pcs: int = 5) -> Dict[str, float]:
        okm = np.isfinite(Z).all(axis=1) & np.isfinite(y)
        res = {}
        for i in range(min(n_pcs, Z.shape[1])):
            res[f"PC{i+1}"] = round(float(abs(np.corrcoef(Z[okm, i], y[okm])[0, 1])), 3)
        return res

    obs = adata.obs
    if "X_pca" in adata.obsm and "total_counts" in obs:
        out["rna_pc_vs_log_total_counts"] = _corr(np.asarray(adata.obsm["X_pca"]), np.log1p(obs["total_counts"].to_numpy(dtype=float)))
        if "pct_counts_mt" in obs:
            out["rna_pc_vs_pct_mt"] = _corr(np.asarray(adata.obsm["X_pca"]), obs["pct_counts_mt"].to_numpy(dtype=float))
    if "X_pca_protein" in adata.obsm and "protein_total_counts" in obs:
        out["protein_pc_vs_log_total_adt"] = _corr(np.asarray(adata.obsm["X_pca_protein"]), np.log1p(obs["protein_total_counts"].to_numpy(dtype=float)))
        if "protein_pct_isotype" in obs:
            out["protein_pc_vs_pct_isotype"] = _corr(np.asarray(adata.obsm["X_pca_protein"]), obs["protein_pct_isotype"].to_numpy(dtype=float))
    if "X_pca" in adata.obsm and "X_pca_protein" in adata.obsm:
        P = np.asarray(adata.obsm["X_pca_protein"])
        ok = np.where(np.isfinite(P).all(axis=1))[0]
        if ok.size > k + 1:
            idx = np.sort(rng.choice(ok, size=min(cfg.multimodal.diagnostic_cells, ok.size), replace=False))
            out["sampled_cells"] = int(idx.size)
            R = np.asarray(adata.obsm["X_pca"])[idx]
            Ps = P[idx]
            nr = NearestNeighbors(n_neighbors=k + 1).fit(R).kneighbors(R, return_distance=False)[:, 1:]
            npn = NearestNeighbors(n_neighbors=k + 1).fit(Ps).kneighbors(Ps, return_distance=False)[:, 1:]
            overlap = np.array([len(set(a) & set(b)) for a, b in zip(nr, npn)]) / k
            out["rna_protein_knn_overlap_mean"] = round(float(overlap.mean()), 4)
            out["rna_protein_knn_overlap_expected_random"] = round(float(k / max(1, idx.size - 1)), 4)
            if "X_multimodal" in adata.obsm:
                M = np.asarray(adata.obsm["X_multimodal"])[idx]
                nm = NearestNeighbors(n_neighbors=k + 1).fit(M).kneighbors(M, return_distance=False)[:, 1:]
                out["multimodal_knn_overlap_with_rna"] = round(float(np.mean([len(set(a) & set(b)) for a, b in zip(nm, nr)]) / k), 4)
                out["multimodal_knn_overlap_with_protein"] = round(float(np.mean([len(set(a) & set(b)) for a, b in zip(nm, npn)]) / k), 4)
    return out
