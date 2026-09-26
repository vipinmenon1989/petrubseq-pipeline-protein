"""Protein representation / PCA numerical reproducibility (float64 working path)."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import base_config, make_dataset  # noqa: E402

from petrubseq_protein.config import Config  # noqa: E402
from petrubseq_protein.preprocessing import embeddings  # noqa: E402


def _multimodal(seed: int = 0, n: int = 600, p: int = 12) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    counts = rng.negative_binomial(3, 0.05, (n, p)).astype(float)
    counts[: n // 3, :4] *= 4  # a structured population
    L = np.log1p(counts)
    clr = (L - L.mean(axis=1, keepdims=True)).astype(np.float32)  # the pipeline's CLR (float32 in obsm)
    names = [f"P{i}" for i in range(p - 2)] + ["Mouse_IgG1", "Mouse_IgG2a"]
    a = ad.AnnData(X=sp.csr_matrix(rng.poisson(1.0, (n, 30)).astype(np.float32)))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obsm["protein"] = pd.DataFrame(clr, index=a.obs_names, columns=names)
    a.uns["protein_features"] = pd.DataFrame({"is_isotype": [nm.startswith("Mouse") for nm in names]}, index=names)
    return a


def cfg(**over) -> Config:
    d = {"inputs": {"rna": {"file": "x"}, "protein": {"file": "y"}}, "protein": {"n_pcs": 5}, "umap": {"enabled": True, "min_features": 5}, "compute": {"seed": 0}}
    for k, v in over.items():
        d[k] = {**d.get(k, {}), **v}
    return Config.from_dict(d)


def test_pca_on_runs_in_float64_and_is_deterministic():
    rng = np.random.default_rng(1)
    X32 = rng.normal(size=(400, 10)).astype(np.float32)
    r1 = embeddings._pca_on(X32, 5, True, 10.0, 0)
    r2 = embeddings._pca_on(X32.copy(), 5, True, 10.0, 0)
    r3 = embeddings._pca_on(X32.astype(np.float64), 5, True, 10.0, 0)
    assert r1["working_dtype"] == "float64" and r1["X_pca"].dtype == np.float32
    for k in ("X_pca", "loadings", "variance", "variance_ratio"):
        np.testing.assert_array_equal(r1[k], r2[k])
        np.testing.assert_array_equal(r1[k], r3[k])  # float32 input is cast to float64 before scaling/PCA
    # sparse input takes the same float64 route
    r4 = embeddings._pca_on(sp.csr_matrix(X32), 5, True, 10.0, 0)
    np.testing.assert_array_equal(r1["X_pca"], r4["X_pca"])


def test_protein_embedding_is_reproducible_and_records_dtype():
    a1, a2 = _multimodal(), _multimodal()
    i1 = embeddings.protein_embedding(a1, cfg(), "protein")
    i2 = embeddings.protein_embedding(a2, cfg(), "protein")
    np.testing.assert_array_equal(a1.obsm["X_pca_protein"], a2.obsm["X_pca_protein"])
    np.testing.assert_array_equal(a1.uns["pca_protein"]["variance"], a2.uns["pca_protein"]["variance"])
    np.testing.assert_array_equal(a1.uns["pca_protein"]["variance_ratio"], a2.uns["pca_protein"]["variance_ratio"])
    pd.testing.assert_frame_equal(a1.uns["pca_protein"]["loadings"], a2.uns["pca_protein"]["loadings"])
    assert a1.uns["pca_protein"]["params"]["working_dtype"] == "float64" and i1["pca"]["working_dtype"] == "float64"
    assert a1.obsm["X_pca_protein"].dtype == np.float32 and a1.uns["pca_protein"]["params"]["stored_dtype"] == "float32"
    assert np.isfinite(a1.obsm["X_pca_protein"]).all() and a1.obsm["X_pca_protein"].shape == (a1.n_obs, 5)
    # isotypes excluded from the embedding
    assert not set(a1.uns["pca_protein"]["loadings"].index) & {"Mouse_IgG1", "Mouse_IgG2a"}
    # neighbour graph and UMAP built on identical PCs are identical (same seed, same process settings)
    d1, d2 = a1.obsp["protein_distances"], a2.obsp["protein_distances"]
    assert (d1 != d2).nnz == 0
    assert np.abs(a1.obsp["protein_connectivities"] - a2.obsp["protein_connectivities"]).max() == 0
    assert i1["umap"] is not None and "X_umap_protein" in a1.obsm
    np.testing.assert_allclose(a1.obsm["X_umap_protein"], a2.obsm["X_umap_protein"], atol=1e-5)


def test_protein_pca_reproducible_through_the_pipeline(tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    d = make_dataset(tmp_path / "synth", n_cells=300, rna_state="raw")
    an = {"perturbation_effects": {"enabled": False}, "clustering": {"enabled": False}}
    r1 = run_pipeline(Config.from_dict(base_config(d, tmp_path / "run1", analysis=an))).adata
    r2 = run_pipeline(Config.from_dict(base_config(d, tmp_path / "run2", analysis=an))).adata
    for key in ("X_pca_protein", "X_pca"):
        np.testing.assert_array_equal(np.nan_to_num(r1.obsm[key]), np.nan_to_num(r2.obsm[key]))
    pd.testing.assert_frame_equal(r1.obsm["protein"], r2.obsm["protein"])
    pd.testing.assert_frame_equal(r1.uns["pca_protein"]["loadings"], r2.uns["pca_protein"]["loadings"])
    assert r1.uns["pca_protein"]["params"]["working_dtype"] == "float64" and r1.uns["pca"]["params"]["working_dtype"] == "float64"
    if "X_umap_protein" in r1.obsm:
        np.testing.assert_allclose(np.nan_to_num(r1.obsm["X_umap_protein"]), np.nan_to_num(r2.obsm["X_umap_protein"]), atol=1e-5)
