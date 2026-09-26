"""Perturbation distance vs control, DistanceTest, distance space, PCoA, phenotype modules,
master tables and distance-protein associations (reference ``distance.py`` / ``meta.py``
semantics) on synthetic data."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import base_config, make_dataset  # noqa: E402
from test_perturbation_effects import cfg as pe_cfg, make  # noqa: E402

from petrubseq_protein import analysis  # noqa: E402
from petrubseq_protein.analysis import distance as D  # noqa: E402
from petrubseq_protein.analysis import master_table as MT  # noqa: E402
from petrubseq_protein.config import Config  # noqa: E402


def dcfg(**over) -> Config:
    base = {"distance": {"enabled": True, "min_cells": 30, "n_permutations": 200}, "distance_space": {"enabled": True, "min_cells": 30, "n_modules": 3, "nearest_neighbors": 3}}
    for k, v in over.items():
        base[k] = {**base.get(k, {}), **v} if isinstance(v, dict) else v
    return pe_cfg(**base)


@pytest.fixture(scope="module")
def data():
    return make(ambiguous=30, rare=5)


@pytest.fixture(scope="module")
def res(data):
    cfg = dcfg()
    return analysis.run_perturbation_effects(data, cfg, analysis.run_perturbation_strength(data, cfg))


# ------------------------------------------------------------------ energy distance
def test_energy_distance_formula_and_properties():
    rng = np.random.default_rng(0)
    X, Y = rng.normal(size=(60, 8)), rng.normal(size=(70, 8)) + 0.5
    e = D.energy_distance(X, Y)
    manual = 2 * cdist(X, Y).mean() - cdist(X, X).mean() - cdist(Y, Y).mean()  # pertpy Edistance convention (diagonal zeros included)
    assert e == pytest.approx(manual, abs=1e-12)
    assert D.energy_distance(Y, X) == pytest.approx(e, abs=1e-12)
    assert D.energy_distance(X, X) == pytest.approx(0.0, abs=1e-12)
    Z = np.vstack([X, Y])
    assert D.energy_distance_from_cdist(cdist(Z, Z), np.arange(60), np.arange(60, 130)) == pytest.approx(e, abs=1e-10)
    far = rng.normal(size=(70, 8)) + 3.0
    assert D.energy_distance(X, far) > e > 0
    near = rng.normal(size=(250, 8))
    assert D.energy_distance(rng.normal(size=(200, 8)), near) < 0.15
    assert D.mmd_rbf(X, far) > D.mmd_rbf(X, rng.normal(size=(70, 8))) >= 0


def test_distance_test_permutation_conventions():
    rng = np.random.default_rng(1)
    X, Y = rng.normal(2, 1, (40, 6)), rng.normal(-2, 1, (90, 6))
    obs, p = D.distance_test_permutation(X, Y, n_permutations=100, seed=7)
    assert obs == pytest.approx(D.energy_distance(X, Y), abs=1e-10)
    assert p == pytest.approx(1 / 101)                       # (1 + 0) / (1 + B): the smallest attainable p
    obs2, p2 = D.distance_test_permutation(X, Y, n_permutations=100, seed=7)
    assert (obs2, p2) == (obs, p)                             # seeded: reproducible
    X0, Y0 = rng.normal(size=(40, 6)), rng.normal(size=(90, 6))
    _, p0 = D.distance_test_permutation(X0, Y0, n_permutations=200, seed=3)
    assert p0 > 0.05
    _, pn = D.distance_test_permutation(X0, Y0, n_permutations=0)
    assert np.isnan(pn)
    assert D.distance_test_permutation(X0[:0], Y0, 10) == (pytest.approx(float("nan"), nan_ok=True), pytest.approx(float("nan"), nan_ok=True))


def test_seed_derivation_and_sampling():
    assert D.derive_seed(123, "0_JAK2") == int(hashlib.sha256(b"123_0_JAK2").hexdigest()[:8], 16) % (2**31 - 1)
    idx = np.arange(5000)
    s1 = D.sample_cell_indices(idx, 1000, np.random.default_rng(123))
    assert len(s1) == 1000 and len(np.unique(s1)) == 1000 and np.all(np.diff(s1) > 0)
    np.testing.assert_array_equal(s1, D.sample_cell_indices(idx, 1000, np.random.default_rng(123)))
    np.testing.assert_array_equal(D.sample_cell_indices(np.arange(50), 100, np.random.default_rng(0)), np.arange(50))
    strata = np.array(["a"] * 800 + ["b"] * 200)
    s2 = D.sample_cell_indices(np.arange(1000), 200, np.random.default_rng(1), strata=strata)
    assert len(s2) == 200 and 140 <= (strata[s2] == "a").sum() <= 180


def test_pcoa_recovers_euclidean_configuration():
    rng = np.random.default_rng(2)
    P = rng.normal(size=(12, 3))
    Dm = cdist(P, P)
    coords, evals = D.compute_pcoa(Dm, n_components=10)
    assert coords.shape == (12, 3) and np.all(evals > 0) and np.all(np.isfinite(coords))
    np.testing.assert_allclose(cdist(coords, coords), Dm, atol=1e-8)
    c2, _ = D.compute_pcoa(Dm, n_components=10)
    np.testing.assert_array_equal(coords, c2)                 # deterministic


# ------------------------------------------------------------------ control distance
def test_control_distance_table(data, res):
    d = res.distance
    assert d is not None and not d.empty
    t = d.table
    assert list(t.columns) == ["target", "n_cells", "n_control", "energy_distance", "pvalue", "mmd_distance", "fdr", "significant"]
    assert list(t["energy_distance"]) == sorted(t["energy_distance"], reverse=True)
    assert set(t["target"]) == {"TA", "TB", "TC", "TD", "TE", "TF"} and list(d.skipped["target"]) == ["TR"]
    assert (t["n_control"] == 200).all() and d.control_used == "ntc"
    st = t.set_index("target")
    assert st.loc[["TA", "TC"], "energy_distance"].min() > st.loc[["TE", "TF"], "energy_distance"].max()
    assert st.loc[["TA", "TB", "TC", "TD"], "significant"].all() and not st.loc[["TE", "TF"], "significant"].any()
    assert (t["fdr"] >= t["pvalue"] - 1e-12).all() and t["fdr"].between(0, 1).all()
    assert d.info["n_permutations"] == 200 and d.info["status"] == "ok"


def test_control_distance_deterministic_and_parallel_invariant(data):
    a = analysis.run_perturbation_effects(data, dcfg(modules={"enabled": False}, ps={"enabled": False}, lochness={"enabled": False}, protein={"enabled": False}, concordance={"enabled": False}))
    cfg2 = dcfg(modules={"enabled": False}, ps={"enabled": False}, lochness={"enabled": False}, protein={"enabled": False}, concordance={"enabled": False})
    cfg2.compute.n_jobs = 2
    b = analysis.run_perturbation_effects(data, cfg2)
    pd.testing.assert_frame_equal(a.distance.table, b.distance.table)
    pd.testing.assert_frame_equal(a.distance_space.distance_matrix, b.distance_space.distance_matrix)


def test_control_distance_edge_cases(data):
    no_perm = analysis.run_perturbation_effects(data, dcfg(distance={"n_permutations": 0}, modules={"enabled": False}, ps={"enabled": False}, lochness={"enabled": False})).distance
    assert no_perm.table["pvalue"].isna().all() and no_perm.table["fdr"].isna().all() and not no_perm.table["significant"].any()
    d2 = data.copy()
    d2.obs["control_class"] = ""
    d2.obs.loc[d2.obs["perturbation_class"] == "single_control", "perturbation_class"] = "ambiguous"
    r = D.compute_perturbation_distance(d2, dcfg())
    assert not r.empty and r.control_used == "other"          # reference fallback: all targeting cells serve as the control
    d2.obs["perturbation_class"] = "ambiguous"
    r = D.compute_perturbation_distance(d2, dcfg())
    assert r.empty and "no control" in r.note
    d3 = data.copy()
    del d3.obsm["X_pca"]
    assert D.compute_perturbation_distance(d3, dcfg()).empty and D.compute_distance_space(d3, dcfg()).empty
    high = dcfg(distance={"min_cells": 100}, distance_space={"min_cells": 100})
    assert D.compute_perturbation_distance(data, high).empty and D.compute_distance_space(data, high).note.startswith("fewer than 2")


# ------------------------------------------------------------------ distance space
def test_distance_space_matrix_neighbors_modules(res):
    s = res.distance_space
    assert s is not None and not s.empty
    M = s.distance_matrix
    assert list(M.index) == list(M.columns) == ["TA", "TB", "TC", "TD", "TE", "TF"]
    np.testing.assert_allclose(M.to_numpy(), M.to_numpy().T, atol=1e-12)
    np.testing.assert_allclose(np.diag(M.to_numpy()), 0.0, atol=1e-12)
    assert (M.to_numpy()[~np.eye(6, dtype=bool)] > 0).all()
    nb = s.neighbors
    assert list(nb.columns) == ["target", "neighbor", "distance", "rank"] and (nb["target"] != nb["neighbor"]).all() and set(nb["rank"]) == {1, 2, 3}
    first = nb[nb["rank"] == 1].set_index("target")["neighbor"]
    assert first["TA"] == "TB" and first["TB"] == "TA" and first["TC"] == "TD" and first["TD"] == "TC"
    for t, g in nb.groupby("target"):
        assert list(g.sort_values("rank")["distance"]) == sorted(g["distance"])
    pm = s.phenotype_modules.set_index("target")["phenotype_module"]
    assert pm["TA"] == pm["TB"] and pm["TC"] == pm["TD"] and pm["TA"] != pm["TC"] and pm.nunique() == 3
    assert s.coordinates.shape[0] == 6 and "PCoA1" in s.coordinates.columns and np.isfinite(s.coordinates.iloc[:, 1:].to_numpy(float)).all()
    assert np.all(np.diff(s.eigenvalues) <= 1e-12)


def test_distance_space_default_module_rule(data):
    r = D.compute_distance_space(data, dcfg(distance_space={"n_modules": None}))
    assert r.info["module_rule"].startswith("reference default") and r.info["n_modules_cut"] == max(2, min(9, 6)) and r.phenotype_modules["phenotype_module"].nunique() == 6
    r2 = D.compute_distance_space(data, dcfg(distance_space={"n_modules": None, "cluster_distance_threshold": 1e9}))
    assert r2.phenotype_modules["phenotype_module"].nunique() == 1 and r2.info["module_rule"] == "cluster_distance_threshold"


# ------------------------------------------------------------------ master tables and protein associations
def test_master_tables(res):
    m = res.master
    assert m is not None and not m.empty
    t = m.perturbation
    for c in ("target", "n_cells", "target_log2fc", "target_fdr", "is_effective_hit", "ps_mean", "ps_median", "ps_responder_fraction", "lochness_mean", "lochness_peak", "energy_distance", "mmd_distance", "distance_pvalue", "distance_fdr", "distance_significant", "cofunctional_module", "phenotype_module"):
        assert c in t.columns, c
    assert not any(("master_score" in c or "combined_score" in c) for c in t.columns)
    ed = t["energy_distance"].to_numpy(float)
    assert list(ed[np.isfinite(ed)]) == sorted(ed[np.isfinite(ed)], reverse=True)
    assert "TR" not in set(t["target"]) and len(t) == 6           # TR (5 cells) enters no target-level table -> absent from the union
    assert "strongest_protein" in t.columns and "protein_effect_P_up" in t.columns and "phenotype_nearest_neighbor" in t.columns
    assert t.set_index("target").loc["TA", "phenotype_nearest_neighbor"] == "TB"
    p = m.protein
    assert not p.empty and set(p["protein_effect_level"]) == {"TARGET_LEVEL"} and set(p["ps_protein_level"]) == {"CELL_LEVEL"} and set(p["target_columns_level"]) == {"TARGET_LEVEL"}
    assert len(p) == p[["target", "protein"]].drop_duplicates().shape[0] and "energy_distance" in p.columns and "phenotype_module" in p.columns
    dp = m.distance_protein
    assert list(dp["protein"]) == ["P_up", "P_down", "P_null"] and {"rho_signed", "fdr_signed", "rho_abs", "fdr_abs", "status_signed", "status_abs"} <= set(dp.columns)
    assert (dp["n_targets"] == 6).all() and (dp["support"] == "few_targets(<10)").all()
    mantel = m.phenotype_protein_mantel
    assert list(mantel["protein"]) == ["all_proteins", "P_up", "P_down", "P_null"] and mantel["p_permutation"].between(0, 1).all()
    assert not m.phenotype_module_means.empty and set(m.phenotype_module_means.columns) >= {"phenotype_module", "protein", "n_targets", "mean_effect"}
    assert m.info["power_note"] and m.info["composite_score"] == "none (reference principle)"


def test_master_table_without_distance(data):
    r = analysis.run_perturbation_effects(data, pe_cfg())
    assert r.distance is None and r.distance_space is None and r.master is not None
    t = r.master.perturbation
    assert "energy_distance" not in t.columns and list(t["target"]) == sorted(t["target"]) and "ps_median" in t.columns
    assert r.master.distance_protein.empty and r.master.phenotype_protein_mantel.empty


def test_attach_and_tables(data, res):
    a = data.copy()
    analysis.attach(a, res)
    assert "perturbation_distance" in a.uns and "perturbation_distance_matrix" in a.uns and "phenotype_modules" in a.uns and "master_perturbation_table" in a.uns
    assert a.uns["perturbation_distance_matrix"].shape == (6, 6) and "perturbation_distance_matrix" not in a.obsm
    tb = analysis.tables(res)
    for k in ("perturbation_distance/perturbation_distance", "perturbation_distance/distance_skipped", "perturbation_distance/perturbation_distance_matrix", "perturbation_distance/perturbation_space_coordinates",
              "perturbation_distance/perturbation_neighbors", "perturbation_distance/phenotype_modules", "master_perturbation_table", "master_perturbation_protein_table", "perturbation_distance/distance_protein_association", "perturbation_distance/phenotype_space_protein_mantel"):
        assert k in tb, k


# ------------------------------------------------------------------ pipeline integration
def test_pipeline_end_to_end_with_distance(tmp_path_factory, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    d = make_dataset(tmp_path_factory.mktemp("dist_synth"), n_cells=400, rna_state="raw")
    pe = {"enabled": True, "modules": {"min_cells_per_perturbation": 10, "min_perturbations": 3, "n_programs": 2}, "ps": {"min_cells_per_target": 10, "compute_lda_umap": False},
          "distance": {"enabled": True, "min_cells": 10, "n_permutations": 50}, "distance_space": {"enabled": True, "min_cells": 10, "nearest_neighbors": 2}, "master_table": {"mantel_permutations": 50}}
    r = run_pipeline(Config.from_dict(base_config(d, tmp_path, analysis={"perturbation_effects": pe, "clustering": {"enabled": False}})))
    tdir = tmp_path / "tables"
    assert (tdir / "perturbation_distance" / "perturbation_distance.csv").is_file() and (tdir / "perturbation_distance" / "perturbation_distance_matrix.csv").is_file()
    assert (tdir / "master_perturbation_table.csv").is_file() and (tdir / "master_perturbation_protein_table.csv").is_file()
    master = pd.read_csv(tdir / "master_perturbation_table.csv")
    assert "energy_distance" in master.columns and master.columns[0] == "target"
    a = ad.read_h5ad(next((tmp_path / "processed").glob("*.h5ad")))
    assert "perturbation_distance" in a.uns and "perturbation_distance_matrix" in a.uns and "master_perturbation_table" in a.uns
    info = a.uns["petrubseq_protein"]["analysis"]["perturbation_effects"]
    assert info["distance"]["n_targets_tested"] >= 2 and "n_modules" in info["distance_space"]
    man = pd.read_csv(tdir / "figure_manifest.csv")
    names = set(man.loc[man["section"] == "perturbation_distance", "name"])
    assert {"perturbation_distance_ranking", "perturbation_atlas", "ps_vs_distance_map", "perturbation_phenotype_space", "perturbation_distance_matrix"} <= names
    html = r.report_html.read_text()
    assert "Perturbation distance vs control" in html and "phenotype modules" in html and "Master perturbation table" in html and "Stage G" not in html
    md = (tmp_path / "report.md").read_text()
    assert "Perturbation distance vs control" in md and "Master perturbation table" in md
