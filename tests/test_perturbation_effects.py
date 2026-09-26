"""Stage E: PS, lochNESS, gene programs / perturbation modules, protein effects,
RNA-protein concordance (synthetic data only)."""

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

from petrubseq_protein import analysis  # noqa: E402
from petrubseq_protein.analysis import lochness as L, modules as M, protein_effects as PE, ps_score as PS  # noqa: E402
from petrubseq_protein.analysis.concordance import compute_concordance  # noqa: E402
from petrubseq_protein.config import Config  # noqa: E402

TARGETS = ["TA", "TB", "TC", "TD", "TE", "TF"]


def make(seed: int = 0, n_ctrl: int = 200, n_per: int = 80, rare: int = 0, ambiguous: int = 0, protein_link: bool = True) -> ad.AnnData:
    """Log-normalized object with known structure.

    TA, TB share response block 0 (genes 10-29) -> same perturbation module; TC, TD
    share block 1 (genes 30-49); TE, TF have no effect. The response of each TA/TC
    cell scales with a per-cell strength w ~ U(0, 1); protein P_up = 2 w + noise in
    TA cells (within-target PS <-> protein link), P_down lower in TC cells, P_null noise.
    """
    rng = np.random.default_rng(seed)
    lab = ["NT"] * n_ctrl + sum([[t] * n_per for t in TARGETS], []) + ["TR"] * rare
    lab = np.array(lab + ["AMB"] * ambiguous)
    n, g = len(lab), 200
    w = rng.uniform(0, 1, n)
    C = rng.poisson(2.0, (n, g)).astype(float)
    for t, block, scale in (("TA", 0, 1.0), ("TB", 0, 1.0), ("TC", 1, 1.0), ("TD", 1, 1.0)):
        m = lab == t
        cols = np.arange(10 + 20 * block, 30 + 20 * block)
        C[np.ix_(m, cols)] += rng.poisson(6.0 * (w[m, None] if t in ("TA", "TC") else 0.8) * scale, (m.sum(), 20))
    amb = lab == "AMB"
    C[np.ix_(amb, np.arange(10, 50))] += 40  # extreme ambiguous cells must not enter any group
    X = np.log1p(C / C.sum(1, keepdims=True) * 1e4).astype(np.float32)
    genes = TARGETS + [f"G{i}" for i in range(len(TARGETS), g)]
    a = ad.AnnData(X=sp.csr_matrix(X), var=pd.DataFrame(index=genes))
    a.obs_names = [f"c{i}" for i in range(n)]
    klass = np.where(lab == "NT", "single_control", np.where(lab == "AMB", "ambiguous", "single_targeting"))
    a.obs["target"] = np.where(lab == "AMB", "", lab)
    a.obs["perturbation_class"] = klass
    a.obs["control_class"] = np.where(lab == "NT", "non_targeting", "")
    a.obs["sample"] = rng.choice(["s1", "s2"], n)
    a.var["highly_variable"] = True
    u, s, _ = np.linalg.svd(X - X.mean(0), full_matrices=False)
    a.obsm["X_pca"] = (u[:, :20] * s[:20]).astype(np.float32)
    up = rng.normal(0, 0.5, n) + np.where(lab == "TA", 2 * w, 0) if protein_link else rng.normal(0, 0.5, n)
    down = rng.normal(0, 0.5, n) - np.where(lab == "TC", 1.5, 0)
    a.obsm["protein"] = pd.DataFrame({"P_up": up, "P_down": down, "P_null": rng.normal(0, 0.5, n)}, index=a.obs_names)
    a.uns["_w"] = w
    return a


def cfg(**pe) -> Config:
    base = {"enabled": True, "lochness": {"n_permutations": 50, "n_neighbors": 50, "max_k_fraction": 1.0}, "modules": {"n_programs": 2, "n_modules": 3, "min_cells_per_perturbation": 20}}
    for k, v in pe.items():
        base[k] = {**base.get(k, {}), **v} if isinstance(v, dict) else v
    return Config.from_dict({"inputs": {"rna": {"file": "x"}}, "analysis": {"perturbation_effects": base}})


@pytest.fixture(scope="module")
def data():
    return make(ambiguous=30, rare=5)


@pytest.fixture(scope="module")
def res(data):
    return analysis.run_perturbation_effects(data, cfg())


# ------------------------------------------------------------------ PS
def test_ps_strong_vs_null(res):
    s = res.ps.summary.set_index("target")
    for t in ("TA", "TB", "TC", "TD"):
        assert s.loc[t, "mean_ps"] > s.loc["TE", "mean_ps"] + 0.1, t
    assert s.loc["TB", "auc_vs_control"] > 0.9 and s.loc["TE", "auc_vs_control"] < s.loc["TB", "auc_vs_control"] - 0.1
    assert s.loc["TB", "d_vs_control"] > 2 * s.loc["TE", "d_vs_control"]
    assert s.loc["TB", "pct_high_ps"] > 5 * max(s.loc["TB", "pct_high_ps_control"], 1)


def test_ps_tracks_per_cell_strength(data, res):
    m = (data.obs["target"] == "TA").to_numpy()
    r = pd.Series(res.ps.own[m].to_numpy()).corr(pd.Series(data.uns["_w"][m]), method="spearman")
    assert r > 0.6


def test_ps_controls_and_exclusions(data, res):
    amb = (data.obs["perturbation_class"] == "ambiguous").to_numpy()
    assert res.ps.scores.loc[amb].isna().all().all() and res.ps.own[amb].isna().all()
    ctrl = (data.obs["perturbation_class"] == "single_control").to_numpy()
    assert abs(res.ps.raw.loc[ctrl, "TA"].mean()) < 0.3  # controls centred near 0
    skipped = set(res.ps.skipped["target"])
    assert "TR" in skipped  # 5 cells < min 10
    # ambiguous extremes do not change the scores: same object without them gives identical PS
    clean = analysis.run_perturbation_effects(data[~amb].copy(), cfg(lochness={"enabled": False}, modules={"enabled": False}))
    pd.testing.assert_series_equal(res.ps.summary.set_index("target")["mean_ps"].sort_index(), clean.ps.summary.set_index("target")["mean_ps"].sort_index(), atol=1e-9)


def test_ps_deterministic(data, res):
    again = PS.compute_ps(data, cfg())
    pd.testing.assert_frame_equal(again.scores, res.ps.scores)


def test_ps_quadrants_use_target_gene(res):
    s = res.ps.summary.set_index("target")
    assert {"pct_successful_kd", "pct_escaper", "net_pct_kd"} <= set(s.columns)
    assert np.isfinite(s.loc["TA", "pct_successful_kd"])


def test_ps_insufficient_controls():
    a = make()
    a.obs["control_class"] = ""
    r = PS.compute_ps(a, cfg())
    assert r.summary.empty and (r.skipped["target"] == "*").any()


# ------------------------------------------------------------ lochNESS
def test_lochness_enrichment_and_null(res):
    s = res.lochness.summary.set_index("target")
    assert s.loc["TB", "mean_lochness_in_own_cells"] > 0.5 and s.loc["TB", "fdr"] < 0.05
    # a null target sits where the controls sit: its own-cell score equals its score in control cells
    assert abs(s.loc["TE", "delta_own_vs_control"]) < 0.2 and s.loc["TB", "delta_own_vs_control"] > 1.0
    assert res.lochness.info["k_used"] == 50 and res.lochness.info["use_rep"] == "X_pca"


def test_lochness_depletion(data, res):
    # TB cells sit in their own region: control cells far from it see TB depleted
    ctrl = (data.obs["perturbation_class"] == "single_control").to_numpy()
    assert np.nanmean(res.lochness.scores.loc[ctrl, "TB"]) < 0
    assert np.nanmean(res.lochness.self_score[ctrl]) > -1  # control label self score defined


def test_lochness_rare_and_deterministic(data, res):
    assert "TR" in set(res.lochness.skipped["target"])
    again = L.compute_lochness(data, cfg())
    pd.testing.assert_frame_equal(again.summary, res.lochness.summary)


def test_lochness_k_cap():
    a = make(n_per=30, n_ctrl=60)
    r = L.compute_lochness(a, cfg(lochness={"n_neighbors": 300, "max_k_fraction": 0.1, "n_permutations": 0}))
    assert r.info["k_capped"] and r.info["k_used"] == max(15, int(0.1 * a.n_obs))


def test_lochness_formula():
    adj = sp.csr_matrix(np.array([[0, 1, 1, 0], [1, 0, 1, 0], [1, 1, 0, 0], [0, 0, 1, 0]], float))
    counts = np.asarray(adj.sum(1)).ravel()
    s = L.lochness_score(adj, counts, np.array([1, 1, 0, 0.0]), 0.5)
    np.testing.assert_allclose(s, [0.0, 0.0, 1.0, -1.0])


# ------------------------------------------------------------- modules
def test_modules_recover_structure(res):
    mods = res.modules.perturbation_modules.set_index("target")["module"]
    assert mods["TA"] == mods["TB"] and mods["TC"] == mods["TD"] and mods["TA"] != mods["TC"]
    gp = res.modules.gene_programs.set_index("gene")["program"]
    b0 = {gp[f"G{i}"] for i in range(10, 30) if f"G{i}" in gp.index}
    b1 = {gp[f"G{i}"] for i in range(30, 50) if f"G{i}" in gp.index}
    assert len(b0) == 1 and len(b1) == 1 and b0 != b1  # each response block is one program
    ppe = res.modules.perturbation_program_effects.set_index(["target", "program"])["mean_log2fc"]
    assert ppe[("TA", b0.pop())] > 0.5
    assert res.modules.program_activity.shape[1] == len(res.modules.program_labels)


def test_modules_deterministic_and_insufficient(data, res):
    again = M.compute_modules(data, cfg())
    pd.testing.assert_frame_equal(again.effect, res.modules.effect)
    assert list(again.perturbation_modules["module"]) == list(res.modules.perturbation_modules["module"])
    few = M.compute_modules(data, cfg(modules={"min_perturbations": 50}))
    assert few.empty and "skipped" in few.info["status"]


# ------------------------------------------------------------- protein
def test_protein_effects_signs(res):
    t = res.protein.table.set_index(["target", "protein"])
    assert t.loc[("TA", "P_up"), "effect"] > 0.5 and t.loc[("TA", "P_up"), "fdr"] < 0.05
    assert t.loc[("TC", "P_down"), "effect"] < -1 and t.loc[("TC", "P_down"), "fdr"] < 0.05
    assert t.loc[("TE", "P_null"), "fdr"] > 0.05
    assert {"n_perturbed", "n_control", "mean_perturbed", "mean_control", "median_perturbed", "cohen_d", "p_value", "n_samples_tested"} <= set(res.protein.table.columns)


def test_protein_missing_and_normalized_only():
    a = make()
    del a.obsm["protein"]
    r = PE.compute_protein_effects(a, cfg())
    assert r.empty and "no normalized protein" in r.info["status"]
    b = make()
    assert "protein_counts" not in b.obsm  # normalized-only input
    r2 = PE.compute_protein_effects(b, cfg())
    assert not r2.empty and r2.info["representation"] == "protein"


# ---------------------------------------------------------- concordance
def test_ps_protein_within_target(res):
    pp = res.concordance.ps_protein
    w = pp[pp["scope"] == "within_target"].set_index(["target", "protein"])
    assert w.loc[("TA", "P_up"), "spearman_rho"] > 0.4 and w.loc[("TA", "P_up"), "status"] == "significant"
    assert w.loc[("TA", "P_null"), "status"] != "significant"
    assert (pp["scope"] == "pooled_all_targets").any() and (pp.loc[pp["scope"] == "pooled_all_targets", "status"] == "descriptive_pooled").all()


def test_ps_protein_null_when_unlinked():
    a = make(seed=3, protein_link=False)
    r = analysis.run_perturbation_effects(a, cfg(lochness={"enabled": False}, modules={"enabled": False}))
    w = r.concordance.ps_protein
    w = w[(w["scope"] == "within_target") & (w["target"] == "TA") & (w["protein"] == "P_up")]
    assert w["status"].iloc[0] != "significant"


def test_program_protein_and_lochness_levels(res):
    ppc = res.concordance.program_protein_cells
    assert {"gene_program", "protein", "n_cells", "spearman_rho", "fdr"} <= set(ppc.columns)
    gp = res.modules.gene_programs.set_index("gene")["program"]
    prog1 = gp["G35"]  # block 1 program: up in TC/TD, P_down lower in TC
    r = ppc.set_index(["gene_program", "protein"]).loc[(prog1, "P_down"), "spearman_rho"]
    assert r < -0.1
    lps = res.concordance.lochness_protein_summary
    assert set(lps["protein"]) == {"P_up", "P_down", "P_null"} and (lps["n_targets"] == 6).all()  # target level, never per cell
    sm = res.concordance.summary.set_index("target")
    assert {"ps_median", "lochness_own_mean", "module", "rna_effect_magnitude", "protein_effect_magnitude", "strongest_programs"} <= set(sm.columns)


def test_attach_and_tables(data, res):
    a = data.copy()
    analysis.attach(a, res)
    assert {"ps_scores", "ps_scores_raw", "lochness", "program_activity"} <= set(a.obsm) and "ps_score" in a.obs
    t = analysis.tables(res)
    assert "perturbation_effects/perturbation_summary" in t and "perturbation_effects/gene_programs" in t


# ------------------------------------------------- pipeline integration
@pytest.fixture(scope="module")
def synth_raw(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("pe_synth"), n_cells=400, rna_state="raw")


def test_pipeline_disabled_adds_nothing(synth_raw, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    a = run_pipeline(Config.from_dict(base_config(synth_raw, tmp_path, umap={"enabled": False}, analysis={"clustering": {"enabled": False}}))).adata
    assert "analysis" not in a.uns["petrubseq_protein"] and "ps_scores" not in a.obsm and "ps_score" not in a.obs
    assert not (tmp_path / "tables" / "perturbation_effects").exists()


def test_pipeline_enabled_end_to_end(synth_raw, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    pe = {"enabled": True, "lochness": {"n_permutations": 20}, "modules": {"min_cells_per_perturbation": 10, "min_perturbations": 3, "n_programs": 2}, "ps": {"min_cells_per_target": 10}}
    r = run_pipeline(Config.from_dict(base_config(synth_raw, tmp_path, analysis={"perturbation_effects": pe, "clustering": {"enabled": False}})))
    a = ad.read_h5ad(next((tmp_path / 'processed').glob('*.h5ad')))
    info = a.uns["petrubseq_protein"]["analysis"]["perturbation_effects"]
    assert info["ps"]["n_targets_scored"] >= 1 and "k_used" in info["lochness"]
    assert "ps_scores" in a.obsm and "lochness" in a.obsm and "ps_score" in a.obs
    tdir = tmp_path / "tables" / "perturbation_effects"
    assert (tdir / "ps_targets.csv").is_file() and (tdir / "protein_effects.csv").is_file() and (tdir / "perturbation_summary.csv").is_file()
    man = pd.read_csv(tmp_path / "tables" / "figure_manifest.csv")
    assert (man["section"] == "perturbation_effects").sum() >= 4
    html = r.report_html.read_text()
    assert "6. Perturbation effects" in html and "7. Outputs and provenance" in html
