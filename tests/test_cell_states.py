"""Stage F: Leiden cell-state clustering and perturbation x cluster enrichment."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import base_config, make_dataset  # noqa: E402

from petrubseq_protein import analysis  # noqa: E402
from petrubseq_protein.analysis._common import bh_fdr  # noqa: E402
from petrubseq_protein.analysis.cluster_enrichment import compute_cluster_enrichment, haldane_log2_or  # noqa: E402
from petrubseq_protein.analysis.clustering import compute_clustering  # noqa: E402
from petrubseq_protein.config import Config, ConfigError  # noqa: E402


def cfg(**cl) -> Config:
    base = {"enabled": True, "resolution": 0.3, "enrichment": {"n_permutations": 50, "min_cells_per_cluster": 20}}
    for k, v in cl.items():
        base[k] = {**base[k], **v} if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return Config.from_dict({"inputs": {"rna": {"file": "x"}}, "analysis": {"clustering": base}})


def labelled(seed: int = 0) -> ad.AnnData:
    """Known cluster labels (0, 1, 2) with a designed perturbation composition.

    controls (300): 1/3 in each cluster.  EN (90): 70 % in cluster 0 (enriched),
    guides EN_1, EN_2 concentrated, EN_3 spread.  DE (90): none in cluster 1
    (depleted).  NU (90): like controls.  60 ambiguous and 40 multi-guide cells
    all sit in cluster 1 and must not change any test.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(300):
        rows.append(("NT", "single_control", "non_targeting", "NT_1", str(i % 3)))
    for i in range(90):
        g = f"EN_{i % 3 + 1}"
        k = "0" if (g != "EN_3" and rng.random() < 0.95) or (g == "EN_3" and rng.random() < 0.2) else str(rng.integers(1, 3))
        rows.append(("EN", "single_targeting", "", g, k))
    for i in range(90):
        rows.append(("DE", "single_targeting", "", f"DE_{i % 3 + 1}", str(rng.choice(["0", "2"]))))
    for i in range(90):
        rows.append(("NU", "single_targeting", "", f"NU_{i % 3 + 1}", str(i % 3)))
    for i in range(60):
        rows.append(("", "ambiguous", "", "", "1"))
    for i in range(40):
        rows.append(("EN;DE", "multi_targeting", "", "", "1"))
    df = pd.DataFrame(rows, columns=["target", "perturbation_class", "control_class", "guide", "leiden"])
    df["sample"] = rng.choice(["s1", "s2"], len(df))
    df.index = [f"c{i}" for i in range(len(df))]
    a = ad.AnnData(obs=df)
    a.obs["leiden"] = pd.Categorical(a.obs["leiden"], categories=["0", "1", "2"])
    return a


# ------------------------------------------------------------- config
def test_config_defaults_and_validation():
    c = Config.from_dict({"inputs": {"rna": {"file": "x"}}})
    cl = c.analysis.clustering
    assert cl.enabled is False and cl.key == "leiden" and cl.resolution == 1.0 and cl.n_iterations == 2
    assert cl.enrichment.control == "non_targeting" and cl.enrichment.stratify_by is None
    for bad in ({"resolution": -1}, {"n_iterations": 0}, {"key": " "}, {"enrichment": {"control": "everything"}}, {"enrichment": {"fdr_alpha": 1.5}}):
        with pytest.raises(ConfigError):
            Config.from_dict({"inputs": {"rna": {"file": "x"}}, "analysis": {"clustering": bad}})


# ---------------------------------------------------------- clustering
@pytest.fixture(scope="module")
def blobs():
    a = sc.datasets.blobs(n_variables=20, n_centers=3, n_observations=450, random_state=0)
    a.obs_names = [f"b{i}" for i in range(a.n_obs)]
    sc.pp.pca(a, n_comps=10, random_state=0)
    sc.pp.neighbors(a, key_added="rna", random_state=0)
    a.obs["perturbation_class"] = "single_control"
    a.obs["control_class"] = "non_targeting"
    a.obs["target"] = "NT"
    a.obs["sample"] = np.where(np.arange(a.n_obs) % 2, "s1", "s2")
    return a


def test_leiden_deterministic_and_categorical(blobs):
    a1, a2 = blobs.copy(), blobs.copy()
    r1, r2 = compute_clustering(a1, cfg()), compute_clustering(a2, cfg())
    assert list(a1.obs["leiden"]) == list(a2.obs["leiden"]) and a1.obs["leiden"].dtype.name == "category"
    assert r1.info["n_clusters"] >= 3 and r1.info["resolution"] == 0.3 and r1.info["random_state"] == 0
    truth = pd.crosstab(a1.obs["blobs"], a1.obs["leiden"]).gt(0).sum(axis=0)
    assert (truth == 1).all()  # no cluster mixes two blobs
    assert int(r1.summary["n_cells"].sum()) == a1.n_obs and {"frac_ambiguous_or_unassigned", "top_sample", "flags"} <= set(r1.summary.columns)
    assert "sample" in r1.by_design


def test_resolution_changes_granularity(blobs):
    lo, hi = blobs.copy(), blobs.copy()
    compute_clustering(lo, cfg(resolution=0.2))
    compute_clustering(hi, cfg(resolution=3.0))
    assert hi.obs["leiden"].nunique() > lo.obs["leiden"].nunique()


def test_key_collision_and_missing_graph(blobs):
    a = blobs.copy()
    a.obs["leiden"] = "user"
    with pytest.raises(ValueError, match="already an obs column"):
        compute_clustering(a, cfg())
    compute_clustering(a, cfg(overwrite=True))
    assert a.obs["leiden"].nunique() > 1
    b = blobs.copy()
    del b.obsp["rna_connectivities"]
    with pytest.raises(ValueError, match="neighbour graph"):
        compute_clustering(b, cfg())


def test_dominated_cluster_flagged(blobs):
    a = blobs.copy()
    a.obs["sample"] = np.where(a.obs["blobs"].astype(str) == "0", "s_only", np.where(np.arange(a.n_obs) % 2, "s1", "s2"))
    r = compute_clustering(a, cfg())
    assert r.flags and any("sample s_only" in f for f in r.flags)


# ---------------------------------------------------------- enrichment
@pytest.fixture(scope="module")
def enr():
    a = labelled()
    return a, compute_cluster_enrichment(a, cfg(), "leiden")


def test_fisher_table_orientation(enr):
    a, r = enr
    row = r.table.set_index(["target", "cluster"]).loc[("EN", "0")]
    obs = a.obs
    tm = (obs["target"] == "EN") & (obs["perturbation_class"] == "single_targeting")
    cm = obs["perturbation_class"] == "single_control"
    ink = obs["leiden"].astype(str) == "0"
    A, B, C, D = int((tm & ink).sum()), int((tm & ~ink).sum()), int((cm & ink).sum()), int((cm & ~ink).sum())
    assert (row["n_target_in_cluster"], row["n_target_total"], row["n_control_in_cluster"], row["n_control_total"]) == (A, A + B, C, C + D)
    orr, p = stats.fisher_exact([[A, B], [C, D]])
    assert np.isclose(row["odds_ratio"], orr) and np.isclose(row["p_value"], p) and np.isclose(row["odds_ratio"], A * D / (B * C))
    assert np.isclose(row["log2_or_haldane"], haldane_log2_or(A, B, C, D))


def test_known_enrichment_depletion_and_null(enr):
    _, r = enr
    t = r.table.set_index(["target", "cluster"])
    assert t.loc[("EN", "0"), "direction"] == "enriched" and t.loc[("EN", "0"), "significant"]
    assert t.loc[("DE", "1"), "direction"] == "depleted" and t.loc[("DE", "1"), "significant"] and t.loc[("DE", "1"), "odds_ratio"] == 0
    assert np.isfinite(t.loc[("DE", "1"), "log2_or_haldane"])
    assert not t.loc["NU"]["significant"].any()


def test_bh_family_and_counts(enr):
    _, r = enr
    assert len(r.table) == 3 * 3 and r.info["n_tests"] == 9
    np.testing.assert_allclose(r.table["fdr"].to_numpy(), bh_fdr(r.table["p_value"].to_numpy()))
    assert r.info["bh_family"].startswith("all tested")
    assert "p_permutation" in r.info["omnibus"]


def test_control_population_and_exclusions(enr):
    a, r = enr
    t = r.table.set_index(["target", "cluster"])
    assert (t["n_control_total"] == 300).all()  # only single-guide non-targeting cells
    assert t.loc[("EN", "1"), "n_target_total"] == 90  # multi-guide 'EN;DE' cells never counted
    b = a[~a.obs["perturbation_class"].isin(["ambiguous", "multi_targeting"])].copy()
    r2 = compute_cluster_enrichment(b, cfg(), "leiden")
    pd.testing.assert_series_equal(r.table.set_index(["target", "cluster"])["p_value"].sort_index(), r2.table.set_index(["target", "cluster"])["p_value"].sort_index())
    other = compute_cluster_enrichment(a, cfg(enrichment={"control": "other"}), "leiden").table.set_index(["target", "cluster"])
    assert other.loc[("EN", "0"), "n_control_total"] == 180 and other.loc[("EN", "0"), "n_control_in_cluster"] < other.loc[("EN", "0"), "n_control_total"]


def test_guide_support_direction(enr):
    a, r = enr
    t = r.table.set_index(["target", "cluster"])
    assert t.loc[("EN", "0"), "n_guides_observed"] == 3 and t.loc[("EN", "0"), "n_guides_supporting_direction"] == 2  # EN_3 is spread
    assert t.loc[("DE", "1"), "n_guides_supporting_direction"] == 3  # depletion supported by guides below the control fraction
    assert np.isclose(t.loc[("EN", "0"), "guide_support_fraction"], 2 / 3)


def test_stratification(enr):
    a, _ = enr
    r = compute_cluster_enrichment(a, cfg(enrichment={"stratify_by": "sample"}), "leiden")
    t = r.table.set_index(["target", "cluster"])
    assert {"cmh_odds_ratio", "cmh_p_value", "cmh_fdr", "cmh_n_strata"} <= set(r.table.columns)
    assert t.loc[("EN", "0"), "cmh_n_strata"] == 2 and t.loc[("EN", "0"), "cmh_significant"]
    assert np.isclose(t.loc[("EN", "0"), "p_value"], compute_cluster_enrichment(a, cfg(), "leiden").table.set_index(["target", "cluster"]).loc[("EN", "0"), "p_value"])  # Fisher result preserved
    with pytest.raises(ValueError, match="not an obs column"):
        compute_cluster_enrichment(a, cfg(enrichment={"stratify_by": "lane_xyz"}), "leiden")


def test_small_targets_and_clusters_skipped(enr):
    a, _ = enr
    r = compute_cluster_enrichment(a, cfg(enrichment={"min_cells_per_target": 100}), "leiden")
    assert r.empty and len(r.skipped) == 3
    r2 = compute_cluster_enrichment(a, cfg(enrichment={"min_cells_per_cluster": 200}), "leiden")
    assert r2.info["n_clusters_tested"] < 3 and r2.info["clusters_not_tested"]


def test_lochness_comparison_is_descriptive(enr):
    a, _ = enr
    loch = pd.DataFrame({"target": ["EN", "DE", "NU"], "mean_lochness_in_own_cells": [1.0, 0.2, 0.0], "delta_own_vs_control": [1.1, 0.3, 0.0], "fdr": [0.01, 0.2, 0.9]})
    r = compute_cluster_enrichment(a, cfg(), "leiden", loch)
    c = r.lochness_comparison.set_index("target")
    assert c.loc["EN", "top_enriched_cluster"] == "0" and c.loc["EN", "lochness_own_mean"] == 1.0
    assert not any("combined" in col for col in c.columns)


# ----------------------------------------------------- pipeline integration
@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("cs_synth"), n_cells=400, rna_state="raw")


def test_pipeline_disabled_adds_nothing(synth, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    r = run_pipeline(Config.from_dict(base_config(synth, tmp_path, umap={"enabled": False})))
    assert "leiden" not in r.adata.obs and "analysis" not in r.adata.uns["petrubseq_protein"]
    assert not (tmp_path / "tables" / "cell_states").exists() and "Cell states and perturbation enrichment" not in r.report_html.read_text()


def test_pipeline_enabled(synth, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    cl = {"enabled": True, "resolution": 0.5, "enrichment": {"n_permutations": 20, "stratify_by": "condition"}}
    r = run_pipeline(Config.from_dict(base_config(synth, tmp_path, analysis={"clustering": cl})))
    a = ad.read_h5ad(next((tmp_path / "processed").glob("*.h5ad")))
    assert a.obs["leiden"].dtype.name == "category" and a.obs["leiden"].nunique() >= 1
    info = a.uns["petrubseq_protein"]["analysis"]["cell_states"]
    assert info["clustering"]["n_clusters"] == a.obs["leiden"].nunique() and "perturbation_effects" not in a.uns["petrubseq_protein"]["analysis"]
    td = tmp_path / "tables" / "cell_states"
    assert (td / "cluster_summary.csv").is_file() and (td / "cluster_composition_by_class.csv").is_file()
    if info["enrichment"].get("status") == "computed":
        assert (td / "perturbation_cluster_enrichment.csv").is_file() and "perturbation_cluster_enrichment" in a.uns
    man = pd.read_csv(tmp_path / "tables" / "figure_manifest.csv")
    assert set(man.loc[man["section"] == "cell_states", "name"]) >= {"umap_leiden", "cluster_composition"}
    html = r.report_html.read_text()
    assert "6. Cell states and perturbation enrichment" in html and "7. Outputs and provenance" in html and html.count('src="figures/') == 0
    md = (tmp_path / "report.md").read_text()
    assert "Cell states and perturbation enrichment" in md


def test_pipeline_both_analyses_numbering(synth, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    an = {"clustering": {"enabled": True, "enrichment": {"n_permutations": 0}}, "perturbation_effects": {"enabled": True, "lochness": {"n_permutations": 5}, "modules": {"min_perturbations": 3, "min_cells_per_perturbation": 10}}}
    r = run_pipeline(Config.from_dict(base_config(synth, tmp_path, umap={"enabled": False}, analysis=an)))
    html = r.report_html.read_text()
    assert "6. Perturbation effects" in html and "7. Cell states and perturbation enrichment" in html and "8. Outputs and provenance" in html
    assert "ps_scores" in r.adata.obsm and "leiden" in r.adata.obs


def test_depth_flag(blobs):
    a = blobs.copy()
    a.obs["total_counts"] = np.where(a.obs["blobs"].astype(str) == "1", 30000.0, 10000.0)
    r = compute_clustering(a, cfg())
    assert any("library size" in f for f in r.flags) and (r.summary["depth_ratio_to_overall"] >= 1.5).any()
