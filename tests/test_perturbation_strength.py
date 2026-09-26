"""Perturbation strength (reference ``perturbation.py`` semantics) on synthetic data."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.stats import ks_2samp

sys.path.insert(0, str(Path(__file__).parent))

from petrubseq_protein.analysis import perturbation_strength as S  # noqa: E402
from petrubseq_protein.analysis._common import bh_fdr  # noqa: E402
from petrubseq_protein.config import Config  # noqa: E402

TARGETS = ["KD1", "KD2", "NULL", "UP", "ABSENT", "SILENT", "RARE"]


def make(seed: int = 0) -> ad.AnnData:
    """KD1 / KD2: strong knockdown of their own gene; NULL: no effect; UP: expression goes up
    (significant but not a hit); ABSENT: gene not in var; SILENT: gene never expressed in
    controls; RARE: 5 cells (< min_cells_per_target). 200 NT controls + 30 ambiguous cells."""
    rng = np.random.default_rng(seed)
    n_per = {"KD1": 80, "KD2": 60, "NULL": 70, "UP": 60, "ABSENT": 40, "SILENT": 40, "RARE": 5}
    lab = ["NT"] * 200 + sum([[t] * n for t, n in n_per.items()], []) + ["AMB"] * 30
    lab = np.array(lab)
    n, g = len(lab), 60
    C = rng.poisson(3.0, (n, g)).astype(float)
    genes = ["KD1", "KD2", "NULL", "UP", "SILENT", "RARE"] + [f"G{i}" for i in range(6, g)]
    gi = {name: i for i, name in enumerate(genes)}
    C[lab == "KD1", gi["KD1"]] = rng.poisson(0.2, (lab == "KD1").sum())
    C[lab == "KD2", gi["KD2"]] = rng.poisson(1.0, (lab == "KD2").sum())
    C[lab == "UP", gi["UP"]] = rng.poisson(9.0, (lab == "UP").sum())
    C[:, gi["SILENT"]] = 0.0
    C[lab == "SILENT", gi["SILENT"]] = rng.poisson(2.0, (lab == "SILENT").sum())
    X = np.log1p(C / C.sum(1, keepdims=True) * 1e4).astype(np.float32)
    a = ad.AnnData(X=sp.csr_matrix(X), var=pd.DataFrame(index=genes))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["target"] = np.where(lab == "AMB", "", lab)
    a.obs["perturbation_class"] = np.where(lab == "NT", "single_control", np.where(lab == "AMB", "ambiguous", "single_targeting"))
    a.obs["control_class"] = np.where(lab == "NT", "non_targeting", "")
    a.obsm["X_umap_rna"] = rng.normal(size=(n, 2)).astype(np.float32)
    return a


def cfg(**st) -> Config:
    base = {"enabled": True, "strength": st}
    return Config.from_dict({"inputs": {"rna": {"file": "x"}}, "analysis": {"perturbation_effects": base}})


@pytest.fixture(scope="module")
def res():
    return S.compute_perturbation_strength(make(), cfg())


def test_defaults_match_reference():
    c = Config.from_dict({"inputs": {"rna": {"file": "x"}}}).analysis.perturbation_effects.strength
    assert c.enabled and c.controls == ["ntc", "other"] and c.primary_control == "ntc" and c.min_cells_per_target == 10 and c.min_control_cells == 10
    assert c.min_pct_expressing_control == 1.0 and c.fdr_alpha == 0.05 and c.max_log2fc_for_hit == 0.0 and c.top_n_report == 12 and c.umap_background_fraction == 0.1
    assert S.PSEUDOCOUNT == 0.01


def test_eligibility_and_skipped(res):
    t = res.table.set_index("target")
    assert set(t.index) == {"KD1", "KD2", "NULL", "UP"}
    sk = res.skipped.set_index("target")["reason"]
    assert "not present" in sk["ABSENT"] and "fewer than 10" in sk["RARE"] and "not detectably expressed" in sk["SILENT"]
    assert res.controls_used == ["ntc", "other"] and res.primary_control == "ntc" and res.n_control_cells == {"ntc": 200, "other": 355}


def test_statistics_match_the_reference_definition(res):
    a = make()
    t = res.table.set_index("target")
    x = a.X[:, a.var_names.get_loc("KD1")].toarray().ravel()
    pert = (a.obs["target"] == "KD1").to_numpy()
    ntc = (a.obs["perturbation_class"] == "single_control").to_numpy()
    other = (a.obs["perturbation_class"] == "single_targeting").to_numpy() & ~pert
    for arm, m in (("ntc", ntc), ("other", other)):
        mp, mc = np.expm1(x[pert]).mean(), np.expm1(x[m]).mean()
        assert np.isclose(t.loc["KD1", f"log2fc_{arm}"], np.log2((mp + 0.01) / (mc + 0.01)))
        assert np.isclose(t.loc["KD1", f"pct_knockdown_{arm}"], 100 * (1 - (mp + 0.01) / (mc + 0.01)))
        ks = ks_2samp(x[pert], x[m])
        assert np.isclose(t.loc["KD1", f"ks_stat_{arm}"], ks.statistic) and np.isclose(t.loc["KD1", f"ks_pval_{arm}"], ks.pvalue)
        assert t.loc["KD1", f"n_control_{arm}"] == int(m.sum())
        np.testing.assert_allclose(res.table[f"ks_fdr_{arm}"].to_numpy(), bh_fdr(res.table[f"ks_pval_{arm}"].to_numpy()))
    assert t.loc["KD1", "n_perturbed"] == 80 and np.isclose(t.loc["KD1", "mean_lognorm_perturbed_ntc"], x[pert].mean())


def test_hit_calls_are_directional_and_ranked(res):
    t = res.table.set_index("target")
    assert t.loc["KD1", "is_hit_ntc"] and t.loc["KD2", "is_hit_ntc"] and not t.loc["NULL", "is_hit_ntc"]
    assert not t.loc["UP", "is_hit_ntc"] and t.loc["UP", "ks_fdr_ntc"] < 0.05 and t.loc["UP", "log2fc_ntc"] > 0  # significant but up: no hit
    assert list(res.table["target"].head(2)) == ["KD1", "KD2"] and list(res.table["rank"]) == [1, 2, 3, 4]
    assert len(res.hits) == 2 and res.display.columns[0] == "Rank" and "log2FC (other)" in res.display.columns


def test_ambiguous_cells_never_enter(res):
    a = make()
    b = a[a.obs["perturbation_class"] != "ambiguous"].copy()
    r2 = S.compute_perturbation_strength(b, cfg())
    pd.testing.assert_frame_equal(res.table.drop(columns=["rank"]).reset_index(drop=True), r2.table.drop(columns=["rank"]).reset_index(drop=True))


def test_single_arm_and_fallback():
    a = make()
    r = S.compute_perturbation_strength(a, cfg(controls=["other"], primary_control="other"))
    assert r.controls_used == ["other"] and "log2fc_ntc" not in r.table.columns and r.hits["target"].tolist()[:2] == ["KD1", "KD2"]
    a.obs["control_class"] = ""  # no ntc cells -> the requested primary falls back to 'other'
    r2 = S.compute_perturbation_strength(a, cfg())
    assert r2.controls_used == ["other"] and r2.primary_control == "other"


def test_target_gene_map():
    a = make()
    a.var_names = ["CD274" if g == "KD1" else g for g in a.var_names]
    r = S.compute_perturbation_strength(a, cfg())
    assert "KD1" in set(r.skipped["target"])
    c = Config.from_dict({"inputs": {"rna": {"file": "x"}}, "analysis": {"perturbation_effects": {"target_gene_map": {"KD1": "CD274"}}}})
    r2 = S.compute_perturbation_strength(a, c)
    assert r2.table.set_index("target").loc["KD1", "gene"] == "CD274" and r2.table.set_index("target").loc["KD1", "is_hit_ntc"]


def test_pipeline_figures_and_tables(tmp_path):
    from petrubseq_protein.reporting.figures import FigureRegistry
    from petrubseq_protein.reporting.strength_plots import perturbation_strength_figures

    a = make()
    r = S.compute_perturbation_strength(a, cfg(top_n_report=2))
    reg = FigureRegistry(tmp_path)
    perturbation_strength_figures(r, a, cfg(top_n_report=2), reg)
    names = {x.name for x in reg.records}
    assert {"perturbation_volcano", "perturbation_waterfall", "perturbation_control_comparison", "perturbation_KD1", "perturbation_NULL"} <= names
    per = [x for x in reg.records if x.stage == "per_target"]
    assert len(per) == 4 and sum(x.in_report for x in per) == 2
