"""Statistical level of every RNA <-> protein analysis (docs/RNA_PROTEIN_LEVELS.md):
pseudoreplication guards for the target-level distance <-> protein association, the
cell-level nature of PS / lochNESS <-> protein, the phenotype-geometry Mantel test and the
module-wise test on synthetic data."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).parent))
from test_perturbation_distance import dcfg, make  # noqa: E402

from petrubseq_protein import analysis  # noqa: E402
from petrubseq_protein.analysis import master_table as MT  # noqa: E402
from petrubseq_protein.analysis.concordance import compute_concordance  # noqa: E402
from petrubseq_protein.analysis.protein_effects import compute_protein_effects  # noqa: E402


class _Dist:
    """Minimal stand-in for DistanceResults: one distance per target."""

    def __init__(self, table):
        self.table = table
        self.empty = table.empty


class _Prot:
    def __init__(self, matrix, fdr=None):
        self.matrix = matrix
        self.fdr_matrix = fdr if fdr is not None else pd.DataFrame(1.0, index=matrix.index, columns=matrix.columns)
        self.d_matrix = matrix
        self.table = pd.DataFrame({"target": np.repeat(matrix.index, matrix.shape[1]), "protein": list(matrix.columns) * len(matrix)})
        self.info = {"fdr_alpha": 0.05}
        self.empty = False


class _Space:
    def __init__(self, M, modules):
        self.distance_matrix = M
        self.phenotype_modules = modules
        self.empty = False
        self.neighbors = pd.DataFrame()
        self.coordinates = pd.DataFrame()


# ------------------------------------------------------------------ Part 12: pseudoreplication
def test_distance_protein_uses_targets_not_cells():
    """3 targets, one distance each, many cells with varying protein: n must be 3 targets."""
    targets = ["T1", "T2", "T3"]
    dist = _Dist(pd.DataFrame({"target": targets, "energy_distance": [3.0, 2.0, 1.0], "fdr": [0.01, 0.2, 0.5], "significant": [True, False, False]}))
    prot = _Prot(pd.DataFrame({"P": [0.9, 0.5, 0.1]}, index=targets))
    summ, long = MT.distance_protein_association(dist, prot, min_targets=3)
    assert list(summ["n_targets"]) == [3] and len(long) == 3 and set(summ["analysis_level"]) == {"TARGET_LEVEL"}
    assert summ.iloc[0]["rho_signed"] == pytest.approx(1.0)
    # duplicating cells cannot change a target-level statistic: the inputs are one row per target
    summ2, _ = MT.distance_protein_association(dist, prot, min_targets=3)
    pd.testing.assert_frame_equal(summ, summ2)
    # too few targets -> no correlation claimed
    summ3, _ = MT.distance_protein_association(dist, prot, min_targets=5)
    assert summ3["rho_signed"].isna().all() and (summ3["status_signed"] == "insufficient_support").all()


@pytest.fixture(scope="module")
def data():
    return make(ambiguous=30, rare=5)


def _dup(a):
    """Every cell duplicated (obs names suffixed); the target-level structure is unchanged."""
    import anndata as ad

    b = ad.concat([a, a], index_unique="-")
    b.obsm["protein"] = pd.concat([a.obsm["protein"], a.obsm["protein"]]).set_axis(b.obs_names)
    b.var["highly_variable"] = True
    b.uns["_w"] = np.concatenate([a.uns["_w"], a.uns["_w"]])
    return b


def test_duplicating_cells_changes_cell_level_but_not_target_level(data):
    """Duplicating every cell: PS <-> protein (cell level) sees 2n cells and its p-values shrink; the
    distance <-> protein table (target level) keeps n = targets and the distance-protein correlation
    is unchanged. Target-level protein effects (difference of means) are unchanged too."""
    cfg = dcfg(modules={"enabled": False}, lochness={"enabled": False}, distance={"n_permutations": 20}, master_table={"mantel_permutations": 20})
    r1 = analysis.run_perturbation_effects(data, cfg)
    r2 = analysis.run_perturbation_effects(_dup(data), cfg)
    # cell level: sample size doubles
    w1 = r1.concordance.ps_protein.query("scope == 'within_target'").set_index(["target", "protein"])
    w2 = r2.concordance.ps_protein.query("scope == 'within_target'").set_index(["target", "protein"])
    assert (w2.loc[w1.index, "n_cells"] == 2 * w1["n_cells"]).all()
    assert np.allclose(w2.loc[w1.index, "spearman_rho"], w1["spearman_rho"], atol=1e-9)
    ok = w1["p_value"].between(1e-12, 0.5)
    assert (w2.loc[w1.index[ok], "p_value"] < w1.loc[ok, "p_value"]).all()          # genuinely cell-level: more cells -> smaller p
    # target level: n stays the number of targets, correlation and p unchanged
    d1, d2 = r1.master.distance_protein.set_index("protein"), r2.master.distance_protein.set_index("protein")
    assert (d1["n_targets"] == 6).all() and (d2["n_targets"] == 6).all()
    assert np.allclose(d1["rho_signed"], d2["rho_signed"], atol=1e-9) and np.allclose(d1["p_signed"], d2["p_signed"], atol=1e-9)
    assert np.allclose(d1["rho_abs"], d2["rho_abs"], atol=1e-9)
    assert np.allclose(r1.protein.matrix.loc[d1.index[0] and r1.protein.matrix.index], r2.protein.matrix.loc[r1.protein.matrix.index], atol=1e-9)
    # energy distances are per target and identical on duplicated cells (means of the same pairwise distances)
    e1, e2 = r1.distance.table.set_index("target")["energy_distance"], r2.distance.table.set_index("target")["energy_distance"]
    assert np.allclose(e1, e2.loc[e1.index], atol=1e-9)
    assert (r2.master.distance_protein_targets.groupby("protein").size() == 6).all()   # one row per target x protein, not per cell


def test_master_protein_table_levels(data):
    cfg = dcfg(distance={"n_permutations": 20}, master_table={"mantel_permutations": 20})
    r = analysis.run_perturbation_effects(data, cfg, analysis.run_perturbation_strength(data, cfg))
    p = r.master.protein
    for col, level in (("protein_effect_level", "TARGET_LEVEL"), ("ps_protein_level", "CELL_LEVEL"), ("lochness_protein_level", "CELL_LEVEL"), ("distance_protein_level", "TARGET_LEVEL"),
                       ("phenotype_module_protein_level", "PHENOTYPE_MODULE_LEVEL"), ("phenotype_geometry_protein_level", "TARGET_PAIR_LEVEL"), ("program_protein_level", "TARGET_LEVEL"), ("target_columns_level", "TARGET_LEVEL")):
        assert col in p.columns and set(p[col].dropna()) == {level}, col
    assert p[["target", "protein"]].drop_duplicates().shape[0] == len(p) == 6 * 3
    # the per-protein target-level statistic is identical on every row of that protein (a repeated summary, labelled as such)
    assert (p.groupby("protein")["distance_protein_rho_signed"].nunique() == 1).all() and (p.groupby("protein")["distance_protein_n_targets"].first() == 6).all()
    # the joint (all-protein) Mantel statistic is not in this table
    assert not any("all_proteins" in str(v) for v in p.columns) and "mantel_rho" not in p.columns
    for lvl_tab, level in ((r.master.distance_protein, "TARGET_LEVEL"), (r.master.phenotype_protein_mantel, "TARGET_PAIR_LEVEL"), (r.master.phenotype_module_protein, "PHENOTYPE_MODULE_LEVEL"),
                           (r.concordance.lochness_protein_cells, "CELL_LEVEL"), (r.concordance.lochness_protein_summary, "TARGET_LEVEL"), (r.concordance.program_protein_cells, "CELL_LEVEL"), (r.concordance.program_protein_targets, "TARGET_LEVEL")):
        assert set(lvl_tab["analysis_level"]) == {level}
    assert set(r.concordance.ps_protein["analysis_level"]) == {"CELL_LEVEL", "CELL_LEVEL_POOLED_DESCRIPTIVE"}
    tb = analysis.tables(r)
    for name in ("ps_protein_associations", "lochness_protein_associations", "distance_protein_associations", "phenotype_module_protein_associations", "rna_protein_geometry_concordance"):
        assert name in tb and "analysis_level" in tb[name].columns, name
    assert "n_guides" not in r.master.perturbation.columns or r.master.perturbation["n_guides"].notna().any()


def test_lochness_protein_is_cell_level(data):
    """Within a target, the per-cell own-target lochNESS varies and is correlated with the per-cell protein."""
    cfg = dcfg(modules={"enabled": False}, ps={"enabled": False}, distance={"enabled": False}, distance_space={"enabled": False})
    r = analysis.run_perturbation_effects(data, cfg)
    lc = r.concordance.lochness_protein_cells
    assert not lc.empty and set(lc["analysis_level"]) == {"CELL_LEVEL"}
    row = lc.set_index(["target", "protein"]).loc[("TA", "P_up")]
    assert row["n_cells"] == 80
    m = ((data.obs["target"] == "TA") & (data.obs["perturbation_class"] == "single_targeting")).to_numpy()
    x = r.lochness.scores.loc[m, "TA"].to_numpy(float)
    assert np.std(x) > 0                                         # the score genuinely varies cell by cell
    from scipy.stats import spearmanr

    assert row["spearman_rho"] == pytest.approx(spearmanr(x, data.obsm["protein"].loc[m, "P_up"]).statistic, abs=1e-12)
    # the target-level summary is a different statistic over n = targets and is labelled so
    assert (r.concordance.lochness_protein_summary["n_targets"] == 6).all() and set(r.concordance.lochness_protein_summary["analysis_level"]) == {"TARGET_LEVEL"}


# ------------------------------------------------------------------ Part 13: phenotype geometry
def _geometry(seed=0, K=8):
    rng = np.random.default_rng(seed)
    targets = [f"T{i}" for i in range(K)]
    pos = rng.normal(size=(K, 3))
    R = pd.DataFrame(cdist(pos, pos), index=targets, columns=targets)          # RNA pairwise geometry
    E = pd.DataFrame(pos[:, :2] + rng.normal(scale=0.05, size=(K, 2)), index=targets, columns=["P1", "P2"])   # protein effects mirror the geometry
    return targets, R, E


def test_mantel_construction_and_permutation():
    targets, R, E = _geometry()
    prot = _Prot(E)
    space = _Space(R, pd.DataFrame({"target": targets, "phenotype_module": ["PM1"] * 4 + ["PM2"] * 4}))
    mantel, kw, means = MT.phenotype_space_protein_association(space, prot, min_targets=5, n_perm=199, seed=1)
    assert list(mantel["protein"]) == ["all_proteins", "P1", "P2"] and (mantel["n_targets"] == 8).all() and (mantel["n_pairs"] == 28).all()
    joint = mantel.iloc[0]
    assert joint["mantel_rho"] > 0.7 and joint["p_permutation"] <= 3 / 200 and joint["analysis_level"] == "TARGET_PAIR_LEVEL"
    assert mantel.loc[mantel["protein"] == "P1", "protein_distance_metric"].iloc[0] == "|E_i,p - E_j,p|"
    # explicit construction: upper triangles in the same target order, Euclidean protein-profile distances
    D1 = R.loc[targets, targets].to_numpy()
    D2 = MT.protein_profile_distances(E.loc[targets].to_numpy())
    np.testing.assert_allclose(D2, cdist(E.loc[targets].to_numpy(), E.loc[targets].to_numpy()))
    rho, p, n_pairs = MT.mantel_test(D1, D2, 199, 1)
    assert rho == pytest.approx(joint["mantel_rho"]) and p == pytest.approx(joint["p_permutation"]) and n_pairs == 28
    assert MT.mantel_test(D1, D2, 199, 1) == MT.mantel_test(D1, D2, 199, 1)              # deterministic
    # a reordered (but consistently reordered) target list gives the same statistic
    perm = np.random.default_rng(3).permutation(8)
    rho2, _, _ = MT.mantel_test(D1[np.ix_(perm, perm)], D2[np.ix_(perm, perm)], 0, 1)
    assert rho2 == pytest.approx(rho)
    # shuffling the protein target labels destroys the concordance
    shuf = np.random.default_rng(5).permutation(8)
    rho_s, p_s, _ = MT.mantel_test(D1, D2[np.ix_(shuf, shuf)], 199, 1)
    assert rho_s < rho - 0.4 and p_s > 0.05
    # the raw pair-level Spearman p-value is never used: with n_perm = 0 the p is NaN
    assert np.isnan(MT.mantel_test(D1, D2, 0, 1)[1])
    with pytest.raises(ValueError):
        MT.mantel_test(D1, D2 + 1.0, 10, 1)                                                # non-zero diagonal rejected
    with pytest.raises(ValueError):
        MT.mantel_test(D1, D2[:5, :5], 10, 1)


def test_module_test_excludes_small_modules():
    targets, R, E = _geometry(K=9)
    prot = _Prot(E)
    mods = pd.DataFrame({"target": targets, "phenotype_module": ["PM1"] * 4 + ["PM2"] * 3 + ["PM3", "PM4"]})
    _, kw, means = MT.phenotype_space_protein_association(_Space(R, mods), prot, min_targets=5, n_perm=0, seed=0, min_targets_per_module=3)
    assert (kw["n_modules_eligible"] == 2).all() and (kw["modules_eligible"] == "PM1;PM2").all() and (kw["n_modules_excluded"] == 2).all() and (kw["n_targets"] == 7).all()
    assert kw["p_value"].between(0, 1).all() and set(kw["analysis_level"]) == {"PHENOTYPE_MODULE_LEVEL"} and kw["test"].str.startswith("kruskal_wallis").all()
    assert set(means.loc[~means["eligible"], "phenotype_module"]) == {"PM3", "PM4"}
    # all singletons -> no test, status insufficient
    mods2 = pd.DataFrame({"target": targets, "phenotype_module": [f"PM{i}" for i in range(9)]})
    _, kw2, _ = MT.phenotype_space_protein_association(_Space(R, mods2), prot, min_targets=5, n_perm=0, seed=0, min_targets_per_module=3)
    assert kw2["p_value"].isna().all() and (kw2["status"] == "insufficient_support").all()
