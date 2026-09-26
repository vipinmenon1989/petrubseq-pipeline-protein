"""Synthetic-data tests for petrubseq-pipeline-protein.

Run with ``pytest`` inside the ``petrubseq-protein`` environment.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import ISO_MAP, ISOTYPES, PROTEINS, TARGETS, base_config, make_dataset  # noqa: E402

from petrubseq_protein.io import readers as pio  # noqa: E402
from petrubseq_protein.preprocessing.align import AlignmentError, align_cells, harmonize_perturbations  # noqa: E402
from petrubseq_protein.config import Config, ConfigError  # noqa: E402
from petrubseq_protein.preprocessing.normalize import normalize_rna, reconstruct_counts  # noqa: E402
from petrubseq_protein.pipeline import run_pipeline, validate_inputs  # noqa: E402

NO_UMAP = {"umap": {"enabled": False}}


@pytest.fixture(scope="session")
def synthetic(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("synthetic"), n_cells=300)


@pytest.fixture(scope="session")
def result(synthetic, tmp_path_factory):
    out = tmp_path_factory.mktemp("results")
    cfg = Config.from_dict(base_config(synthetic, out))
    return run_pipeline(cfg)


# ---------------------------------------------------------------- config ---


def test_default_config_values():
    cfg = Config.from_dict({"inputs": {"rna": {"file": "x.csv"}}})
    assert cfg.rna.hvg.n_top_genes == 3000 and cfg.rna.hvg.flavor == "seurat" and cfg.rna.n_pcs == 50 and cfg.rna.scale
    assert cfg.rna.normalize == "auto" and cfg.rna.target_sum is None and cfg.rna.reconstruct_counts == "auto"  # reference: median library size
    assert cfg.protein.normalization == "auto" and cfg.protein.embedding_representation == "auto" and cfg.protein.n_pcs == 10
    assert cfg.protein.exclude_isotypes_from_embedding and cfg.protein.use_isotypes_for_qc and cfg.protein.clr_axis == "cells"
    assert cfg.neighbors.n_neighbors == 15 and cfg.umap.min_dist == 0.5 and cfg.compute.seed == 0
    assert cfg.qc.prefilter.enabled and cfg.qc.prefilter.min_genes_per_cell == 200 and cfg.qc.prefilter.min_cells_per_gene == 3
    assert cfg.qc.filter.enabled and cfg.qc.filter.rna.min_genes == 1000 and cfg.qc.filter.rna.max_pct_mt == 20.0 and cfg.qc.filter.rna.min_counts is None  # reference qc defaults
    assert cfg.qc.filter.rna.min_cells_per_gene == 3 and cfg.qc.filter.rna.max_pct_hb is None and cfg.rna.hb_pattern == "^HB[^(P)]"
    assert cfg.analysis.perturbation_effects.enabled and cfg.analysis.perturbation_effects.strength.enabled and cfg.analysis.clustering.enabled
    assert cfg.qc.filter.protein.min_total_counts is None and not cfg.qc.filter.protein.remove_extreme_counts and cfg.qc.filter.perturbation.cells == "all"
    assert cfg.qc.flags.rna_n_mads == 5.0 and cfg.qc.flags.extreme_fold == 10.0
    assert cfg.report.embed_figures and cfg.report.write_markdown and cfg.report.figure_dpi == 120 and cfg.report.max_table_rows == 100
    assert not cfg.output.write_prefilter_h5ad and not cfg.output.archive
    assert cfg.perturbation.min_cells_per_guide == 10 and cfg.perturbation.min_cells_per_target == 10
    assert cfg.perturbation.preserve_multiguide and cfg.perturbation.preserve_unassigned
    assert not cfg.multimodal.enabled and cfg.multimodal.method == "concat_pcs"


def test_config_override_and_roundtrip(synthetic, tmp_path):
    cfg_path = tmp_path / "c.yaml"
    d = base_config(synthetic, tmp_path / "r", rna={"n_pcs": 7, "hvg": {"flavor": "cell_ranger"}}, neighbors={"n_neighbors": 9}, umap={"min_dist": 0.1}, compute={"seed": 42})
    yaml.safe_dump(d, open(cfg_path, "w"))
    cfg = Config.from_yaml(cfg_path)
    assert cfg.rna.n_pcs == 7 and cfg.rna.hvg.flavor == "cell_ranger" and cfg.rna.hvg.n_top_genes == 50
    assert cfg.neighbors.n_neighbors == 9 and cfg.umap.min_dist == 0.1 and cfg.compute.seed == 42
    assert cfg.inputs.rna.paths() == ["rna_a.csv.gz", "rna_b.csv.gz"]
    rd = cfg.to_dict()
    assert Config.from_dict(rd).to_dict() == rd


def test_config_unknown_key_rejected(synthetic, tmp_path):
    with pytest.raises(ConfigError, match="Unknown key"):
        Config.from_dict(base_config(synthetic, tmp_path, rna={"normalisation": "auto"}))
    with pytest.raises(ConfigError, match="Unknown key"):
        Config.from_dict(base_config(synthetic, tmp_path, embeddings={"rna": {}}))


def test_config_invalid_value_rejected(synthetic, tmp_path):
    with pytest.raises(ConfigError, match="rna.normalize"):
        Config.from_dict(base_config(synthetic, tmp_path, rna={"normalize": "bogus"}))
    with pytest.raises(ConfigError, match="protein.clr_axis"):
        Config.from_dict(base_config(synthetic, tmp_path, protein={"clr_axis": "rows"}))


def test_missing_file_detected(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path, inputs={"protein": {"file": "nope.csv.gz"}}))
    with pytest.raises(FileNotFoundError, match="nope.csv.gz"):
        validate_inputs(cfg)


# ------------------------------------------------------------------- io ---


def test_dense_reader_matches_ground_truth(synthetic):
    m = pio.read_dense_matrix([synthetic["dir"] / "protein_counts.csv.gz"], chunk_rows=3)
    assert list(m.features) == synthetic["protein_cols"]
    assert list(m.cells) == synthetic["cells"]
    np.testing.assert_array_equal(m.X.toarray(), synthetic["protein_counts"])


def test_multifile_reader_concatenates_cells(synthetic):
    m = pio.read_dense_matrix([synthetic["dir"] / "rna_a.csv.gz", synthetic["dir"] / "rna_b.csv.gz"], chunk_rows=50, n_jobs=2)
    assert m.shape == (300, 120)
    assert list(m.cells) == synthetic["cells"]
    assert m.n_cells_per_source == [150, 150]


def test_scp_metadata_type_row_skipped(synthetic):
    meta = pio.read_table(synthetic["dir"] / "metadata.csv", fmt="scp_metadata", index_col="NAME")
    assert "TYPE" not in meta.index and meta.shape[0] == 300


def test_value_state_detection(synthetic, tmp_path):
    raw = make_dataset(tmp_path / "raw", n_cells=60, rna_state="raw")
    m_raw = pio.read_dense_matrix([raw["dir"] / "rna_a.csv.gz", raw["dir"] / "rna_b.csv.gz"])
    assert pio.detect_value_state(m_raw.X).state == "raw_counts"
    vs = pio.detect_value_state(pio.read_dense_matrix([synthetic["dir"] / "rna_a.csv.gz"]).X)
    assert vs.state == "log_normalized" and vs.log_scale == pytest.approx(1e6, rel=1e-3)


# ------------------------------------------------------------ alignment ---


def test_duplicate_barcode_detected(tmp_path):
    d = make_dataset(tmp_path / "dup", n_cells=50, duplicate_cell=True)
    cfg = Config.from_dict(base_config(d, tmp_path / "r"))
    ids = {"rna": d["cells"], "protein_counts": d["cells"] + [d["cells"][0]], "metadata": d["cells"]}
    cfg.alignment.required = ["rna", "protein_counts", "metadata"]
    with pytest.raises(AlignmentError, match="duplicate"):
        align_cells(ids, cfg)


def test_mismatched_modalities(tmp_path):
    d = make_dataset(tmp_path / "mm", n_cells=50)
    cfg = Config.from_dict(base_config(d, tmp_path / "r"))
    ids = {"rna": d["cells"], "protein": d["cells"][:40], "metadata": d["cells"]}
    cfg.alignment.min_overlap_fraction = 0.9
    with pytest.raises(AlignmentError, match="present in all required"):
        align_cells(ids, cfg)
    cfg.alignment.min_overlap_fraction = 0.5
    res = align_cells(ids, cfg)
    assert res.n_kept == 40 and res.cells == d["cells"][:40]


def test_optional_modality_missing_cells_are_kept(tmp_path):
    d = make_dataset(tmp_path / "opt", n_cells=80, drop_protein_cells=3)
    cfg = Config.from_dict(base_config(d, tmp_path / "r", **NO_UMAP))
    a = run_pipeline(cfg).adata
    assert a.n_obs == 80
    assert int((~a.obs["has_protein_counts"]).sum()) == 3
    assert a.obsm["protein_counts"].isna().any(axis=1).sum() == 3
    assert np.isnan(a.obsm["X_pca_protein"]).any(axis=1).sum() == 3
    assert a.uns["petrubseq_protein"]["alignment"]["optional_modality_missing_cells"]["protein_counts"] == 3


# -------------------------------------------------------- perturbations ---


def test_perturbation_harmonization_classes(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path))
    guides = pio.read_guide_assignments(synthetic["dir"] / "guides.txt")
    meta = pio.read_table(synthetic["dir"] / "metadata.csv", fmt="scp_metadata", index_col="NAME")
    df = harmonize_perturbations(synthetic["cells"], cfg, guides, meta)
    lists = synthetic["lists"]
    assert (df["n_guides"].to_numpy() == [len(l) for l in lists]).all()
    for i, l in enumerate(lists):
        row = df.iloc[i]
        ctrl = [g.startswith("NO_SITE") for g in l]
        if len(l) == 0:
            assert row["perturbation_class"] == "unassigned" and row["perturbation"] == "unassigned" and row["guides"] == ""
        elif len(l) == 1:
            assert row["guide"] == l[0] and row["target"] == l[0].rsplit("_", 1)[0]
            assert row["perturbation_class"] == ("single_control" if ctrl[0] else "single_targeting")
            assert row["control_class"] == ("non_targeting" if ctrl[0] else "")
            assert bool(row["is_control"]) == ctrl[0] and bool(row["is_targeting"]) == (not ctrl[0])
        else:
            assert row["perturbation"] == "multi" and row["guide"] == "" and set(row["guides"].split(";")) == set(l)
            expected = "multi_control" if all(ctrl) else ("multi_targeting" if not any(ctrl) else "mixed_control_targeting")
            assert row["perturbation_class"] == expected
            assert row["n_targets"] == len({g.rsplit("_", 1)[0] for g in l})
    assert df.attrs["metadata_guide_mismatches"] == 0
    assert df.attrs["control_classes"] == ["non_targeting"]


def test_multiguide_and_unassigned_preserved_by_default(result, synthetic):
    a = result.adata
    lists = synthetic["lists"]
    assert a.n_obs == len(lists)
    assert (a.obs["n_guides"] == 0).sum() == sum(len(l) == 0 for l in lists)
    assert (a.obs["n_guides"] > 1).sum() == sum(len(l) > 1 for l in lists)
    assert a.uns["petrubseq_protein"]["perturbations"]["cells_dropped_by_preserve_filters"] == 0


def test_preserve_filters_drop_when_disabled(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", perturbation={"preserve_multiguide": False, "preserve_unassigned": False}, output={"write_h5ad": False}, **NO_UMAP))
    a = run_pipeline(cfg).adata
    assert (a.obs["n_guides"] == 1).all()
    assert a.n_obs == sum(len(l) == 1 for l in synthetic["lists"])


def test_coverage_tables(result, synthetic):
    r = result.run_dir / "tables"
    gc = pd.read_csv(r / "guide_coverage.csv", index_col=0)
    tc = pd.read_csv(r / "target_coverage.csv", index_col=0)
    cg = pd.read_csv(r / "condition_guide_coverage.csv")
    ct = pd.read_csv(r / "condition_target_coverage.csv")
    lists = synthetic["lists"]
    single = [l[0] for l in lists if len(l) == 1]
    # every guide seen anywhere is listed; single-guide counts match ground truth
    all_guides = {g for l in lists for g in l}
    assert set(gc.index) == all_guides
    for g in all_guides:
        assert gc.loc[g, "n_cells_single_guide"] == single.count(g)
        assert gc.loc[g, "n_cells_any"] == sum(g in l for l in lists)
    assert set(gc.loc[gc["class"] == "non_targeting"].index) == {"NO_SITE_1", "NO_SITE_2"}
    assert (gc["low_coverage"] == (gc["n_cells_single_guide"] < 10)).all()
    # targets
    assert "NO_SITE" in tc.index and tc.loc["NO_SITE", "class"] == "non_targeting"
    assert set(tc.loc[tc["class"] == "targeting"].index) <= set(TARGETS)
    assert (tc["low_coverage"] == (tc["n_cells_single_guide"] < 20)).all()
    assert tc.loc[TARGETS[0], "n_guides"] == 2
    # condition x guide/target long tables
    cond = synthetic["condition"]
    assert set(cg.columns) >= {"guide", "condition", "n_cells", "class", "low_coverage"}
    for g in ["NO_SITE_1", TARGETS[0] + "_1"]:
        for c in ("Control", "Treated"):
            truth = sum(1 for l, cc in zip(lists, cond) if len(l) == 1 and l[0] == g and cc == c)
            assert int(cg.loc[(cg["guide"] == g) & (cg["condition"] == c), "n_cells"].iloc[0]) == truth
    assert ct.groupby("target")["n_cells"].sum().reindex(tc.index).fillna(0).astype(int).eq(tc["n_cells_single_guide"]).all()
    pqc = result.adata.uns["petrubseq_protein"]["perturbations"]["qc"]
    assert pqc["n_targeting_targets"] == len(set(tc.loc[tc["class"] == "targeting"].index))
    assert pqc["class_counts"]["unassigned"] == sum(len(l) == 0 for l in lists)


# ---------------------------------------------------------- normalization ---


def test_reconstruct_counts_exact(synthetic):
    m = pio.read_dense_matrix([synthetic["dir"] / "rna_a.csv.gz", synthetic["dir"] / "rna_b.csv.gz"])
    tot = synthetic["counts"].sum(axis=1).astype(float)
    tot[tot == 0] = 1
    counts, info = reconstruct_counts(m.X, tot, 1e6)
    assert info["accepted"]
    np.testing.assert_array_equal(counts.toarray(), synthetic["counts"])


def test_reconstructed_counts_are_named_and_documented(result):
    a = ad.read_h5ad(result.h5ad)
    assert "reconstructed_counts" in a.layers and "counts" not in a.layers
    rna = a.uns["petrubseq_protein"]["rna"]
    assert rna["counts_layer"] == "reconstructed_counts"
    assert "NOT the directly observed" in rna["counts_source"]
    assert rna["counts_reconstruction"]["accepted"]
    assert "reconstructed" in a.uns["petrubseq_protein"]["schema"]["layers['reconstructed_counts']"]["status"]
    assert a.uns["petrubseq_protein"]["qc"]["rna"]["source"] == "layers['reconstructed_counts']"


def test_no_double_normalization(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path))
    m = pio.read_dense_matrix([synthetic["dir"] / "rna_a.csv.gz"])
    X, counts, layer, info = normalize_rna(m.X, cfg.rna, None, "auto")
    assert info["input_state_used"] == "log_normalized" and info["method"].startswith("none") and layer is None
    np.testing.assert_array_equal(X.toarray(), m.X.toarray())
    cfg.rna.normalize = "always"
    with pytest.raises(ValueError, match="normalize twice"):
        normalize_rna(m.X, cfg.rna, None, "auto")


def test_raw_counts_are_normalized_and_kept(tmp_path):
    d = make_dataset(tmp_path / "raw", n_cells=80, rna_state="raw")
    cfg = Config.from_dict(base_config(d, tmp_path / "r", **NO_UMAP))
    a = run_pipeline(cfg).adata
    assert "counts" in a.layers and "reconstructed_counts" not in a.layers
    np.testing.assert_array_equal(a.layers["counts"].toarray(), d["counts"])
    assert a.uns["petrubseq_protein"]["rna"]["method"] == "normalize_total+log1p"
    rs = np.expm1(a.X.toarray()).sum(axis=1)
    # reference default: normalize_total(target_sum=None) scales every cell to the median library size
    lib = np.asarray(a.layers["counts"].sum(axis=1)).ravel()
    assert np.allclose(rs[rs > 0], np.median(lib), rtol=1e-3) and not np.allclose(np.median(lib), 1e4)


# ----------------------------------------------------------------- protein ---


def test_protein_representations_and_isotypes(result, synthetic):
    a = ad.read_h5ad(result.h5ad)
    # three separate representations, none overwritten
    assert list(a.obsm["protein"].columns) == PROTEINS
    assert list(a.obsm["protein_counts"].columns) == synthetic["protein_cols"]
    assert list(a.obsm["protein_clr"].columns) == PROTEINS
    np.testing.assert_array_equal(a.obsm["protein_counts"].to_numpy(), synthetic["protein_counts"])
    np.testing.assert_allclose(a.obsm["protein"].to_numpy(), synthetic["protein_norm"], atol=1e-5)
    C = synthetic["protein_counts"][:, [synthetic["protein_cols"].index(p) for p in PROTEINS]]
    np.testing.assert_allclose(a.obsm["protein_clr"].to_numpy(), np.log1p(C) - np.log1p(C).mean(axis=0, keepdims=True), atol=1e-5)
    feats = a.uns["protein_features"]
    assert set(feats.index) == set(synthetic["protein_cols"])
    assert feats.loc[ISOTYPES, "is_isotype"].all() and not feats.loc[PROTEINS, "is_isotype"].any()
    assert dict(feats.loc[PROTEINS, "isotype_control"]) == ISO_MAP
    # isotypes: excluded from embedding, retained for QC
    assert not feats.loc[ISOTYPES, "used_in_embedding"].any() and feats.loc[PROTEINS, "used_in_embedding"].all()
    assert feats.loc[ISOTYPES, "used_in_qc"].all()
    assert set(a.uns["pca_protein"]["params"]["features"]) == set(PROTEINS)
    assert set(a.uns["pca_protein"]["loadings"].index) == set(PROTEINS)
    assert a.uns["petrubseq_protein"]["protein"]["embedding_key"] == "protein_clr"
    assert a.uns["petrubseq_protein"]["protein_check"]["formula_reproduced"]
    iso_total = synthetic["protein_counts"][:, [synthetic["protein_cols"].index(i) for i in ISOTYPES]].sum(axis=1)
    np.testing.assert_array_equal(a.obs["protein_isotype_counts"].to_numpy(), iso_total)
    assert "protein_pct_isotype" in a.obs and "protein_extreme_counts" in a.obs
    tab = pd.read_csv(result.tables["protein_qc"], index_col=0)
    assert "background_dominated" in tab.columns and "frac_cells_above_isotype" in tab.columns


def test_protein_embedding_on_provided_when_requested(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", protein={"embedding_representation": "provided"}, output={"write_h5ad": False}, **NO_UMAP))
    a = run_pipeline(cfg).adata
    assert a.uns["petrubseq_protein"]["protein"]["embedding_key"] == "protein"
    assert "protein_clr" in a.obsm


# ---------------------------------------------------------- end to end ---


def test_outputs_created(result):
    r = result.run_dir
    assert result.h5ad.exists()
    for f in ("report.html", "report.md", "logs/run.log", "logs/resolved_config.yaml", "logs/run_manifest.json"):
        assert (r / f).exists(), f
    for t in ("perturbation_counts", "guide_coverage", "target_coverage", "condition_guide_coverage", "condition_target_coverage", "condition_summary", "protein_qc", "protein_features", "protein_pca_loadings", "rna_gene_qc", "run_summary", "moi_summary", "qc_filtering_steps", "qc_summary", "perturbation_qc", "figure_manifest"):
        assert (r / "tables" / f"{t}.csv").exists(), t
    for t in ("cell_qc", "cell_qc_prefilter"):
        assert (r / "tables" / f"{t}.csv.gz").exists(), t
    figs = list((r / "figures").rglob("*.png"))
    assert len(figs) >= 20
    for name in ("rna_umap_condition", "protein_umap_condition", "coverage_histograms", "condition_target_coverage", "protein_pca_loadings"):
        assert any(name in f.name for f in figs), name
    assert "X_umap_multimodal" not in result.adata.obsm  # off by default


def test_processed_object_schema(result, synthetic):
    a = ad.read_h5ad(result.h5ad)
    assert a.n_obs == 300 and a.n_vars == 120
    assert list(a.obs_names) == synthetic["cells"] and list(a.var_names) == synthetic["genes"]
    for col in ("condition", "guides", "n_guides", "targets", "n_targets", "guide", "target", "control_class", "perturbation", "perturbation_class", "is_control", "is_targeting", "is_single_guide", "moi_provided", "rna_total_counts_provided"):
        assert col in a.obs.columns, col
    for col in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "rna_qc_fail", "rna_outlier_low_counts", "protein_total_counts", "protein_total_counts_targeting", "protein_isotype_counts", "protein_n_detected", "protein_pct_isotype", "protein_qc_fail"):
        assert col in a.obs.columns, col
    for k in ("X_pca", "X_umap_rna", "X_pca_protein", "X_umap_protein", "X_umap_provided", "protein", "protein_counts", "protein_clr"):
        assert k in a.obsm, k
    assert a.obsm["X_pca"].shape == (300, 10) and a.obsm["X_pca_protein"].shape == (300, 4)
    assert "highly_variable" in a.var and a.var["highly_variable"].sum() == 50
    assert "PCs" in a.varm
    u = a.uns["petrubseq_protein"]
    for k in ("version", "rna", "protein", "qc", "embeddings", "alignment", "perturbations", "provenance", "config", "schema"):
        assert k in u, k
    assert u["qc"]["filter"]["cells_after"] == 300
    assert u["embeddings"]["rna"]["pca"]["scaled"] and u["embeddings"]["rna"]["hvg"]["flavor"] == "seurat"
    assert u["embeddings"]["rna"]["umap"]["random_state"] == 0
    assert u["config"]["dataset"]["name"] == "synthetic"
    prov = u["provenance"]
    for k in ("python", "scanpy", "anndata", "numpy"):
        assert k in prov["packages"], k
    assert "run_timestamp" in prov
    for slot in ("X", "layers['reconstructed_counts']", "obsm['protein']", "obsm['protein_counts']", "obsm['protein_clr']", "obsm['X_pca']", "obsm['X_pca_protein']"):
        assert slot in u["schema"], slot
    d = u["embeddings"]["diagnostics"]
    assert "rna_protein_knn_overlap_mean" in d and "protein_pc_vs_log_total_adt" in d


def test_reproducible_with_same_seed(synthetic, tmp_path):
    runs = []
    for i in range(2):
        cfg = Config.from_dict(base_config(synthetic, tmp_path / f"r{i}", output={"write_h5ad": False}))
        runs.append(run_pipeline(cfg).adata)
    a, b = runs
    assert list(a.obs_names) == list(b.obs_names) and list(a.var_names) == list(b.var_names)
    np.testing.assert_allclose(a.obsm["X_pca"], b.obsm["X_pca"], atol=1e-5)
    np.testing.assert_allclose(a.obsm["X_pca_protein"], b.obsm["X_pca_protein"], atol=1e-5)
    np.testing.assert_allclose(a.obsm["X_umap_rna"], b.obsm["X_umap_rna"], atol=1e-4)
    np.testing.assert_allclose(a.obsm["X_umap_protein"], b.obsm["X_umap_protein"], atol=1e-4)
    pd.testing.assert_frame_equal(a.obs, b.obs)
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r2", output={"write_h5ad": False}, compute={"seed": 7}))
    c = run_pipeline(cfg).adata
    assert not np.allclose(a.obsm["X_umap_rna"], c.obsm["X_umap_rna"], atol=1e-4)


def test_multimodal_optional(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", multimodal={"enabled": True, "protein_weight": 0.5}, output={"write_h5ad": False}))
    a = run_pipeline(cfg).adata
    assert a.obsm["X_multimodal"].shape == (300, 14) and a.obsm["X_umap_multimodal"].shape == (300, 2)
    d = a.uns["petrubseq_protein"]["embeddings"]["diagnostics"]
    assert "multimodal_knn_overlap_with_rna" in d
    assert a.uns["petrubseq_protein"]["embeddings"]["multimodal"]["params"]["protein_weight"] == 0.5


def test_qc_filter_when_enabled(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", qc={"filter": {"enabled": True, "rna": {"min_genes": 200}}}, output={"write_h5ad": False}, **NO_UMAP))
    with pytest.raises(ValueError, match="removed every cell"):
        run_pipeline(cfg)
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r2", qc={"filter": {"enabled": True, "rna": {"min_genes": 104, "max_pct_mt": 100}}}, output={"write_h5ad": False}, **NO_UMAP))
    res = run_pipeline(cfg)
    assert 0 < res.n_cells < 300
    assert res.adata.uns["petrubseq_protein"]["qc"]["filter"]["cells_after"] == res.n_cells


def test_summary_tables_consistent(result):
    pc = pd.read_csv(result.tables["perturbation_counts"], index_col=0)
    assert pc["n_cells"].sum() == 300
    assert set(pc.loc[pc["perturbation_class"] == "single_targeting"].index) <= set(TARGETS)
    assert pc.loc["NO_SITE", "perturbation_class"] == "single_control"
    cs = pd.read_csv(result.tables["condition_summary"], index_col=0)
    assert cs["n_cells"].sum() == 300
    rs = pd.read_csv(result.tables["run_summary"], index_col=0)["value"]
    assert rs["rna_counts_layer"] == "reconstructed_counts" and rs["protein_embedding_representation"] == "protein_clr"
