"""Stage C: generic input adapters, guide calling, multi-lane inputs, protein
feature annotation and the cross-format equivalence of the processed object."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import PROTEINS, base_config, make_dataset, write_formats, write_lanes, write_mtx  # noqa: E402

from petrubseq_protein import audit, validation  # noqa: E402
from petrubseq_protein.config import Config, ConfigError  # noqa: E402
from petrubseq_protein.io import adapters  # noqa: E402
from petrubseq_protein.io.adapters import InputError  # noqa: E402
from petrubseq_protein.io.tenx import make_unique, read_10x_h5, read_mtx_dir  # noqa: E402
from petrubseq_protein.pipeline import run_pipeline  # noqa: E402
from petrubseq_protein.preprocessing import guides as guidecall  # noqa: E402

NO_PROTEIN_NORM = {"file": None, "required": False}
NO_DENSE = {"rna": {"files": [], "file": None}, "protein_counts": {"file": None}, "protein": NO_PROTEIN_NORM}
COMMON = {"embedding": {"file": None}, "metadata": {"file": "metadata.csv", "format": "scp_metadata"}}


@pytest.fixture(scope="session")
def raw(tmp_path_factory):
    """Raw-count experiment (no provided normalized protein) in every layout."""
    root = tmp_path_factory.mktemp("crossformat")
    d = make_dataset(root / "dense", n_cells=200, rna_state="raw", protein_norm=False)
    d["formats"] = write_formats(d, root / "formats")
    d["lanes"] = write_lanes(d, root / "lanes")
    return d


def _cfg(d, out, inputs, pert=None, columns=None, **kw):
    pert = {"assignment": pert or {"source": "guide_counts", "method": "threshold"}}
    return Config.from_dict(base_config(d, out, inputs=inputs, columns={"guide": None, **(columns or {})}, perturbation=pert, **kw))


def variants(d):
    f = d["formats"]
    return {
        "dense": {"protein": NO_PROTEIN_NORM, "guide_counts": {"file": "guide_counts.csv.gz"}, "guide_assignments": {"file": None}, **COMMON},
        "mtx_combined": {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"]), "var_names": "name"}, "guide_assignments": {"file": None}, **NO_DENSE, **COMMON},
        "mtx_separate": {"rna": {"format": "mtx", "file": str(f["mtx_rna"]), "files": [], "var_names": "name"}, "protein": NO_PROTEIN_NORM, "protein_counts": {"format": "mtx", "file": str(f["mtx_protein"])}, "guide_counts": {"format": "mtx", "file": str(f["mtx_guide"])}, "guide_assignments": {"file": None}, **COMMON},
        "10x_h5": {"multiplexed": {"format": "10x_h5", "path": str(f["h5_combined"]), "var_names": "name"}, "guide_assignments": {"file": None}, **NO_DENSE, **COMMON},
        "h5ad_combined": {"multiplexed": {"format": "h5ad", "path": str(f["h5ad_combined"])}, "guide_assignments": {"file": None}, **NO_DENSE, **COMMON},
        "h5ad_slots": {"multiplexed": {"format": "h5ad", "path": str(f["h5ad_slots"]), "slots": {"rna": {"slot": "layers", "key": "counts"}, "protein_counts": {"slot": "obsm", "key": "protein_counts"}, "guide_counts": {"slot": "obsm", "key": "guide_counts"}}}, "guide_assignments": {"file": None}, **NO_DENSE, **COMMON},
    }


@pytest.fixture(scope="session")
def processed(raw, tmp_path_factory):
    out = tmp_path_factory.mktemp("cf_out")
    res = {}
    for name, inputs in variants(raw).items():
        res[name] = run_pipeline(_cfg(raw, out / name, inputs)).adata
    return res


# ------------------------------------------------------------- readers ---


def test_readers_agree(raw):
    f = raw["formats"]
    m = read_mtx_dir(f["mtx_combined"])
    h = read_10x_h5(f["h5_combined"])
    assert list(m.features) == list(h.features) and list(m.cells) == list(h.cells)
    assert (m.X != h.X).nnz == 0
    assert sorted(m.var["feature_types"].unique()) == ["Antibody Capture", "CRISPR Guide Capture", "Gene Expression"]
    assert m.X.dtype == np.float32


def test_make_unique():
    names, n = make_unique(np.array(["A", "B", "A", "A", "B-1"]))
    assert n == 2 and len(set(names)) == 5 and names[0] == "A" and names[2] == "A-1" and names[3] == "A-2"


def test_canonical_inputs_identical_across_formats(raw, tmp_path):
    ref = None
    for name, inputs in variants(raw).items():
        ci = adapters.load_inputs(_cfg(raw, tmp_path, inputs))
        errs, _ = validation.validate_canonical(ci, _cfg(raw, tmp_path, inputs))
        assert errs == [], (name, errs)
        mats = ci.matrices()
        assert set(mats) == {"rna", "protein_counts", "guide_counts"}, name
        if ref is None:
            ref = mats
            continue
        for k in ref:
            assert list(mats[k].features) == list(ref[k].features), (name, k)
            assert list(mats[k].cells) == list(ref[k].cells), (name, k)
            assert (mats[k].X != ref[k].X).nnz == 0, (name, k)
        assert ci.provenance["guide_counts_available"] and not ci.provenance["provided_assignments_available"]


# ------------------------------------------------------- equivalence ---


def _same_frame(a: pd.DataFrame, b: pd.DataFrame, atol: float):
    assert list(a.columns) == list(b.columns)
    np.testing.assert_allclose(a.to_numpy(float), b.to_numpy(float), atol=atol, equal_nan=True)


def test_processed_object_equivalent_across_formats(processed):
    ref = processed["dense"]
    assert "guide_counts" in ref.obsm and "guide_features" in ref.uns
    for name, a in processed.items():
        if name == "dense":
            continue
        assert list(a.obs_names) == list(ref.obs_names) and list(a.var_names) == list(ref.var_names), name
        assert (a.X != ref.X).nnz == 0 and (a.layers["counts"] != ref.layers["counts"]).nnz == 0, name
        for col in ("perturbation", "perturbation_class", "guides", "guide", "target", "n_guides", "control_class", "condition", "has_guide_counts", "guide_dominant_call", "guides_detected"):
            assert (a.obs[col].astype(str).to_numpy() == ref.obs[col].astype(str).to_numpy()).all(), (name, col)
        for col in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "protein_total_counts", "guide_total_counts", "guide_top_count"):
            np.testing.assert_allclose(a.obs[col].to_numpy(float), ref.obs[col].to_numpy(float), atol=1e-5, err_msg=f"{name}:{col}")
        for k in ("protein", "protein_counts"):
            _same_frame(a.obsm[k], ref.obsm[k], 1e-5)
        assert (a.obsm["guide_counts"] != ref.obsm["guide_counts"]).nnz == 0, name
        np.testing.assert_allclose(a.obsm["X_pca"], ref.obsm["X_pca"], atol=1e-4, err_msg=name)
        np.testing.assert_allclose(a.obsm["X_pca_protein"], ref.obsm["X_pca_protein"], atol=1e-4, err_msg=name)
        np.testing.assert_allclose(a.obsm["X_umap_rna"], ref.obsm["X_umap_rna"], atol=1e-3, err_msg=name)
        gf, gr = a.uns["guide_features"], ref.uns["guide_features"]
        assert list(gf.index) == list(gr.index) and (gf["target"] == gr["target"]).all() and (gf["n_cells_detected"] == gr["n_cells_detected"]).all(), name
        pf, pr = a.uns["protein_features"], ref.uns["protein_features"]
        assert list(pf.index) == list(pr.index), name
        for col in ("role", "is_isotype", "isotype_control", "in_counts", "used_in_embedding", "used_in_qc"):
            assert (pf[col].astype(str) == pr[col].astype(str)).all(), (name, col)
        prov = a.uns["petrubseq_protein"]["inputs"]
        assert "provenance" in prov and prov["provenance"]["guide_counts_available"], name


def test_threshold_rule_reproduces_provided_lists(processed, raw):
    a = processed["dense"]
    truth = [";".join(l) for l in raw["lists"]]
    assert list(a.obs["guides"].astype(str)) == truth
    assert a.uns["petrubseq_protein"]["perturbations"]["source"] == "guide_counts:threshold"
    asg = a.uns["petrubseq_protein"]["perturbations"]["assignment"]
    assert asg["effective_source"] == "guide_counts" and asg["guide_calling"]["method"] == "threshold"


def test_input_provenance_recorded(processed):
    for name, a in processed.items():
        p = a.uns["petrubseq_protein"]["inputs"]["provenance"]
        assert p["modalities"]["rna"]["format"] == ("dense_csv" if name == "dense" else ("mtx" if name.startswith("mtx") else ("10x_h5" if name == "10x_h5" else "h5ad"))), name
        assert "guide_counts" in p["modalities"]
        if name in ("mtx_combined", "10x_h5", "h5ad_combined", "h5ad_slots"):
            assert p["multiplexed"] is not None and "split" in p["multiplexed"]


# ------------------------------------------------ guide sources / rules ---


def test_guide_source_provided_keeps_counts_and_compares(raw, tmp_path):
    inputs = {"protein": NO_PROTEIN_NORM, "guide_counts": {"file": "guide_counts.csv.gz"}, **COMMON}
    a = run_pipeline(_cfg(raw, tmp_path, inputs, {"source": "provided"}, umap={"enabled": False})).adata
    pi = a.uns["petrubseq_protein"]["perturbations"]
    assert pi["source"] == "guide_assignments" and pi["assignment"]["effective_source"] == "provided"
    assert "guide_counts" in a.obsm and a.obsm["guide_counts"].shape[1] == len(raw["guides"])
    agr = pi["assignment"]["agreement"]
    assert agr["set_agreement"]["disagree"] == 0 and agr["dominant_vs_provided_single"]["disagree"] == 0
    assert set(guidecall.DIAG_COLUMNS) <= set(a.obs.columns)
    assert list(a.obs["guides"].astype(str)) == [";".join(l) for l in raw["lists"]]


def test_guide_source_auto_with_both_sources_fails(raw, tmp_path):
    inputs = {"protein": NO_PROTEIN_NORM, "guide_counts": {"file": "guide_counts.csv.gz"}, **COMMON}
    with pytest.raises(InputError, match="perturbation.assignment.source"):
        run_pipeline(_cfg(raw, tmp_path, inputs, {"source": "auto"}, umap={"enabled": False}))


def test_source_guide_counts_without_counts_rejected(raw, tmp_path):
    inputs = {"protein": NO_PROTEIN_NORM, **COMMON}
    with pytest.raises(FileNotFoundError, match="guide_counts"):
        validation.validate_input_files(_cfg(raw, tmp_path, inputs, {"source": "guide_counts"}))


def test_dominant_rule_marks_ambiguous_and_keeps_detected_lists(raw, tmp_path):
    inputs = variants(raw)["mtx_combined"]
    a = run_pipeline(_cfg(raw, tmp_path, inputs, {"source": "guide_counts", "method": "dominant"}, umap={"enabled": False})).adata
    obs = a.obs
    amb = obs["perturbation_class"] == "ambiguous"
    assert amb.any() and (obs.loc[amb, "perturbation"] == "ambiguous").all() and (obs.loc[amb, "guide_dominant_call"] == "ambiguous").all()
    assert (obs.loc[amb, "n_guides_detected"] >= 2).all() and (obs.loc[amb, "guides"] == "").all()
    single = obs["perturbation_class"].isin(["single_targeting", "single_control"])
    assert (obs.loc[single, "guide"].astype(str) == obs.loc[single, "guide_dominant_call"].astype(str)).all()
    truth = [";".join(l) for l in raw["lists"]]
    assert list(obs["guides_detected"].astype(str)) == truth  # detected lists are not collapsed
    unassigned = obs["perturbation_class"] == "unassigned"
    assert (obs.loc[unassigned, "guide_top_count"] < 3).all()
    info = a.uns["petrubseq_protein"]["perturbations"]["assignment"]["guide_calling"]
    assert info["dominant_rule"]["ambiguous"] == int(amb.sum())
    assert a.uns["petrubseq_protein"]["perturbations"]["qc"]["class_counts"]["ambiguous"] == int(amb.sum())


def test_call_guides_rules():
    from petrubseq_protein.io.readers import Matrix
    import scipy.sparse as sp

    X = sp.csr_matrix(np.array([[10, 1, 0], [10, 6, 0], [2, 1, 0], [0, 0, 0], [5, 5, 5]], dtype=np.float32))
    gm = Matrix(np.array(["c1", "c2", "c3", "c4", "c5"], dtype=object), np.array(["A_1", "B_1", "NO_SITE_1"], dtype=object), X)
    cfg = Config.from_dict({"inputs": {"rna": {"file": "x"}}})
    calls = guidecall.call_guides(gm, ["c1", "c2", "c3", "c4", "c5", "c6"], cfg)
    d = calls.diagnostics
    assert list(d["guide_dominant_call"]) == ["A_1", "ambiguous", "unassigned", "unassigned", "ambiguous", "unassigned"]
    assert list(d["n_guides_detected"]) == [1, 2, 0, 0, 3, 0]
    assert list(d["has_guide_counts"]) == [True] * 5 + [False]
    assert calls.X.shape == (6, 3) and calls.X[5].nnz == 0
    assert list(calls.features["target"]) == ["A", "B", "NO_SITE"] and list(calls.features["control_class"]) == ["", "", "non_targeting"]
    cfg2 = Config.from_dict({"inputs": {"rna": {"file": "x"}}, "perturbation": {"assignment": {"method": "threshold", "detection_min_umi": 5}}})
    calls2 = guidecall.call_guides(gm, ["c1", "c2", "c5"], cfg2)
    assert [list(l) for l in calls2.lists] == [["A_1"], ["A_1", "B_1"], ["A_1", "B_1", "NO_SITE_1"]]


# ---------------------------------------------------------- multi-lane ---


def test_multilane_ids_and_lane_metadata(raw, tmp_path):
    L = raw["lanes"]
    inputs = {"multiplexed": {"format": "mtx", "lanes": [{"id": "lane1", "path": str(L["dir"] / "lane1")}, {"id": "lane2", "path": str(L["dir"] / "lane2")}], "var_names": "name"},
              "lane_metadata": {"file": str(L["dir"] / "lane_metadata.csv")}, "metadata": {"file": str(L["dir"] / "metadata_lanes.csv"), "format": "auto"},
              "guide_assignments": {"file": str(L["dir"] / "guides_lanes.txt")}, "embedding": {"file": None}, **NO_DENSE}
    cfg = _cfg(raw, tmp_path, inputs, {"source": "provided"}, umap={"enabled": False})
    a = run_pipeline(cfg).adata
    assert a.n_obs == 200 and list(a.obs_names) == L["global_ids"]
    assert a.obs["lane_id"].value_counts().to_dict() == {"lane1": 100, "lane2": 100}
    assert a.obs["barcode_original"].iloc[0] == "BC0000-1" and a.obs["barcode_original"].nunique() == 100
    assert (a.obs.loc[a.obs["lane_id"] == "lane1", "condition"] == "Control").all() and (a.obs.loc[a.obs["lane_id"] == "lane2", "condition"] == "Treated").all()
    assert list(a.obs["guides"].astype(str)) == [";".join(l) for l in raw["lists"]]
    prov = a.uns["petrubseq_protein"]["inputs"]["provenance"]
    assert prov["lanes"]["n_lanes"] == 2 and "__" in prov["cell_id_transformation"]
    assert "library" in a.obs.columns or "library" not in cfg.columns.keep  # lane sheet columns are available to columns.*


def test_single_lane_not_suffixed(processed):
    assert "lane_id" not in processed["mtx_combined"].obs.columns
    assert processed["mtx_combined"].obs_names[0] == "CELL_1"


def test_lane_errors(raw, tmp_path):
    L = raw["lanes"]
    lanes = [{"id": "lane1", "path": str(L["dir"] / "lane1")}, {"id": "lane2", "path": str(L["dir"] / "lane2")}]
    with pytest.raises(ConfigError, match="duplicate lane id"):
        _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "lanes": [lanes[0], {"id": "lane1", "path": lanes[1]["path"]}]}, **NO_DENSE})
    bad = tmp_path / "lm_dup.csv"
    bad.write_text("lane_id,condition\nlane1,A\nlane1,B\nlane2,C\n")
    cfg = _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "lanes": lanes}, "lane_metadata": {"file": str(bad)}, "metadata": {"file": None}, "guide_assignments": {"file": None}, "embedding": {"file": None}, **NO_DENSE}, columns={"condition": None, "moi": None, "rna_total_counts": None, "guide": None})
    with pytest.raises(InputError, match="duplicate lane rows"):
        adapters.load_inputs(cfg)
    bad.write_text("lane_id,condition\nlane1,A\n")
    with pytest.raises(InputError, match="no row for lane"):
        adapters.load_inputs(cfg)
    with pytest.raises(ConfigError, match="lane_metadata needs multi-lane"):
        _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": lanes[0]["path"]}, "lane_metadata": {"file": str(bad)}, **NO_DENSE})


# ------------------------------------------------------------ validation ---


def test_conflicting_declarations_rejected(raw, tmp_path):
    f = raw["formats"]
    with pytest.raises(ConfigError, match="declared both"):
        _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"])}, "protein": NO_PROTEIN_NORM})
    with pytest.raises(ConfigError, match="declared both"):
        _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"])}, "guide_counts": {"format": "mtx", "file": str(f["mtx_guide"])}, **NO_DENSE})


def test_protein_availability_rule(raw, tmp_path):
    f = raw["formats"]
    # counts only, protein required (default): fine
    validation.validate_input_files(_cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"])}, **NO_DENSE}))
    # no protein at all but required: error, independent of the format
    with pytest.raises(FileNotFoundError, match="protein is required"):
        validation.validate_input_files(_cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"]), "feature_types": {"protein": []}}, **NO_DENSE, "protein": {"file": None, "required": True}}))
    with pytest.raises(FileNotFoundError, match="protein is required"):
        validation.validate_input_files(_cfg(raw, tmp_path, {"protein": {"file": None, "required": True}, "protein_counts": {"file": None}}))
    # "protein" in alignment.required is satisfied by counts alone
    inputs = {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"])}, **NO_DENSE, "guide_assignments": {"file": None}, **COMMON}
    cfg = _cfg(raw, tmp_path, inputs, alignment={"required": ["rna", "protein", "metadata"]})
    ci = adapters.load_inputs(cfg)
    from petrubseq_protein.preprocessing.align import align_cells, effective_required
    req = effective_required(cfg.alignment.required, ci.id_sets())
    assert req == ["rna", "protein_counts", "metadata"]
    assert align_cells(ci.id_sets(), cfg, required=req).n_kept == 200


def test_unclaimed_feature_type_rejected(raw, tmp_path):
    f = raw["formats"]
    cfg = _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"]), "feature_types": {"guide": []}}, **NO_DENSE, **COMMON, "guide_assignments": {"file": None}}, {"source": "auto"})
    with pytest.raises(InputError, match="assigned to no modality"):
        adapters.load_inputs(cfg)
    cfg = _cfg(raw, tmp_path, {"multiplexed": {"format": "mtx", "path": str(f["mtx_combined"]), "feature_types": {"protein": ["Antibody Capture"], "guide": ["CRISPR Guide Capture"], "rna": ["Gene Expressionn"]}}, **NO_DENSE})
    with pytest.raises(InputError, match="no features of type"):
        adapters.load_inputs(cfg)


def test_h5ad_missing_slot_rejected(raw, tmp_path):
    f = raw["formats"]
    cfg = _cfg(raw, tmp_path, {"multiplexed": {"format": "h5ad", "path": str(f["h5ad_slots"]), "slots": {"rna": {"slot": "layers", "key": "nope"}}}, **NO_DENSE})
    with pytest.raises(InputError, match="layers\\['nope'\\]"):
        adapters.load_inputs(cfg)
    cfg = _cfg(raw, tmp_path, {"rna": {"format": "h5ad", "file": str(f["h5ad_slots"]), "files": [], "slot": "obsm", "key": "guide_counts"}, "protein": NO_PROTEIN_NORM, "protein_counts": {"file": None}})
    ci_rna = adapters.read_matrix_input(cfg, "rna", cfg.inputs.rna, 40, 1)[0]
    assert ci_rna.shape == (200, len(raw["guides"]))  # unlabelled obsm array + uns['<key>_features'] works


def test_negative_and_nonint_counts_rejected(raw, tmp_path):
    f = raw["formats"]
    a = ad.read_h5ad(f["h5ad_slots"])
    a.layers["counts"] = a.layers["counts"].copy()
    a.layers["counts"][0, 0] = -1
    a.obsm["protein_counts"].iloc[0, 0] = 0.5
    p = tmp_path / "bad.h5ad"
    a.write_h5ad(p)
    cfg = _cfg(raw, tmp_path, {"multiplexed": {"format": "h5ad", "path": str(p), "rna_state": "raw_counts", "slots": {"rna": {"slot": "layers", "key": "counts"}, "protein_counts": {"slot": "obsm", "key": "protein_counts"}}}, **NO_DENSE, **COMMON, "guide_assignments": {"file": None}})
    errs, _ = validation.validate_canonical(adapters.load_inputs(cfg), cfg)
    assert any("rna: declared raw counts but contains negative" in e for e in errs)
    assert any("protein_counts: declared raw counts but contains non-integer" in e for e in errs)


def test_metadata_mismatch_and_bad_guide_table(raw, tmp_path):
    L = raw["lanes"]
    inputs = {"multiplexed": {"format": "mtx", "lanes": [{"id": "lane1", "path": str(L["dir"] / "lane1")}, {"id": "lane2", "path": str(L["dir"] / "lane2")}]}, "metadata": {"file": "metadata.csv", "format": "scp_metadata"}, "guide_assignments": {"file": None}, "embedding": {"file": None}, **NO_DENSE}
    cfg = _cfg(raw, tmp_path, inputs, {"source": "auto"})
    errs, _ = validation.validate_canonical(adapters.load_inputs(cfg), cfg)
    assert any("shares no cell ID" in e for e in errs)
    gt = tmp_path / "gt.csv"
    gt.write_text("guide,target\nGENEA_1,GENEA\nGENEA_1,GENEB\n")
    cfg = _cfg(raw, tmp_path, variants(raw)["dense"], perturbation={"guide_target_table": str(gt), "assignment": {"source": "guide_counts"}}) if False else Config.from_dict(base_config(raw, tmp_path, inputs=variants(raw)["dense"], columns={"guide": None}, perturbation={"guide_target_table": str(gt), "assignment": {"source": "guide_counts"}}))
    errs, _ = validation.validate_canonical(adapters.load_inputs(cfg), cfg)
    assert any("duplicate guide IDs" in e for e in errs)


def test_guide_name_mismatch_warns(raw, tmp_path):
    gfile = tmp_path / "guides_odd.txt"
    gfile.write_text("Cell,sgRNAs\n" + "\n".join(f"{c},{'FOO_9' if i == 0 else ''}" for i, c in enumerate(raw["cells"])) + "\n")
    inputs = {"protein": NO_PROTEIN_NORM, "guide_counts": {"file": "guide_counts.csv.gz"}, "guide_assignments": {"file": str(gfile)}, **COMMON}
    cfg = _cfg(raw, tmp_path, inputs, {"source": "provided"})
    errs, warns = validation.validate_canonical(adapters.load_inputs(cfg), cfg)
    assert errs == [] and any("absent from the guide-count features" in w for w in warns)


# ------------------------------------------------- protein annotation ---


def test_protein_feature_annotation(raw, tmp_path):
    table = tmp_path / "antibodies.csv"
    table.write_text("feature_id,antibody_name,protein_name,gene_symbol,clone,isotype\nCD1,anti-CD1,CD1 antigen,CD1A,HI149,Mouse IgG1\nCD2,anti-CD2,,CD2,RPA-2.10,\nNOPE,x,x,x,x,x\n")
    inputs = variants(raw)["mtx_combined"]
    a = run_pipeline(_cfg(raw, tmp_path, inputs, umap={"enabled": False}, protein={"feature_table": str(table)})).adata
    pf = a.uns["protein_features"]
    for c in ("feature_id", "antibody_name", "protein_name", "gene_symbol", "clone", "feature_type", "isotype", "isotype_control", "annotation_source"):
        assert c in pf.columns
    assert pf.loc["CD1", "clone"] == "HI149" and pf.loc["CD1", "gene_symbol"] == "CD1A" and pf.loc["CD1", "isotype"] == "Mouse IgG1"
    assert pf.loc["CD1", "annotation_source"] == "input_features+feature_table"
    assert pf.loc["CD2", "protein_name"] == "" and pf.loc["CD2", "isotype"] == ""  # never invented
    assert pf.loc["CD3", "annotation_source"] == "input_features" and pf.loc["CD3", "clone"] == ""
    assert (pf["feature_type"] == "Antibody Capture").all() and (pf["feature_id"] == pf.index).all()
    assert pf.loc["Mouse_IgG1", "role"] == "isotype_control" and not pf.loc["Mouse_IgG1", "used_in_embedding"]
    assert set(a.obsm["protein"].columns) == set(PROTEINS)


def test_dense_annotation_is_na_without_sources(processed):
    pf = processed["dense"].uns["protein_features"]
    assert (pf["annotation_source"] == "none").all() and (pf["antibody_name"] == "").all() and (pf["feature_type"] == "").all()
    assert (pf["feature_id"] == pf.index).all()


# --------------------------------------------------------------- audit ---


def test_audit_is_format_aware(raw, tmp_path):
    cfg = _cfg(raw, tmp_path, variants(raw)["mtx_combined"])
    rec = audit.audit_config(cfg)
    assert rec["roles"]["multiplexed"]["format"] == "mtx" and rec["input_layout"].startswith("feature-barcode")
    assert rec["cell_ids"]["n_common"] == 200 if "n_common" in rec["cell_ids"] else True
    assert rec["guide_counts"]["n_guides"] == len(raw["guides"]) and rec["guide_counts"]["method"] == "threshold"
    assert rec["rna"]["state"] == "raw_counts" and rec["rna"]["format"] == "mtx"
    assert "counts_matrix" in rec["protein"]
    md = audit.render_markdown(rec)
    assert "guide" in md.lower()


def test_scp1064_configs_still_valid():
    root = Path(__file__).resolve().parents[1]
    for name in ("scp1064.yaml", "scp1064_smoke.yaml"):
        cfg = Config.from_yaml(root / "config" / name)
        assert cfg.perturbation.assignment.source == "auto" and not cfg.inputs.multiplexed.is_set()
        assert cfg.inputs.rna.format == "dense_csv" and cfg.inputs.guide_counts.file is None
