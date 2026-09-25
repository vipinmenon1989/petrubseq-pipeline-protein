"""Stage D: Papalexi ECCITE-seq demo preparation (synthetic GEO-like fixture) and,
when the local source data exist, a skippable real-data integration check."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PREP = REPO / "demo" / "prepare_papalexi_demo.py"
CHECK = REPO / "demo" / "check_demo_reproducibility.py"
sys.path.insert(0, str(REPO / "demo"))
import prepare_papalexi_demo as prep  # noqa: E402

from petrubseq_protein.config import Config  # noqa: E402
from petrubseq_protein.io import adapters  # noqa: E402
from petrubseq_protein import validation  # noqa: E402

LOCAL_SOURCE = Path(os.environ.get("PETRUBSEQ_DATA_ROOT", REPO / ".." / "data")) / "ECCITE-seq"


def _write_tsv(path: Path, df: pd.DataFrame) -> None:
    with gzip.open(path, "wt") as fh:
        df.to_csv(fh, sep="\t", quoting=2)


@pytest.fixture(scope="module")
def geo_like(tmp_path_factory):
    """Tiny pooled-screen look-alike: 8 lanes, 12 HTOs, 4 ADTs, 3 targets x 2 guides + NT."""
    rng = np.random.default_rng(1)
    src = tmp_path_factory.mktemp("geo_src")
    n = 400
    cells = [f"l{1 + i % 8}_{''.join(rng.choice(list('ACGT'), 16))}" for i in range(n)]
    genes = [f"G{i:03d}" for i in range(60)] + ["MT-CO1", "MT-ND1"]
    rna = pd.DataFrame(rng.poisson(2, (len(genes), n)), index=genes, columns=cells)
    adt = pd.DataFrame(rng.poisson(50, (4, n)), index=["CD86", "PDL1", "PDL2", "CD366"], columns=cells)
    htos = ["rep1-tx", "rep1-ctrl", "rep2-tx", "rep2-ctrl", "PDL1g1-tx", "PDL1g1-ctrl", "PDL1g2-tx", "PDL1g2-ctrl", "rep3-tx", "rep3-ctrl", "rep4-tx", "rep4-ctrl"]
    hto = pd.DataFrame(rng.poisson(1, (12, n)), index=htos, columns=cells)
    for j in range(n):
        hto.iloc[[0, 8, 10][j % 3], j] += 60  # three real samples, all 'tx'
    hto.iloc[1, :10] += 60  # first 10 cells: doublets (two strong HTOs)
    guides = [f"{t}g{k}" for t in ("AAA", "BBB", "CCC") for k in (1, 2)] + ["NTg1", "NTg2", "eGFPg1"]
    gdo = pd.DataFrame(rng.poisson(0.2, (len(guides), n)), index=guides, columns=cells)
    for j in range(n):
        if j % 10 == 9:
            gdo.iloc[[0, 2], j] += 20  # ambiguous: two comparable guides
        else:
            gdo.iloc[j % 8, j] += 20
    for k, df in (("rna", rna), ("adt", adt), ("hto", hto), ("gdo", gdo)):
        _write_tsv(src / prep.SOURCES[k], df)
    with gzip.open(src / prep.SOURCES["adt_barcodes"], "wt") as fh:
        fh.write("".join(f"SEQ{i},{a}\n" for i, a in enumerate(adt.index)))
    with gzip.open(src / prep.SOURCES["gdo_barcodes"], "wt") as fh:
        fh.write("".join(f"SEQ{i},{g}\n" for i, g in enumerate(guides)))
    return {"src": src, "cells": cells, "genes": genes, "rna": rna, "adt": adt, "gdo": gdo, "guides": guides}


def _run_prep(src: Path, out: Path, n_cells: int = 120, seed: int = 0) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(PREP), "--source", str(src), "--output", str(out), "--n-cells", str(n_cells), "--seed", str(seed), "--no-source-hash", "--min-per-stratum", "3"], capture_output=True, text=True)


@pytest.fixture(scope="module")
def prepared(geo_like, tmp_path_factory):
    out = tmp_path_factory.mktemp("demo_out") / "papalexi_eccite"
    r = _run_prep(geo_like["src"], out)
    assert r.returncode == 0, r.stderr[-2000:]
    return out


def test_outputs_and_manifest(prepared, geo_like):
    for f in ("filtered_feature_bc_matrix/barcodes.tsv.gz", "filtered_feature_bc_matrix/features.tsv.gz", "filtered_feature_bc_matrix/matrix.mtx.gz", "cells.csv", "antibodies.csv", "guide_targets.csv", "selected_cells.txt", "demo_manifest.json"):
        assert (prepared / f).is_file(), f
    m = json.load(open(prepared / "demo_manifest.json"))
    assert m["selected"]["n_cells"] == 120 and m["parameters"]["seed"] == 0
    assert m["derived"]["filtered_feature_bc_matrix"] == {"n_features": 62 + 4 + 9, "n_cells": 120, "gene_expression": 62, "antibody_capture": 4, "crispr_guide_capture": 9, "nnz": m["derived"]["filtered_feature_bc_matrix"]["nnz"]}
    assert set(m["sources"]) == set(prep.SOURCES) and all("size_bytes" in v for v in m["sources"].values())
    assert m["source_summary"]["hto_classes"]["doublet_or_negative"] >= 10
    sel = (prepared / "selected_cells.txt").read_text().split()
    assert len(sel) == 120 and len(set(sel)) == 120 and all(c in set(geo_like["cells"]) for c in sel)


def test_values_unchanged_and_strata_preserved(prepared, geo_like):
    from petrubseq_protein.io.tenx import read_mtx_dir

    fm = read_mtx_dir(prepared / "filtered_feature_bc_matrix")
    cells = list(fm.cells)
    rna = geo_like["rna"][cells].to_numpy().T
    adt = geo_like["adt"][cells].to_numpy().T
    gdo = geo_like["gdo"][cells].to_numpy().T
    X = fm.X.toarray()
    t = fm.var["feature_types"].to_numpy()
    np.testing.assert_array_equal(X[:, t == "Gene Expression"], rna)
    np.testing.assert_array_equal(X[:, t == "Antibody Capture"], adt)
    np.testing.assert_array_equal(X[:, t == "CRISPR Guide Capture"], gdo)
    assert list(fm.features[t == "Gene Expression"]) == geo_like["genes"]
    meta = pd.read_csv(prepared / "cells.csv", index_col=0)
    assert list(meta.index) == cells and (meta["hto_class"] == "singlet").all()
    assert set(meta["hto_sample"]) == {"rep1-tx", "rep3-tx", "rep4-tx"} and (meta["stimulation"] == "IFNg").all()
    targets = set(meta["sampling_stratum_target"])
    assert {"AAA", "BBB", "CCC", "NT", "ambiguous"} <= targets  # every target, controls and ambiguous cells survive
    strata = meta.groupby(["hto_sample", "sampling_stratum_target"]).size()
    assert (strata >= 3).all()
    gt = pd.read_csv(prepared / "guide_targets.csv")
    assert dict(zip(gt["guide"], gt["target"]))["NTg2"] == "NT" and dict(zip(gt["guide"], gt["target"]))["AAAg1"] == "AAA"
    ab = pd.read_csv(prepared / "antibodies.csv").fillna("")
    assert list(ab["feature_id"]) == ["CD86", "PDL1", "PDL2", "CD366"] and (ab["clone"] == "").all() and ab.loc[1, "gene_symbol"] == "CD274"


def test_deterministic_and_source_untouched(prepared, geo_like, tmp_path):
    before = {p.name: (p.stat().st_size, p.stat().st_mtime) for p in geo_like["src"].iterdir()}
    out2 = tmp_path / "again"
    assert _run_prep(geo_like["src"], out2).returncode == 0
    r = subprocess.run([sys.executable, str(CHECK), str(prepared), str(out2)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout
    out3 = tmp_path / "seed1"
    assert _run_prep(geo_like["src"], out3, seed=1).returncode == 0
    assert (out3 / "selected_cells.txt").read_text() != (prepared / "selected_cells.txt").read_text()
    after = {p.name: (p.stat().st_size, p.stat().st_mtime) for p in geo_like["src"].iterdir()}
    assert before == after
    assert _run_prep(geo_like["src"], geo_like["src"] / "inside").returncode != 0  # output inside source refused


def test_demo_config_loads_through_stage_c_adapter(prepared):
    cfg = Config.from_yaml(REPO / "config" / "demo_papalexi.yaml")
    d = cfg.to_dict()
    d["dataset"]["input_dir"] = str(prepared)
    d["output"]["dir"] = str(prepared.parent / "out")
    cfg = Config.from_dict(d)
    assert cfg.inputs.multiplexed.format == "mtx" and cfg.inputs.multiplexed.modalities() == ["rna", "protein_counts", "guide_counts"]
    validation.validate_input_files(cfg)
    ci = adapters.load_inputs(cfg)
    errs, warns = validation.validate_canonical(ci, cfg)
    assert errs == []
    assert ci.rna.shape == (120, 62) and ci.protein_counts.shape == (120, 4) and ci.guide_counts.shape == (120, 9)
    assert ci.rna.state == "raw_counts" and ci.protein_normalized is None and ci.guide_assignments is None
    assert ci.provenance["multiplexed"]["format"] == "mtx" and set(ci.metadata.columns) >= {"lane", "hto_sample", "stimulation"}


def test_demo_pipeline_on_fixture(prepared, tmp_path):
    from petrubseq_protein.pipeline import run_pipeline

    cfg = Config.from_yaml(REPO / "config" / "demo_papalexi.yaml")
    d = cfg.to_dict()
    d["dataset"]["input_dir"] = str(prepared)
    d["output"]["dir"] = str(tmp_path / "out")
    d["umap"]["enabled"] = False
    d["rna"]["n_pcs"] = 5
    d["rna"]["hvg"]["n_top_genes"] = 30
    d["qc"]["prefilter"]["min_genes_per_cell"] = 5
    d["qc"]["filter"]["rna"]["min_genes"] = 5
    d["compute"]["n_jobs"] = 1
    a = run_pipeline(Config.from_dict(d)).adata
    pi = a.uns["petrubseq_protein"]["perturbations"]
    assert pi["assignment"]["effective_source"] == "guide_counts" and pi["source"] == "guide_counts:dominant"
    assert "guide_counts" in a.obsm and a.obsm["guide_counts"].shape[1] == 9
    assert set(a.obs["perturbation_class"]) >= {"single_targeting", "single_control", "ambiguous"}
    assert (a.obs.loc[a.obs["perturbation_class"] == "single_control", "control_class"] == "non_targeting").all()
    assert a.uns["protein_features"]["is_isotype"].sum() == 0 and (a.uns["protein_features"]["gene_symbol"] != "").all()
    assert a.obs["condition"].nunique() == 1 and a.obs["sample"].nunique() == 3 and a.obs["lane"].nunique() == 8


@pytest.mark.skipif(not (LOCAL_SOURCE / prep.SOURCES["gdo"]).is_file(), reason="local ECCITE-seq source data not available")
def test_local_source_inventory_matches_audit():
    """Lightweight: file sizes only (no matrix is loaded on the login node)."""
    sys.path.insert(0, str(REPO / "demo"))
    import fetch_papalexi_data as fetch

    for gsm, name, size in fetch.POOLED:
        p = LOCAL_SOURCE / name
        assert p.is_file() and p.stat().st_size == size, name


# ------------------------------------------------------- bundled demo ---

BUNDLE = REPO / "demo" / "data" / "papalexi_eccite"


def test_bundled_demo_is_the_default_input():
    cfg = Config.from_yaml(REPO / "config" / "demo_papalexi.yaml")
    assert Path(cfg.dataset.input_dir).resolve() == BUNDLE.resolve()
    assert "results" in Path(cfg.output.dir).parts and REPO not in Path(cfg.output.dir).resolve().parents
    validation.validate_input_files(cfg)  # every bundled input exists
    m = json.load(open(BUNDLE / "demo_manifest.json"))
    assert m["parameters"]["seed"] == 0 and m["selected"]["n_cells"] == 1800 and m["accession"] == "GSE153056"
    assert m["derived"]["filtered_feature_bc_matrix"] == {"n_features": 18764, "n_cells": 1800, "gene_expression": 18649, "antibody_capture": 4, "crispr_guide_capture": 111, "nnz": 6347995}
    sel = (BUNDLE / "selected_cells.txt").read_text().split()
    with gzip.open(BUNDLE / "filtered_feature_bc_matrix" / "barcodes.tsv.gz", "rt") as fh:
        bcs = fh.read().split()
    assert sel == bcs and len(bcs) == 1800
    with gzip.open(BUNDLE / "filtered_feature_bc_matrix" / "features.tsv.gz", "rt") as fh:
        types = [line.rstrip("\n").split("\t")[2] for line in fh]
    assert types.count("Gene Expression") == 18649 and types.count("Antibody Capture") == 4 and types.count("CRISPR Guide Capture") == 111
    assert len(pd.read_csv(BUNDLE / "guide_targets.csv")) == 111 and len(pd.read_csv(BUNDLE / "cells.csv")) == 1800
    total = sum(p.stat().st_size for p in BUNDLE.rglob("*") if p.is_file())
    assert total < 25_000_000  # the bundle must stay small
