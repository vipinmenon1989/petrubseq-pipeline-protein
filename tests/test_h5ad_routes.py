"""End-to-end routes from an existing AnnData: RNA + guides (no protein) and
RNA + guides + protein, with cell metadata and guide calls taken from obs;
clustering on by default; the documented example configs validate."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import base_config, make_dataset, write_h5ad_slots  # noqa: E402

from petrubseq_protein.config import Config  # noqa: E402
from petrubseq_protein.io import adapters  # noqa: E402
from petrubseq_protein.io.h5ad import read_h5ad_obs  # noqa: E402
from petrubseq_protein.pipeline import run_pipeline  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
NO_DENSE = {"rna": {"files": [], "file": None}, "protein_counts": {"file": None}, "protein": {"file": None, "required": False}, "guide_assignments": {"file": None}, "embedding": {"file": None}}


@pytest.fixture(scope="module")
def h5(tmp_path_factory):
    root = tmp_path_factory.mktemp("h5ad_routes")
    d = make_dataset(root / "dense", n_cells=300, rna_state="raw", protein_norm=False)
    d["h5ad"] = write_h5ad_slots(d, root / "screen.h5ad")
    return d


def _h5_inputs(path, protein: bool, guide_obs: bool):
    inp = {**NO_DENSE,
           "rna": {"format": "h5ad", "file": str(path), "files": [], "slot": "layers", "key": "counts", "state": "raw_counts"},
           "guide_counts": {"format": "h5ad", "file": str(path), "slot": "obsm", "key": "guide_counts"},
           "metadata": {"file": str(path), "format": "h5ad"}}
    if protein:
        inp["protein_counts"] = {"format": "h5ad", "file": str(path), "slot": "obsm", "key": "protein_counts"}
        inp["protein"] = {"file": None, "required": True}
    if guide_obs:
        inp["guide_assignments"] = {"file": str(path), "format": "h5ad", "guides_column": "guides", "list_separator": ";"}
    return inp


def test_read_h5ad_obs(h5):
    obs = read_h5ad_obs(h5["h5ad"])
    assert list(obs.index) == h5["cells"] and {"condition", "UMI_count", "guides"} <= set(obs.columns)
    obs2 = read_h5ad_obs(h5["h5ad"], "condition")
    assert obs2.index.name == "condition" and set(obs2.index) == {"Control", "Treated"}


def test_h5ad_rna_guides_no_protein(h5, tmp_path):
    cfg = Config.from_dict(base_config(h5, tmp_path, inputs=_h5_inputs(h5["h5ad"], protein=False, guide_obs=False),
                                       columns={"cell_id": "obs_names", "guide": None, "moi": None, "rna_total_counts": "UMI_count"},
                                       alignment={"required": ["rna", "metadata"]}, perturbation={"assignment": {"source": "auto"}}, umap={"enabled": False},
                                       analysis={"clustering": {"enrichment": {"n_permutations": 20}}, "perturbation_effects": {"enabled": False}}))
    ci = adapters.load_inputs(cfg)
    assert ci.protein_counts is None and ci.protein_normalized is None and ci.guide_assignments is None
    assert ci.rna.state == "raw_counts" and ci.metadata.shape[0] == 300 and "condition" in ci.metadata.columns
    a = run_pipeline(cfg).adata
    assert a.obs["condition"].nunique() == 2 and "protein" not in a.obsm and "guide_counts" in a.obsm
    assert a.uns["petrubseq_protein"]["perturbations"]["assignment"]["effective_source"] == "guide_counts"
    assert "leiden" in a.obs and "cell_states" in a.uns["petrubseq_protein"]["analysis"]  # clustering on by default
    assert a.obs["leiden"].dtype.name == "category"


def test_h5ad_rna_guides_protein_obs_calls(h5, tmp_path):
    cfg = Config.from_dict(base_config(h5, tmp_path, inputs=_h5_inputs(h5["h5ad"], protein=True, guide_obs=True),
                                       columns={"cell_id": "obs_names", "guide": None, "moi": None, "rna_total_counts": "UMI_count"},
                                       perturbation={"assignment": {"source": "provided"}}, umap={"enabled": False},
                                       analysis={"clustering": {"enrichment": {"n_permutations": 20}}, "perturbation_effects": {"enabled": True, "lochness": {"n_permutations": 5}, "modules": {"min_perturbations": 3, "min_cells_per_perturbation": 10}}}))
    ci = adapters.load_inputs(cfg)
    assert ci.guide_assignments is not None and list(ci.guide_assignments.index) == h5["cells"]
    assert [list(l) for l in ci.guide_assignments] == h5["lists"]
    assert ci.provenance["modalities"]["guide_assignments"]["source"].endswith("obs['guides']")
    a = run_pipeline(cfg).adata
    pi = a.uns["petrubseq_protein"]["perturbations"]
    assert pi["source"] == "guide_assignments" and pi["assignment"]["effective_source"] == "provided" and "agreement" in pi["assignment"]
    assert list(a.obs["guides"].astype(str)) == [";".join(l) for l in h5["lists"]]
    assert "protein" in a.obsm and a.uns["petrubseq_protein"]["protein"]["primary_method"].startswith("clr")
    assert "protein_effects" in a.uns and "ps_scores" in a.obsm and "leiden" in a.obs


def test_clustering_default_on_and_explicit_off(h5, tmp_path):
    on = Config.from_dict(base_config(h5, tmp_path))
    off = Config.from_dict(base_config(h5, tmp_path, analysis={"clustering": {"enabled": False}}))
    assert on.analysis.clustering.enabled and not off.analysis.clustering.enabled
    for name in ("scp1064.yaml", "scp1064_smoke.yaml"):
        c = Config.from_yaml(REPO / "config" / name)
        assert not c.analysis.clustering.enabled and not c.analysis.perturbation_effects.enabled  # frozen-baseline semantics


@pytest.mark.parametrize("name", ["h5ad_perturbseq.yaml", "h5ad_perturb_cite_seq.yaml", "h5ad_slots.yaml"])
def test_h5ad_examples_validate(name):
    c = Config.from_yaml(REPO / "config" / "examples" / name)
    assert c.analysis.clustering.enabled
    if name != "h5ad_slots.yaml":
        assert c.inputs.rna.format == "h5ad" and c.inputs.metadata.format == "h5ad" and c.inputs.guide_counts.slot == "obsm"


def test_notebook_expects_existing_demo_outputs():
    """The notebook reads only files the CLI writes; check the names it uses are produced."""
    import json

    nb = json.load(open(REPO / "notebooks" / "demo_papalexi_eccite.ipynb"))
    src = "\n".join("".join(c["source"]) if isinstance(c["source"], list) else c["source"] for c in nb["cells"] if c["cell_type"] == "code")
    for token in ("papalexi_eccite_demo_processed.h5ad", "perturbation_cluster_enrichment.csv", "cluster_summary.csv", "ps_targets.csv", "lochness_targets.csv", "cluster_enrichment_vs_lochness.csv",
                  "perturbation_modules.csv", "gene_programs.csv", "protein_effects.csv", "ps_protein_association.csv", "perturbation_summary.csv", "report.html"):
        assert token in src, token
    for stale in ("Stage C", "preprocessing/QC report", "Stage E", "Stage F"):
        assert stale not in src
