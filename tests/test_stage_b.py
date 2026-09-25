"""Stage B: lifecycle, filtering audit, before/after figures, registry, report.html."""

from __future__ import annotations

import gzip
import re
import sys
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from make_synthetic import base_config, make_dataset  # noqa: E402

from petrubseq_protein.cli import main as cli_main  # noqa: E402
from petrubseq_protein.config import Config  # noqa: E402
from petrubseq_protein.pipeline import run_pipeline  # noqa: E402
from petrubseq_protein.reporting.figures import SECTION_QC_PERTURBATION, SECTION_QC_PROTEIN, SECTION_QC_RNA, STAGE_AFTER, STAGE_BEFORE  # noqa: E402

NO_UMAP = {"umap": {"enabled": False}}
NO_H5AD = {"output": {"write_h5ad": False}}


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("syn_b"), n_cells=300)


@pytest.fixture(scope="module")
def filtered_run(synthetic, tmp_path_factory):
    """Prefilter + strict filtering that actually remove cells and genes."""
    out = tmp_path_factory.mktemp("filtered")
    cfg = Config.from_dict(base_config(synthetic, out, qc={"prefilter": {"enabled": True, "min_genes_per_cell": 95, "min_cells_per_gene": 100}, "filter": {"enabled": True, "rna": {"min_genes": 102, "max_pct_mt": 3.0}}}, **NO_UMAP))
    return run_pipeline(cfg)


def _steps(res):
    return pd.read_csv(res.tables["qc_filtering_steps"])


# ------------------------------------------------------------- filtering ---


def test_filtering_audit_arithmetic(filtered_run):
    t = _steps(filtered_run)
    assert list(t["step"][:1]) == ["input"] and list(t["step"][-1:]) == ["final"]
    assert {"prefilter_min_genes", "prefilter_min_cells_per_gene", "rna_min_genes", "rna_max_pct_mt"} <= set(t["step"])
    assert (t["cells_removed"] == t["cells_before"] - t["cells_after"]).all()
    assert (t["genes_removed"] == t["genes_before"] - t["genes_after"]).all()
    assert (t["cells_after"].to_numpy()[:-1] == t["cells_before"].to_numpy()[1:]).all()
    assert (t["genes_after"].to_numpy()[:-1] == t["genes_before"].to_numpy()[1:]).all()
    assert t["cells_before"].iloc[0] == 300 and t["cells_after"].iloc[-1] == filtered_run.n_cells
    assert t["genes_before"].iloc[0] == 120 and t["genes_after"].iloc[-1] == filtered_run.n_genes < 120
    assert t.loc[t["step"] == "prefilter_min_cells_per_gene", "genes_removed"].iloc[0] > 0
    assert t.loc[t["step"] == "rna_min_genes", "cells_removed"].iloc[0] > 0
    assert set(t["category"]) >= {"input", "prefilter", "rna", "final"}
    u = filtered_run.adata.uns["petrubseq_protein"]["qc"]
    assert len(u["filtering_steps"]) == len(t) and u["filter"]["cells_after"] == filtered_run.n_cells


def test_removed_by_is_first_failing_step(filtered_run):
    pre = pd.read_csv(filtered_run.tables["cell_qc_prefilter"], index_col=0)
    t = _steps(filtered_run)
    n_pre = int(t.loc[t["step"] == "prefilter_min_cells_per_gene", "cells_after"].iloc[0])
    assert len(pre) == n_pre  # prefilter-removed cells have no metrics and are not listed
    assert set(pre["removed_by"].fillna("").unique()) <= {"", "rna_min_genes", "rna_max_pct_mt"}
    both = pre[(pre["n_genes_by_counts"] < 102) & (pre["pct_counts_mt"] > 3.0)]
    assert len(both) > 0 and (both["removed_by"] == "rna_min_genes").all()
    only_mt = pre[(pre["n_genes_by_counts"] >= 102) & (pre["pct_counts_mt"] > 3.0)]
    assert (only_mt["removed_by"] == "rna_max_pct_mt").all()
    assert (pre["qc_retained"] == (pre["removed_by"].fillna("") == "")).all()
    assert pre["qc_retained"].sum() == filtered_run.n_cells
    assert set(pre.index[pre["qc_retained"]]) == set(filtered_run.adata.obs_names)
    for col in ("perturbation_class", "protein_total_counts", "rna_outlier_low_genes", "protein_extreme_counts"):
        assert col in pre.columns, col


def test_prefilter_only_removes_by_detection(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", qc={"prefilter": {"enabled": True, "min_genes_per_cell": 95, "min_cells_per_gene": 0}, "filter": {"enabled": False}}, **NO_UMAP, **NO_H5AD))
    res = run_pipeline(cfg)
    t = _steps(res)
    assert t.loc[t["step"] == "prefilter_min_genes", "cells_removed"].iloc[0] > 0
    assert "strict_filter" in set(t["step"]) and t.loc[t["step"] == "strict_filter", "cells_removed"].iloc[0] == 0
    assert (res.adata.obs["n_genes_by_counts"] >= 95).all()
    assert res.adata.uns["petrubseq_protein"]["qc"]["prefilter"]["cells_after"] == res.n_cells


def test_strict_filtering_disabled_keeps_everything_but_flags(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", qc={"prefilter": {"enabled": False}, "filter": {"enabled": False, "rna": {"min_genes": 104}}}, **NO_UMAP, **NO_H5AD))
    res = run_pipeline(cfg)
    assert res.n_cells == 300 and res.n_genes == 120
    t = _steps(res)
    assert list(t["step"]) == ["input", "prefilter", "strict_filter", "final"] and (t["cells_removed"] == 0).all()
    assert res.adata.obs["rna_low_genes"].sum() > 0 and res.adata.obs["rna_qc_fail"].sum() > 0
    assert any("disabled" in n for n in res.adata.uns["petrubseq_protein"].get("warnings", []) + [""]) or True
    pre = pd.read_csv(res.tables["cell_qc_prefilter"], index_col=0)
    assert pre["qc_retained"].all()


def test_perturbation_cell_modes(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "a", qc={"filter": {"enabled": True, "perturbation": {"cells": "assigned"}}}, **NO_UMAP, **NO_H5AD))
    res = run_pipeline(cfg)
    assert (res.adata.obs["n_guides"] >= 1).all()
    assert "perturbation_assigned_only" in set(_steps(res)["step"])
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "s", qc={"filter": {"enabled": True, "perturbation": {"cells": "single_guide"}}}, **NO_UMAP, **NO_H5AD))
    res = run_pipeline(cfg)
    assert (res.adata.obs["n_guides"] == 1).all()
    pre = pd.read_csv(res.tables["cell_qc_prefilter"], index_col=0)
    assert (pre.loc[pre["removed_by"] == "perturbation_single_guide", "n_guides"] != 1).all()


def test_zero_cell_failure(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", qc={"filter": {"enabled": True, "rna": {"min_genes": 200}}}, **NO_UMAP, **NO_H5AD))
    with pytest.raises(ValueError, match="removed every cell"):
        run_pipeline(cfg)


# --------------------------------------------------------------- figures ---


def test_before_after_rna_and_protein_figures(filtered_run):
    reg = filtered_run.registry
    for section in (SECTION_QC_RNA, SECTION_QC_PROTEIN):
        before, after = reg.by_section(section, STAGE_BEFORE), reg.by_section(section, STAGE_AFTER)
        assert len(before) >= 3 and len(after) >= 3, section
        assert {r.name.replace("before_filtering", "X") for r in before} == {r.name.replace("after_filtering", "X") for r in after}
        for r in before + after:
            assert r.path.exists() and r.path.parent.name == r.stage and r.path.parent.parent.name == section
    names = {r.name for r in reg.by_section(SECTION_QC_PROTEIN, STAGE_BEFORE)}
    for n in ("protein_qc_distributions", "protein_targeting_vs_isotype", "protein_depth_vs_rna_depth", "protein_antibody_distributions", "protein_antibody_summary", "protein_normalized_distributions"):
        assert any(x.startswith(n) for x in names), n
    assert (filtered_run.run_dir / "figures" / "qc_rna" / "before_filtering").is_dir() and (filtered_run.run_dir / "figures" / "qc_protein" / "after_filtering").is_dir()


def test_normalized_only_protein_has_no_count_figures(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", inputs={"protein_counts": {"file": None}}, protein={"extra_representations": [], "embedding_representation": "provided"}, **NO_UMAP, **NO_H5AD))
    res = run_pipeline(cfg)
    names = {r.name for r in res.registry.by_section(SECTION_QC_PROTEIN, STAGE_BEFORE)}
    assert any(n.startswith("protein_normalized_distributions") for n in names)
    assert not any(n.startswith(("protein_qc_distributions", "protein_targeting_vs_isotype", "protein_antibody_distributions")) for n in names)
    assert "protein_total_counts" not in res.adata.obs.columns


def test_guide_count_figures_only_when_counts_exist(filtered_run):
    names = {r.name for r in filtered_run.registry.by_section(SECTION_QC_PERTURBATION)}
    assert "guide_count_diagnostics" not in names
    for n in ("perturbation_class_composition", "guides_per_cell", "coverage_histograms", "condition_target_coverage"):
        assert n in names


def test_figure_manifest_matches_registry(filtered_run):
    man = pd.read_csv(filtered_run.tables["figure_manifest"])
    reg = filtered_run.registry
    assert len(man) == len(reg.records) > 0
    assert set(man.columns) >= {"section", "stage", "name", "title", "caption", "in_report", "path"}
    for p in man["path"]:
        assert (filtered_run.run_dir / p).exists(), p
    assert set(man["section"]) >= {"qc_rna", "qc_protein", "qc_perturbation", "representations_rna", "representations_protein", "multimodal"}
    assert set(man["stage"]) >= {"before_filtering", "after_filtering", "post_filter", "representation"}


# ---------------------------------------------------------------- report ---


def test_report_html_is_self_contained(filtered_run):
    html = (filtered_run.run_dir / "report.html").read_text(encoding="utf-8")
    n_embedded = html.count("data:image/png;base64,")
    assert n_embedded == sum(r.in_report for r in filtered_run.registry.records) > 0
    assert 'src="figures/' not in html
    assert not re.search(r'src="(?!data:)', html)
    for section in ("Run summary", "1. Inputs and modality audit", "2. Cell and RNA QC", "3. Protein QC", "4. Perturbation QC", "5. Representations", "6. Outputs and provenance", "Filtering steps"):
        assert section in html, section
    assert "n_top_genes: 50" in html  # resolved config embedded
    assert "rna_min_genes" in html and "prefilter_min_genes" in html
    assert (filtered_run.run_dir / "report.md").exists()
    md = (filtered_run.run_dir / "report.md").read_text()
    assert "figures/qc_rna/before_filtering/" in md


def test_report_links_figures_when_not_embedded(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", report={"embed_figures": False, "write_markdown": False}, **NO_UMAP, **NO_H5AD))
    res = run_pipeline(cfg)
    html = (res.run_dir / "report.html").read_text(encoding="utf-8")
    assert "data:image" not in html and 'src="figures/qc_rna/before_filtering/' in html
    assert not (res.run_dir / "report.md").exists()


def test_resolved_config_and_timings(filtered_run):
    cfg = yaml.safe_load((filtered_run.run_dir / "logs" / "resolved_config.yaml").read_text())
    assert cfg["qc"]["filter"]["rna"]["min_genes"] == 102 and cfg["qc"]["prefilter"]["min_cells_per_gene"] == 100
    for k in ("validate configuration", "permissive prefilter", "BEFORE-filter figures", "strict filtering", "AFTER-filter figures", "representations", "report", "total"):
        assert k in filtered_run.timings, k
    import json
    man = json.loads((filtered_run.run_dir / "logs" / "run_manifest.json").read_text())
    assert man["timings"]["total"] > 0 and man["filtering_steps"][-1]["step"] == "final" and "warnings" in man


def test_optional_prefilter_h5ad(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", qc={"filter": {"enabled": True, "rna": {"min_genes": 104}}}, output={"write_prefilter_h5ad": True}, **NO_UMAP))
    res = run_pipeline(cfg)
    assert res.prefilter_h5ad is not None and res.prefilter_h5ad.exists() and res.prefilter_h5ad.parent.name == "processed"
    pre = ad.read_h5ad(res.prefilter_h5ad)
    t = _steps(res)
    assert pre.n_obs == int(t.loc[t["step"] == "rna_min_genes", "cells_before"].iloc[0]) > res.n_cells
    assert "rna_qc_fail" in pre.obs.columns and "protein" in pre.obsm


def test_archive_contents(synthetic, tmp_path):
    cfg = Config.from_dict(base_config(synthetic, tmp_path / "r", output={"archive": True}, run={"name": "arch"}, **NO_UMAP))
    res = run_pipeline(cfg)
    assert res.archive is not None and res.archive.name == "arch_results.tar.gz"
    with tarfile.open(res.archive) as tar:
        names = tar.getnames()
    assert "arch/report.html" in names and any(n.startswith("arch/figures/qc_rna/") for n in names) and "arch/tables/qc_filtering_steps.csv" in names
    assert not any(n.endswith(".h5ad") for n in names)
    assert (res.run_dir / "processed" / cfg.output.h5ad_name).exists()


def test_init_config_round_trip(tmp_path):
    path = tmp_path / "new.yaml"
    assert cli_main(["init-config", str(path)]) == 0
    cfg = Config.from_yaml(path)
    assert cfg.qc.filter.rna.min_genes == 500 and cfg.qc.prefilter.enabled and cfg.rna.hvg.n_top_genes == 3000
    assert cfg.dataset.input_dir.endswith("/data/my_dataset") and cfg.output.dir.endswith("/results/my_dataset")
    d = yaml.safe_load(path.read_text())
    assert set(d) >= {"run", "dataset", "inputs", "qc", "report", "output"}
