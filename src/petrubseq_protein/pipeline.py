"""The processing/QC pipeline as explicit, logged, timed stages.

``run_pipeline(cfg)`` is the single entry point for the CLI, the tests and
(later) the demo notebook.

Lifecycle
---------
 1 validate configuration        13 filtering audit table
 2 create output tree            14 post-filter QC summaries + perturbation QC
 3 logging / resolved config     15 AFTER-filter figures
 4 input audit                   16 RNA / protein representations
 5 load data                     17 cross-modality diagnostics + figures
 6 align cells                   18 processed .h5ad
 7 harmonize perturbations       19 tables
 8 RNA input state / count layer 20 figure manifest
 9 permissive prefilter          21 report.html (+ report.md)
10 pre-filter QC metrics         22 run manifest / provenance
11 BEFORE-filter figures         23 optional archive
12 strict filtering
"""

from __future__ import annotations

import datetime as _dt
import logging
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from . import __version__
from .config import Config
from . import validation
from .io import adapters
from .io import readers as pio
from .preprocessing import align, embeddings, normalize
from .preprocessing import guides as guidecall
from .qc import filtering
from .qc import metrics as qcm
from .qc import perturbation as pqc
from .reporting import plots, provenance
from .reporting.archive import archive_results
from .reporting.context import ReportContext
from .reporting.figures import FigureRegistry
from .reporting.html import write_html
from .reporting.markdown import write_markdown

logger = logging.getLogger("petrubseq_protein")
N_STAGES = 23


@dataclass
class PipelineResult:
    run_dir: Path
    h5ad: Optional[Path]
    report_html: Path
    report_md: Optional[Path] = None
    prefilter_h5ad: Optional[Path] = None
    archive: Optional[Path] = None
    tables: Dict[str, Path] = field(default_factory=dict)
    registry: Optional[FigureRegistry] = None
    adata: Optional[ad.AnnData] = None
    filtering_steps: Optional[pd.DataFrame] = None
    warnings: List[str] = field(default_factory=list)
    timings: Dict[str, float] = field(default_factory=dict)
    n_cells: int = 0
    n_genes: int = 0
    n_proteins: int = 0
    runtime_seconds: float = 0.0

    @property
    def report(self) -> Path:  # backwards-compatible alias
        return self.report_html

    def summary(self) -> str:
        return f"{self.n_cells:,} cells x {self.n_genes:,} genes x {self.n_proteins} proteins in {self.runtime_seconds:.0f}s | report: {self.report_html}"


def setup_logging(run_dir: Path, level: int = logging.INFO) -> Path:
    log_path = run_dir / "logs" / "run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger.setLevel(level)
    logger.propagate = False
    logger.handlers = [h for h in logger.handlers if not isinstance(h, logging.FileHandler)]
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler) for h in logger.handlers):
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        logger.addHandler(sh)
    return log_path


def validate_inputs(cfg: Config) -> Dict[str, Any]:
    """Check that every configured input exists (format-aware); see :mod:`validation`."""
    return validation.validate_input_files(cfg)


class _Stages:
    """Stage banners and timings."""

    def __init__(self, n_stages: int = N_STAGES) -> None:
        self.timings: Dict[str, float] = {}
        self.n = 0
        self.n_stages = n_stages

    @contextmanager
    def stage(self, name: str):
        self.n += 1
        logger.info("=== Stage %d/%d: %s ===", self.n, self.n_stages, name)
        t = time.time()
        yield
        self.timings[name] = time.time() - t
        logger.info("--- %s done in %.1f s", name, self.timings[name])


def run_pipeline(cfg: Config) -> PipelineResult:
    t0 = time.time()
    pe_enabled = cfg.analysis.perturbation_effects.enabled
    cs_enabled = cfg.analysis.clustering.enabled
    st = _Stages(N_STAGES + (2 if pe_enabled else 0) + (1 if cs_enabled else 0) + (1 if (cs_enabled and cfg.analysis.clustering.enrichment.enabled) else 0))
    warnings: List[str] = []
    notes: List[str] = []
    tables: Dict[str, pd.DataFrame] = {}

    # 1 ------------------------------------------------------------------
    with st.stage("validate configuration"):
        cfg.validate()
        run_dir = cfg.output_dir()
        run_name = cfg.run_name()
    # 2 ------------------------------------------------------------------
    with st.stage("create output tree"):
        for sub in ("processed", "tables", "figures", "logs"):
            (run_dir / sub).mkdir(parents=True, exist_ok=True)
        registry = FigureRegistry(run_dir, cfg.report.figure_format, cfg.report.figure_dpi)
    # 3 ------------------------------------------------------------------
    with st.stage("logging / resolved config"):
        log_path = setup_logging(run_dir)
        logger.info("petrubseq-pipeline-protein %s: run %r (dataset %s) -> %s", __version__, run_name, cfg.dataset.name, run_dir)
        cfg.to_yaml(run_dir / "logs" / "resolved_config.yaml")
    # 4 ------------------------------------------------------------------
    with st.stage("input audit"):
        input_rec = validate_inputs(cfg)
        input_files = pd.DataFrame([{"role": k, **(v or {})} for k, v in input_rec.items()])
    # 5 ------------------------------------------------------------------
    with st.stage("load data"):
        max_cells = cfg.compute.max_cells
        ci = adapters.load_inputs(cfg, max_cells=max_cells)
        errs, warns = validation.validate_canonical(ci, cfg)
        if errs:
            raise adapters.InputError("Input validation failed:\n  - " + "\n  - ".join(errs))
        warnings.extend(warns)
        meta, guides, emb = ci.metadata, ci.guide_assignments, ci.embedding
        rna, prot, counts, gcounts = ci.rna, ci.protein_normalized, ci.protein_counts, ci.guide_counts
    # 6 ------------------------------------------------------------------
    with st.stage("align cells"):
        id_sets: Dict[str, Any] = ci.id_sets()
        if max_cells is not None:
            keep = set(id_sets["rna"])
            id_sets = {k: [c for c in v if c in keep] for k, v in id_sets.items()}
        ares = align.align_cells(id_sets, cfg, reference="rna", required=align.effective_required(cfg.alignment.required, id_sets))
        cells = ares.cells
        if emb is not None:
            ares.optional_missing["embedding"] = int(len(set(cells) - set(emb.index)))
        for k, n in ares.optional_missing.items():
            if n:
                warnings.append(f"{n} cells lack the optional modality '{k}' (kept, values NaN / flagged).")
        for k, n in ares.only_in.items():
            if n:
                notes.append(f"{n} cell IDs occur only in '{k}' and were not retained.")
    # 7 ------------------------------------------------------------------
    with st.stage("harmonize perturbations"):
        gt_table = None
        if cfg.perturbation.guide_target_table:
            gt = pio.read_table(cfg.resolve(cfg.perturbation.guide_target_table))
            gt_table = dict(zip(gt.iloc[:, 0].astype(str), gt.iloc[:, 1].astype(str)))
        # -- which assignment feeds the perturbation columns ------------------
        asg = cfg.perturbation.assignment
        provided_available = guides is not None or bool(meta is not None and ((cfg.columns.guide and cfg.columns.guide in meta.columns) or (cfg.columns.perturbation and cfg.columns.perturbation in meta.columns)))
        counts_available = gcounts is not None
        if asg.source == "auto":
            if provided_available and counts_available:
                raise adapters.InputError("both provided guide assignments and a guide-count matrix are available; set perturbation.assignment.source to 'provided' or 'guide_counts'")
            effective = "guide_counts" if counts_available else "provided"
        else:
            effective = asg.source
        if effective == "guide_counts" and not counts_available:
            raise adapters.InputError("perturbation.assignment.source = guide_counts but no guide-count matrix was loaded")
        calls = None
        if counts_available:
            calls = guidecall.call_guides(gcounts, cells, cfg, gt_table)
        lists_for_harmonize = guides
        if effective == "guide_counts":
            lists_for_harmonize = calls.lists
        pert = align.harmonize_perturbations(cells, cfg, lists_for_harmonize, meta, gt_table)
        if effective == "guide_counts":
            pert.attrs["source"] = f"guide_counts:{asg.method}"
            if asg.method == "dominant" and calls.ambiguous.any():
                amb = calls.ambiguous.to_numpy()
                pert.loc[amb, "perturbation"] = asg.ambiguous_label
                pert.loc[amb, "perturbation_class"] = "ambiguous"
        pert_info = {"source": pert.attrs.get("source"), "control_classes": pert.attrs.get("control_classes"), "metadata_guide_mismatches": pert.attrs.get("metadata_guide_mismatches"), "metadata_guide_compared": pert.attrs.get("metadata_guide_compared"), "cells_before_preserve_filters": len(cells), "guide_counts_available": counts_available,
                     "assignment": {"provided_assignments_available": provided_available, "guide_counts_available": counts_available, "configured_source": asg.source, "effective_source": effective, "guide_calling": calls.info if calls is not None else None}}
        if calls is not None and provided_available:
            prov_lists = guides if guides is not None else align.harmonize_perturbations(cells, cfg, None, meta, gt_table)["guides"].map(lambda s: s.split(";") if s else [])
            agreement = guidecall.compare_assignments(prov_lists, calls, cells)
            pert_info["assignment"]["agreement"] = agreement
            n_dis = agreement["set_agreement"]["disagree"]
            if n_dis:
                warnings.append(f"Provided guide assignments and count-derived detected guides disagree for {n_dis} of {agreement['cells_compared']} cells with guide counts (effective source: {effective}); see uns['petrubseq_protein']['perturbations']['assignment']['agreement'].")
        if calls is not None and calls.info["cells_without_guide_counts"]:
            warnings.append(f"{calls.info['cells_without_guide_counts']} cells have no guide counts (has_guide_counts=False; zero rows in obsm['guide_counts']).")
        if effective == "guide_counts" and asg.method == "dominant":
            notes.append(f"Guide calling: dominant rule (min_umi {asg.min_umi}, dominance_ratio {asg.dominance_ratio}); {calls.info['dominant_rule']['ambiguous']} cells ambiguous, {calls.info['dominant_rule']['unassigned']} unassigned; detected guide lists kept in obs['guides_detected'].")
        keep = np.ones(len(cells), dtype=bool)
        if not cfg.perturbation.preserve_multiguide:
            keep &= (pert["n_guides"] <= 1).to_numpy()
        if not cfg.perturbation.preserve_unassigned:
            keep &= (pert["n_guides"] >= 1).to_numpy()
        if not keep.all():
            warnings.append(f"perturbation.preserve_* dropped {int((~keep).sum())} multi-guide/unassigned cells before QC.")
            cells = [c for c, k in zip(cells, keep) if k]
            pert = pert.loc[cells]
        pert_info["cells_dropped_by_preserve_filters"] = int((~keep).sum())
        if pert.attrs.get("metadata_guide_mismatches"):
            warnings.append(f"{pert.attrs['metadata_guide_mismatches']} single-guide cells disagree between the guide assignments and the metadata guide column.")
        obs = align.build_obs(cells, cfg, meta, pert)
        if calls is not None:
            diag = calls.diagnostics.loc[cells]
            for c in diag.columns:
                obs[c] = diag[c].to_numpy()
        if ci.lane_table is not None:
            lt = ci.lane_table.reindex(obs.index)
            obs["lane_id"] = lt["lane_id"].fillna("").astype(str).to_numpy()
            obs["barcode_original"] = lt["barcode_original"].fillna("").astype(str).to_numpy()
    # 8 ------------------------------------------------------------------
    with st.stage("RNA input state / count layer"):
        if len(cells) == len(rna.cells) and np.array_equal(np.asarray(cells, dtype=object), rna.cells):
            rna_sub = rna
        else:
            rna_sub = rna.subset_cells(cells)
        del rna
        tot = obs["rna_total_counts_provided"].to_numpy(dtype=float) if "rna_total_counts_provided" in obs else None
        X, counts_mat, counts_layer, rna_info = normalize.normalize_rna(rna_sub.X, cfg.rna, tot, cfg.inputs.rna.state, apply_normalization=False)
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index([str(f) for f in rna_sub.features], name="gene")))
        adata.obs_names = [str(c) for c in cells]
        if counts_mat is not None and counts_layer:
            adata.layers[counts_layer] = counts_mat
        del rna_sub
        if calls is not None:
            gpos = pd.Index(calls.diagnostics.index).get_indexer(cells)
            adata.obsm["guide_counts"] = sp.csr_matrix(calls.X[gpos])
            adata.uns["guide_features"] = calls.features.copy()
        if rna_info.get("input_state_used") == "log_normalized" and counts_layer is None:
            warnings.append("RNA input is log-normalized and no per-cell library size is available: count-based QC uses expm1(X) fractions and no count layer is stored.")
        if rna_info.get("input_state_used") == "uncertain":
            warnings.append("RNA value state could not be determined; values were preserved unchanged. Set rna.input_state explicitly.")
        # protein representations are per-cell and independent of filtering
        prot_info: Dict[str, Any] = {}
        prot_emb_key = None
        if prot is not None or counts is not None:
            ftab = pio.read_table(cfg.resolve(cfg.protein.feature_table)) if cfg.protein.feature_table else None
            prot_df, counts_df, extras, feats_df, prot_info = normalize.normalize_protein(prot, counts, cells, cfg.protein, feature_table=ftab)
            if prot_df is not None:
                adata.obsm["protein"] = prot_df
            if counts_df is not None:
                adata.obsm["protein_counts"] = counts_df
            for name, df in extras.items():
                adata.obsm[f"protein_{name}"] = df
            adata.uns["protein_features"] = feats_df
            prot_emb_key = prot_info.get("embedding_key")
            chk = prot_info.get("provided_vs_counts_check", {})
            if chk.get("checked") and not chk.get("formula_reproduced"):
                warnings.append("The provided normalized protein matrix could not be reproduced from the raw counts with the isotype-ratio formula; both are kept, provenance of the provided values is unknown.")
        if emb is not None:
            adata.obsm[cfg.inputs.embedding.key] = emb.reindex(adata.obs_names).to_numpy(dtype=np.float32)
        n_input = adata.n_obs
    # 9 ------------------------------------------------------------------
    with st.stage("permissive prefilter"):
        audit = filtering.FilterAudit()
        filtering.start(audit, adata)
        adata = filtering.prefilter(adata, cfg, audit)
        n_prefilter = adata.n_obs
    # 10 -----------------------------------------------------------------
    with st.stage("pre-filter QC metrics"):
        rna_qc_info = qcm.rna_qc(adata, cfg)
        prot_qc_info: Dict[str, Any] = {}
        prot_table = None
        if "protein" in adata.obsm or "protein_counts" in adata.obsm:
            prot_qc_info, prot_table = qcm.protein_qc(adata, cfg)
            for n in prot_qc_info.get("background_dominated", []):
                pass
            if prot_qc_info.get("background_dominated"):
                warnings.append(f"{len(prot_qc_info['background_dominated'])} antibodies are background-dominated (signal above isotype in < {100*cfg.qc.flags.background_min_fraction_above_isotype:.0f}% of cells): {', '.join(prot_qc_info['background_dominated'])}.")
            ext = prot_qc_info.get("extreme_counts_threshold", {}).get("n_cells", 0)
            if ext:
                warnings.append(f"{ext} cells have extreme total ADT counts (> {cfg.qc.flags.extreme_fold} x the 99th percentile); flagged as protein_extreme_counts{' and removed' if cfg.qc.filter.enabled and cfg.qc.filter.protein.remove_extreme_counts else ', not removed'}.")
        if rna_qc_info.get("n_mito_genes", 0) == 0:
            warnings.append(f"No mitochondrial genes matched prefix {cfg.rna.mito_prefix!r}; pct_counts_mt is zero everywhere.")
    # 11 -----------------------------------------------------------------
    with st.stage("BEFORE-filter figures"):
        _qc_figures(adata, prot_table, cfg, registry, "before_filtering", prot_info)
    # 12 -----------------------------------------------------------------
    with st.stage("strict filtering"):
        prefilter_obs = adata.obs.copy()
        prefilter_h5ad = None
        if cfg.output.write_prefilter_h5ad:
            prefilter_h5ad = run_dir / "processed" / (Path(cfg.output.h5ad_name).stem + "_prefilter.h5ad")
            tmp = adata.copy()
            _sanitize_for_h5ad(tmp)
            tmp.write_h5ad(prefilter_h5ad, compression=cfg.output.compression)
            del tmp
            logger.info("wrote pre-filter object %s", prefilter_h5ad)
        adata, removed_by = filtering.strict_filter(adata, cfg, audit)
        n_final = adata.n_obs
        if not cfg.qc.filter.enabled:
            notes.append("Strict filtering is disabled for this run (qc.filter.enabled: false); QC flags are recorded but no cell was removed, so before- and after-filter figures describe the same cells.")
    # 13 -----------------------------------------------------------------
    with st.stage("filtering audit table"):
        steps = filtering.finalize(audit, adata)
        tables["qc_filtering_steps"] = steps
        prefilter_obs["qc_retained"] = (removed_by.reindex(prefilter_obs.index) == "").to_numpy()
        prefilter_obs["removed_by"] = removed_by.reindex(prefilter_obs.index).to_numpy()
        tables["cell_qc_prefilter"] = prefilter_obs
        logger.info("filtering audit: input %d -> prefilter %d -> strict %d cells", n_input, n_prefilter, n_final)
    # 14 -----------------------------------------------------------------
    with st.stage("post-filter QC summaries"):
        if n_final < n_prefilter:
            rna_qc_info = qcm.rna_qc(adata, cfg)  # refresh per-gene metrics and flag counts on retained cells
            if prot_table is not None:
                prot_qc_info, prot_table = qcm.protein_qc(adata, cfg)
        obs = adata.obs
        group_col = plots.choose_group(obs)
        tables["qc_summary"] = qcm.qc_summary_table(obs, group_col)
        pc = pqc.perturbation_counts(obs)
        gc_ = pqc.guide_counts(obs, cfg)
        tc = pqc.target_counts(obs, cfg)
        cond_guide = pqc.condition_coverage(obs, "guide", cfg)
        cond_target = pqc.condition_coverage(obs, "target", cfg)
        moi = pqc.moi_summary(obs)
        pert_qc = pqc.perturbation_qc(obs, cfg, gc_, tc, cond_guide, cond_target)
        asg_rec = pert_info.get("assignment", {})
        pert_qc["assignment_source"] = f"{asg_rec.get('effective_source')} (configured: {asg_rec.get('configured_source')})"
        if asg_rec.get("guide_calling"):
            gcall = asg_rec["guide_calling"]
            pert_qc["guide_calling"] = f"{gcall['method']} rule; min_umi {gcall['min_umi']}, dominance_ratio {gcall['dominance_ratio']}, detection_min_umi {gcall['detection_min_umi']}"
            pert_qc["cells_with_guide_counts"] = gcall["cells_with_guide_counts"]
            pert_qc["dominant_rule_assigned_ambiguous_unassigned"] = "{assigned} / {ambiguous} / {unassigned}".format(**gcall["dominant_rule"])
            pert_qc["n_cells_ambiguous"] = int((obs["perturbation_class"] == "ambiguous").sum())
        if asg_rec.get("agreement"):
            ag = asg_rec["agreement"]
            pert_qc["provided_vs_counts_set_agreement"] = f"{ag['set_agreement']['agree']} / {ag['cells_compared']} cells ({100 * ag['set_agreement']['fraction_agree']:.1f}%)"
            pert_qc["provided_single_vs_dominant_call"] = f"{ag['dominant_vs_provided_single']['agree']} agree / {ag['dominant_vs_provided_single']['disagree']} disagree"
        class_by_cond = pd.crosstab(obs["perturbation_class"], obs["condition"]) if "condition" in obs else obs["perturbation_class"].value_counts().to_frame()
        group_sums = {k: s for k in ("condition", "sample", "donor", "batch", "replicate", "lane") if (s := pqc.group_summary(obs, k)) is not None}
        tables.update({"perturbation_counts": pc, "guide_coverage": gc_, "target_coverage": tc, "moi_summary": moi, "perturbation_qc": pd.Series({k: v for k, v in pert_qc.items() if not isinstance(v, dict)}, name="value").to_frame()})
        if cond_guide is not None:
            tables["condition_guide_coverage"] = cond_guide
        if cond_target is not None:
            tables["condition_target_coverage"] = cond_target
        for k, df in group_sums.items():
            tables[f"{k}_summary"] = df
        if pert_qc["n_low_coverage_guides"]:
            warnings.append(f"{pert_qc['n_low_coverage_guides']} targeting guides have fewer than {cfg.perturbation.min_cells_per_guide} single-guide cells (flagged low_coverage, not removed).")
        if pert_qc["n_low_coverage_targets"]:
            warnings.append(f"{pert_qc['n_low_coverage_targets']} targets have fewer than {cfg.perturbation.min_cells_per_target} single-guide cells (flagged low_coverage, not removed).")
        if pert_qc.get("n_low_coverage_target_condition_pairs"):
            warnings.append(f"{pert_qc['n_low_coverage_target_condition_pairs']} of {pert_qc.get('n_target_condition_pairs')} target x condition pairs are below {cfg.perturbation.min_cells_per_target} cells.")
        if pert_qc["frac_multi_guide"] > 0.3:
            warnings.append(f"{100*pert_qc['frac_multi_guide']:.1f}% of cells carry more than one guide; single-guide analyses will use a subset.")
        if pert_qc["frac_unassigned"] > 0.3:
            warnings.append(f"{100*pert_qc['frac_unassigned']:.1f}% of cells have no guide assignment.")
        if not pert_qc["control_classes"]:
            warnings.append("No control guides were recognised; check perturbation.control_classes against the guide naming.")
    # 15 -----------------------------------------------------------------
    with st.stage("AFTER-filter figures"):
        _qc_figures(adata, prot_table, cfg, registry, "after_filtering", prot_info)
        plots.perturbation_qc_figures(obs, pc, gc_, tc, cond_target, cfg, registry)
    # 16 -----------------------------------------------------------------
    with st.stage("representations"):
        normalize.log_normalize_rna(adata, cfg.rna, rna_info)
        emb_info = {"rna": embeddings.rna_embedding(adata, cfg), "protein": embeddings.protein_embedding(adata, cfg, prot_emb_key)}
        emb_info["multimodal"] = embeddings.multimodal_embedding(adata, cfg)
    # 17 -----------------------------------------------------------------
    with st.stage("cross-modality diagnostics"):
        emb_info["diagnostics"] = embeddings.modality_diagnostics(adata, cfg)
        color_by = [c for c in cfg.umap.color_by if c in obs.columns] + [c for c in ("n_guides", "pct_counts_mt", "total_counts", "protein_total_counts") if c in obs.columns]
        seen: set = set()
        color_by = [c for c in color_by if not (c in seen or seen.add(c))]
        plots.representation_figures(adata, cfg, registry, color_by, cfg.inputs.embedding.key if emb is not None else None, emb_info["diagnostics"])
        d = emb_info["diagnostics"]
        r = d.get("protein_pc_vs_log_total_adt", {}).get("PC1")
        if r is not None and r > 0.7:
            warnings.append(f"Protein PC1 correlates with ADT depth (|r| = {r:.2f}); the protein embedding is depth-driven (documented CITE-seq effect, not corrected).")
    # 17b-17e: the reference analysis order (docs/reference/REFERENCE_PIPELINE_COMPLETE_AUDIT.md):
    #   Leiden clustering -> perturbation strength -> perturbation x cluster enrichment
    #   -> modules / programs -> PS -> lochNESS -> protein effects -> concordance
    pe_res = None
    cs_res = None
    if cs_enabled:
        with st.stage("cell states (Leiden clustering)"):
            from . import analysis as cs_analysis
            from .reporting.cell_state_plots import clustering_figures
            cs_res = cs_analysis.run_clustering(adata, cfg)
            clustering_figures(cs_res, adata, cfg, registry)
            for f in cs_res.clustering.flags:
                warnings.append(f"Leiden {f} (see tables/cell_states/cluster_summary.csv).")
    if pe_enabled:
        with st.stage("perturbation strength"):
            from . import analysis as pe_analysis
            from .reporting.strength_plots import perturbation_strength_figures
            strength = pe_analysis.run_perturbation_strength(adata, cfg)
            if strength is not None:
                perturbation_strength_figures(strength, adata, cfg, registry)
                if strength.empty:
                    warnings.append(f"Perturbation strength not computed: {strength.info.get('status')}.")
                elif len(strength.skipped):
                    notes.append(f"Perturbation strength: {len(strength.skipped)} target(s) not testable (tables/perturbation_strength/skipped.csv).")
                if strength.primary_control != cfg.analysis.perturbation_effects.strength.primary_control:
                    warnings.append(f"Perturbation strength: requested primary control {cfg.analysis.perturbation_effects.strength.primary_control!r} unavailable; {strength.primary_control!r} drives ranking and hit calls.")
    else:
        strength = None
    if cs_enabled and cfg.analysis.clustering.enrichment.enabled:
        with st.stage("perturbation x cluster enrichment"):
            from .reporting.enrichment_plots import enrichment_figures
            cs_res = cs_analysis.run_enrichment(adata, cfg, cs_res)
            enrichment_figures(cs_res.enrichment, adata, cfg, registry)
            if cs_res.enrichment is not None and cs_res.enrichment.empty:
                notes.append(f"Perturbation x cluster enrichment not computed: {cs_res.enrichment.info.get('status')}.")
    if pe_enabled:
        with st.stage("perturbation effects (modules, PS, lochNESS, distance, protein, concordance, master tables)"):
            from .reporting.perturbation_plots import perturbation_effect_figures
            pe_res = pe_analysis.run_perturbation_effects(adata, cfg, strength, cs_res)
            integrated = pe_analysis.integrate_target_summary(pe_res, cs_res)
            if integrated is not None:
                pe_res.concordance.summary = integrated
            pe_analysis.attach(adata, pe_res)
            perturbation_effect_figures(pe_res, adata, cfg, registry)
            for name, r in (("Perturbation distance", pe_res.distance), ("Perturbation distance space", pe_res.distance_space)):
                if r is not None and r.empty:
                    warnings.append(f"{name} not computed: {r.note}.")
                elif r is not None and len(r.skipped):
                    notes.append(f"{name}: {len(r.skipped)} target(s) below min_cells = {r.info.get('min_cells')} (tables/perturbation_distance/*_skipped.csv).")
            if pe_res.master is not None and pe_res.master.info.get("power_note"):
                notes.append(f"Distance-protein associations: {pe_res.master.info['power_note']}.")
            for name, r in (("PS", pe_res.ps), ("lochNESS", pe_res.lochness), ("protein effects", pe_res.protein)):
                if r is not None and getattr(r, "skipped", None) is not None and len(r.skipped):
                    notes.append(f"{name}: {len(r.skipped)} target(s) not analysed (see tables/perturbation_effects/*_skipped.csv or the {name} table).")
            if pe_res.modules is not None and pe_res.modules.empty:
                warnings.append(f"Gene programs / perturbation modules not computed: {pe_res.modules.info.get('status')}.")
            if pe_res.ps is not None and pe_res.ps.lda_note:
                notes.append(f"PS: the supervised LDA embedding was not built: {pe_res.ps.lda_note}.")
            if pe_res.lochness is not None and pe_res.lochness.info.get("k_capped"):
                notes.append(f"lochNESS k reduced to {pe_res.lochness.info['k_used']} (max_k_fraction {cfg.analysis.perturbation_effects.lochness.max_k_fraction:g} x {adata.n_obs} cells) instead of {cfg.analysis.perturbation_effects.lochness.n_neighbors}.")
            if cs_res is not None and cs_res.enrichment is not None and not cs_res.enrichment.empty and pe_res.lochness is not None:
                from .analysis.cluster_enrichment import compare_with_lochness
                cs_res.enrichment.lochness_comparison = compare_with_lochness(cs_res.enrichment, pe_res.lochness.summary)
    # 18 -----------------------------------------------------------------
    with st.stage("processed h5ad"):
        prov = provenance.collect(cfg.source_path, input_rec, extra={"alignment": ares.to_dict(), "run_name": run_name})
        adata.uns["petrubseq_protein"] = {
            "version": __version__,
            "schema_version": "0.2",
            "run_name": run_name,
            "dataset": cfg.dataset.name,
            "inputs": _jsonable(_input_audit(cfg, rna_info, prot_info, pert_info, input_rec, ci.provenance)),
            "rna": _jsonable(rna_info),
            "protein": _jsonable({k: v for k, v in prot_info.items() if k != "provided_vs_counts_check"}),
            "protein_check": _jsonable(prot_info.get("provided_vs_counts_check", {})),
            "qc": _jsonable({"rna": rna_qc_info, "protein": prot_qc_info, "prefilter": {"enabled": cfg.qc.prefilter.enabled, "cells_before": n_input, "cells_after": n_prefilter}, "filter": {"enabled": cfg.qc.filter.enabled, "cells_before": n_prefilter, "cells_after": n_final, "cells_removed": n_prefilter - n_final}, "filtering_steps": audit.table()}),
            "embeddings": _jsonable(emb_info),
            "alignment": _jsonable(ares.to_dict()),
            "perturbations": _jsonable({**pert_info, "qc": pert_qc}),
            "warnings": list(warnings),
            "provenance": _jsonable({k: v for k, v in prov.items() if k != "config"}),
            "config": _jsonable(cfg.to_dict()),
        }
        if pe_res is not None or cs_res is not None:
            adata.uns["petrubseq_protein"]["analysis"] = {}
        if pe_res is not None:
            adata.uns["petrubseq_protein"]["analysis"]["perturbation_effects"] = _jsonable(pe_res.info())
        if cs_res is not None:
            adata.uns["petrubseq_protein"]["analysis"]["cell_states"] = _jsonable(cs_res.info())
            if cs_res.enrichment is not None and not cs_res.enrichment.empty:
                from .analysis import _h5
                adata.uns["perturbation_cluster_enrichment"] = _h5(cs_res.enrichment.table)
        adata.uns["petrubseq_protein"]["schema"] = _jsonable(object_schema(adata, rna_info, prot_info, cfg))
        h5ad_path = None
        if cfg.output.write_h5ad:
            h5ad_path = run_dir / "processed" / cfg.output.h5ad_name
            _sanitize_for_h5ad(adata)
            adata.write_h5ad(h5ad_path, compression=cfg.output.compression)
            logger.info("wrote %s (%.1f MB)", h5ad_path, h5ad_path.stat().st_size / 1e6)
    # 19 -----------------------------------------------------------------
    with st.stage("tables"):
        tdir = run_dir / "tables"
        qc_cols = [c for c in obs.columns if c.startswith(("total_counts", "n_genes", "pct_counts", "rna_", "protein_", "has_protein"))]
        base_cols = [c for c in ("condition", "sample", "donor", "batch", "perturbation", "perturbation_class", "control_class", "guide", "target", "n_guides", "n_targets") if c in obs.columns]
        tables["cell_qc"] = obs[base_cols + qc_cols]
        if prot_table is not None:
            tables["protein_qc"] = prot_table
        if "protein_features" in adata.uns:
            tables["protein_features"] = adata.uns["protein_features"]
        if "guide_features" in adata.uns:
            tables["guide_features"] = adata.uns["guide_features"]
        tables["rna_gene_qc"] = adata.var
        if "pca_protein" in adata.uns:
            tables["protein_pca_loadings"] = adata.uns["pca_protein"]["loadings"]
        if pe_res is not None:
            tables.update(pe_analysis.tables(pe_res))
        if cs_res is not None:
            tables.update(cs_analysis.cell_state_tables(cs_res))
        summary = _run_summary(cfg, adata, rna_info, prot_info, pert_qc, n_input, n_prefilter, n_final)
        tables["run_summary"] = pd.Series(summary, name="value").to_frame()
        table_paths: Dict[str, Path] = {}
        for name, df in tables.items():
            if df is None:
                continue
            long_form = name in ("condition_guide_coverage", "condition_target_coverage", "qc_filtering_steps", "qc_summary", "master_perturbation_table", "master_perturbation_protein_table",
                                 "ps_protein_associations", "lochness_protein_associations", "distance_protein_associations", "phenotype_module_protein_associations", "rna_protein_geometry_concordance")
            gz = name in ("cell_qc", "cell_qc_prefilter")
            p = tdir / f"{name}.csv{'.gz' if gz else ''}"
            p.parent.mkdir(parents=True, exist_ok=True)
            long_form = long_form or (name.startswith(("perturbation_effects/", "cell_states/", "perturbation_strength/", "perturbation_distance/")) and isinstance(df.index, pd.RangeIndex))
            df.to_csv(p, index=not long_form)
            table_paths[name] = p
    # 20 -----------------------------------------------------------------
    with st.stage("figure manifest"):
        man = registry.manifest()
        p = tdir / "figure_manifest.csv"
        man.to_csv(p, index=False)
        table_paths["figure_manifest"] = p
    # 21 -----------------------------------------------------------------
    with st.stage("report"):
        report_html = run_dir / "report.html"
        report_md = run_dir / "report.md" if cfg.report.write_markdown else None
        outputs = {"report.html": str(report_html), **({"report.md": str(report_md)} if report_md else {}), "processed h5ad": str(h5ad_path) if h5ad_path else "(not written)", **({"pre-filter h5ad": str(prefilter_h5ad)} if prefilter_h5ad else {}), "figures": str(run_dir / "figures"), "tables": str(tdir), "log": str(log_path), "resolved config": str(run_dir / "logs" / "resolved_config.yaml"), "run manifest": str(run_dir / "logs" / "run_manifest.json")}
        if cfg.output.archive:
            outputs["archive"] = str(run_dir / (cfg.output.archive_name or f"{run_name}_results.tar.gz"))
        ctx = ReportContext(
            run_dir=run_dir, run_name=run_name, dataset=cfg.dataset.name,
            title=cfg.report.title or f"{run_name} - Perturb-CITE-seq preprocessing report", version=__version__,
            generated_at=_dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            cards=_cards(adata, summary, n_input, n_prefilter, n_final, pert_qc, prot_info),
            warnings=warnings, notes=notes,
            inputs=adata.uns["petrubseq_protein"]["inputs"]["modalities"], input_files=input_files,
            alignment=ares.to_dict(), filtering_steps=steps, qc_summary=tables["qc_summary"], rna_qc=rna_qc_info,
            protein_features=prot_table if prot_table is not None else adata.uns.get("protein_features"),
            protein_qc={k: v for k, v in prot_qc_info.items()}, protein_info={k: v for k, v in prot_info.items() if k not in ("provided_vs_counts_check", "counts_state", "provided_state")} | {"provided_vs_counts_check": {k: v for k, v in prot_info.get("provided_vs_counts_check", {}).items() if k != "max_abs_diff"}},
            perturbation_qc={**{k: v for k, v in pert_info.items()}, **{k: v for k, v in pert_qc.items() if k != "class_counts"}},
            class_by_condition=class_by_cond, target_coverage=tc, guide_coverage=gc_, moi_summary=moi, group_summaries=group_sums,
            representations={k: emb_info.get(k, {}) for k in ("rna", "protein", "multimodal")}, diagnostics=emb_info["diagnostics"],
            registry=registry, tables=table_paths, outputs=outputs, timings=st.timings, config_yaml=yaml.safe_dump(cfg.to_dict(), sort_keys=False, allow_unicode=True),
            versions=prov["packages"], provenance={k: v for k, v in prov.items() if k not in ("packages", "inputs")}, schema=adata.uns["petrubseq_protein"]["schema"],
            adata_repr=repr(adata), filter_enabled=cfg.qc.filter.enabled,
        )
        if pe_res is not None:
            ctx.extra["perturbation_effects"] = pe_res
        if cs_res is not None:
            ctx.extra["cell_states"] = cs_res
        write_html(ctx, report_html, embed=cfg.report.embed_figures, max_rows=cfg.report.max_table_rows)
        if report_md:
            write_markdown(ctx, report_md)
    # 22 -----------------------------------------------------------------
    with st.stage("run manifest"):
        st.timings["total"] = time.time() - t0
        provenance.write_json({**prov, "summary": summary, "warnings": warnings, "timings": st.timings, "outputs": outputs, "filtering_steps": audit.to_records(), "tables": {k: str(v) for k, v in table_paths.items()}, "n_figures": len(registry.records)}, run_dir / "logs" / "run_manifest.json")
    # 23 -----------------------------------------------------------------
    archive = None
    if cfg.output.archive:
        with st.stage("archive"):
            archive = archive_results(run_dir, run_name, cfg.output.archive_exclude, cfg.output.archive_name)
    logger.info("done in %.1f s; %d warnings; report: %s", st.timings["total"], len(warnings), report_html)
    return PipelineResult(run_dir, h5ad_path, report_html, report_md, prefilter_h5ad, archive, table_paths, registry, adata, steps, warnings, st.timings, int(adata.n_obs), int(adata.n_vars), summary["n_proteins"], st.timings["total"])


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _qc_figures(adata: ad.AnnData, prot_table: Optional[pd.DataFrame], cfg: Config, registry: FigureRegistry, stage: str, prot_info: Dict[str, Any]) -> None:
    plots.rna_qc_figures(adata.obs, cfg, registry, stage)
    counts = adata.obsm.get("protein_counts")
    prot = adata.obsm.get("protein")
    if counts is not None or prot is not None:
        plots.protein_qc_figures(adata.obs, counts, prot, prot_table, cfg, registry, stage, prot_info.get("primary_method", ""))


def _input_audit(cfg: Config, rna_info: Dict[str, Any], prot_info: Dict[str, Any], pert_info: Dict[str, Any], input_rec: Dict[str, Any], in_prov: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    in_prov = in_prov or {}
    mods = in_prov.get("modalities", {})

    def fmt(name: str) -> str:
        return str(mods.get(name, {}).get("format", getattr(getattr(cfg.inputs, name, None), "format", "-")))

    def src(name: str) -> str:
        m = getattr(cfg.inputs, name, None)
        fallback = ", ".join(Path(p).name for p in m.paths()) if hasattr(m, "paths") else (Path(m.file).name if getattr(m, "file", None) else "-")
        return str(mods.get(name, {}).get("source", fallback))

    rows = []
    rows.append({"modality": "RNA", "format": fmt("rna"), "source": src("rna"), "declared_state": mods.get("rna", {}).get("state", cfg.inputs.rna.state), "detected_state": rna_info.get("input_state_used"), "raw_counts": ("observed" if rna_info.get("counts_layer") == "counts" else "reconstructed" if rna_info.get("counts_layer") == "reconstructed_counts" else "unavailable"), "treatment": rna_info.get("method")})
    if "protein" in mods:
        rows.append({"modality": "protein (normalized)", "format": fmt("protein"), "source": src("protein"), "declared_state": mods["protein"].get("state", cfg.inputs.protein.state), "detected_state": prot_info.get("provided_state", {}).get("state"), "raw_counts": "-", "treatment": prot_info.get("primary_method")})
    if "protein_counts" in mods:
        rows.append({"modality": "protein (raw ADT)", "format": fmt("protein_counts"), "source": src("protein_counts"), "declared_state": mods["protein_counts"].get("state", "raw_counts"), "detected_state": prot_info.get("counts_state", {}).get("state"), "raw_counts": "available", "treatment": "kept in obsm['protein_counts']; extras: " + ", ".join(prot_info.get("extra", {}).keys())})
    asg = pert_info.get("assignment", {})
    if "guide_counts" in mods:
        gcall = asg.get("guide_calling") or {}
        rows.append({"modality": "guide counts", "format": fmt("guide_counts"), "source": src("guide_counts"), "declared_state": "raw_counts", "detected_state": f"{gcall.get('n_guides')} guides; median {gcall.get('median_guide_umis_per_cell')} UMIs/cell", "raw_counts": "available", "treatment": f"obsm['guide_counts'] + diagnostics; {gcall.get('method')} rule " + ("drives the assignment" if asg.get("effective_source") == "guide_counts" else "for comparison only")})
    pert_src = pert_info.get("source")
    rows.append({"modality": "perturbation", "format": "guide list" if pert_src == "guide_assignments" else pert_src, "source": mods.get("guide_assignments", {}).get("source", cfg.columns.guide or "-") if asg.get("effective_source") != "guide_counts" else src("guide_counts"), "declared_state": f"assignments (source: {asg.get('effective_source', 'provided')})", "detected_state": f"controls: {pert_info.get('control_classes')}", "raw_counts": "available (guide-count matrix)" if pert_info.get("guide_counts_available") else "unavailable (no guide-count matrix)", "treatment": "harmonized into obs; multi-guide/unassigned kept" if cfg.perturbation.preserve_multiguide and cfg.perturbation.preserve_unassigned else "harmonized into obs"})
    if cfg.inputs.metadata.file:
        rows.append({"modality": "metadata", "format": cfg.inputs.metadata.format, "source": Path(cfg.inputs.metadata.file).name, "declared_state": "-", "detected_state": "-", "raw_counts": "-", "treatment": f"columns mapped: condition={cfg.columns.condition}, sample={cfg.columns.sample}, donor={cfg.columns.donor}, batch={cfg.columns.batch}"})
    if cfg.inputs.embedding.file:
        rows.append({"modality": "embedding (provided)", "format": cfg.inputs.embedding.format, "source": Path(cfg.inputs.embedding.file).name, "declared_state": "-", "detected_state": "-", "raw_counts": "-", "treatment": f"kept as obsm['{cfg.inputs.embedding.key}'] for comparison"})
    if "lane_metadata" in mods:
        rows.append({"modality": "lane metadata", "format": fmt("lane_metadata"), "source": src("lane_metadata"), "declared_state": "-", "detected_state": "-", "raw_counts": "-", "treatment": f"joined onto cells by lane_id; columns: {mods['lane_metadata'].get('columns')}"})
    modalities = pd.DataFrame(rows).fillna("").astype(str)
    return {"modalities": modalities, "files": {k: {kk: str(vv) for kk, vv in (v or {}).items()} for k, v in input_rec.items()}, "protein_embedding_representation": prot_info.get("embedding_key") or "", "isotype_map": prot_info.get("isotype_map", {}), "provenance": in_prov}


def _run_summary(cfg: Config, adata: ad.AnnData, rna_info, prot_info, pert_qc, n_input, n_prefilter, n_final) -> Dict[str, Any]:
    obs = adata.obs
    return {
        "run_name": cfg.run_name(), "dataset": cfg.dataset.name, "version": __version__,
        "n_cells_input": int(n_input), "n_cells_after_prefilter": int(n_prefilter), "n_cells": int(n_final), "n_genes": int(adata.n_vars),
        "n_proteins": int(adata.obsm["protein"].shape[1]) if "protein" in adata.obsm else 0,
        "n_protein_counts_features": int(adata.obsm["protein_counts"].shape[1]) if "protein_counts" in adata.obsm else 0,
        "rna_input_state": rna_info.get("input_state_used"), "rna_method": rna_info.get("method"), "rna_counts_layer": rna_info.get("counts_layer"), "rna_counts_source": rna_info.get("counts_source"),
        "protein_primary": prot_info.get("primary_method"), "protein_extra_representations": list(prot_info.get("extra", {}).keys()), "protein_embedding_representation": prot_info.get("embedding_key"), "n_isotype_controls": prot_info.get("n_isotype_features"),
        "prefilter_enabled": cfg.qc.prefilter.enabled, "qc_filter_enabled": cfg.qc.filter.enabled, "multimodal_enabled": cfg.multimodal.enabled, "seed": cfg.compute.seed,
        **{f"perturbation_{k}": v for k, v in pert_qc["class_counts"].items()},
        "n_targeting_guides": pert_qc["n_targeting_guides"], "n_targeting_targets": pert_qc["n_targeting_targets"], "n_low_coverage_guides": pert_qc["n_low_coverage_guides"], "n_low_coverage_targets": pert_qc["n_low_coverage_targets"],
        "frac_unassigned": pert_qc["frac_unassigned"], "frac_multi_guide": pert_qc["frac_multi_guide"],
        "conditions": {str(k): int(v) for k, v in obs["condition"].value_counts().items()} if "condition" in obs else None,
    }


def _cards(adata, summary, n_input, n_prefilter, n_final, pert_qc, prot_info):
    cards = [("input cells", f"{n_input:,}"), ("after prefilter", f"{n_prefilter:,}"), ("after strict filtering", f"{n_final:,}"), ("genes", f"{adata.n_vars:,}"), ("proteins (targeting)", f"{summary['n_proteins']}"), ("antibodies incl. isotypes", f"{summary['n_protein_counts_features']}"), ("targeting guides", f"{pert_qc['n_targeting_guides']:,}"), ("targets", f"{pert_qc['n_targeting_targets']:,}")]
    if summary.get("conditions"):
        cards.append(("conditions", f"{len(summary['conditions'])}"))
    return cards


def _save_tsv(df: pd.DataFrame, path: Path) -> Path:  # kept for compatibility
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t")
    return path


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items() if v is not None}
    if isinstance(obj, (list, tuple)):
        if obj and all(isinstance(v, dict) for v in obj):
            return pd.DataFrame([{str(k): _jsonable(x) for k, x in v.items()} for v in obj]).astype(str)  # h5ad cannot hold lists of dicts
        return [_jsonable(v) for v in obj if v is not None]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def _sanitize_for_h5ad(adata: ad.AnnData) -> None:
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str)
    for c in adata.var.columns:
        if adata.var[c].dtype == object:
            adata.var[c] = adata.var[c].astype(str)
    for k in list(adata.uns):
        if isinstance(adata.uns[k], dict) and not adata.uns[k]:
            del adata.uns[k]


def object_schema(adata: ad.AnnData, rna_info: Dict[str, Any], prot_info: Dict[str, Any], cfg: Config) -> Dict[str, Dict[str, str]]:
    """Human-readable description of every non-standard slot (also written to the report)."""
    sch: Dict[str, Dict[str, str]] = {}
    sch["X"] = {"shape": f"{adata.n_obs} x {adata.n_vars}", "content": "normalized RNA expression", "status": rna_info.get("method", ""), "state": rna_info.get("input_state_used", "")}
    for layer in adata.layers:
        if layer == "counts":
            sch[f"layers['{layer}']"] = {"shape": f"{adata.n_obs} x {adata.n_vars}", "content": "observed raw RNA UMI counts (input)", "status": "raw"}
        elif layer == "reconstructed_counts":
            sch[f"layers['{layer}']"] = {"shape": f"{adata.n_obs} x {adata.n_vars}", "content": "integer counts reconstructed from log-normalized X and the per-cell library size", "status": "reconstructed, NOT observed raw counts", "provenance": rna_info.get("counts_source", "")}
    for key in adata.obsm:
        m = adata.obsm[key]
        shape = f"{m.shape[0]} x {m.shape[1]}"
        if key == "protein":
            sch["obsm['protein']"] = {"shape": shape, "content": "primary normalized protein values, targeting antibodies (DataFrame)", "status": prot_info.get("primary_method", ""), "provenance": prot_info.get("primary_source", "")}
        elif key == "protein_counts":
            sch["obsm['protein_counts']"] = {"shape": shape, "content": "raw ADT UMI counts incl. isotype controls (DataFrame; NaN = cell without counts)", "status": "raw"}
        elif key.startswith("protein_"):
            sch[f"obsm['{key}']"] = {"shape": shape, "content": f"protein representation '{key[8:]}' (DataFrame, targeting antibodies)", "status": prot_info.get("extra", {}).get(key[8:], "")}
        elif key in ("ps_scores", "ps_scores_raw"):
            sch[f"obsm['{key}']"] = {"shape": shape, "content": "per-target PS (cells x targets; NaN outside the target's perturbed + control cells)" + (", clipped to [0, scale_factor], not max-normalized" if key.endswith("raw") else ", max-normalized to [0, 1]"), "status": "derived (Stage E)"}
        elif key == "lochness":
            sch["obsm['lochness']"] = {"shape": shape, "content": "lochNESS per target for every cell (local/overall fraction - 1, kNN in PCA space)", "status": "derived (Stage E)"}
        elif key == "program_activity":
            sch["obsm['program_activity']"] = {"shape": shape, "content": "per-cell gene-program activity (mean z-scored log expression of program genes)", "status": "derived (Stage E)"}
        elif key == "guide_counts":
            sch["obsm['guide_counts']"] = {"shape": shape, "content": "guide (sgRNA) UMI counts per cell, sparse; columns = uns['guide_features'].index; zero rows where obs['has_guide_counts'] is False", "status": "raw"}
        elif key == "X_pca":
            sch["obsm['X_pca']"] = {"shape": shape, "content": "RNA PCA on scaled HVGs", "status": "derived"}
        elif key == "X_pca_protein":
            sch["obsm['X_pca_protein']"] = {"shape": shape, "content": f"protein PCA on obsm['{prot_info.get('embedding_key')}'] (isotypes excluded, scaled and decomposed in float64, stored float32)", "status": "derived"}
        elif key == "X_multimodal":
            sch["obsm['X_multimodal']"] = {"shape": shape, "content": f"block-normalized RNA PCs + {cfg.multimodal.protein_weight} x protein PCs", "status": "derived"}
        elif key == "X_lda_umap":
            continue
        elif key.startswith("X_umap_"):
            sch[f"obsm['{key}']"] = {"shape": shape, "content": {"X_umap_rna": "UMAP of RNA neighbors", "X_umap_protein": "UMAP of protein neighbors", "X_umap_multimodal": "UMAP of multimodal neighbors"}.get(key, "provided embedding (input)"), "status": "derived" if key != cfg.inputs.embedding.key else "input"}
    ckey = cfg.analysis.clustering.key
    if cfg.analysis.clustering.enabled and ckey in adata.obs.columns:
        sch[f"obs['{ckey}']"] = {"shape": f"{adata.n_obs} cells, {adata.obs[ckey].nunique()} clusters", "content": "Leiden cluster on the RNA neighbour graph (numbered states, not cell types)", "status": "derived (Stage F)"}
    if "perturbation_cluster_enrichment" in adata.uns:
        sch["uns['perturbation_cluster_enrichment']"] = {"shape": f"{len(adata.uns['perturbation_cluster_enrichment'])} rows", "content": "target x cluster Fisher/CMH enrichment under both control arms (Haldane odds ratio, p, BH-FDR per arm, direction, guide concordance)", "status": "derived"}
    if "perturbation_distance" in adata.uns:
        sch["uns['perturbation_distance']"] = {"shape": f"{len(adata.uns['perturbation_distance'])} rows", "content": "per-target energy distance vs control in X_pca, DistanceTest permutation p, BH-FDR, call (reference distance stage)", "status": "derived"}
    if "perturbation_distance_matrix" in adata.uns:
        sch["uns['perturbation_distance_matrix']"] = {"shape": f"{adata.uns['perturbation_distance_matrix'].shape[0]} x {adata.uns['perturbation_distance_matrix'].shape[1]} targets", "content": "pairwise target x target energy distance (symmetric, zero diagonal); PCoA / phenotype modules derive from it", "status": "derived"}
    if "phenotype_modules" in adata.uns:
        sch["uns['phenotype_modules']"] = {"shape": f"{len(adata.uns['phenotype_modules'])} rows", "content": "phenotype module (average-linkage cluster of the pairwise distance matrix) per target", "status": "derived"}
    if "master_perturbation_table" in adata.uns:
        sch["uns['master_perturbation_table']"] = {"shape": f"{len(adata.uns['master_perturbation_table'])} rows", "content": "master perturbation table: efficacy, PS, lochNESS, distance, modules per target (reference meta) + protein / cluster extension columns", "status": "derived"}
    if "perturbation_strength" in adata.uns:
        sch["uns['perturbation_strength']"] = {"shape": f"{len(adata.uns['perturbation_strength'])} rows", "content": "per-target knockdown of the target's own expression vs ntc / other controls (log2FC, KS, MWU, BH-FDR, hit call, rank)", "status": "derived"}
    if "X_lda_umap" in adata.obsm:
        sch["obsm['X_lda_umap']"] = {"shape": f"{adata.n_obs} x 2", "content": "supervised LDA-UMAP of the scored targets + controls (PS_python compute_lda_umap); NaN for cells outside the trained classes", "status": "derived"}
    if "guide_features" in adata.uns:
        sch["uns['guide_features']"] = {"shape": f"{len(adata.uns['guide_features'])} rows", "content": "one row per guide: target, control_class, total_umis, n_cells_detected, n_cells_dominant (+ 10x feature columns)"}
    sch["uns['protein_features']"] = {"shape": f"{len(adata.uns['protein_features'])} rows" if "protein_features" in adata.uns else "", "content": "one row per antibody: role, isotype control, presence in each matrix, used_in_embedding/used_in_qc"}
    sch["uns['petrubseq_protein']"] = {"content": "version, inputs audit, normalization decisions, QC settings/results incl. filtering steps, alignment, perturbation QC, embeddings, warnings, provenance, config, this schema"}
    return sch
