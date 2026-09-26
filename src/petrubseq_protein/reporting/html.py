"""Self-contained HTML report (jinja2 template + optional base64-embedded figures).

Section order follows the reference Perturb-seq report (QC, clustering,
perturbation strength, cluster enrichment, PS, lochNESS, modules/programs) with
the protein extension appended (protein QC, protein effects, RNA-protein
concordance, multimodal summaries). Sections whose analysis did not run are
omitted and the numbering closes up.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup

from .context import ReportContext
from .figures import SECTION_MULTIMODAL, SECTION_QC_PERTURBATION, SECTION_QC_PROTEIN, SECTION_QC_RNA, SECTION_REP_PROTEIN, SECTION_REP_RNA, STAGE_AFTER, STAGE_BEFORE, FigureRecord, FigureRegistry

logger = logging.getLogger("petrubseq_protein")
TEMPLATE_DIR = Path(__file__).parent / "templates"


def df_to_html(df: Optional[pd.DataFrame], max_rows: int = 100, index: bool = True) -> Markup:
    if df is None or len(df) == 0:
        return Markup('<p class="sub">Not available for this run.</p>')
    shown = df.head(max_rows)
    html = shown.to_html(index=index, escape=True, border=0, na_rep="", float_format=lambda v: f"{v:.4g}")
    if len(df) > max_rows:
        html += f'<p class="sub">Showing {max_rows} of {len(df)} rows; the full table is in <code>tables/</code>.</p>'
    return Markup(html)


def kv_table(d: Dict[str, Any]) -> Markup:
    rows = []
    for k, v in (d or {}).items():
        if isinstance(v, (dict, list)):
            v = json.dumps(v, default=str)
            if len(v) > 400:
                v = v[:400] + " ..."
        rows.append({"key": k, "value": v})
    return df_to_html(pd.DataFrame(rows), max_rows=500, index=False)


def render_figure(fig: FigureRecord, embed: bool, run_dir: Path) -> Markup:
    src = fig.data_uri() if embed else fig.rel(run_dir)
    return Markup(f'<figure><img src="{src}" alt="{Markup.escape(fig.title)}"><figcaption><b>{Markup.escape(fig.title)}.</b> {Markup.escape(fig.caption)}</figcaption></figure>')


def _t(df, cols=None, n=100, index=False):
    if df is None or getattr(df, "empty", True):
        return Markup("<p class='sub'>(not computed)</p>")
    d = df[[c for c in cols if c in df.columns]] if cols else df
    return df_to_html(d.round(4) if hasattr(d, "round") else d, max_rows=n, index=index)


def _extras(reg: FigureRegistry, section: str, stage: str) -> List[str]:
    return [f"{r.name}.{r.path.suffix.lstrip('.')}" for r in reg.by_section(section, stage, only_in_report=False) if not r.in_report]


def section_views(ctx: ReportContext, max_rows: int) -> Dict[str, Any]:
    """Everything the template needs beyond the QC tables: per-section contexts, figure lists, numbering."""
    from ..analysis.perturbation_strength import CONTROL_LABELS
    from .cell_state_plots import SECTION as CS_SECTION, ST_CLUSTERS
    from .enrichment_plots import ST_ENRICH, ST_PER_TARGET as ST_ENRICH_TARGET
    from .perturbation_plots import SECTION as PE_SECTION, ST_CONC, ST_LDA, ST_LDA_TARGET, ST_LOCH, ST_LOCH_TARGET, ST_PROG, ST_PROG_UMAP, ST_PROT, ST_PS, ST_PS_TARGET
    from .strength_plots import SECTION as STR_SECTION, ST_OVERVIEW, ST_PER_TARGET as ST_STR_TARGET

    reg = ctx.registry
    figs = {
        "qc_rna_before": reg.by_section(SECTION_QC_RNA, STAGE_BEFORE), "qc_rna_after": reg.by_section(SECTION_QC_RNA, STAGE_AFTER),
        "qc_protein_before": reg.by_section(SECTION_QC_PROTEIN, STAGE_BEFORE), "qc_protein_after": reg.by_section(SECTION_QC_PROTEIN, STAGE_AFTER),
        "qc_perturbation": reg.by_section(SECTION_QC_PERTURBATION), "rep_rna": reg.by_section(SECTION_REP_RNA), "rep_protein": reg.by_section(SECTION_REP_PROTEIN), "multimodal": reg.by_section(SECTION_MULTIMODAL),
        "cs_clusters": reg.by_section(CS_SECTION, ST_CLUSTERS), "cs_enrichment": reg.by_section(CS_SECTION, ST_ENRICH), "cs_enrichment_per_target": reg.by_section(CS_SECTION, ST_ENRICH_TARGET),
        "strength_overview": reg.by_section(STR_SECTION, ST_OVERVIEW), "strength_per_target": reg.by_section(STR_SECTION, ST_STR_TARGET),
        "pe_ps": reg.by_section(PE_SECTION, ST_PS), "pe_ps_per_target": reg.by_section(PE_SECTION, ST_PS_TARGET), "pe_ps_lda": reg.by_section(PE_SECTION, ST_LDA), "pe_ps_lda_per_target": reg.by_section(PE_SECTION, ST_LDA_TARGET),
        "pe_lochness": reg.by_section(PE_SECTION, ST_LOCH), "pe_lochness_per_target": reg.by_section(PE_SECTION, ST_LOCH_TARGET),
        "pe_gene_programs": reg.by_section(PE_SECTION, ST_PROG), "pe_gene_programs_umap": reg.by_section(PE_SECTION, ST_PROG_UMAP),
        "pe_protein_effects": reg.by_section(PE_SECTION, ST_PROT), "pe_concordance": reg.by_section(PE_SECTION, ST_CONC),
    }
    extras = {"strength": _extras(reg, STR_SECTION, ST_STR_TARGET), "enrichment": _extras(reg, CS_SECTION, ST_ENRICH_TARGET), "ps": _extras(reg, PE_SECTION, ST_PS_TARGET), "lda": _extras(reg, PE_SECTION, ST_LDA_TARGET),
              "lochness": _extras(reg, PE_SECTION, ST_LOCH_TARGET), "programs": _extras(reg, PE_SECTION, ST_PROG_UMAP)}
    pe = ctx.extra.get("perturbation_effects")
    cs = ctx.extra.get("cell_states")
    rep = ctx.representations.get("rna", {})
    # --- clustering ----------------------------------------------------------------------
    cs_view = None
    if cs is not None and cs.clustering is not None:
        info = cs.clustering.info
        cs_view = {"summary": _t(cs.clustering.summary, n=max_rows), "flags": cs.clustering.flags, "resolution": info.get("resolution"), "n_clusters": info.get("n_clusters"),
                   "n_hvg": (rep.get("hvg") or {}).get("n_top_genes", ""), "n_pcs": (rep.get("pca") or {}).get("n_comps", ""), "n_neighbors": (rep.get("neighbors") or {}).get("n_neighbors", "")}
    # --- perturbation strength ------------------------------------------------------------
    strength_view = None
    if pe is not None and pe.strength is not None and not pe.strength.empty:
        s = pe.strength
        strength_view = {"controls_described": " and ".join(CONTROL_LABELS[c] for c in s.controls_used), "primary_label": CONTROL_LABELS[s.primary_control], "fallback": s.primary_control != s.info.get("primary_control", s.primary_control) or s.info.get("requested_primary_unavailable", False),
                         "fdr_alpha": s.info.get("fdr_alpha"), "cards": [("Targets tested", f"{len(s.table):,}"), ("Effective knockdowns", f"{len(s.hits):,}"), ("Control cells", f"{s.n_control_cells.get(s.primary_control, 0):,}"), ("Targets not testable", f"{len(s.skipped):,}")],
                         "table": df_to_html(s.display, max_rows=max_rows, index=False), "skipped": df_to_html(s.skipped, max_rows=max_rows, index=False) if len(s.skipped) else ""}
    # --- enrichment ------------------------------------------------------------------------
    enr_view = None
    if cs is not None and cs.enrichment is not None and not cs.enrichment.empty:
        from ..analysis.cluster_enrichment import format_enrichment_table
        e = cs.enrichment
        om = e.omnibus or {}
        enr_view = {"n_hits": int(e.table["significant"].sum()), "n_targets_with_hits": len(e.targets_with_hits()), "n_targets": int(e.composition.shape[0]), "n_clusters": int(e.composition.shape[1]), "n_tests": int(e.composition.shape[0] * e.composition.shape[1]),
                    "control_label": CONTROL_LABELS[e.primary_control], "controls_described": " and ".join(CONTROL_LABELS[c] for c in e.controls_used), "chi2": f"{om.get('chi2', float('nan')):.0f}", "dof": om.get("dof", 0), "p_perm": f"{om.get('p_permutation', float('nan')):.3g}",
                    "pct_small": f"{om.get('pct_expected_below_5', 0):.0f}", "stratified": e.stratified, "stratify_by": e.stratify_by, "n_low_power": int(e.info.get("n_low_power", 0)), "fdr_alpha": e.info.get("fdr_alpha"),
                    "top_shift": e.effect_magnitude.iloc[0].to_dict() if len(e.effect_magnitude) else {}, "table": df_to_html(format_enrichment_table(e), max_rows=max_rows, index=False),
                    "comparison": df_to_html(e.lochness_comparison.round(4), max_rows=max_rows, index=False) if not e.lochness_comparison.empty else ""}
    # --- PS / lochNESS / modules ---------------------------------------------------------------
    ps_view = loch_view = mods_view = None
    if pe is not None and pe.ps is not None and not pe.ps.summary.empty:
        summ = pe.ps.summary
        ps_view = {"n_targets": int(len(summ)), "threshold": pe.ps.ps_threshold, "cut_method": pe.ps.info.get("expression_cut"), "median_kd": f"{summ['pct_successful_kd'].median():.0f}", "median_escaper": f"{summ['pct_escaper'].median():.0f}",
                   "best": summ.iloc[0]["target"], "best_kd": f"{summ.iloc[0]['pct_successful_kd']:.0f}", "worst_escaper": summ.sort_values("pct_escaper").iloc[-1]["target"], "worst_escaper_pct": f"{summ['pct_escaper'].max():.0f}",
                   "n_skipped": int(len(pe.ps.skipped)), "has_lda": pe.ps.lda_umap is not None, "lda_note": pe.ps.lda_note,
                   "table": _t(summ, ["target", "gene", "n_perturbed_cells", "n_control_cells", "mean_ps", "median_ps", "pct_high_ps", "pct_successful_kd", "pct_escaper", "pct_non_responder", "pct_low_signal", "pct_controls_called_kd", "net_pct_kd", "expression_cut", "auc_vs_control", "d_vs_control", "warnings"], max_rows),
                   "skipped": _t(pe.ps.skipped, n=max_rows) if len(pe.ps.skipped) else ""}
    if pe is not None and pe.lochness is not None and not pe.lochness.summary.empty:
        s0 = pe.lochness.summary
        cut = pe.lochness.info.get("enrichment_cut")
        loch_view = {"n_targets": int(len(s0)), "k": pe.lochness.n_neighbors, "cut": cut, "best": s0.iloc[0]["target"], "best_score": f"{s0.iloc[0]['mean_lochness_in_own_cells']:.2f}", "n_positive": int((s0["mean_lochness_in_own_cells"] > cut).sum()),
                     "table": _t(s0, ["target", "n_cells", "overall_fraction_pct", "mean_lochness_all_cells", "mean_lochness_in_own_cells", "max_lochness", "pct_cells_enriched", "top_cluster", "top_cluster_mean", "mean_lochness_in_control_cells", "delta_own_vs_control", "null_mean", "null_sd", "z_score", "p_empirical", "fdr", "max_sample_share"], max_rows)}
    if pe is not None and pe.modules is not None and not pe.modules.empty:
        m = pe.modules
        mods_view = {"n_modules": m.n_modules, "n_programs": m.n_programs, "n_perturbations": int(m.effect.shape[0]), "n_genes": int(m.effect.shape[1]), "control_label": CONTROL_LABELS.get(m.control, m.control), "gene_selection": m.info.get("gene_selection"),
                     "module_correlation": m.info.get("module_correlation"), "program_correlation": m.info.get("program_correlation"), "linkage": m.info.get("linkage_method"), "hub_lfc": m.info.get("de_lfc_threshold"), "de_fdr": m.info.get("de_fdr_alpha"),
                     "top_hub": m.info.get("top_hub", ""), "top_hub_n": m.info.get("top_hub_n_de", 0), "n_tf_edges": m.info.get("n_tf_edges", 0), "program_genes": {p: m.program_genes.get(p, [])[:10] for p in m.program_labels},
                     "modules": _t(m.perturbation_modules, n=max_rows), "programs": _t(m.gene_programs, ["gene", "program", "program_size", "mean_log2fc", "mean_abs_log2fc", "n_targets_de", "n_targets_up", "n_targets_down"], max_rows),
                     "module_program": _t(m.module_program.astype(float), n=max_rows, index=True), "hubs": _t(m.hubs, n=max_rows)}
    pe_view = None
    if pe is not None:
        info = pe.info()
        pe_view = {"info": {k: kv_table({kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict))}) for k, v in info.items() if isinstance(v, dict)},
                   "protein": _t(pe.protein.table if pe.protein else None, ["target", "protein", "n_perturbed", "n_control", "mean_perturbed", "mean_control", "effect", "cohen_d", "p_value", "fdr", "n_samples_tested", "n_samples_same_sign", "n_guides_tested", "n_guides_same_sign", "status"], 60),
                   "ps_protein": _t(pe.concordance.ps_protein if pe.concordance else None, ["scope", "target", "protein", "n_cells", "spearman_rho", "p_value", "fdr", "status"], 40),
                   "lochness_protein": _t(pe.concordance.lochness_protein_summary if pe.concordance else None), "program_protein": _t(pe.concordance.program_protein_cells if pe.concordance else None, None, 40),
                   "program_protein_targets": _t(pe.concordance.program_protein_targets if pe.concordance else None, None, 40), "summary": _t(pe.concordance.summary if pe.concordance else None, None, 60)}
    # --- numbering ---------------------------------------------------------------------------
    nav: List[tuple] = [("inputs", "Inputs"), ("qc", "Quality control"), ("clustering", "Clustering" if cs_view else "Embedding")]
    if strength_view: nav.append(("strength", "Perturbation strength"))
    if enr_view: nav.append(("enrichment", "Cluster enrichment"))
    if ps_view: nav.append(("ps", "Per-cell response"))
    if loch_view: nav.append(("lochness", "lochNESS"))
    if mods_view: nav.append(("modules", "Modules & programs"))
    nav.append(("protein-qc", "Protein QC"))
    if pe_view:
        nav += [("protein-effects", "Protein effects"), ("concordance", "RNA-protein concordance")]
    nav += [("multimodal", "Multimodal summaries"), ("outputs", "Outputs & provenance")]
    secn = {key: i + 1 for i, (key, _) in enumerate(nav)}
    secn.update({"protein_qc": secn["protein-qc"], "protein_effects": secn.get("protein-effects", 0), "outputs": secn["outputs"]})
    return {"figs": figs, "extras": extras, "cs": cs_view, "strength": strength_view, "enr": enr_view, "ps": ps_view, "loch": loch_view, "mods": mods_view, "pe": pe_view, "nav": nav, "secn": secn}


def write_html(ctx: ReportContext, path: Path, embed: bool = True, max_rows: int = 100) -> Path:
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), autoescape=select_autoescape(["html"]))
    tpl = env.get_template("report.html")
    reg = ctx.registry
    tcov = ctx.target_coverage.sort_values("n_cells_single_guide").head(25) if ctx.target_coverage is not None else None
    tables = {
        "inputs": df_to_html(ctx.inputs, index=False), "input_files": df_to_html(ctx.input_files, index=False), "alignment": kv_table(ctx.alignment),
        "filtering_steps": df_to_html(ctx.filtering_steps, index=False), "qc_summary": df_to_html(ctx.qc_summary, index=False),
        "protein_info": kv_table(ctx.protein_info), "protein_features": df_to_html(ctx.protein_features, max_rows=max_rows), "protein_qc": kv_table(ctx.protein_qc),
        "perturbation_qc": kv_table(ctx.perturbation_qc), "class_by_condition": df_to_html(ctx.class_by_condition), "moi_summary": df_to_html(ctx.moi_summary), "target_coverage": df_to_html(tcov, max_rows=max_rows),
        "rep_rna": kv_table(ctx.representations.get("rna", {})), "rep_protein": kv_table(ctx.representations.get("protein", {})),
        "diagnostics": kv_table({**ctx.representations.get("multimodal", {}), **{f"diag:{k}": v for k, v in (ctx.diagnostics or {}).items()}}),
        "outputs": df_to_html(pd.DataFrame([{"deliverable": k, "path": v} for k, v in ctx.outputs.items()]), index=False),
        "schema": df_to_html(pd.DataFrame.from_dict(ctx.schema, orient="index").reset_index().rename(columns={"index": "slot"}), max_rows=100, index=False),
        "timings": df_to_html(pd.DataFrame([{"stage": k, "seconds": round(v, 1)} for k, v in ctx.timings.items()]), max_rows=100, index=False),
        "table_manifest": df_to_html(pd.DataFrame([{"table": k, "path": str(Path(v).relative_to(ctx.run_dir)) if str(v).startswith(str(ctx.run_dir)) else str(v)} for k, v in ctx.tables.items()]), max_rows=200, index=False),
        "figure_manifest": df_to_html(reg.manifest().drop(columns=["caption"], errors="ignore"), max_rows=500, index=False),
        "versions": kv_table(ctx.versions), "provenance": kv_table({k: v for k, v in ctx.provenance.items() if k not in ("packages", "inputs", "config")}),
    }
    views = section_views(ctx, max_rows)
    html = tpl.render(ctx=ctx, tables=tables, render=lambda f: render_figure(f, embed, ctx.run_dir), group_summaries=[(k, df_to_html(v)) for k, v in ctx.group_summaries.items()], n_tables=len(ctx.tables), n_figures=len(reg.records), **views)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    logger.info("wrote %s (%.1f MB, figures %s)", path.name, path.stat().st_size / 1e6, "embedded" if embed else "linked")
    return path
