"""Self-contained HTML report (jinja2 template + optional base64-embedded figures)."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup

from .context import ReportContext
from .figures import (
    SECTION_MULTIMODAL,
    SECTION_QC_PERTURBATION,
    SECTION_QC_PROTEIN,
    SECTION_QC_RNA,
    SECTION_REP_PROTEIN,
    SECTION_REP_RNA,
    STAGE_AFTER,
    STAGE_BEFORE,
    FigureRecord,
)

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


def write_html(ctx: ReportContext, path: Path, embed: bool = True, max_rows: int = 100) -> Path:
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), autoescape=select_autoescape(["html"]))
    tpl = env.get_template("report.html")
    reg = ctx.registry
    figs = {
        "qc_rna_before": reg.by_section(SECTION_QC_RNA, STAGE_BEFORE),
        "qc_rna_after": reg.by_section(SECTION_QC_RNA, STAGE_AFTER),
        "qc_protein_before": reg.by_section(SECTION_QC_PROTEIN, STAGE_BEFORE),
        "qc_protein_after": reg.by_section(SECTION_QC_PROTEIN, STAGE_AFTER),
        "qc_perturbation": reg.by_section(SECTION_QC_PERTURBATION),
        "rep_rna": reg.by_section(SECTION_REP_RNA),
        "rep_protein": reg.by_section(SECTION_REP_PROTEIN),
        "multimodal": reg.by_section(SECTION_MULTIMODAL),
    }
    tcov = ctx.target_coverage.sort_values("n_cells_single_guide").head(25) if ctx.target_coverage is not None else None
    tables = {
        "inputs": df_to_html(ctx.inputs, index=False),
        "input_files": df_to_html(ctx.input_files, index=False),
        "alignment": kv_table(ctx.alignment),
        "filtering_steps": df_to_html(ctx.filtering_steps, index=False),
        "qc_summary": df_to_html(ctx.qc_summary, index=False),
        "protein_info": kv_table(ctx.protein_info),
        "protein_features": df_to_html(ctx.protein_features, max_rows=max_rows),
        "protein_qc": kv_table(ctx.protein_qc),
        "perturbation_qc": kv_table(ctx.perturbation_qc),
        "class_by_condition": df_to_html(ctx.class_by_condition),
        "moi_summary": df_to_html(ctx.moi_summary),
        "target_coverage": df_to_html(tcov, max_rows=max_rows),
        "rep_rna": kv_table(ctx.representations.get("rna", {})),
        "rep_protein": kv_table(ctx.representations.get("protein", {})),
        "diagnostics": kv_table({**ctx.representations.get("multimodal", {}), **{f"diag:{k}": v for k, v in (ctx.diagnostics or {}).items()}}),
        "outputs": df_to_html(pd.DataFrame([{"deliverable": k, "path": v} for k, v in ctx.outputs.items()]), index=False),
        "schema": df_to_html(pd.DataFrame.from_dict(ctx.schema, orient="index").reset_index().rename(columns={"index": "slot"}), max_rows=100, index=False),
        "timings": df_to_html(pd.DataFrame([{"stage": k, "seconds": round(v, 1)} for k, v in ctx.timings.items()]), max_rows=100, index=False),
        "table_manifest": df_to_html(pd.DataFrame([{"table": k, "path": str(Path(v).relative_to(ctx.run_dir)) if str(v).startswith(str(ctx.run_dir)) else str(v)} for k, v in ctx.tables.items()]), max_rows=100, index=False),
        "figure_manifest": df_to_html(reg.manifest().drop(columns=["caption"], errors="ignore"), max_rows=500, index=False),
        "versions": kv_table(ctx.versions),
        "provenance": kv_table({k: v for k, v in ctx.provenance.items() if k not in ("packages", "inputs", "config")}),
    }
    pe = ctx.extra.get("perturbation_effects")
    pe_view = None
    if pe is not None:
        from .perturbation_plots import SECTION as PE_SECTION, ST_CONC, ST_LOCH, ST_PROG, ST_PROT, ST_PS

        def _t(df, cols=None, n=max_rows):
            if df is None or getattr(df, "empty", True):
                return Markup("<p class='sub'>(not computed)</p>")
            d = df[[c for c in cols if c in df.columns]] if cols else df
            return df_to_html(d.round(4) if hasattr(d, "round") else d, max_rows=n, index=False)

        info = pe.info()
        pe_view = {
            "info": {k: kv_table({kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict))}) for k, v in info.items() if isinstance(v, dict)},
            "ps": _t(pe.ps.summary if pe.ps else None, ["target", "n_perturbed_cells", "n_control_cells", "median_ps", "mean_ps", "q25_ps", "q75_ps", "pct_high_ps", "pct_high_ps_control", "auc_vs_control", "d_vs_control", "net_pct_kd", "status", "warnings"]),
            "ps_skipped": _t(pe.ps.skipped if pe.ps else None),
            "lochness": _t(pe.lochness.summary if pe.lochness else None, ["target", "n_cells", "mean_lochness_in_own_cells", "median_lochness_in_own_cells", "mean_lochness_in_control_cells", "delta_own_vs_control", "pct_own_cells_enriched", "null_mean", "null_sd", "z_score", "p_empirical", "fdr", "max_sample_share"]),
            "modules": _t(pe.modules.perturbation_modules if pe.modules else None),
            "programs": _t(pe.modules.gene_programs.groupby("program").head(10) if pe.modules is not None and not pe.modules.empty else None, ["program", "program_size", "gene", "mean_log2fc", "mean_abs_log2fc", "n_targets_de", "n_targets_up", "n_targets_down"]),
            "protein": _t(pe.protein.table if pe.protein else None, ["target", "protein", "n_perturbed", "n_control", "mean_perturbed", "mean_control", "effect", "cohen_d", "p_value", "fdr", "n_samples_tested", "n_samples_same_sign", "n_guides_tested", "n_guides_same_sign", "status"], 60),
            "ps_protein": _t(pe.concordance.ps_protein if pe.concordance else None, ["scope", "target", "protein", "n_cells", "spearman_rho", "p_value", "fdr", "status"], 40),
            "lochness_protein": _t(pe.concordance.lochness_protein_summary if pe.concordance else None),
            "program_protein": _t(pe.concordance.program_protein_cells if pe.concordance else None, None, 40),
            "program_protein_targets": _t(pe.concordance.program_protein_targets if pe.concordance else None, None, 40),
            "summary": _t(pe.concordance.summary if pe.concordance else None, None, 60),
        }
        for st_ in (ST_PS, ST_LOCH, ST_PROG, ST_PROT, ST_CONC):
            figs[f"pe_{st_}"] = reg.by_section(PE_SECTION, st_)
    html = tpl.render(ctx=ctx, tables=tables, figs=figs, pe=pe_view, render=lambda f: render_figure(f, embed, ctx.run_dir), group_summaries=[(k, df_to_html(v)) for k, v in ctx.group_summaries.items()], n_tables=len(ctx.tables), n_figures=len(reg.records))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    logger.info("wrote %s (%.1f MB, figures %s)", path.name, path.stat().st_size / 1e6, "embedded" if embed else "linked")
    return path
