"""Markdown mirror of the HTML report, rendered from the same ReportContext and section views."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from .context import ReportContext
from .figures import SECTION_MULTIMODAL, SECTION_QC_PERTURBATION, SECTION_QC_PROTEIN, SECTION_QC_RNA, SECTION_REP_PROTEIN, SECTION_REP_RNA, STAGE_AFTER, STAGE_BEFORE


def md_table(df: Optional[pd.DataFrame], max_rows: int = 40, max_cols: int = 14, index: bool = True) -> str:
    if df is None or len(df) == 0:
        return "_(not available)_"
    d = df.reset_index() if index else df.copy()
    d = d.iloc[:max_rows, :max_cols]
    cols = [str(c) for c in d.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in d.iterrows():
        lines.append("| " + " | ".join(f"{v:.4g}" if isinstance(v, float) else str(v).replace("|", "\\|") for v in row.tolist()) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_({len(df)} rows; first {max_rows} shown)_")
    return "\n".join(lines)


def kv(d: Dict[str, Any]) -> str:
    lines = ["| key | value |", "|---|---|"]
    for k, v in (d or {}).items():
        if isinstance(v, (dict, list)):
            v = json.dumps(v, default=str)
            if len(v) > 300:
                v = v[:300] + " ..."
        lines.append(f"| {k} | {v} |")
    return "\n".join(lines)


def _figs(ctx: ReportContext, recs) -> str:
    return "\n\n".join(f"![{r.title}]({r.rel(ctx.run_dir)})\n\n_{r.title}. {r.caption}_" for r in recs) if recs else "_(none)_"


def _strip(html_markup) -> str:
    """Plain text of a small HTML snippet (used for the section views' tables)."""
    return re.sub(r"<[^>]+>", " ", str(html_markup)).strip()


def write_markdown(ctx: ReportContext, path: Path) -> Path:
    from .html import section_views

    v = section_views(ctx, max_rows=40)
    figs, secn = v["figs"], v["secn"]
    pe = ctx.extra.get("perturbation_effects")
    cs = ctx.extra.get("cell_states")
    s: List[str] = [f"# {ctx.title}", "", f"Run {ctx.run_name} (dataset {ctx.dataset}), generated {ctx.generated_at}, petrubseq-pipeline-protein v{ctx.version}. Figures are linked by relative path (report.html embeds them).", ""]
    s += ["## Run summary", "", kv(dict(ctx.cards)), ""]
    s += ["## Warnings", ""] + ([f"- {w}" for w in ctx.warnings] or ["_(none)_"]) + [""]
    s += [f"> {n}" for n in ctx.notes] + [""]
    s += [f"## {secn['inputs']}. Inputs and modality audit", "", md_table(ctx.inputs, index=False), "", md_table(ctx.input_files, index=False), "", kv(ctx.alignment), ""]
    s += [f"## {secn['qc']}. Quality control", "", "### Filtering steps", "", md_table(ctx.filtering_steps, index=False), "", "### QC summary", "", md_table(ctx.qc_summary, index=False), "", "### Before filtering", "", _figs(ctx, figs["qc_rna_before"]), "", "### After filtering", "", _figs(ctx, figs["qc_rna_after"]), "",
          "### Perturb-seq guide QC", "", kv(ctx.perturbation_qc), "", md_table(ctx.class_by_condition), "", md_table(ctx.moi_summary), "", md_table(ctx.target_coverage.sort_values("n_cells_single_guide").head(20) if ctx.target_coverage is not None else None, max_rows=20), "", _figs(ctx, figs["qc_perturbation"]), ""]
    s += [f"## {secn['clustering']}. " + ("Clustering analysis" if v["cs"] else "Embedding"), "", kv(ctx.representations.get("rna", {})), ""]
    if cs is not None and cs.clustering is not None:
        s += [kv({k: x for k, x in cs.clustering.info.items() if not isinstance(x, (list, dict))}), ""] + [f"- warning: {f}" for f in cs.clustering.flags] + ["", md_table(cs.clustering.summary.round(3), index=False, max_rows=60), ""]
    s += [_figs(ctx, figs["rep_rna"]), "", _figs(ctx, figs["cs_clusters"]), ""]
    if v["strength"] and pe is not None and pe.strength is not None:
        st = pe.strength
        s += [f"## {secn['strength']}. Perturbation strength", "", f"Controls: {v['strength']['controls_described']}; primary {v['strength']['primary_label']}.", "", kv(dict(v["strength"]["cards"])), "", _figs(ctx, figs["strength_overview"]), "", md_table(st.display, index=False, max_rows=60), "", "### Targets that could not be tested", "", md_table(st.skipped, index=False), "", "### Strongest perturbation effects", "", _figs(ctx, figs["strength_per_target"]), ""]
    if v["enr"] and cs is not None and cs.enrichment is not None:
        from ..analysis.cluster_enrichment import format_enrichment_table
        e = cs.enrichment
        s += [f"## {secn['enrichment']}. Perturbation enrichment across clusters", "", kv({k: x for k, x in e.info.items() if not isinstance(x, (list, dict))}), "", kv(e.omnibus), "", _figs(ctx, figs["cs_enrichment"]), "", md_table(format_enrichment_table(e), index=False, max_rows=40), "", _figs(ctx, figs["cs_enrichment_per_target"]), ""]
        if not e.lochness_comparison.empty:
            s += ["### Cluster enrichment and lochNESS (descriptive)", "", md_table(e.lochness_comparison.round(4), index=False, max_rows=60), ""]
    if v["ps"] and pe is not None and pe.ps is not None:
        s += [f"## {secn['ps']}. Per-cell perturbation response", "", kv({k: x for k, x in pe.ps.info.items() if not isinstance(x, (list, dict))}), "", _figs(ctx, figs["pe_ps"]), "", md_table(pe.ps.summary.round(4), index=False, max_rows=40), "", _figs(ctx, figs["pe_ps_per_target"]), "", "### Supervised LDA embedding", "", _figs(ctx, figs["pe_ps_lda"]), "", _figs(ctx, figs["pe_ps_lda_per_target"]), ""]
    if v["loch"] and pe is not None and pe.lochness is not None:
        s += [f"## {secn['lochness']}. lochNESS neighbourhood enrichment", "", kv({k: x for k, x in pe.lochness.info.items() if not isinstance(x, (list, dict))}), "", _figs(ctx, figs["pe_lochness"]), "", md_table(pe.lochness.summary.round(4), index=False, max_rows=40), "", _figs(ctx, figs["pe_lochness_per_target"]), ""]
    if v["mods"] and pe is not None and pe.modules is not None:
        m = pe.modules
        s += [f"## {secn['modules']}. Co-functional modules and gene programs", "", kv({k: x for k, x in m.info.items() if not isinstance(x, (list, dict))}), "", _figs(ctx, figs["pe_gene_programs"]), "", _figs(ctx, figs["pe_gene_programs_umap"]), "", md_table(m.perturbation_modules.round(4), index=False, max_rows=60), "", md_table(m.gene_programs.round(4), index=False, max_rows=40), "", md_table(m.module_program.astype(float).round(4)), "", md_table(m.hubs, index=False, max_rows=40), ""]
    s += [f"## {secn['protein_qc']}. Protein QC", "", kv(ctx.protein_info), "", md_table(ctx.protein_features, max_rows=60), "", kv(ctx.protein_qc), "", "### Before filtering", "", _figs(ctx, figs["qc_protein_before"]), "", "### After filtering", "", _figs(ctx, figs["qc_protein_after"]), "", "### Protein representation", "", kv(ctx.representations.get("protein", {})), "", _figs(ctx, figs["rep_protein"]), ""]
    if pe is not None:
        s += [f"## {secn['protein_effects']}. Protein effects", "", kv({k: x for k, x in (pe.protein.info if pe.protein is not None else {"status": "disabled"}).items() if not isinstance(x, (list, dict))}), "", md_table(pe.protein.table.round(4) if pe.protein is not None and not pe.protein.empty else None, index=False, max_rows=40), "", _figs(ctx, figs["pe_protein_effects"]), ""]
        s += [f"## {secn['concordance']}. RNA-protein concordance", "", kv({k: x for k, x in (pe.concordance.info if pe.concordance is not None else {"status": "disabled"}).items() if not isinstance(x, (list, dict))}), "", md_table(pe.concordance.ps_protein.round(4) if pe.concordance is not None else None, index=False, max_rows=40), "", _figs(ctx, figs["pe_concordance"]), ""]
    s += [f"## {secn['multimodal']}. Multimodal summaries", "", kv({**ctx.representations.get("multimodal", {}), **ctx.diagnostics}), "", _figs(ctx, figs["multimodal"]), ""]
    for name, df in ctx.group_summaries.items():
        s += [f"### {name} summary", "", md_table(df), ""]
    if pe is not None and pe.concordance is not None:
        s += ["### Integrated perturbation summary", "", md_table(pe.concordance.summary.round(4), index=False, max_rows=60), ""]
    s += [f"## {secn['outputs']}. Outputs and provenance", "", kv(ctx.outputs), "", "```", ctx.adata_repr, "```", "", "| slot | shape | content | status | provenance |", "|---|---|---|---|---|"]
    s += [f"| `{k}` | {d.get('shape','')} | {d.get('content','')} | {d.get('status','')} | {d.get('provenance','')} |" for k, d in ctx.schema.items()]
    s += ["", "### Stage timings (seconds)", "", kv({k: round(x, 1) for k, x in ctx.timings.items()}), "", "### Software versions", "", kv(ctx.versions), "", "### Provenance", "", kv({k: x for k, x in ctx.provenance.items() if k not in ("packages", "inputs", "config")}), ""]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(s))
    return path
