"""Markdown mirror of the HTML report, rendered from the same ReportContext."""

from __future__ import annotations

import json
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


def _figs(ctx: ReportContext, section: str, stage: Optional[str] = None) -> str:
    recs = ctx.registry.by_section(section, stage)
    return "\n\n".join(f"![{r.title}]({r.rel(ctx.run_dir)})\n\n_{r.title}. {r.caption}_" for r in recs) if recs else "_(none)_"


def write_markdown(ctx: ReportContext, path: Path) -> Path:
    s: List[str] = [f"# {ctx.title}", "", f"Run {ctx.run_name} (dataset {ctx.dataset}), generated {ctx.generated_at}, petrubseq-pipeline-protein v{ctx.version}. Preprocessing/QC report; figures are linked by relative path (report.html embeds them).", ""]
    s += ["## Run summary", "", kv(dict(ctx.cards)), ""]
    s += ["## Warnings", ""] + ([f"- {w}" for w in ctx.warnings] or ["_(none)_"]) + [""]
    s += [f"> {n}" for n in ctx.notes] + [""]
    s += ["## 1. Inputs and modality audit", "", md_table(ctx.inputs, index=False), "", md_table(ctx.input_files, index=False), "", kv(ctx.alignment), ""]
    s += ["## 2. Cell and RNA QC", "", "### Filtering steps", "", md_table(ctx.filtering_steps, index=False), "", "### QC summary", "", md_table(ctx.qc_summary, index=False), "", "### Before filtering", "", _figs(ctx, SECTION_QC_RNA, STAGE_BEFORE), "", "### After filtering", "", _figs(ctx, SECTION_QC_RNA, STAGE_AFTER), ""]
    s += ["## 3. Protein QC", "", kv(ctx.protein_info), "", md_table(ctx.protein_features, max_rows=60), "", kv(ctx.protein_qc), "", "### Before filtering", "", _figs(ctx, SECTION_QC_PROTEIN, STAGE_BEFORE), "", "### After filtering", "", _figs(ctx, SECTION_QC_PROTEIN, STAGE_AFTER), ""]
    s += ["## 4. Perturbation QC", "", kv(ctx.perturbation_qc), "", md_table(ctx.class_by_condition), "", md_table(ctx.moi_summary), "", md_table(ctx.target_coverage.sort_values("n_cells_single_guide").head(20) if ctx.target_coverage is not None else None, max_rows=20), "", _figs(ctx, SECTION_QC_PERTURBATION), ""]
    s += ["## 5. Representations", "", "### RNA", "", kv(ctx.representations.get("rna", {})), "", _figs(ctx, SECTION_REP_RNA), "", "### Protein", "", kv(ctx.representations.get("protein", {})), "", _figs(ctx, SECTION_REP_PROTEIN), "", "### Cross-modality diagnostics", "", kv({**ctx.representations.get("multimodal", {}), **ctx.diagnostics}), "", _figs(ctx, SECTION_MULTIMODAL), ""]
    for name, df in ctx.group_summaries.items():
        s += [f"### {name} summary", "", md_table(df), ""]
    pe = ctx.extra.get("perturbation_effects")
    if pe is not None:
        from .perturbation_plots import SECTION as PE, ST_CONC, ST_LOCH, ST_PROG, ST_PROT, ST_PS
        s += ["## 6. Perturbation effects", "", "Perturbed = single-guide targeting cells of one target; controls = single-guide non-targeting cells; ambiguous / multi-guide cells excluded. Associations are correlations, not mediation or causality (docs/PERTURBATION_EFFECTS.md).", ""]
        for title, r, df, st_ in (("PS", pe.ps, pe.ps.summary if pe.ps else None, ST_PS), ("lochNESS", pe.lochness, pe.lochness.summary if pe.lochness else None, ST_LOCH),
                                  ("Gene programs and perturbation modules", pe.modules, pe.modules.perturbation_modules if pe.modules is not None and not pe.modules.empty else None, ST_PROG),
                                  ("Protein effects", pe.protein, pe.protein.table if pe.protein is not None and not pe.protein.empty else None, ST_PROT),
                                  ("RNA-protein concordance", pe.concordance, pe.concordance.ps_protein if pe.concordance else None, ST_CONC)):
            s += [f"### {title}", "", kv({k: v for k, v in (r.info if r is not None else {"status": "disabled"}).items() if not isinstance(v, (list, dict))}), "", md_table(df.round(4) if df is not None else None, index=False, max_rows=40), "", _figs(ctx, PE, st_), ""]
        if pe.concordance is not None:
            s += ["### Integrated perturbation summary", "", md_table(pe.concordance.summary.round(4), index=False, max_rows=60), ""]
    s += [f"## {7 if pe is not None else 6}. Outputs and provenance", "", kv(ctx.outputs), "", "```", ctx.adata_repr, "```", "", "| slot | shape | content | status | provenance |", "|---|---|---|---|---|"]
    s += [f"| `{k}` | {d.get('shape','')} | {d.get('content','')} | {d.get('status','')} | {d.get('provenance','')} |" for k, d in ctx.schema.items()]
    s += ["", "### Stage timings (seconds)", "", kv({k: round(v, 1) for k, v in ctx.timings.items()}), "", "### Software versions", "", kv(ctx.versions), "", "### Provenance", "", kv({k: v for k, v in ctx.provenance.items() if k not in ("packages", "inputs", "config")}), ""]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(s))
    return path
