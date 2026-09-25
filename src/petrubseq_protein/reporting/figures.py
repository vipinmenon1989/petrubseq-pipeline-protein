"""Figure registry: every figure is saved once under ``figures/<section>[/<stage>]/``
and described by a record (path, section, stage, title, caption) that drives the
report and ``tables/figure_manifest.csv``. Modelled on the Wei Li pipeline's
``FigureRegistry``; report code never touches filesystem paths directly."""

from __future__ import annotations

import base64
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

SECTION_QC_RNA = "qc_rna"
SECTION_QC_PROTEIN = "qc_protein"
SECTION_QC_PERTURBATION = "qc_perturbation"
SECTION_REP_RNA = "representations_rna"
SECTION_REP_PROTEIN = "representations_protein"
SECTION_MULTIMODAL = "multimodal"
STAGE_BEFORE = "before_filtering"
STAGE_AFTER = "after_filtering"
STAGE_POST = "post_filter"
STAGE_REP = "representation"
STAGED = (STAGE_BEFORE, STAGE_AFTER)


def _slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


@dataclass
class FigureRecord:
    path: Path
    name: str
    section: str
    stage: str
    title: str
    caption: str = ""
    in_report: bool = True

    def data_uri(self) -> str:
        mime = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}.get(self.path.suffix.lstrip("."), "image/png")
        return f"data:{mime};base64," + base64.b64encode(self.path.read_bytes()).decode("ascii")

    def rel(self, run_dir: Path) -> str:
        return self.path.relative_to(run_dir).as_posix()


@dataclass
class FigureRegistry:
    run_dir: Path
    fmt: str = "png"
    dpi: int = 120
    records: List[FigureRecord] = field(default_factory=list)

    @property
    def figdir(self) -> Path:
        return Path(self.run_dir) / "figures"

    def path_for(self, name: str, section: str, stage: str) -> Path:
        d = self.figdir / section / (stage if stage in STAGED else "")
        return d / f"{_slug(name)}.{self.fmt}"

    def save(self, fig, name: str, section: str, stage: str, title: str, caption: str = "", in_report: bool = True) -> Optional[FigureRecord]:
        """Write ``fig`` (or skip when None) and register it. Always closes the figure."""
        if fig is None:
            return None
        path = self.path_for(name, section, stage)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        rec = FigureRecord(path, _slug(name), section, stage, title, caption, in_report)
        self.records.append(rec)
        return rec

    def by_section(self, section: str, stage: Optional[str] = None, only_in_report: bool = True) -> List[FigureRecord]:
        return [r for r in self.records if r.section == section and (stage is None or r.stage == stage) and (r.in_report or not only_in_report)]

    def manifest(self) -> pd.DataFrame:
        return pd.DataFrame([{"section": r.section, "stage": r.stage, "name": r.name, "title": r.title, "caption": r.caption, "in_report": r.in_report, "path": r.rel(Path(self.run_dir))} for r in self.records])
