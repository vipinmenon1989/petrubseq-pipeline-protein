"""One report context consumed by both the HTML and the Markdown renderer."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from .figures import FigureRegistry


@dataclass
class ReportContext:
    run_dir: Path
    run_name: str
    dataset: str
    title: str
    version: str
    generated_at: str
    cards: List[Tuple[str, str]]
    warnings: List[str]
    notes: List[str]
    inputs: pd.DataFrame                 # modality audit rows
    input_files: pd.DataFrame            # file records
    alignment: Dict[str, Any]
    filtering_steps: pd.DataFrame
    qc_summary: Optional[pd.DataFrame]
    rna_qc: Dict[str, Any]
    protein_features: Optional[pd.DataFrame]
    protein_qc: Dict[str, Any]
    protein_info: Dict[str, Any]
    perturbation_qc: Dict[str, Any]
    class_by_condition: Optional[pd.DataFrame]
    target_coverage: Optional[pd.DataFrame]
    guide_coverage: Optional[pd.DataFrame]
    moi_summary: Optional[pd.DataFrame]
    group_summaries: Dict[str, pd.DataFrame]
    representations: Dict[str, Any]
    diagnostics: Dict[str, Any]
    registry: FigureRegistry
    tables: Dict[str, Path]
    outputs: Dict[str, str]
    timings: Dict[str, float]
    config_yaml: str
    versions: Dict[str, str]
    provenance: Dict[str, Any]
    schema: Dict[str, Dict[str, str]]
    adata_repr: str
    filter_enabled: bool = True
    extra: Dict[str, Any] = field(default_factory=dict)
