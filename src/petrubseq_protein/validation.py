"""Centralized input validation.

Two layers:

* :func:`validate_input_files` runs before anything is read: every configured
  path exists and has the shape its format needs (a 10x directory has its three
  files, ...), and the modality declarations are possible (RNA present; a
  protein representation exists when protein is required; raw ADT counts exist
  when they are required).
* :func:`validate_canonical` runs on the :class:`~petrubseq_protein.io.adapters.CanonicalInput`
  and returns ``(errors, warnings)``: duplicate feature IDs, negative or
  non-integer values in matrices declared as counts, guide names that do not
  match the guide-count features or the target regex, malformed guide-target
  tables, metadata that shares no cell with the RNA matrix, ...

Errors are raised when the processed object would be ambiguous or wrong;
everything else is a warning that ends up in the report.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import Config, ConfigError
from .io import readers as pio
from .io.adapters import CanonicalInput, MATRIX_MODALITIES
from .reporting import provenance

MTX_FILES = ("barcodes.tsv", "features.tsv", "matrix.mtx")


class InputValidationError(FileNotFoundError):
    pass


def _check_path(fmt: str, p: Path, what: str, problems: List[str]) -> None:
    if not p.exists():
        problems.append(f"{what}: not found: {p}")
        return
    if fmt == "mtx":
        if not p.is_dir():
            problems.append(f"{what}: format mtx expects a 10x directory, got a file: {p}")
            return
        missing = [f for f in MTX_FILES if not (p / f).is_file() and not (p / f"{f}.gz").is_file()]
        if missing:
            problems.append(f"{what}: {p} is not a 10x MTX directory; missing {missing} (with or without .gz)")
    elif fmt in ("10x_h5", "h5ad", "dense_csv") and not p.is_file():
        problems.append(f"{what}: format {fmt} expects a file, got a directory: {p}")


def validate_input_files(cfg: Config) -> Dict[str, Any]:
    """Check that every configured input exists and the modality declarations are satisfiable.

    Returns the file records (path, size, mtime) keyed by role, for provenance.
    """
    problems: List[str] = []
    rec: Dict[str, Any] = {}
    mp = cfg.inputs.multiplexed
    mp_mods = mp.modalities() if mp.is_set() else []
    if mp.is_set():
        for i, rel in enumerate(mp.paths()):
            p = cfg.resolve(rel)
            key = f"multiplexed:{mp.lanes[i]['id']}" if mp.lanes else "multiplexed"
            rec[key] = provenance.file_record(p)
            _check_path(mp.format, p, "inputs.multiplexed", problems)
    available = set(mp_mods)
    for name in MATRIX_MODALITIES:
        m = getattr(cfg.inputs, name)
        paths = [cfg.resolve(p) for p in m.paths()]
        if paths:
            available.add(name)
        for i, p in enumerate(paths):
            key = f"{name}:{m.lanes[i]['id']}" if m.lanes else f"{name}:{p.name}"
            rec[key] = provenance.file_record(p)
            _check_path(m.format, p, f"inputs.{name}", problems)
    # modality availability (format-independent rules) ---------------------
    if "rna" not in available:
        problems.append("inputs.rna: required but no RNA source is configured")
    has_protein = "protein" in available or "protein_counts" in available
    if cfg.inputs.protein.required and not has_protein:
        problems.append("protein is required (inputs.protein.required) but neither normalized protein values nor raw ADT counts are configured; set inputs.protein / inputs.protein_counts / inputs.multiplexed.feature_types.protein, or inputs.protein.required: false")
    if cfg.inputs.protein_counts.required and "protein_counts" not in available:
        problems.append("raw ADT counts are required (inputs.protein_counts.required) but not configured")
    if cfg.inputs.guide_counts.required and "guide_counts" not in available:
        problems.append("guide counts are required (inputs.guide_counts.required) but not configured")
    for mod in cfg.alignment.required:
        if mod == "protein":
            if not has_protein:
                problems.append("alignment.required lists 'protein' but no protein input is configured")
        elif mod in ("protein_counts", "guide_counts") and mod not in available:
            problems.append(f"alignment.required lists '{mod}' but it is not configured")
        elif mod == "metadata" and not cfg.inputs.metadata.file and not cfg.inputs.lane_metadata.file:
            problems.append("alignment.required lists 'metadata' but inputs.metadata.file is not set")
        elif mod == "guide_assignments" and not cfg.inputs.guide_assignments.file and "guide_assignments" not in mp_mods:
            problems.append("alignment.required lists 'guide_assignments' but none is configured")
    # tables -----------------------------------------------------------------
    for name in ("metadata", "lane_metadata", "guide_assignments", "embedding"):
        t = getattr(cfg.inputs, name)
        p = cfg.resolve(t.file)
        if p is None:
            if getattr(t, "required", False):
                problems.append(f"inputs.{name}: required but no file given")
            continue
        rec[name] = provenance.file_record(p)
        if not p.exists():
            problems.append(f"inputs.{name}: file not found: {p}")
    for key, rel in (("perturbation.guide_target_table", cfg.perturbation.guide_target_table), ("protein.feature_table", cfg.protein.feature_table)):
        p = cfg.resolve(rel)
        if p is not None:
            rec[key.split(".")[-1]] = provenance.file_record(p)
            if not p.exists():
                problems.append(f"{key}: file not found: {p}")
    asg = cfg.perturbation.assignment
    has_provided = bool(cfg.inputs.guide_assignments.file) or "guide_assignments" in mp_mods or bool(cfg.columns.guide) or bool(cfg.columns.perturbation)
    if asg.source == "guide_counts" and "guide_counts" not in available:
        problems.append("perturbation.assignment.source is 'guide_counts' but no guide-count matrix is configured")
    if asg.source == "provided" and not has_provided:
        problems.append("perturbation.assignment.source is 'provided' but no guide assignment file / metadata guide column is configured")
    if problems:
        raise InputValidationError("Input validation failed:\n  - " + "\n  - ".join(problems))
    return rec


# ---------------------------------------------------------------------------
# canonical-input checks
# ---------------------------------------------------------------------------


def _count_like(X: sp.spmatrix, max_check: int = 5_000_000) -> Tuple[bool, bool]:
    """(any negative, any non-integer) on (a prefix of) the stored values."""
    data = X.data[:max_check] if X.data.size > max_check else X.data
    if data.size == 0:
        return False, False
    neg = bool(data.min() < 0)
    nonint = bool(np.any(np.abs(data - np.round(data)) > 1e-6))
    return neg, nonint


def _dup_examples(idx: pd.Index, n: int = 3) -> List[str]:
    return idx[idx.duplicated()].unique()[:n].astype(str).tolist()


def validate_canonical(ci: CanonicalInput, cfg: Config) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []
    for name, m in ci.matrices().items():
        feats = pd.Index(m.features.astype(str))
        if feats.has_duplicates:
            hint = " (10x: use var_names: id, which is unique)" if m.format in ("mtx", "10x_h5") else ""
            errors.append(f"{name}: {int(feats.duplicated().sum())} duplicate feature IDs (e.g. {_dup_examples(feats)}){hint}")
        neg, nonint = _count_like(m.X)
        if m.state == "raw_counts" or name == "guide_counts":
            if neg:
                errors.append(f"{name}: declared raw counts but contains negative values")
            if nonint:
                errors.append(f"{name}: declared raw counts but contains non-integer values")
        elif m.state == "normalized" and not neg and not nonint and m.X.nnz:
            warnings.append(f"{name}: declared 'normalized' but every stored value is a non-negative integer (looks like raw counts).")
        if m.X.nnz == 0:
            errors.append(f"{name}: the matrix is empty (no non-zero values)")
    # metadata ----------------------------------------------------------------
    rna_cells = pd.Index(ci.rna.cells.astype(str))
    if ci.metadata is not None:
        if ci.metadata.index.has_duplicates:
            errors.append(f"metadata: duplicate cell IDs (e.g. {_dup_examples(ci.metadata.index)})")
        overlap = rna_cells.isin(ci.metadata.index).sum()
        if overlap == 0:
            errors.append(f"metadata shares no cell ID with the RNA matrix (RNA e.g. {rna_cells[:2].tolist()}, metadata e.g. {ci.metadata.index[:2].astype(str).tolist()}); check columns.cell_id and lane suffixes ('<barcode>__<lane>')")
        for key in ("condition", "guide", "moi", "rna_total_counts", "sample", "donor", "batch", "replicate", "lane", "perturbation", "target"):
            col = getattr(cfg.columns, key)
            if col and col not in ci.metadata.columns:
                errors.append(f"columns.{key} = {col!r} is not a metadata column (available: {list(ci.metadata.columns)[:12]} ...)")
        for col in cfg.columns.keep:
            if col not in ci.metadata.columns:
                errors.append(f"columns.keep: {col!r} is not a metadata column")
    elif any(getattr(cfg.columns, k) for k in ("condition", "guide", "moi", "rna_total_counts", "sample", "donor", "batch")):
        errors.append("columns.* map metadata columns but no metadata table is configured")
    # guides ----------------------------------------------------------------
    regex = cfg.perturbation.guide_target_regex
    try:
        re.compile(regex)
    except re.error as exc:
        errors.append(f"perturbation.guide_target_regex is not a valid regex: {exc}")
        regex = None
    gt_table: Optional[Dict[str, str]] = None
    if cfg.perturbation.guide_target_table:
        p = cfg.resolve(cfg.perturbation.guide_target_table)
        try:
            gt = pio.read_table(p)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"perturbation.guide_target_table: cannot read {p}: {exc}")
            gt = None
        if gt is not None:
            if gt.shape[1] < 2:
                errors.append("perturbation.guide_target_table needs at least two columns (guide, target)")
            else:
                g = gt.iloc[:, 0].astype(str)
                if g.duplicated().any():
                    errors.append(f"perturbation.guide_target_table: duplicate guide IDs {g[g.duplicated()].unique()[:3].tolist()}")
                if gt.iloc[:, 1].isna().any() or (gt.iloc[:, 1].astype(str).str.strip() == "").any():
                    errors.append("perturbation.guide_target_table: empty target values")
                gt_table = dict(zip(g, gt.iloc[:, 1].astype(str)))
    guide_ids: Optional[pd.Index] = None
    if ci.guide_counts is not None:
        guide_ids = pd.Index(ci.guide_counts.features.astype(str))
        if regex:
            unparsed = [g for g in guide_ids if not (gt_table and g in gt_table) and not (re.match(regex, g) and re.match(regex, g).groupdict().get("target"))]
            if unparsed:
                warnings.append(f"{len(unparsed)} of {len(guide_ids)} guide-count features do not match perturbation.guide_target_regex (e.g. {unparsed[:3]}); their target defaults to the guide ID.")
        if ci.guide_counts.shape[1] < 2 and cfg.perturbation.assignment.source != "provided":
            warnings.append("guide-count matrix has a single guide; the dominant-guide rule degenerates to a UMI threshold.")
    if ci.guide_assignments is not None:
        if ci.guide_assignments.index.has_duplicates:
            errors.append(f"guide_assignments: duplicate cell IDs (e.g. {_dup_examples(ci.guide_assignments.index)})")
        used = pd.Index(sorted({g for l in ci.guide_assignments for g in l}))
        if guide_ids is not None and len(used):
            missing = used[~used.isin(guide_ids)]
            if len(missing):
                warnings.append(f"{len(missing)} of {len(used)} provided guide IDs are absent from the guide-count features (e.g. {missing[:3].tolist()}); provided-vs-count comparison covers the shared guides only.")
        if regex and len(used):
            unparsed = [g for g in used if not (gt_table and g in gt_table) and not (re.match(regex, g) and re.match(regex, g).groupdict().get("target"))]
            if unparsed:
                warnings.append(f"{len(unparsed)} provided guide IDs do not match perturbation.guide_target_regex (e.g. {unparsed[:3]}); their target defaults to the guide ID.")
    asg = cfg.perturbation.assignment
    if asg.source == "auto" and ci.guide_counts is not None and ci.guide_assignments is not None:
        errors.append("both a guide-count matrix and provided guide assignments are available; set perturbation.assignment.source to 'provided' (use the supplied assignment, counts feed diagnostics) or 'guide_counts' (call guides from counts, compare with the provided assignment)")
    # protein -------------------------------------------------------------------
    if ci.protein_normalized is not None and ci.protein_counts is not None:
        pass  # name reconciliation and formula verification live in preprocessing.normalize
    if cfg.protein.feature_table:
        p = cfg.resolve(cfg.protein.feature_table)
        try:
            ft = pio.read_table(p)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"protein.feature_table: cannot read {p}: {exc}")
        else:
            if "feature_id" not in ft.columns:
                errors.append(f"protein.feature_table: a 'feature_id' column is required (columns: {list(ft.columns)})")
            elif ft["feature_id"].astype(str).duplicated().any():
                errors.append("protein.feature_table: duplicate feature_id rows")
    return errors, warnings
