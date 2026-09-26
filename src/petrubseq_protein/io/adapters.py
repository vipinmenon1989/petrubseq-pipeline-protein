"""Format adapters: every supported input layout is turned into one
:class:`CanonicalInput`, so nothing after the "load data" stage knows whether
the data came from dense text, 10x MTX / H5 or an h5ad file.

Supported layouts

* separate per-modality inputs (``inputs.rna`` / ``protein`` / ``protein_counts``
  / ``guide_counts``), each in ``dense_csv`` (existing v0.1 reader), ``mtx``,
  ``10x_h5`` or ``h5ad`` (``slot``/``key`` addressing);
* one ``inputs.multiplexed`` matrix (``mtx`` / ``10x_h5`` / ``h5ad``) read once
  and split by feature type, or an h5ad split by explicit ``slots``;
* several lanes of either (``lanes: [{id, path}]``): cell IDs become
  ``<barcode>__<lane id>`` and a lane table (``barcode_original``, ``lane_id``)
  is kept; a per-lane sample sheet (``inputs.lane_metadata``) is joined onto
  the cells by lane ID.

Value-state vocabulary: matrix inputs declare ``auto`` | ``raw_counts`` |
``normalized`` (``Matrix.state``); ``preprocessing.normalize.normalize_rna``
maps ``normalized`` to its internal ``log_normalized`` hint. Feature-barcode
matrices (mtx / 10x_h5) are always ``raw_counts``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import Config, ConfigError, MatrixInput
from . import readers as pio
from .readers import Matrix

logger = logging.getLogger(__name__)

LANE_SEP = "__"
MATRIX_MODALITIES = ("rna", "protein", "protein_counts", "guide_counts")


class InputError(ValueError):
    """An input cannot be interpreted unambiguously."""


@dataclass
class FeatureMatrix:
    """A raw feature-barcode matrix as read from disk (cells x features)."""

    cells: np.ndarray
    features: np.ndarray
    X: sp.csr_matrix
    #: per-feature metadata, index aligned with ``features``
    var: pd.DataFrame
    source: str


@dataclass
class CanonicalInput:
    """Everything the pipeline needs, independent of the source format."""

    rna: Matrix
    protein_counts: Optional[Matrix] = None
    protein_normalized: Optional[Matrix] = None
    guide_counts: Optional[Matrix] = None
    #: cell -> list of guide IDs (an upstream / provided assignment)
    guide_assignments: Optional[pd.Series] = None
    metadata: Optional[pd.DataFrame] = None
    embedding: Optional[pd.DataFrame] = None
    #: index = global cell ID; columns ``barcode_original``, ``lane_id`` (multi-lane only)
    lane_table: Optional[pd.DataFrame] = None
    provenance: Dict[str, Any] = field(default_factory=dict)

    def matrices(self) -> Dict[str, Matrix]:
        out = {"rna": self.rna}
        if self.protein_normalized is not None:
            out["protein"] = self.protein_normalized
        if self.protein_counts is not None:
            out["protein_counts"] = self.protein_counts
        if self.guide_counts is not None:
            out["guide_counts"] = self.guide_counts
        return out

    def id_sets(self) -> Dict[str, List[str]]:
        sets = {k: m.cells.tolist() for k, m in self.matrices().items()}
        if self.metadata is not None:
            sets["metadata"] = self.metadata.index.tolist()
        if self.guide_assignments is not None:
            sets["guide_assignments"] = self.guide_assignments.index.tolist()
        return sets

    @property
    def has_protein(self) -> bool:
        return self.protein_counts is not None or self.protein_normalized is not None


# ---------------------------------------------------------------------------
# lanes
# ---------------------------------------------------------------------------


def lane_specs(lanes: Sequence[Dict[str, str]]) -> List[Tuple[str, str]]:
    return [(str(l["id"]), str(l["path"])) for l in lanes]


def concat_lanes(parts: List[Tuple[str, FeatureMatrix]], what: str) -> Tuple[FeatureMatrix, pd.DataFrame]:
    """Concatenate per-lane matrices along cells with ``<barcode>__<lane>`` IDs.

    Returns the merged matrix and the lane table (index = global ID).
    """
    ids = [lid for lid, _ in parts]
    if len(set(ids)) != len(ids):
        raise InputError(f"{what}: lane IDs must be unique, got {ids}")
    ref = parts[0][1]
    for lid, fm in parts[1:]:
        if len(fm.features) != len(ref.features) or not np.array_equal(fm.features, ref.features):
            raise InputError(f"{what}: lane {lid!r} has a different feature set than lane {ids[0]!r} ({len(fm.features)} vs {len(ref.features)} features); all lanes must share one reference")
    cells, orig, lane, blocks = [], [], [], []
    for lid, fm in parts:
        if LANE_SEP in lid:
            raise InputError(f"{what}: lane id {lid!r} must not contain {LANE_SEP!r}")
        cells.extend(f"{c}{LANE_SEP}{lid}" for c in fm.cells)
        orig.extend(str(c) for c in fm.cells)
        lane.extend([lid] * len(fm.cells))
        blocks.append(fm.X)
    cells_arr = np.asarray(cells, dtype=object)
    if pd.Index(cells_arr).has_duplicates:
        dup = pd.Index(cells_arr)[pd.Index(cells_arr).duplicated()][:3].tolist()
        raise InputError(f"{what}: global cell IDs are not unique after lane suffixing (e.g. {dup}); barcodes repeat inside one lane")
    table = pd.DataFrame({"barcode_original": orig, "lane_id": lane}, index=pd.Index(cells_arr, name="cell"))
    merged = FeatureMatrix(cells_arr, ref.features, sp.vstack(blocks, format="csr"), ref.var.copy(), "; ".join(fm.source for _, fm in parts))
    return merged, table


def _read_feature_matrix(fmt: str, path: Path, feature_type_column: str) -> FeatureMatrix:
    if fmt == "mtx":
        from .tenx import read_mtx_dir
        return read_mtx_dir(path)
    if fmt == "10x_h5":
        from .tenx import read_10x_h5
        return read_10x_h5(path)
    if fmt == "h5ad":
        from .h5ad import read_h5ad_matrix
        return read_h5ad_matrix(path, slot="X", key=None, feature_type_column=feature_type_column)
    raise ConfigError(f"unsupported feature-matrix format {fmt!r}")


def read_lanes_or_path(fmt: str, path: Optional[str], lanes: Sequence[Dict[str, str]], resolve: Callable[[str], Path], what: str, reader: Callable[[Path], FeatureMatrix]) -> Tuple[FeatureMatrix, Optional[pd.DataFrame]]:
    """Read a single path (no suffixing) or several lanes (suffixed IDs)."""
    if lanes:
        parts = []
        for lid, lp in lane_specs(lanes):
            fm = reader(resolve(lp))
            logger.info("%s lane %s: %d cells x %d features (%s)", what, lid, len(fm.cells), len(fm.features), fm.source)
            parts.append((lid, fm))
        return concat_lanes(parts, what)
    assert path is not None
    return reader(resolve(path)), None


# ---------------------------------------------------------------------------
# feature-type splitting
# ---------------------------------------------------------------------------


def split_feature_types(fm: FeatureMatrix, type_column: str, wanted: Dict[str, List[str]], what: str) -> Dict[str, FeatureMatrix]:
    """Split a combined matrix into modalities by ``var[type_column]``.

    ``wanted`` maps modality -> accepted type labels (empty list = modality not
    expected). Every non-empty modality must match at least one feature, and
    every feature type present in the file must be claimed by some modality
    (otherwise the data would be silently dropped).
    """
    if type_column not in fm.var.columns:
        raise InputError(f"{what}: no feature-type column {type_column!r} in the feature table (columns: {list(fm.var.columns)}); cannot split modalities")
    types = fm.var[type_column].astype(str).to_numpy()
    present = sorted(set(types))
    out: Dict[str, FeatureMatrix] = {}
    claimed: set = set()
    for mod, labels in wanted.items():
        if not labels:
            continue
        mask = np.isin(types, labels)
        if not mask.any():
            raise InputError(f"{what}: no features of type {labels} for modality '{mod}'; types present: {present}")
        claimed.update(l for l in labels if l in present)
        out[mod] = FeatureMatrix(fm.cells, fm.features[mask], sp.csr_matrix(fm.X[:, mask]), fm.var.loc[mask].copy(), fm.source)
    unclaimed = [t for t in present if t not in claimed]
    if unclaimed:
        raise InputError(f"{what}: feature type(s) {unclaimed} are present in the file but assigned to no modality in feature_types; list them under a modality or remove them upstream")
    return out


def _to_matrix(fm: FeatureMatrix, state: str, fmt: str, var_names: str = "id") -> Matrix:
    if fmt in ("mtx", "10x_h5") and var_names != "id":
        from .tenx import rename_features
        fm = rename_features(fm, var_names)
    var = fm.var.copy()
    var.index = pd.Index(fm.features, name="feature")
    return Matrix(fm.cells, fm.features, sp.csr_matrix(fm.X), [fm.source], [len(fm.cells)], var, state, fmt)


# ---------------------------------------------------------------------------
# per-modality readers
# ---------------------------------------------------------------------------


def read_matrix_input(cfg: Config, name: str, m: MatrixInput, chunk_rows: int, n_jobs: int, max_cells: Optional[int] = None) -> Tuple[Matrix, Optional[pd.DataFrame]]:
    """Read one separately configured modality. Returns (matrix, lane table)."""
    if m.format == "dense_csv":
        paths = [cfg.resolve(p) for p in m.paths()]
        mat = pio.read_dense_matrix(paths, m.orientation, m.sep, chunk_rows, n_jobs, max_cells=max_cells)
        mat.state, mat.format = m.state, "dense_csv"
        return mat, None
    if m.format in ("mtx", "10x_h5"):
        reader = lambda p: _read_feature_matrix(m.format, p, m.feature_type_column)  # noqa: E731
        fm, lane_table = read_lanes_or_path(m.format, m.file, m.lanes, cfg.resolve, f"inputs.{name}", reader)
        state = "raw_counts" if m.state == "auto" else m.state
    else:  # h5ad
        from .h5ad import read_h5ad_matrix
        reader = lambda p: read_h5ad_matrix(p, slot=m.slot, key=m.key, feature_type_column=m.feature_type_column)  # noqa: E731
        fm, lane_table = read_lanes_or_path("h5ad", m.file, m.lanes, cfg.resolve, f"inputs.{name}", reader)
        state = m.state
    if m.feature_types:
        fm = _select_types(fm, m.feature_type_column, m.feature_types, f"inputs.{name}")
    mat = _to_matrix(fm, state, m.format, m.var_names)
    if max_cells is not None and mat.shape[0] > max_cells:
        mat = mat.subset_cells(mat.cells[:max_cells].tolist())
    return mat, lane_table


def _select_types(fm: FeatureMatrix, type_column: str, labels: Sequence[str], what: str) -> FeatureMatrix:
    if type_column not in fm.var.columns:
        raise InputError(f"{what}: feature_types given but the feature table has no column {type_column!r}")
    mask = fm.var[type_column].astype(str).isin(list(labels)).to_numpy()
    if not mask.any():
        raise InputError(f"{what}: no features of type {list(labels)}; present: {sorted(fm.var[type_column].astype(str).unique())}")
    return FeatureMatrix(fm.cells, fm.features[mask], sp.csr_matrix(fm.X[:, mask]), fm.var.loc[mask].copy(), fm.source)


def read_multiplexed(cfg: Config) -> Tuple[Dict[str, Any], Optional[pd.DataFrame], Dict[str, Any]]:
    """Read ``inputs.multiplexed`` once and split it.

    Returns (parts, lane table, provenance) where ``parts`` maps ``rna`` /
    ``protein_counts`` / ``protein`` / ``guide_counts`` to :class:`Matrix` and
    ``guide_assignments`` to a Series of guide lists.
    """
    mp = cfg.inputs.multiplexed
    what = "inputs.multiplexed"
    prov: Dict[str, Any] = {"format": mp.format, "paths": [str(cfg.resolve(p)) for p in mp.paths()]}
    parts: Dict[str, Any] = {}
    if mp.format == "h5ad" and mp.slots.declared():
        from .h5ad import read_h5ad_slots
        lane_table = None
        if mp.lanes:
            per_lane = [(lid, read_h5ad_slots(cfg.resolve(lp), mp.slots, mp.feature_type_column)) for lid, lp in lane_specs(mp.lanes)]
            merged: Dict[str, Any] = {}
            for mod in per_lane[0][1]:
                objs = [(lid, d[mod]) for lid, d in per_lane]
                if mod == "guide_assignments":
                    merged[mod] = pd.concat([s.rename(index=lambda c, lid=lid: f"{c}{LANE_SEP}{lid}") for lid, s in objs])
                else:
                    fm, lt = concat_lanes(objs, f"{what}[{mod}]")
                    merged[mod] = fm
                    lane_table = lt if lane_table is None else lane_table
            raw = merged
        else:
            raw = read_h5ad_slots(cfg.resolve(mp.path), mp.slots, mp.feature_type_column)
        prov["split"] = "slots: " + ", ".join(f"{k}={getattr(mp.slots, k).slot}[{getattr(mp.slots, k).key or ''}]" for k in mp.slots.declared())
        for mod, obj in raw.items():
            if mod == "guide_assignments":
                parts[mod] = obj
            else:
                state = mp.rna_state if mod == "rna" else ("raw_counts" if mod in ("protein_counts", "guide_counts") else "normalized")
                parts[mod] = _to_matrix(obj, state, "h5ad")
        return parts, lane_table, prov
    reader = lambda p: _read_feature_matrix(mp.format, p, mp.feature_type_column)  # noqa: E731
    fm, lane_table = read_lanes_or_path(mp.format, mp.path, mp.lanes, cfg.resolve, what, reader)
    wanted = {"rna": list(mp.feature_types.rna), "protein_counts": list(mp.feature_types.protein), "guide_counts": list(mp.feature_types.guide)}
    split = split_feature_types(fm, mp.feature_type_column, wanted, what)
    prov["split"] = {k: int(len(v.features)) for k, v in split.items()}
    prov["feature_types"] = {k: v for k, v in wanted.items() if v}
    rna_state = "raw_counts" if mp.format in ("mtx", "10x_h5") else mp.rna_state
    for mod, sub in split.items():
        # var_names (id | name) applies to RNA only; antibodies and guides keep their unique feature IDs
        parts[mod] = _to_matrix(sub, rna_state if mod == "rna" else "raw_counts", mp.format, mp.var_names if mod == "rna" else "id")
    return parts, lane_table, prov


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------


def load_inputs(cfg: Config, max_cells: Optional[int] = None) -> CanonicalInput:
    """Read every configured input and return the canonical model."""
    prov: Dict[str, Any] = {"modalities": {}, "multiplexed": None, "lanes": None, "cell_id_transformation": "none"}
    mats: Dict[str, Matrix] = {}
    guides: Optional[pd.Series] = None
    lane_tables: Dict[str, pd.DataFrame] = {}
    mp = cfg.inputs.multiplexed
    if mp.is_set():
        parts, lt, mprov = read_multiplexed(cfg)
        prov["multiplexed"] = mprov
        if lt is not None:
            lane_tables["multiplexed"] = lt
        for mod, obj in parts.items():
            if mod == "guide_assignments":
                guides = obj
                prov["modalities"]["guide_assignments"] = {"format": mp.format, "source": "multiplexed obs column", "state": "assignments"}
            else:
                mats[mod] = obj
                prov["modalities"][mod] = {"format": mp.format, "source": "multiplexed", "state": obj.state, "n_features": int(obj.shape[1])}
                logger.info("%s (multiplexed %s): %d cells x %d features", mod, mp.format, *obj.shape)
    for name in MATRIX_MODALITIES:
        m: MatrixInput = getattr(cfg.inputs, name)
        if not m.is_set():
            continue
        if name in mats:
            raise ConfigError(f"inputs.{name} is provided both by inputs.multiplexed and separately")
        is_rna = name == "rna"
        mat, lt = read_matrix_input(cfg, name, m, cfg.compute.chunk_rows if is_rna else 64, cfg.compute.n_jobs if is_rna else 1, max_cells if is_rna else None)
        mats[name] = mat
        if lt is not None:
            lane_tables[name] = lt
        prov["modalities"][name] = {"format": m.format, "source": ", ".join(Path(p).name for p in m.paths()), "state": mat.state, "n_features": int(mat.shape[1])}
        logger.info("%s (%s): %d cells x %d features from %d source(s)", name, m.format, mat.shape[0], mat.shape[1], len(mat.sources))
    if "rna" not in mats:
        raise ConfigError("no RNA matrix was loaded")
    # lane bookkeeping ------------------------------------------------------
    lane_table = None
    if lane_tables:
        lane_table = pd.concat(lane_tables.values())
        lane_table = lane_table[~lane_table.index.duplicated()]
        ids = sorted(lane_table["lane_id"].unique())
        prov["lanes"] = {"ids": ids, "n_lanes": len(ids), "cells_per_lane": {k: int(v) for k, v in lane_table["lane_id"].value_counts().items()}}
        prov["cell_id_transformation"] = f"<barcode>{LANE_SEP}<lane_id> (original barcode kept in obs['barcode_original'])"
    # tables ----------------------------------------------------------------
    meta = None
    if cfg.inputs.metadata.file:
        t = cfg.inputs.metadata
        if t.format == "h5ad":
            from .h5ad import read_h5ad_obs
            meta = read_h5ad_obs(cfg.resolve(t.file), cfg.columns.cell_id)
        else:
            meta = pio.read_table(cfg.resolve(t.file), t.format, t.sep, index_col=cfg.columns.cell_id)
        meta.index = meta.index.astype(str)
        logger.info("metadata: %d cells x %d columns", *meta.shape)
        prov["modalities"]["metadata"] = {"format": t.format, "source": Path(t.file).name}
    if cfg.inputs.lane_metadata.file:
        lm_cfg = cfg.inputs.lane_metadata
        if lane_table is None:
            raise InputError("inputs.lane_metadata given but the inputs are not multi-lane")
        lm = pio.read_table(cfg.resolve(lm_cfg.file), lm_cfg.format, lm_cfg.sep)
        if lm_cfg.lane_column not in lm.columns:
            raise InputError(f"inputs.lane_metadata: column {lm_cfg.lane_column!r} not found (columns: {list(lm.columns)})")
        lm = lm.set_index(lm.pop(lm_cfg.lane_column).astype(str))
        if lm.index.has_duplicates:
            raise InputError(f"inputs.lane_metadata: duplicate lane rows {lm.index[lm.index.duplicated()].unique().tolist()}; the join would be ambiguous")
        missing = [l for l in prov["lanes"]["ids"] if l not in lm.index]
        if missing:
            raise InputError(f"inputs.lane_metadata: no row for lane(s) {missing}")
        per_cell = lane_table[["lane_id"]].join(lm, on="lane_id").drop(columns="lane_id")
        clash = [c for c in per_cell.columns if meta is not None and c in meta.columns]
        if clash:
            raise InputError(f"inputs.lane_metadata columns {clash} also exist in inputs.metadata; rename one side")
        meta = per_cell if meta is None else meta.join(per_cell, how="left")
        prov["modalities"]["lane_metadata"] = {"format": lm_cfg.format, "source": Path(lm_cfg.file).name, "columns": list(lm.columns)}
    if cfg.inputs.guide_assignments.file:
        g = cfg.inputs.guide_assignments
        if guides is not None:
            raise ConfigError("guide assignments are given both in inputs.multiplexed.slots and inputs.guide_assignments")
        if g.format == "h5ad":
            from .h5ad import guide_lists_from_obs, read_h5ad_obs
            guides = guide_lists_from_obs(read_h5ad_obs(cfg.resolve(g.file)), g.guides_column, g.list_separator, str(cfg.resolve(g.file)))
        else:
            guides = pio.read_guide_assignments(cfg.resolve(g.file), g.cell_column, g.guides_column, g.list_separator, g.format, g.sep)
        prov["modalities"]["guide_assignments"] = {"format": g.format, "source": Path(g.file).name + (f" obs[{g.guides_column!r}]" if g.format == "h5ad" else ""), "state": "assignments"}
    if guides is not None:
        logger.info("guide assignments: %d cells", len(guides))
    emb = None
    if cfg.inputs.embedding.file:
        e = cfg.inputs.embedding
        emb = pio.read_embedding(cfg.resolve(e.file), e.columns, cfg.columns.cell_id, e.format, e.sep)
        prov["modalities"]["embedding"] = {"format": e.format, "source": Path(e.file).name}
    prov["guide_counts_available"] = "guide_counts" in mats
    prov["provided_assignments_available"] = guides is not None
    return CanonicalInput(
        rna=mats["rna"],
        protein_counts=mats.get("protein_counts"),
        protein_normalized=mats.get("protein"),
        guide_counts=mats.get("guide_counts"),
        guide_assignments=guides,
        metadata=meta,
        embedding=emb,
        lane_table=lane_table,
        provenance=prov,
    )
