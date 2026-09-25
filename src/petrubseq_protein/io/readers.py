"""Input readers: dense text matrices (streamed to sparse), SCP-style tables,
guide-assignment lists and precomputed embeddings.

Everything here is memory-aware: a features-by-cells CSV is parsed in row
chunks, each chunk converted to a sparse block, and the blocks stacked, so the
dense matrix never exists in memory. Multi-file matrices (cells partitioned
across files) are read file by file, optionally in parallel processes.
"""

from __future__ import annotations

import gzip
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger("petrubseq_protein")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text(path: Path):
    return gzip.open(path, "rt") if is_gzip(path) else open(path, "rt")


def guess_sep(path: Path, sep: str = "auto") -> str:
    if sep != "auto":
        return {"\\t": "\t", "tab": "\t", "comma": ","}.get(sep, sep)
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".csv"):
        return ","
    if name.endswith((".tsv", ".tab")):
        return "\t"
    # Sniff the header line.
    with open_text(path) as fh:
        head = fh.readline()
    return "\t" if head.count("\t") > head.count(",") else ","


def read_header(path: Path, sep: str = "auto") -> List[str]:
    """Return the header tokens of a text table without reading the body."""
    s = guess_sep(path, sep)
    with open_text(path) as fh:
        line = fh.readline().rstrip("\r\n")
    return [t.strip('"') for t in line.split(s)]


def count_lines(path: Path) -> int:
    """Count newline-terminated lines by streaming (works for gzip)."""
    n = 0
    opener = gzip.open if is_gzip(path) else open
    with opener(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            n += block.count(b"\n")
    return n


# ---------------------------------------------------------------------------
# Dense text matrix -> sparse
# ---------------------------------------------------------------------------


@dataclass
class Matrix:
    """A cells-by-features sparse matrix with labelled axes."""

    cells: np.ndarray
    features: np.ndarray
    X: sp.csr_matrix
    sources: List[str] = field(default_factory=list)
    #: Per-source (file) cell counts, in reading order.
    n_cells_per_source: List[int] = field(default_factory=list)
    #: Optional per-feature metadata (index = ``features``; e.g. 10x ``feature_id``,
    #: ``feature_name``, ``feature_type``, ``genome``). ``None`` for dense text inputs.
    features_meta: Optional[pd.DataFrame] = None
    #: Declared/canonical value state (``auto`` | ``raw_counts`` | ``normalized``).
    state: str = "auto"
    #: Input format the matrix came from (``dense_csv`` | ``mtx`` | ``10x_h5`` | ``h5ad``).
    format: str = "dense_csv"

    @property
    def shape(self) -> Tuple[int, int]:
        return self.X.shape

    def subset_cells(self, cells: Sequence[str]) -> "Matrix":
        idx = pd.Index(self.cells).get_indexer(cells)
        if (idx < 0).any():
            missing = [c for c, i in zip(cells, idx) if i < 0][:5]
            raise KeyError(f"Cells not present in matrix: {missing} ...")
        return Matrix(np.asarray(cells, dtype=object), self.features, self.X[idx], list(self.sources), [], self.features_meta, self.state, self.format)

    def subset_features(self, mask: np.ndarray) -> "Matrix":
        mask = np.asarray(mask, dtype=bool)
        fm = self.features_meta.loc[mask] if self.features_meta is not None else None
        return Matrix(self.cells, self.features[mask], sp.csr_matrix(self.X[:, mask]), list(self.sources), list(self.n_cells_per_source), fm, self.state, self.format)


def _read_dense_file(
    path: str,
    sep: str,
    orientation: str,
    chunk_rows: int,
    usecols: Optional[Sequence[str]],
    max_cells: Optional[int],
    dtype=np.float32,
) -> Tuple[np.ndarray, np.ndarray, sp.csr_matrix]:
    """Read one dense text matrix and return (row_names, col_names, csr rows x cols).

    ``usecols``/``max_cells`` restrict *cells* (columns for features_x_cells).
    """
    p = Path(path)
    s = guess_sep(p, sep)
    header = read_header(p, s)
    col_names = header[1:]
    if orientation == "features_x_cells":
        if usecols is not None:
            keep = [c for c in col_names if c in set(usecols)]
        elif max_cells is not None:
            keep = col_names[:max_cells]
        else:
            keep = None
        cols_arg = None if keep is None else [header[0]] + keep
    else:
        cols_arg = None  # row subsetting is handled after reading

    blocks: List[sp.csr_matrix] = []
    row_names: List[str] = []
    value_cols = cols_arg[1:] if cols_arg is not None else col_names
    dtypes = {c: dtype for c in value_cols}
    dtypes[header[0]] = str
    reader = pd.read_csv(
        p,
        sep=s,
        index_col=0,
        usecols=cols_arg,
        chunksize=chunk_rows,
        dtype=dtypes,
        engine="c",
        low_memory=False,
        compression="infer",
        memory_map=False,
    )
    n_rows_read = 0
    for chunk in reader:
        if orientation == "cells_x_features" and max_cells is not None and n_rows_read >= max_cells:
            break
        if orientation == "cells_x_features" and usecols is not None:
            chunk = chunk.loc[chunk.index.isin(set(usecols))]
        if orientation == "cells_x_features" and max_cells is not None:
            chunk = chunk.iloc[: max_cells - n_rows_read]
        n_rows_read += chunk.shape[0]
        row_names.extend(chunk.index.astype(str).tolist())
        blocks.append(sp.csr_matrix(np.ascontiguousarray(chunk.to_numpy(dtype=dtype))))
    if not blocks:
        raise ValueError(f"No rows read from {p}")
    mat = sp.vstack(blocks, format="csr")
    if cols_arg is not None:
        col_names = list(keep)  # type: ignore[arg-type]
    return np.asarray(row_names, dtype=object), np.asarray(col_names, dtype=object), mat


def _read_dense_file_star(args):
    return _read_dense_file(*args)


def read_dense_matrix(
    paths: Sequence[str | Path],
    orientation: str = "features_x_cells",
    sep: str = "auto",
    chunk_rows: int = 500,
    n_jobs: int = 1,
    cells: Optional[Sequence[str]] = None,
    max_cells: Optional[int] = None,
) -> Matrix:
    """Read one or more dense text matrices into a cells-by-features CSR matrix.

    Parameters
    ----------
    paths
        Files that partition the cells. Features must be identical (same order)
        across files.
    orientation
        Layout of each file.
    cells
        Optional subset of cell IDs to load (order of the file is kept).
    max_cells
        Optional cap on cells per file (smoke tests).
    """
    paths = [str(p) for p in paths]
    args = [(p, sep, orientation, chunk_rows, cells, max_cells) for p in paths]
    if n_jobs > 1 and len(paths) > 1:
        with ProcessPoolExecutor(max_workers=min(n_jobs, len(paths))) as ex:
            results = list(ex.map(_read_dense_file_star, args))
    else:
        results = [_read_dense_file_star(a) for a in args]

    mats: List[sp.csr_matrix] = []
    all_cells: List[str] = []
    features: Optional[np.ndarray] = None
    per_source: List[int] = []
    for p, (rows, cols, mat) in zip(paths, results):
        if orientation == "features_x_cells":
            feats, cell_ids, block = rows, cols, mat.T.tocsr()
        else:
            feats, cell_ids, block = cols, rows, mat
        if features is None:
            features = feats
        elif len(feats) != len(features) or not np.array_equal(feats, features):
            raise ValueError(
                f"Feature axis of {p} differs from the first file "
                f"({len(feats)} vs {len(features)} features, or different order)."
            )
        logger.info("read %s: %d cells x %d features (nnz=%d)", Path(p).name, block.shape[0], block.shape[1], block.nnz)
        mats.append(block)
        all_cells.extend(cell_ids.tolist())
        per_source.append(block.shape[0])
    X = sp.vstack(mats, format="csr") if len(mats) > 1 else mats[0]
    X.sort_indices()
    return Matrix(np.asarray(all_cells, dtype=object), np.asarray(features, dtype=object), X, paths, per_source)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def _has_scp_type_row(path: Path, sep: str) -> bool:
    with open_text(path) as fh:
        fh.readline()
        second = fh.readline()
    return second.split(sep)[0].strip('"').upper() == "TYPE"


def read_table(path: str | Path, fmt: str = "auto", sep: str = "auto", index_col: Optional[str] = None) -> pd.DataFrame:
    """Read a metadata / cluster table. Handles the SCP ``TYPE`` second row."""
    p = Path(path)
    s = guess_sep(p, sep)
    skip = None
    if fmt in ("scp_metadata", "scp_cluster") or (fmt == "auto" and _has_scp_type_row(p, s)):
        skip = [1]
    df = pd.read_csv(p, sep=s, skiprows=skip, dtype=str, keep_default_na=False, compression="infer")
    if index_col is not None:
        if index_col not in df.columns:
            raise KeyError(f"Column '{index_col}' not found in {p} (columns: {list(df.columns)[:10]})")
        df = df.set_index(index_col)
    df.index = df.index.astype(str)
    return df


def read_guide_assignments(
    path: str | Path,
    cell_column: str = "Cell",
    guides_column: str = "sgRNAs",
    list_separator: str = ",",
    fmt: str = "auto",
    sep: str = "auto",
) -> pd.Series:
    """Return a Series cell -> list[str] of assigned guides (empty list = none)."""
    df = read_table(path, fmt=fmt, sep=sep)
    for col in (cell_column, guides_column):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in guide assignment file {path}")
    ser = df.set_index(cell_column)[guides_column].fillna("").astype(str)
    lists = ser.map(lambda v: [g.strip() for g in v.split(list_separator) if g.strip()] if v else [])
    lists.index = lists.index.astype(str)
    return lists


def read_embedding(path: str | Path, columns: Sequence[str], index_col: str = "NAME", fmt: str = "auto", sep: str = "auto") -> pd.DataFrame:
    df = read_table(path, fmt=fmt, sep=sep, index_col=index_col)
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"Embedding file {path} lacks columns {missing}")
    return df[list(columns)].astype(np.float32)


# ---------------------------------------------------------------------------
# Value-state detection
# ---------------------------------------------------------------------------


@dataclass
class ValueState:
    state: str  # raw_counts | log_normalized | normalized | uncertain
    is_integer: bool
    min: float
    max: float
    #: For log-normalized data: the per-cell sum of expm1(x) if constant (e.g. 1e6).
    log_scale: Optional[float] = None
    log_scale_cv: Optional[float] = None
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return dict(state=self.state, is_integer=self.is_integer, min=self.min, max=self.max,
                    log_scale=self.log_scale, log_scale_cv=self.log_scale_cv, notes=list(self.notes))


def detect_value_state(X: sp.spmatrix, sample_cells: int = 2000, seed: int = 0) -> ValueState:
    """Infer whether a cells-by-features matrix holds raw counts or (log-)normalized values."""
    X = sp.csr_matrix(X)
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=min(sample_cells, n), replace=False)) if n > sample_cells else np.arange(n)
    S = X[idx]
    data = S.data
    if data.size == 0:
        return ValueState("uncertain", True, 0.0, 0.0, notes=["matrix is empty"])
    vmin, vmax = float(data.min()), float(data.max())
    is_int = bool(np.allclose(data, np.round(data), atol=1e-6))
    notes: List[str] = []
    if vmin < 0:
        return ValueState("normalized", False, vmin, vmax, notes=["negative values present (scaled/centred?)"])
    if is_int:
        return ValueState("raw_counts", True, vmin, vmax, notes=["all sampled values are integers"])
    # Log-normalized data: expm1 row sums are constant across cells.
    rs = np.asarray(np.expm1(S).sum(axis=1)).ravel()
    rs = rs[rs > 0]
    if rs.size:
        cv = float(rs.std() / rs.mean())
        if cv < 0.01:
            scale = float(np.median(rs))
            notes.append(f"expm1 row sums constant (~{scale:.6g}) -> log1p of size-normalized values")
            return ValueState("log_normalized", False, vmin, vmax, scale, cv, notes)
        notes.append(f"expm1 row sums not constant (cv={cv:.3f})")
    rs2 = np.asarray(S.sum(axis=1)).ravel()
    rs2 = rs2[rs2 > 0]
    if rs2.size and float(rs2.std() / rs2.mean()) < 0.01:
        notes.append("row sums constant -> size-normalized without log")
        return ValueState("normalized", False, vmin, vmax, notes=notes)
    if vmax < 30:
        notes.append("non-integer, bounded values -> treated as transformed/normalized")
        return ValueState("normalized", False, vmin, vmax, notes=notes)
    return ValueState("uncertain", False, vmin, vmax, notes=notes)


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


TABULAR_SUFFIXES = (".csv", ".tsv", ".txt", ".tab", ".csv.gz", ".tsv.gz", ".txt.gz", ".h5", ".h5ad", ".mtx", ".mtx.gz")


def discover_files(data_dir: str | Path) -> List[Path]:
    """List candidate data files under a directory (skips .git and hidden files)."""
    root = Path(data_dir)
    out: List[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith(".") and d not in ("results", "__pycache__", "src", "tests")]
        for fn in sorted(filenames):
            if fn.startswith("."):
                continue
            if fn.lower().endswith(TABULAR_SUFFIXES):
                out.append(Path(dirpath) / fn)
    return sorted(out)
