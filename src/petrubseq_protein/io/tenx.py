"""10x Genomics feature-barcode matrices: MTX directories and HDF5 files.

Both readers return a :class:`~petrubseq_protein.io.adapters.FeatureMatrix`
(cells x features, float32 CSR) whose ``var`` table carries ``feature_id``,
``feature_name``, ``feature_types`` (10x column name) and any further feature
columns (``genome``, ...). Matrices are indexed by the 10x feature ID; :func:`rename_features` re-indexes
the RNA part by ``name`` when ``var_names: name`` is configured (gene symbols) (gene symbols; duplicates are made
unique with ``-1``, ``-2`` ... suffixes, Scanpy-style, and the count is stored
in ``var.attrs['n_names_made_unique']``).

MTX: ``barcodes.tsv``, ``features.tsv`` (or Cell Ranger v2 ``genes.tsv``) and
``matrix.mtx``, gzipped or not; the third features column is the feature type
(absent in v2 -> ``Gene Expression``).

H5: Cell Ranger v3+ layout (``/matrix/{data,indices,indptr,shape,barcodes}``,
``/matrix/features/{id,name,feature_type,...}``); the v2 layout (one group per
genome with ``genes`` / ``gene_names``) is read as gene expression only.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import List, Optional, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from .adapters import FeatureMatrix, InputError

logger = logging.getLogger(__name__)

FEATURE_TYPE_COL = "feature_types"


def _find(path: Path, stem: str, alternatives: Tuple[str, ...] = ()) -> Optional[Path]:
    for s in (stem, *alternatives):
        for cand in (path / s, path / f"{s}.gz"):
            if cand.is_file():
                return cand
    return None


def _read_tsv(p: Path, names: Optional[List[str]] = None) -> pd.DataFrame:
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as fh:
        return pd.read_csv(fh, sep="\t", header=None, dtype=str, keep_default_na=False, names=names)


def make_unique(names: np.ndarray) -> Tuple[np.ndarray, int]:
    """Scanpy-style ``var_names_make_unique``: repeated names get ``-1``, ``-2`` ..."""
    s = pd.Index(names.astype(str))
    if not s.has_duplicates:
        return s.to_numpy(dtype=object), 0
    seen: dict = {}
    out = []
    n_changed = 0
    for n in s:
        if n in seen:
            seen[n] += 1
            new = f"{n}-{seen[n]}"
            while new in seen:
                seen[n] += 1
                new = f"{n}-{seen[n]}"
            out.append(new)
            seen[new] = 0
            n_changed += 1
        else:
            seen[n] = 0
            out.append(n)
    return np.asarray(out, dtype=object), n_changed


def _build(cells: np.ndarray, X: sp.csr_matrix, var: pd.DataFrame, source: str) -> FeatureMatrix:
    var = var.copy()
    var.index = pd.Index(var["feature_id"].to_numpy(dtype=object), name="feature")
    var.attrs["var_names"] = "id"
    return FeatureMatrix(np.asarray(cells, dtype=object), var.index.to_numpy(dtype=object), X, var, source)


def rename_features(fm: FeatureMatrix, var_names: str) -> FeatureMatrix:
    """Re-index a feature matrix by ``feature_id`` or ``feature_name`` (made unique)."""
    if var_names not in ("id", "name"):
        raise InputError(f"var_names must be 'id' or 'name', got {var_names!r}")
    col = "feature_id" if var_names == "id" else "feature_name"
    if col not in fm.var.columns:
        return fm
    feats, n_changed = make_unique(fm.var[col].to_numpy())
    if n_changed:
        logger.warning("%s: %d duplicate feature %ss made unique with -1/-2 suffixes", fm.source, n_changed, var_names)
    var = fm.var.copy()
    var.index = pd.Index(feats, name="feature")
    var.attrs["n_names_made_unique"] = n_changed
    var.attrs["var_names"] = var_names
    return FeatureMatrix(fm.cells, np.asarray(feats, dtype=object), fm.X, var, fm.source)


# ---------------------------------------------------------------------------
# MTX directory
# ---------------------------------------------------------------------------


def read_mtx_dir(path: str | Path) -> FeatureMatrix:
    path = Path(path)
    if not path.is_dir():
        raise InputError(f"not a directory: {path}")
    bc, ft, mm = _find(path, "barcodes.tsv"), _find(path, "features.tsv", ("genes.tsv",)), _find(path, "matrix.mtx")
    missing = [n for n, p in (("barcodes.tsv", bc), ("features.tsv/genes.tsv", ft), ("matrix.mtx", mm)) if p is None]
    if missing:
        raise InputError(f"{path} is not a 10x MTX directory; missing {missing} (with or without .gz)")
    barcodes = _read_tsv(bc).iloc[:, 0].to_numpy(dtype=object)
    feats = _read_tsv(ft)
    if feats.shape[1] < 2:
        raise InputError(f"{ft}: expected at least two columns (id, name)")
    var = pd.DataFrame({"feature_id": feats.iloc[:, 0].to_numpy(), "feature_name": feats.iloc[:, 1].to_numpy()})
    var[FEATURE_TYPE_COL] = feats.iloc[:, 2].to_numpy() if feats.shape[1] >= 3 else "Gene Expression"
    for j in range(3, feats.shape[1]):
        var[f"feature_col{j}"] = feats.iloc[:, j].to_numpy()
    opener = gzip.open if mm.suffix == ".gz" else open
    with opener(mm, "rb") as fh:
        M = scipy.io.mmread(fh)  # features x cells
    if M.shape != (len(var), len(barcodes)):
        raise InputError(f"{mm}: matrix is {M.shape} but features.tsv has {len(var)} rows and barcodes.tsv {len(barcodes)} rows")
    X = sp.csr_matrix(M.T, dtype=np.float32)
    X.sum_duplicates()
    return _build(barcodes, X, var, str(path))


# ---------------------------------------------------------------------------
# 10x HDF5
# ---------------------------------------------------------------------------


def _str(ds) -> np.ndarray:
    arr = ds[()]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)


def read_10x_h5(path: str | Path) -> FeatureMatrix:
    path = Path(path)
    if not path.is_file():
        raise InputError(f"not a file: {path}")
    with h5py.File(path, "r") as f:
        if "matrix" in f:
            g = f["matrix"]
            if "features" not in g:
                raise InputError(f"{path}: /matrix has no 'features' group (not a Cell Ranger v3+ feature-barcode H5)")
            fg = g["features"]
            for k in ("id", "name", "feature_type"):
                if k not in fg:
                    raise InputError(f"{path}: /matrix/features lacks '{k}'")
            var = pd.DataFrame({"feature_id": _str(fg["id"]), "feature_name": _str(fg["name"]), FEATURE_TYPE_COL: _str(fg["feature_type"])})
            for k in fg:
                if k in ("id", "name", "feature_type", "_all_tag_keys") or not isinstance(fg[k], h5py.Dataset):
                    continue
                v = fg[k][()]
                if v.shape == (len(var),):
                    var[k] = _str(fg[k]) if v.dtype.kind in ("S", "O", "U") else v
        else:
            genomes = [k for k in f if isinstance(f[k], h5py.Group) and "genes" in f[k]]
            if len(genomes) != 1:
                raise InputError(f"{path}: neither /matrix (v3) nor a single genome group (v2) found; groups: {list(f)}")
            g = f[genomes[0]]
            var = pd.DataFrame({"feature_id": _str(g["genes"]), "feature_name": _str(g["gene_names"]), FEATURE_TYPE_COL: "Gene Expression", "genome": genomes[0]})
        for k in ("data", "indices", "indptr", "shape", "barcodes"):
            if k not in g:
                raise InputError(f"{path}: matrix group lacks '{k}'")
        shape = tuple(int(x) for x in g["shape"][()])  # (features, cells)
        data = g["data"][()].astype(np.float32)
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _str(g["barcodes"])
    if shape[0] != len(var) or shape[1] != len(barcodes):
        raise InputError(f"{path}: shape {shape} disagrees with {len(var)} features / {len(barcodes)} barcodes")
    # 10x stores features x cells in CSC; the same arrays are cells x features in CSR.
    X = sp.csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))
    X.sum_duplicates()
    return _build(barcodes, X, var, str(path))
