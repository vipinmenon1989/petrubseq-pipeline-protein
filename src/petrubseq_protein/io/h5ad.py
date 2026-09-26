"""AnnData (``.h5ad``) inputs with slot addressing.

A modality lives in ``X``, ``layers[key]`` or ``obsm[key]``. ``obsm`` matrices
need feature names: a DataFrame carries them as columns; a plain array needs
``uns['<key>_features']`` (a list) or ``uns['<key>']['features']``. Guide
assignments can come from an ``obs`` column (guides separated by
``list_separator``). Feature metadata comes from ``var`` (X / layers) or from
``uns['<key>_var']`` when present (obsm).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .adapters import FeatureMatrix, InputError

logger = logging.getLogger(__name__)


def _csr32(X) -> sp.csr_matrix:
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    return sp.csr_matrix(X, dtype=np.float32)


def _obsm_names(a: ad.AnnData, key: str, n: int) -> Optional[np.ndarray]:
    for cand in (a.uns.get(f"{key}_features"), (a.uns.get(key) or {}).get("features") if isinstance(a.uns.get(key), dict) else None):
        if cand is not None:
            arr = np.asarray(list(cand), dtype=object)
            if len(arr) == n:
                return arr
    return None


def extract(a: ad.AnnData, slot: str, key: Optional[str], source: str, feature_type_column: str = "feature_types") -> FeatureMatrix:
    cells = np.asarray(a.obs_names, dtype=object)
    if slot == "X":
        if a.X is None:
            raise InputError(f"{source}: X is empty")
        return FeatureMatrix(cells, np.asarray(a.var_names, dtype=object), _csr32(a.X), a.var.copy(), f"{source}[X]")
    if slot == "layers":
        if key not in a.layers:
            raise InputError(f"{source}: layers[{key!r}] not found; layers present: {list(a.layers)}")
        return FeatureMatrix(cells, np.asarray(a.var_names, dtype=object), _csr32(a.layers[key]), a.var.copy(), f"{source}[layers/{key}]")
    if slot == "obsm":
        if key not in a.obsm:
            raise InputError(f"{source}: obsm[{key!r}] not found; obsm present: {list(a.obsm)}")
        obj = a.obsm[key]
        if isinstance(obj, pd.DataFrame):
            names = np.asarray(obj.columns, dtype=object)
            X = _csr32(obj.to_numpy(dtype=np.float32))
        else:
            names = _obsm_names(a, key, obj.shape[1])
            if names is None:
                raise InputError(f"{source}: obsm[{key!r}] is an unlabelled array with {obj.shape[1]} columns; feature names are needed in uns['{key}_features'] (or store a DataFrame)")
            X = _csr32(obj)
        var = a.uns.get(f"{key}_var")
        var = var.copy() if isinstance(var, pd.DataFrame) and len(var) == len(names) else pd.DataFrame(index=pd.Index(names))
        var.index = pd.Index(names, name="feature")
        return FeatureMatrix(cells, names, X, var, f"{source}[obsm/{key}]")
    raise InputError(f"{source}: unsupported slot {slot!r}")


def read_h5ad_matrix(path: str | Path, slot: str = "X", key: Optional[str] = None, feature_type_column: str = "feature_types") -> FeatureMatrix:
    path = Path(path)
    if not path.is_file():
        raise InputError(f"not a file: {path}")
    a = ad.read_h5ad(path)
    return extract(a, slot, key, str(path), feature_type_column)


def read_h5ad_obs(path: str | Path, cell_id: Optional[str] = None) -> pd.DataFrame:
    """``obs`` of an h5ad as a cell table (only obs is read). ``cell_id`` may name an
    obs column to use as the index; otherwise ``obs_names`` are the cell IDs."""
    import h5py
    from anndata.io import read_elem

    path = Path(path)
    if not path.is_file():
        raise InputError(f"not a file: {path}")
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
    obs = obs.copy()
    for c in obs.columns:
        if obs[c].dtype.name == "category":
            obs[c] = obs[c].astype(str)
    if cell_id and cell_id in obs.columns:
        obs = obs.set_index(obs[cell_id].astype(str))
    obs.index = obs.index.astype(str)
    obs.index.name = cell_id or "obs_names"
    return obs


def guide_lists_from_obs(a, key: str, sep: str, source: str) -> pd.Series:
    obs = a if isinstance(a, pd.DataFrame) else a.obs
    if key not in obs.columns:
        raise InputError(f"{source}: obs[{key!r}] not found; columns: {list(obs.columns)[:12]} ...")
    col = obs[key].astype(str).replace({"nan": "", "None": ""}).fillna("")
    lists = col.map(lambda v: [g.strip() for g in v.split(sep) if g.strip()] if v else [])
    lists.index = lists.index.astype(str)
    lists.name = key
    return lists


def read_h5ad_slots(path: str | Path, slots: Any, feature_type_column: str = "feature_types") -> Dict[str, Any]:
    """Read one h5ad and extract every declared slot (``MultiplexedSlots``)."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"not a file: {path}")
    a = ad.read_h5ad(path)
    out: Dict[str, Any] = {}
    for mod in slots.declared():
        sa = getattr(slots, mod)
        if mod == "guide_assignments":
            out[mod] = guide_lists_from_obs(a, sa.key, sa.list_separator, str(path))
        else:
            out[mod] = extract(a, sa.slot, sa.key, str(path), feature_type_column)
    return out
