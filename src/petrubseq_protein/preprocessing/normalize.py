"""Modality-specific normalization with input-state detection.

RNA
---
* raw counts          -> ``layers['counts']`` kept; ``X`` = normalize_total + log1p
* log-normalized      -> ``X`` kept untouched (never normalized twice); if a
                         per-cell library size is available, integer counts are
                         reconstructed into ``layers['counts']`` and validated.

Protein (ADT)
-------------
* provided normalized -> kept as the primary ``obsm['protein']``
* raw counts          -> primary normalization = CLR (default), log1p, or the
                         isotype-ratio transform of Frangieh et al. 2021
                         (``max(0, log((x+1)/(isotype+1)))``); counts are always
                         kept in ``obsm['protein_counts']``.
A secondary raw-derived transform (default CLR) can be stored alongside a
provided normalization so both are available downstream, clearly labelled.
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import ProteinConfig
from ..io.readers import Matrix, ValueState, detect_value_state

logger = logging.getLogger("petrubseq_protein")


# ---------------------------------------------------------------------------
# RNA
# ---------------------------------------------------------------------------


def reconstruct_counts(X: sp.csr_matrix, total_counts: np.ndarray, scale: float, tolerance: float = 0.05, block_cells: int = 20000) -> Tuple[Optional[sp.csr_matrix], Dict[str, Any]]:
    """Recover integer counts from ``log1p(scale * count / total)`` values.

    Returns (counts or None, info). Counts are accepted only if every
    non-zero entry is within ``tolerance`` of an integer (and >= 1);
    otherwise ``None`` is returned and the info dict says why. Work is done
    in blocks of ``block_cells`` rows so memory stays close to the size of
    the result.
    """
    X = sp.csr_matrix(X)
    tot = np.asarray(total_counts, dtype=np.float64)
    info: Dict[str, Any] = {"attempted": True, "scale": float(scale), "tolerance": tolerance}
    if np.isnan(tot).any() or (tot <= 0).any():
        info.update(accepted=False, reason="missing or non-positive per-cell total counts")
        return None, info
    data = np.empty(X.nnz, dtype=np.int32)
    max_dev = 0.0
    n_bad = 0
    min_nz = np.iinfo(np.int32).max
    n_rows_ok = 0
    for r0 in range(0, X.shape[0], block_cells):
        r1 = min(r0 + block_cells, X.shape[0])
        a, b = X.indptr[r0], X.indptr[r1]
        rows = np.repeat(np.arange(r0, r1), np.diff(X.indptr[r0 : r1 + 1]))
        raw = np.expm1(X.data[a:b].astype(np.float64)) * tot[rows] / float(scale)
        rounded = np.rint(raw)
        dev = np.abs(raw - rounded)
        if dev.size:
            max_dev = max(max_dev, float(dev.max()))
            n_bad += int((dev > tolerance).sum())
            min_nz = min(min_nz, int(rounded.min()))
        data[a:b] = rounded.astype(np.int32)
        rs = np.bincount(rows - r0, weights=rounded, minlength=r1 - r0)
        n_rows_ok += int((np.abs(rs - tot[r0:r1]) <= 0.5).sum())
        del raw, rounded, dev, rows
    info["max_abs_deviation"] = max_dev
    info["frac_entries_within_tolerance"] = 1.0 - n_bad / max(1, X.nnz)
    info["row_sum_matches_total"] = n_rows_ok / max(1, X.shape[0])
    info["min_nonzero_count"] = int(min_nz) if X.nnz else 0
    ok = info["frac_entries_within_tolerance"] >= 0.999 and info["min_nonzero_count"] >= 1
    info["accepted"] = bool(ok)
    if not ok:
        info["reason"] = "values do not map to integers (data may have been filtered/rescaled after normalization)"
        return None, info
    if info["row_sum_matches_total"] < 0.99:
        info["note"] = "row sums differ from provided totals: some genes were removed after normalization; counts are exact per gene but totals are partial"
    counts = sp.csr_matrix((data, X.indices.copy(), X.indptr.copy()), shape=X.shape)
    return counts, info


def normalize_rna(X: sp.csr_matrix, cfg_rna, total_counts: Optional[np.ndarray] = None, state_hint: str = "auto", apply_normalization: bool = True) -> Tuple[sp.csr_matrix, Optional[sp.csr_matrix], Optional[str], Dict[str, Any]]:
    """Return (X_normalized, counts_matrix_or_None, counts_layer_name_or_None, info).

    ``counts_layer_name`` is ``'counts'`` only for directly observed raw
    counts and ``'reconstructed_counts'`` for counts recovered from
    log-normalized values, so downstream users can never mistake one for the
    other.
    """
    import scanpy as sc
    import anndata as ad

    vs: ValueState = detect_value_state(X)
    forced = cfg_rna.input_state if cfg_rna.input_state != "auto" else (
        "raw_counts" if state_hint == "raw_counts" else ("log_normalized" if state_hint == "normalized" else None)
    )
    detected = forced or vs.state
    if forced and forced != vs.state:
        logger.warning("RNA input state forced to %s but values look like %s", forced, vs.state)
    info: Dict[str, Any] = {"detected_state": vs.to_dict(), "input_state_used": detected, "method": None, "counts_layer": None, "counts_source": None}
    mode = cfg_rna.normalize
    if detected == "raw_counts" and mode in ("auto", "always") and not apply_normalization:
        # Lifecycle: normalization is applied after filtering (see log_normalize_rna).
        counts = sp.csr_matrix(X, dtype=np.int32)
        info.update(method="pending (normalize_total" + ("+log1p" if cfg_rna.log1p else "") + " after filtering)", target_sum=cfg_rna.target_sum, counts_layer="counts", counts_source="input (observed raw counts)", pending_normalization=True)
        return counts.astype(np.float32), counts, "counts", info
    if detected == "raw_counts" and mode in ("auto", "always"):
        counts = sp.csr_matrix(X, dtype=np.float32)
        tmp = ad.AnnData(X=counts.copy())
        sc.pp.normalize_total(tmp, target_sum=cfg_rna.target_sum)
        if cfg_rna.log1p:
            sc.pp.log1p(tmp)
        info.update(method="normalize_total" + ("+log1p" if cfg_rna.log1p else ""), target_sum=cfg_rna.target_sum, counts_layer="counts", counts_source="input (observed raw counts)")
        return sp.csr_matrix(tmp.X, dtype=np.float32), counts.astype(np.int32), "counts", info
    if detected == "raw_counts" and mode == "never":
        info.update(method="none (raw counts kept in X)", counts_layer="counts", counts_source="input (observed raw counts)")
        return sp.csr_matrix(X, dtype=np.float32), sp.csr_matrix(X, dtype=np.int32), "counts", info
    if mode == "always" and detected != "raw_counts":
        raise ValueError("rna.normalize=always requested but the input is not raw counts; refusing to normalize twice.")
    # Already normalized: keep as is.
    info.update(method=f"none (input already {detected}; preserved as X)")
    counts = None
    layer = None
    if detected == "log_normalized" and cfg_rna.reconstruct_counts == "auto" and total_counts is not None and vs.log_scale:
        counts, rinfo = reconstruct_counts(X, total_counts, vs.log_scale, cfg_rna.reconstruct_tolerance)
        info["counts_reconstruction"] = rinfo
        if counts is not None:
            layer = "reconstructed_counts"
            info["counts_layer"] = layer
            info["counts_source"] = (
                "reconstructed: round(expm1(X) * total_counts / scale) from the provided log-normalized "
                "matrix and the per-cell library size; NOT the directly observed UMI matrix (genes removed "
                "upstream are absent)"
            )
            logger.info("RNA counts reconstructed from log-normalized values (scale=%.6g) -> layers['%s']", vs.log_scale, layer)
        else:
            logger.warning("RNA counts could not be reconstructed: %s", rinfo.get("reason"))
    return sp.csr_matrix(X, dtype=np.float32), counts, layer, info


def log_normalize_rna(adata, cfg_rna, info: Dict[str, Any]) -> None:
    """Apply the pending normalize_total(+log1p) to ``adata.X`` from ``layers['counts']``."""
    import scanpy as sc

    if not info.get("pending_normalization"):
        return
    adata.X = sp.csr_matrix(adata.layers["counts"], dtype=np.float32).copy()
    sc.pp.normalize_total(adata, target_sum=cfg_rna.target_sum)
    if cfg_rna.log1p:
        sc.pp.log1p(adata)
        adata.uns.pop("log1p", None)
    info["method"] = "normalize_total" + ("+log1p" if cfg_rna.log1p else "")
    info["pending_normalization"] = False
    logger.info("RNA normalized: %s (target_sum=%s)", info["method"], cfg_rna.target_sum)


# ---------------------------------------------------------------------------
# Protein
# ---------------------------------------------------------------------------


def is_isotype(name: str, patterns: Sequence[str]) -> bool:
    return any(re.search(p, name, flags=re.IGNORECASE) for p in patterns)


def canonical_protein_name(name: str, pcfg: ProteinConfig) -> str:
    n = str(name)
    if pcfg.strip_feature_suffix and n.endswith(pcfg.strip_feature_suffix):
        n = n[: -len(pcfg.strip_feature_suffix)]
    n = n.strip()
    return pcfg.rename.get(n, n)


def clr(X: np.ndarray, axis: str = "features", pseudocount: float = 1.0) -> np.ndarray:
    """Centred log-ratio: log1p(x) minus the mean log1p over the chosen axis."""
    L = np.log1p(np.asarray(X, dtype=np.float64) / 1.0 if pseudocount == 1.0 else np.asarray(X, dtype=np.float64))
    if pseudocount != 1.0:
        L = np.log(np.asarray(X, dtype=np.float64) + pseudocount)
    if axis == "features":
        return (L - L.mean(axis=1, keepdims=True)).astype(np.float32)
    return (L - L.mean(axis=0, keepdims=True)).astype(np.float32)


def isotype_ratio(X: np.ndarray, names: Sequence[str], iso_map: Dict[str, str]) -> Tuple[np.ndarray, List[str]]:
    """Frangieh et al. transform: max(0, log((x+1)/(isotype+1))) per antibody."""
    names = list(names)
    targets = [n for n in names if n in iso_map]
    out = np.zeros((X.shape[0], len(targets)), dtype=np.float32)
    for j, n in enumerate(targets):
        i, c = names.index(n), names.index(iso_map[n])
        out[:, j] = np.maximum(0.0, np.log((X[:, i] + 1.0) / (X[:, c] + 1.0)))
    return out, targets


def verify_isotype_ratio(norm: Matrix, counts: Matrix, pcfg: ProteinConfig, max_cells: int = 20000) -> Dict[str, Any]:
    """Check whether the provided normalized protein matrix equals the isotype-ratio
    transform of the raw counts for *some* isotype control per antibody.
    Returns the inferred antibody -> isotype map and residuals."""
    common = pd.Index(norm.cells).intersection(pd.Index(counts.cells))[:max_cells]
    if len(common) == 0:
        return {"checked": False, "reason": "no shared cells"}
    N = norm.subset_cells(common).X.toarray().astype(np.float64)
    C = counts.subset_cells(common).X.toarray().astype(np.float64)
    nnames = [canonical_protein_name(n, pcfg) for n in norm.features]
    cnames = [str(n) for n in counts.features]
    iso = [n for n in cnames if is_isotype(n, pcfg.isotype_patterns)]
    result: Dict[str, Any] = {"checked": True, "n_cells_checked": int(len(common)), "isotype_controls": iso, "inferred_isotype_map": {}, "max_abs_diff": {}, "unmatched": []}
    for j, n in enumerate(nnames):
        if n not in cnames:
            result["unmatched"].append(n)
            continue
        x = C[:, cnames.index(n)]
        best, best_err = None, np.inf
        for c in iso:
            pred = np.maximum(0.0, np.log((x + 1.0) / (C[:, cnames.index(c)] + 1.0)))
            err = float(np.abs(pred - N[:, j]).max())
            if err < best_err:
                best, best_err = c, err
        # Also consider plain log1p and CLR-free transforms
        result["inferred_isotype_map"][n] = best
        result["max_abs_diff"][n] = round(best_err, 6)
    result["formula_reproduced"] = bool(result["max_abs_diff"]) and max(result["max_abs_diff"].values()) < 1e-4
    result["formula"] = "max(0, log((count+1)/(isotype_count+1)))"
    return result


def normalize_protein(
    norm: Optional[Matrix],
    counts: Optional[Matrix],
    cells: Sequence[str],
    pcfg: ProteinConfig,
    feature_table: Optional[pd.DataFrame] = None,
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Dict[str, pd.DataFrame], pd.DataFrame, Dict[str, Any]]:
    """Return (protein_df, protein_counts_df, extra_dfs, features_df, info), all indexed by ``cells``.

    ``feature_table`` (``protein.feature_table``) adds antibody annotations; see
    :func:`annotate_protein_features`.

    * ``protein_df``: primary normalized values, targeting antibodies only (``obsm['protein']``).
    * ``protein_counts_df``: raw counts for every antibody incl. isotype controls
      (``obsm['protein_counts']``; NaN if a cell lacks counts).
    * ``extra_dfs``: additional raw-derived matrices keyed by name (``obsm['protein_<name>']``).
    * ``features_df``: one row per antibody with role/isotype/summary info.
    """
    idx = pd.Index(cells)
    info: Dict[str, Any] = {"primary_method": None, "primary_source": None, "extra": {}, "isotypes": []}
    counts_df = None
    feats: Dict[str, Dict[str, Any]] = {}
    iso_map: Dict[str, str] = dict(pcfg.isotype_map)
    if counts is not None:
        vs = detect_value_state(counts.X)
        info["counts_state"] = vs.to_dict()
        cnames = [str(n) for n in counts.features]
        sub_idx = pd.Index(counts.cells).get_indexer(idx)
        arr = np.full((len(idx), len(cnames)), np.nan, dtype=np.float32)
        present = sub_idx >= 0
        arr[present] = counts.X[sub_idx[present]].toarray()
        counts_df = pd.DataFrame(arr, index=idx, columns=cnames)
        info["cells_without_counts"] = int((~present).sum())
        for n in cnames:
            feats[n] = {"in_counts": True, "is_isotype": is_isotype(n, pcfg.isotype_patterns)}
        info["isotypes"] = [n for n in cnames if feats[n]["is_isotype"]]
    provided_df = None
    if norm is not None:
        vs = detect_value_state(norm.X)
        info["provided_state"] = vs.to_dict()
        nnames = [canonical_protein_name(n, pcfg) for n in norm.features]
        sub_idx = pd.Index(norm.cells).get_indexer(idx)
        if (sub_idx < 0).any():
            raise KeyError(f"{int((sub_idx < 0).sum())} cells missing from the normalized protein matrix")
        provided_df = pd.DataFrame(norm.X[sub_idx].toarray().astype(np.float32), index=idx, columns=nnames)
        provided_df.columns.name = None
        for n, raw_n in zip(nnames, norm.features):
            feats.setdefault(n, {"in_counts": False, "is_isotype": is_isotype(n, pcfg.isotype_patterns)})
            feats[n]["in_normalized"] = True
            feats[n]["name_in_normalized_file"] = str(raw_n)
        iso_in_provided = [n for n in nnames if feats[n]["is_isotype"]]
        if iso_in_provided and pcfg.exclude_isotypes_from_embedding:
            provided_df = provided_df.drop(columns=iso_in_provided)
            info["isotypes_dropped_from_provided"] = iso_in_provided
    # Verify / infer the isotype map when both matrices exist.
    if counts is not None and norm is not None:
        chk = verify_isotype_ratio(norm, counts, pcfg)
        info["provided_vs_counts_check"] = chk
        if chk.get("formula_reproduced"):
            for n, c in chk["inferred_isotype_map"].items():
                feats.setdefault(n, {})["isotype_control"] = c
                iso_map.setdefault(n, c)
    for n, c in pcfg.isotype_map.items():
        feats.setdefault(n, {})["isotype_control"] = c
    # Primary matrix
    method = pcfg.normalization
    if method == "auto":
        method = "provided" if provided_df is not None else ("clr" if counts is not None else "none")
    if method == "provided":
        if provided_df is None:
            raise ValueError("protein.normalization=provided but no normalized protein matrix was given")
        prot_df, used, src = provided_df, "provided (kept as-is)", "inputs.protein"
    elif method == "none":
        prot_df, used, src = provided_df, "none", "inputs.protein" if provided_df is not None else None
    else:
        if counts_df is None:
            raise ValueError(f"protein.normalization={method} needs raw counts (inputs.protein_counts)")
        prot_df, used = _from_counts(counts_df, method, pcfg, iso_map)
        src = "inputs.protein_counts"
    info["primary_method"], info["primary_source"] = used, src
    if prot_df is None and counts_df is None:
        raise ValueError("No protein data loaded")
    # Extra representations from raw counts
    extras: Dict[str, pd.DataFrame] = {}
    if counts_df is not None:
        for name in pcfg.extra_representations:
            if name == method:
                continue
            try:
                extras[name], desc = _from_counts(counts_df, name, pcfg, iso_map)
                info["extra"][name] = desc
            except ValueError as exc:
                logger.warning("protein extra representation '%s' skipped: %s", name, exc)
        if provided_df is not None and method != "provided":
            extras["provided"] = provided_df
            info["extra"]["provided"] = "provided normalized matrix (kept as-is)"
    # Embedding representation
    emb = pcfg.embedding_representation
    if emb == "auto":
        emb = "clr" if counts_df is not None else ("provided" if provided_df is not None else method)
    if emb == method:
        emb_key = "protein"
    elif emb in extras:
        emb_key = f"protein_{emb}"
    elif counts_df is not None and emb in ("clr", "log1p"):
        extras[emb], desc = _from_counts(counts_df, emb, pcfg, iso_map)
        info["extra"][emb] = desc
        emb_key = f"protein_{emb}"
    elif emb == "provided" and provided_df is not None:
        extras["provided"] = provided_df
        emb_key = "protein_provided"
    else:
        raise ValueError(f"protein.embedding_representation={emb} is not available for this dataset")
    info["embedding_representation"], info["embedding_key"] = emb, emb_key
    features_df = pd.DataFrame.from_dict(feats, orient="index")
    features_df.index.name = "protein"
    for col in ("in_counts", "in_normalized", "is_isotype"):
        if col not in features_df.columns:
            features_df[col] = False
        features_df[col] = features_df[col].map(lambda v: bool(v) if v == v and v is not None else False).astype(bool)
    for col in ("name_in_normalized_file", "isotype_control"):
        if col in features_df.columns:
            features_df[col] = features_df[col].fillna("").astype(str)
    features_df["role"] = np.where(features_df["is_isotype"], "isotype_control", "target")
    emb_df = prot_df if emb_key == "protein" else extras[emb_key.removeprefix("protein_")]
    features_df["in_primary_matrix"] = features_df.index.isin(prot_df.columns if prot_df is not None else [])
    features_df["used_in_embedding"] = features_df.index.isin(emb_df.columns if emb_df is not None else [])
    features_df["used_in_qc"] = features_df["in_counts"] & (pcfg.use_isotypes_for_qc | ~features_df["is_isotype"])
    info["isotype_map"] = iso_map
    info["n_targeting_features"] = int((~features_df["is_isotype"]).sum())
    info["n_isotype_features"] = int(features_df["is_isotype"].sum())
    features_df = annotate_protein_features(features_df, counts, norm, pcfg, feature_table)
    info["feature_annotation_sources"] = sorted({s for v in features_df["annotation_source"] for s in v.split("+") if s and s != "none"})
    return prot_df, counts_df, extras, features_df, info


ANNOTATION_COLUMNS = ["feature_id", "antibody_name", "protein_name", "gene_symbol", "clone", "feature_type", "isotype", "isotype_control", "annotation_source"]


def annotate_protein_features(features_df: pd.DataFrame, counts: Optional[Matrix], norm: Optional[Matrix], pcfg: ProteinConfig, table: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Add structured antibody metadata (``ANNOTATION_COLUMNS``) to the feature table.

    Sources, in increasing precedence: nothing (empty string = not available,
    never invented), the input's own feature metadata (10x ``features.tsv``:
    ``feature_name`` -> ``antibody_name``, ``feature_types`` -> ``feature_type``),
    then ``protein.feature_table`` (``feature_id`` + any of the columns).
    ``feature_id`` is the canonical antibody name (the row index);
    ``isotype_control`` keeps the matched control inferred or configured
    earlier unless the table overrides it. ``annotation_source`` lists what
    contributed per row.
    """
    df = features_df.copy()
    n = len(df)
    src = pd.Series([[] for _ in range(n)], index=df.index, dtype=object)
    for c in ANNOTATION_COLUMNS:
        if c == "isotype_control":
            df[c] = df[c].fillna("").astype(str) if c in df.columns else ""
        elif c != "annotation_source":
            df[c] = ""
    df["feature_id"] = df.index.astype(str)
    # input feature metadata (10x) ------------------------------------------
    for m in (norm, counts):
        if m is None or m.features_meta is None:
            continue
        fm = m.features_meta.copy()
        fm.index = pd.Index([canonical_protein_name(str(f), pcfg) for f in fm.index])
        fm = fm[~fm.index.duplicated()]
        hit = df.index.isin(fm.index)
        if not hit.any():
            continue
        if "feature_name" in fm.columns:
            df.loc[hit, "antibody_name"] = fm["feature_name"].astype(object).reindex(df.index[hit]).fillna("").astype(str).to_numpy()
        if "feature_types" in fm.columns:
            df.loc[hit, "feature_type"] = fm["feature_types"].astype(object).reindex(df.index[hit]).fillna("").astype(str).to_numpy()
        if "feature_id" in fm.columns:
            df.loc[hit, "feature_id"] = fm["feature_id"].astype(object).reindex(df.index[hit]).fillna("").astype(str).to_numpy()
        for i in df.index[hit]:
            src[i].append("input_features")
    # user table ---------------------------------------------------------------
    if table is not None:
        t = table.copy()
        if "feature_id" not in t.columns:
            raise ValueError("protein.feature_table needs a 'feature_id' column")
        t["feature_id"] = t["feature_id"].astype(str)
        t = t.set_index(t["feature_id"].map(lambda f: canonical_protein_name(f, pcfg)))
        unknown = [f for f in t.index if f not in df.index]
        if unknown:
            logger.warning("protein.feature_table: %d rows match no antibody in the data (e.g. %s)", len(unknown), unknown[:3])
        cols = [c for c in ANNOTATION_COLUMNS if c in t.columns and c not in ("annotation_source", "feature_id")]
        hit = df.index.isin(t.index)
        for c in cols:
            vals = t[c].reindex(df.index[hit])
            ok = vals.notna() & (vals.astype(str).str.strip() != "")
            df.loc[df.index[hit][ok.to_numpy()], c] = vals[ok].astype(str).to_numpy()
        for i in df.index[hit]:
            src[i].append("feature_table")
    df["annotation_source"] = src.map(lambda l: "+".join(l) if l else "none")
    for c in ANNOTATION_COLUMNS:
        df[c] = df[c].fillna("").astype(str)
    return df


def _from_counts(counts_df: pd.DataFrame, method: str, pcfg: ProteinConfig, iso_map: Dict[str, str]) -> Tuple[pd.DataFrame, str]:
    names = list(counts_df.columns)
    targets = [n for n in names if not is_isotype(n, pcfg.isotype_patterns)]
    X = counts_df.to_numpy(dtype=np.float64)
    if method == "clr":
        sub = counts_df[targets].to_numpy(dtype=np.float64)
        out = clr(np.nan_to_num(sub, nan=0.0), axis=pcfg.clr_axis)
        out[np.isnan(sub).any(axis=1)] = np.nan
        return pd.DataFrame(out, index=counts_df.index, columns=targets), f"clr (axis={pcfg.clr_axis}, log1p, targeting antibodies only)"
    if method == "log1p":
        out = np.log1p(counts_df[targets].to_numpy(dtype=np.float64)).astype(np.float32)
        return pd.DataFrame(out, index=counts_df.index, columns=targets), "log1p"
    if method == "isotype_ratio":
        if not iso_map:
            raise ValueError("protein.normalization=isotype_ratio requires protein.isotype_map")
        out, cols = isotype_ratio(np.nan_to_num(X, nan=0.0), names, iso_map)
        out[np.isnan(X).any(axis=1)] = np.nan
        return pd.DataFrame(out, index=counts_df.index, columns=cols), "isotype_ratio max(0, log((x+1)/(iso+1)))"
    raise ValueError(f"Unknown protein normalization '{method}'")
