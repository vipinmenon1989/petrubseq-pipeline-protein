#!/usr/bin/env python
"""Compare two processed objects produced by petrubseq-pipeline-protein.

    python scripts/compare_processed.py REFERENCE.h5ad NEW.h5ad [--atol 1e-5] [--allow-uns KEY ...]

Regression convention (docs/REGRESSION.md): the immutable references live under
``results/baselines/v0.1/{SCP1064,SCP1064_smoke}``, the current outputs under
``results/{SCP1064,SCP1064_smoke}``, e.g.

    python scripts/compare_processed.py ../results/baselines/v0.1/SCP1064_smoke/processed/scp1064_smoke_processed.h5ad \
        ../results/SCP1064_smoke/processed/scp1064_smoke_processed.h5ad

Reports, per slot, whether the two objects are identical (labels, matrices),
numerically close (embeddings), or different. Differences outside the
whitelisted ``uns`` keys make the script exit non-zero, so it can be used as a
regression gate after refactors.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# uns keys that legitimately differ between runs of the same config (provenance)
DEFAULT_ALLOW_UNS = ["petrubseq_protein/provenance", "petrubseq_protein/config", "petrubseq_protein/schema", "petrubseq_protein/version", "petrubseq_protein/embeddings/diagnostics/sampled_cells"]
# Slots that a newer schema may ADD (present only in NEW) without being a regression.
# Everything else that exists in only one object is a DIFF; anything present in
# both is compared strictly regardless of these lists.
ADDITIVE_OBS = {"has_guide_counts", "guide_total_counts", "n_guides_detected", "guides_detected", "guide_top", "guide_top_count", "guide_second_count", "guide_dominant_call", "lane_id", "barcode_original",
                # reference-parity QC additions (haemoglobin gene class, as in perturbseq-pipeline qc.py)
                "pct_counts_hb", "total_counts_hb", "log1p_total_counts_hb"}
ADDITIVE_VAR = {"hb"}
ADDITIVE_OBSM = {"guide_counts"}
ADDITIVE_UNS_TOP = {"guide_features"}
ADDITIVE_PROTEIN_FEATURE_COLS = {"feature_id", "antibody_name", "protein_name", "gene_symbol", "clone", "feature_type", "isotype", "annotation_source"}


def _dense(m):
    return m.toarray() if sp.issparse(m) else np.asarray(m)


def _same_sparse(a, b, atol: float) -> Tuple[bool, str]:
    if sp.issparse(a) and sp.issparse(b):
        a, b = sp.csr_matrix(a), sp.csr_matrix(b)
        a.sort_indices(); b.sort_indices()
        if a.shape != b.shape:
            return False, f"shape {a.shape} vs {b.shape}"
        if a.nnz != b.nnz or not np.array_equal(a.indptr, b.indptr) or not np.array_equal(a.indices, b.indices):
            return False, f"sparsity pattern differs (nnz {a.nnz} vs {b.nnz})"
        d = np.abs(a.data.astype(np.float64) - b.data.astype(np.float64))
        return bool(d.max() <= atol if d.size else True), f"max|diff|={d.max() if d.size else 0:.3g}"
    a, b = _dense(a).astype(np.float64), _dense(b).astype(np.float64)
    if a.shape != b.shape:
        return False, f"shape {a.shape} vs {b.shape}"
    both_nan = np.isnan(a) & np.isnan(b)
    d = np.abs(np.where(both_nan, 0, a - b))
    d = np.where(np.isnan(d), np.inf, d)
    return bool(d.max() <= atol), f"max|diff|={d.max():.3g}"


def _flatten(d: Any, prefix: str = "") -> dict:
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flatten(v, f"{prefix}/{k}" if prefix else str(k)))
    else:
        out[prefix] = d
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("new")
    ap.add_argument("--atol", type=float, default=1e-5, help="tolerance for matrices and embeddings")
    ap.add_argument("--umap-atol", type=float, default=1e-3, help="tolerance for UMAP coordinates")
    ap.add_argument("--allow-uns", nargs="*", default=[], help="extra uns key prefixes allowed to differ")
    a = ap.parse_args()
    A, B = ad.read_h5ad(a.reference), ad.read_h5ad(a.new)
    rows: List[Tuple[str, str, str]] = []
    fail = False

    def rec(slot, ok, detail="", soft=False, added=False):
        nonlocal fail
        rows.append((slot, "SAME" if ok else ("ADD" if added else ("NOTE" if soft else "DIFF")), detail))
        if not ok and not soft and not added:
            fail = True

    rec("obs_names", list(A.obs_names) == list(B.obs_names), f"{A.n_obs} vs {B.n_obs} cells")
    rec("var_names", list(A.var_names) == list(B.var_names), f"{A.n_vars} vs {B.n_vars} genes")
    # obs / var columns
    for name, da, db in (("obs", A.obs, B.obs), ("var", A.var, B.var)):
        only_a, only_b = sorted(set(da.columns) - set(db.columns)), sorted(set(db.columns) - set(da.columns))
        additive = not only_a and all(c in (ADDITIVE_OBS if name == "obs" else ADDITIVE_VAR) for c in only_b)
        rec(f"{name} columns", not only_a and not only_b, f"only in reference: {only_a}; only in new: {only_b}", added=additive)
        for c in sorted(set(da.columns) & set(db.columns)):
            x, y = da[c], db[c]
            if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y) and not pd.api.types.is_bool_dtype(x):
                xv, yv = x.to_numpy(float), y.to_numpy(float)
                ok = np.allclose(xv, yv, atol=a.atol, equal_nan=True)
                rec(f"{name}[{c}]", ok, "" if ok else f"max|diff|={np.nanmax(np.abs(xv - yv)):.3g}")
            else:
                ok = (x.astype(str).to_numpy() == y.astype(str).to_numpy()).all()
                rec(f"{name}[{c}]", ok, "" if ok else f"{int((x.astype(str).to_numpy() != y.astype(str).to_numpy()).sum())} cells differ")
    ok, d = _same_sparse(A.X, B.X, a.atol); rec("X", ok, d)
    for L in sorted(set(A.layers) | set(B.layers)):
        if L in A.layers and L in B.layers:
            ok, d = _same_sparse(A.layers[L], B.layers[L], 0); rec(f"layers[{L}]", ok, d)
        else:
            rec(f"layers[{L}]", False, "present in only one object")
    for k in sorted(set(A.obsm) | set(B.obsm)):
        if k not in A.obsm or k not in B.obsm:
            rec(f"obsm[{k}]", False, "present only in new (expected additive slot)" if k in B.obsm and k in ADDITIVE_OBSM else "present in only one object", added=k in B.obsm and k in ADDITIVE_OBSM); continue
        ma, mb = A.obsm[k], B.obsm[k]
        if isinstance(ma, pd.DataFrame) or isinstance(mb, pd.DataFrame):
            if list(ma.columns) != list(mb.columns):
                rec(f"obsm[{k}]", False, f"columns differ: {list(ma.columns)[:5]} vs {list(mb.columns)[:5]}"); continue
            ok, d = _same_sparse(ma.to_numpy(float), mb.to_numpy(float), a.atol)
        else:
            tol = a.umap_atol if "umap" in k else a.atol
            ok, d = _same_sparse(ma, mb, tol)
        rec(f"obsm[{k}]", ok, d)
    for k in sorted(set(A.varm) | set(B.varm)):
        if k in A.varm and k in B.varm:
            ok, d = _same_sparse(A.varm[k], B.varm[k], a.atol); rec(f"varm[{k}]", ok, d)
        else:
            rec(f"varm[{k}]", False, "present in only one object")
    # uns: protein_features + flattened petrubseq_protein
    if "protein_features" in A.uns and "protein_features" in B.uns:
        fa, fb = A.uns["protein_features"], B.uns["protein_features"]
        common = sorted(set(fa.columns) & set(fb.columns))
        ok = list(fa.index) == list(fb.index) and all((fa[c].astype(str) == fb[c].astype(str)).all() for c in common)
        rec("uns[protein_features]", ok, f"common cols {len(common)}; only ref {sorted(set(fa.columns) - set(fb.columns))}; only new {sorted(set(fb.columns) - set(fa.columns))}")
        new_cols = sorted(set(fb.columns) - set(fa.columns))
        if new_cols or sorted(set(fa.columns) - set(fb.columns)):
            rec("uns[protein_features] columns", False, f"only new {new_cols}; only ref {sorted(set(fa.columns) - set(fb.columns))}", added=not (set(fa.columns) - set(fb.columns)) and all(c in ADDITIVE_PROTEIN_FEATURE_COLS for c in new_cols))
    for k in sorted((set(A.uns) | set(B.uns)) - {"petrubseq_protein", "protein_features"}):
        if (k in A.uns) != (k in B.uns):
            rec(f"uns[{k}]", False, "present only in new (expected additive slot)" if k in B.uns and k in ADDITIVE_UNS_TOP else "present in only one object", added=k in B.uns and k in ADDITIVE_UNS_TOP)
    allow = DEFAULT_ALLOW_UNS + a.allow_uns
    ua, ub = _flatten({"petrubseq_protein": A.uns.get("petrubseq_protein", {})}), _flatten({"petrubseq_protein": B.uns.get("petrubseq_protein", {})})
    diffs = []
    for k in sorted(set(ua) | set(ub)):
        if any(k.startswith(p) for p in allow):
            continue
        va, vb = ua.get(k, "<absent>"), ub.get(k, "<absent>")
        same = (json.dumps(va, sort_keys=True, default=str) == json.dumps(vb, sort_keys=True, default=str)) if not isinstance(va, np.ndarray) else np.array_equal(va, vb)
        if not same:
            diffs.append(k)
    rec("uns[petrubseq_protein] (non-whitelisted keys)", not diffs, f"{len(diffs)} differing: {diffs[:8]}", soft=True)
    w = max(len(r[0]) for r in rows)
    for slot, st, det in rows:
        if st != "SAME" or "--all" in sys.argv:
            print(f"{st:4s}  {slot:{w}s}  {det}")
    n_same = sum(r[1] == "SAME" for r in rows)
    print(f"\n{len(rows)} slots compared: {n_same} identical/close, {sum(r[1]=='DIFF' for r in rows)} DIFF, {sum(r[1]=='ADD' for r in rows)} ADD (expected additive slots), {sum(r[1]=='NOTE' for r in rows)} NOTE (whitelisted provenance differences)")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
