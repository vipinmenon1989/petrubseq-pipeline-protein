"""Dataset audit: inventory files, inspect headers/dimensions/value states and
cross-check cell identifiers, without loading the full RNA matrix.

Two entry points:

* :func:`audit_directory` -- generic, config-free inventory of a data folder.
* :func:`audit_config`    -- role-aware audit driven by a run config (RNA,
  protein, metadata, guides, embedding): cell overlap, perturbation design,
  protein panel, value states.

Both return plain dicts (JSON-serialisable) and :func:`render_markdown` turns
them into a readable summary.
"""

from __future__ import annotations

import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import Config
from .io import readers as pio

logger = logging.getLogger("petrubseq_protein")


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def inspect_file(path: Path, count_rows_max_mb: float = 1000.0, deep: bool = False, sample_values: int = 5) -> Dict[str, Any]:
    """Cheap per-file inspection: size, compression, header, row count (if small)."""
    st = path.stat()
    rec: Dict[str, Any] = {
        "path": str(path),
        "name": path.name,
        "size_bytes": st.st_size,
        "size_human": _human(st.st_size),
        "gzip": pio.is_gzip(path),
    }
    if path.suffix in (".h5", ".h5ad"):
        rec["type"] = "hdf5"
        return rec
    try:
        header = pio.read_header(path)
    except Exception as exc:  # pragma: no cover
        rec["error"] = str(exc)
        return rec
    rec["type"] = "text_table"
    rec["sep"] = "," if pio.guess_sep(path) == "," else "tab"
    rec["n_columns"] = len(header)
    rec["first_column"] = header[0]
    rec["header_preview"] = header[1:6]
    rec["last_column"] = header[-1]
    # Second line preview (row names + first values)
    with pio.open_text(path) as fh:
        fh.readline()
        second = fh.readline().rstrip("\n")
    toks = second.split("," if rec["sep"] == "," else "\t")
    rec["row2_first_column"] = toks[0] if toks else None
    rec["row2_values_preview"] = toks[1 : 1 + sample_values]
    rec["scp_type_row"] = toks[0].strip('"').upper() == "TYPE" if toks else False
    mb = st.st_size / 1e6
    if deep or mb <= count_rows_max_mb:
        n = pio.count_lines(path)
        rec["n_lines"] = n
        rec["n_data_rows"] = n - 1 - (1 if rec["scp_type_row"] else 0)
    else:
        rec["n_lines"] = None
        rec["note"] = f"row count skipped (> {count_rows_max_mb:.0f} MB; use --deep)"
    if rec["gzip"]:
        try:
            with open(path, "rb") as fh:
                fh.seek(-4, 2)
                isize = int.from_bytes(fh.read(4), "little")
            rec["uncompressed_size_mod4GB"] = isize
        except Exception:  # pragma: no cover
            pass
    return rec


def audit_directory(data_dir: str | Path, deep: bool = False) -> Dict[str, Any]:
    root = Path(data_dir)
    files = pio.discover_files(root)
    records = [inspect_file(p, deep=deep) for p in files]
    total = sum(r["size_bytes"] for r in records)
    return {"data_dir": str(root), "n_files": len(records), "total_size_bytes": total, "total_size_human": _human(total), "files": records}


# ---------------------------------------------------------------------------
# Role-aware audit
# ---------------------------------------------------------------------------


def _cells_from_matrix_header(paths: Sequence[Path], orientation: str, sep: str) -> List[str]:
    cells: List[str] = []
    if orientation == "features_x_cells":
        for p in paths:
            cells.extend(pio.read_header(p, sep)[1:])
    else:
        for p in paths:
            with pio.open_text(p) as fh:
                fh.readline()
                s = pio.guess_sep(p, sep)
                for line in fh:
                    cells.append(line.split(s, 1)[0].strip('"'))
    return cells


def id_overlap(sets: Dict[str, Sequence[str]]) -> Dict[str, Any]:
    """Pairwise and global overlap statistics between named cell-ID collections."""
    out: Dict[str, Any] = {"per_set": {}, "pairwise": {}, "all": {}}
    uniq: Dict[str, set] = {}
    for name, ids in sets.items():
        ids = list(ids)
        c = Counter(ids)
        dups = [k for k, v in c.items() if v > 1]
        uniq[name] = set(c)
        out["per_set"][name] = {"total": len(ids), "unique": len(c), "duplicates": len(dups), "duplicate_examples": dups[:5], "examples": ids[:3]}
    names = list(uniq)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            inter = uniq[a] & uniq[b]
            out["pairwise"][f"{a}~{b}"] = {
                "intersection": len(inter),
                f"only_{a}": len(uniq[a] - uniq[b]),
                f"only_{b}": len(uniq[b] - uniq[a]),
                f"pct_of_{a}": round(100 * len(inter) / max(1, len(uniq[a])), 3),
                f"pct_of_{b}": round(100 * len(inter) / max(1, len(uniq[b])), 3),
                f"only_{a}_examples": sorted(uniq[a] - uniq[b])[:5],
                f"only_{b}_examples": sorted(uniq[b] - uniq[a])[:5],
            }
    if names:
        common = set.intersection(*uniq.values())
        union = set.union(*uniq.values())
        out["all"] = {"intersection": len(common), "union": len(union), "pct_of_union": round(100 * len(common) / max(1, len(union)), 3)}
    # Suffix / prefix diagnostics
    pat: Dict[str, Any] = {}
    for name, s in uniq.items():
        sample = list(s)[:2000]
        pat[name] = {
            "ends_with_-1": round(100 * np.mean([x.endswith("-1") for x in sample]), 1),
            "contains_underscore": round(100 * np.mean(["_" in x for x in sample]), 1),
            "looks_like_10x_barcode": round(100 * np.mean([bool(re.match(r"^[ACGT]{14,18}(-\d+)?$", x.split("_")[-1])) for x in sample]), 1),
        }
    out["id_patterns"] = pat
    return out


def _matrix_value_state(paths: Sequence[Path], orientation: str, sep: str, max_cells: int = 2000, chunk_rows: int = 500) -> Tuple[Dict[str, Any], pio.Matrix]:
    """Load a small cell subset of a matrix and report value state + basic stats."""
    m = pio.read_dense_matrix([paths[0]], orientation=orientation, sep=sep, chunk_rows=chunk_rows, max_cells=max_cells)
    vs = pio.detect_value_state(m.X)
    X = m.X
    nnz_per_cell = np.diff(X.indptr)
    rec = vs.to_dict()
    rec.update(
        {
            "sampled_file": str(paths[0]),
            "sampled_cells": int(X.shape[0]),
            "n_features": int(X.shape[1]),
            "zero_fraction": round(1 - X.nnz / (X.shape[0] * X.shape[1]), 4),
            "features_detected_per_cell": {"median": float(np.median(nnz_per_cell)), "min": int(nnz_per_cell.min()), "max": int(nnz_per_cell.max())},
            "dtype_text": "float" if not vs.is_integer else "integer",
        }
    )
    return rec, m


def _guide_stats(guides: pd.Series, cfg: Config, meta: Optional[pd.DataFrame]) -> Dict[str, Any]:
    from .preprocessing.align import classify_target, parse_target

    n_guides = guides.map(len)
    single = guides[n_guides == 1].map(lambda l: l[0])
    all_guides = Counter(g for l in guides for g in l)
    targets_all = Counter(parse_target(g, cfg.perturbation.guide_target_regex) for g in all_guides)
    guides_per_target = Counter(parse_target(g, cfg.perturbation.guide_target_regex) for g in all_guides)
    single_targets = single.map(lambda g: parse_target(g, cfg.perturbation.guide_target_regex))
    cls = single_targets.map(lambda t: classify_target(t, cfg.perturbation.control_classes))
    cells_per_target = single_targets.value_counts()
    cells_per_guide = single.value_counts()
    rec: Dict[str, Any] = {
        "n_cells": int(len(guides)),
        "guides_per_cell_distribution": {str(k): int(v) for k, v in sorted(Counter(n_guides).items())},
        "n_unassigned_cells": int((n_guides == 0).sum()),
        "n_single_guide_cells": int((n_guides == 1).sum()),
        "n_multi_guide_cells": int((n_guides > 1).sum()),
        "n_unique_guides": len(all_guides),
        "n_unique_targets": len(targets_all),
        "guides_per_target": {"median": float(np.median(list(guides_per_target.values()))), "min": int(min(guides_per_target.values())), "max": int(max(guides_per_target.values()))},
        "single_guide_cells_by_class": {k: int(v) for k, v in cls.value_counts().items()},
        "control_targets": sorted({t for t in single_targets.unique() if classify_target(t, cfg.perturbation.control_classes) != "targeting"}),
        "cells_per_target_single_guide": {
            "median": float(cells_per_target.median()),
            "min": int(cells_per_target.min()),
            "max": int(cells_per_target.max()),
            "n_lt_10": int((cells_per_target < 10).sum()),
            "n_lt_20": int((cells_per_target < 20).sum()),
            "n_lt_30": int((cells_per_target < 30).sum()),
            "n_lt_50": int((cells_per_target < 50).sum()),
            "top": {k: int(v) for k, v in cells_per_target.head(10).items()},
        },
        "cells_per_guide_single_guide": {
            "median": float(cells_per_guide.median()),
            "min": int(cells_per_guide.min()),
            "max": int(cells_per_guide.max()),
            "n_lt_10": int((cells_per_guide < 10).sum()),
            "n_lt_20": int((cells_per_guide < 20).sum()),
            "n_lt_30": int((cells_per_guide < 30).sum()),
            "n_lt_50": int((cells_per_guide < 50).sum()),
        },
    }
    cond = cfg.columns.condition
    if meta is not None and cond and cond in meta.columns:
        c = meta[cond].reindex(guides.index)
        ct = pd.crosstab(cls.reindex(guides.index).fillna("none"), c)
        rec["class_by_condition"] = {str(k): {str(kk): int(vv) for kk, vv in row.items()} for k, row in ct.iterrows()}
        tc = pd.crosstab(single_targets, c.reindex(single_targets.index))
        rec["targets_by_condition_min_cells"] = {str(k): int(v) for k, v in tc.min(axis=1).describe().round(1).items()}
    return rec


def _protein_audit(cfg: Config, prot: Optional[pio.Matrix] = None, counts: Optional[pio.Matrix] = None) -> Dict[str, Any]:
    """Protein audit from already-loaded matrices, or from dense files named in ``cfg``."""
    from .preprocessing.normalize import is_isotype, verify_isotype_ratio

    out: Dict[str, Any] = {}
    pin = cfg.inputs.protein
    pcin = cfg.inputs.protein_counts
    if prot is None and pin.paths() and pin.format == "dense_csv":
        prot = pio.read_dense_matrix([cfg.resolve(p) for p in pin.paths()], orientation=pin.orientation, sep=pin.sep, chunk_rows=64)
    if counts is None and pcin.paths() and pcin.format == "dense_csv":
        counts = pio.read_dense_matrix([cfg.resolve(p) for p in pcin.paths()], orientation=pcin.orientation, sep=pcin.sep, chunk_rows=64)
    if prot is not None:
        paths = prot.sources
        vs = pio.detect_value_state(prot.X)
        Xd = prot.X.toarray()
        names = [str(f) for f in prot.features]
        out["normalized_matrix"] = {
            "files": [str(p) for p in paths],
            "n_cells": int(prot.shape[0]),
            "n_features": int(prot.shape[1]),
            "features": names,
            "value_state": vs.to_dict(),
            "per_feature": {
                n: {"zero_frac": round(float((Xd[:, i] == 0).mean()), 4), "median": round(float(np.median(Xd[:, i])), 4), "q90": round(float(np.quantile(Xd[:, i], 0.9)), 4), "max": round(float(Xd[:, i].max()), 4)}
                for i, n in enumerate(names)
            },
        }
    if counts is not None:
        paths = counts.sources
        vs = pio.detect_value_state(counts.X)
        Xd = counts.X.toarray()
        names = [str(f) for f in counts.features]
        iso = [n for n in names if is_isotype(n, cfg.protein.isotype_patterns)]
        tot = Xd.sum(axis=1)
        det = (Xd > 0).sum(axis=1)
        out["counts_matrix"] = {
            "files": [str(p) for p in paths],
            "n_cells": int(counts.shape[0]),
            "n_features": int(counts.shape[1]),
            "features": names,
            "isotype_controls": iso,
            "value_state": vs.to_dict(),
            "total_counts_per_cell": {"median": float(np.median(tot)), "min": float(tot.min()), "max": float(tot.max()), "q05": float(np.quantile(tot, 0.05)), "q95": float(np.quantile(tot, 0.95))},
            "proteins_detected_per_cell": {"median": float(np.median(det)), "min": int(det.min()), "max": int(det.max())},
            "per_feature": {
                n: {"zero_frac": round(float((Xd[:, i] == 0).mean()), 4), "median": float(np.median(Xd[:, i])), "q90": float(np.quantile(Xd[:, i], 0.9)), "max": float(Xd[:, i].max()), "mean": round(float(Xd[:, i].mean()), 2)}
                for i, n in enumerate(names)
            },
        }
        iso_median = {n: float(np.median(Xd[:, names.index(n)])) for n in iso}
        if iso_median:
            ref = max(iso_median.values())
            out["counts_matrix"]["background_dominated_candidates"] = [n for n in names if n not in iso and float(np.median(Xd[:, names.index(n)])) <= ref]
    if prot is not None and counts is not None:
        out["normalization_check"] = verify_isotype_ratio(prot, counts, cfg.protein)
    return out


def audit_config(cfg: Config, deep: bool = False, rna_sample_cells: int = 2000) -> Dict[str, Any]:
    """Role-aware audit of the inputs named in ``cfg``."""
    rec: Dict[str, Any] = {"dataset": cfg.dataset.name, "input_dir": str(Path(cfg.dataset.input_dir).resolve())}
    rec["inventory"] = audit_directory(cfg.dataset.input_dir, deep=deep)
    mp = cfg.inputs.multiplexed
    dense_only = not mp.is_set() and all(getattr(cfg.inputs, r).format == "dense_csv" for r in ("rna", "protein", "protein_counts", "guide_counts"))
    # --- roles / files ---------------------------------------------------
    roles: Dict[str, Any] = {}
    if mp.is_set():
        paths = [cfg.resolve(p) for p in mp.paths()]
        roles["multiplexed"] = {"format": mp.format, "files": [str(p) for p in paths], "exists": [p.exists() for p in paths], "modalities": mp.modalities(), "lanes": [l["id"] for l in mp.lanes]}
    for role in ("rna", "protein", "protein_counts", "guide_counts"):
        m = getattr(cfg.inputs, role)
        paths = [cfg.resolve(p) for p in m.paths()]
        if paths:
            roles[role] = {"format": m.format, "files": [str(p) for p in paths], "exists": [p.exists() for p in paths], "orientation": m.orientation if m.format == "dense_csv" else None, "required": m.required, "lanes": [l["id"] for l in m.lanes]}
    for role in ("metadata", "lane_metadata", "guide_assignments", "embedding"):
        t = getattr(cfg.inputs, role)
        p = cfg.resolve(t.file)
        roles[role] = {"file": str(p) if p else None, "exists": bool(p and p.exists()), "required": getattr(t, "required", False)}
    rec["roles"] = roles
    rec["input_layout"] = "dense text (header-based audit)" if dense_only else "feature-barcode / h5ad (matrices loaded through the format adapters)"
    ci = None
    if not dense_only:
        from .io import adapters

        ci = adapters.load_inputs(cfg)
        rec["input_provenance"] = ci.provenance
    # --- cell ids --------------------------------------------------------
    id_sets: Dict[str, List[str]] = {}
    if ci is not None:
        id_sets.update({k: list(v) for k, v in ci.id_sets().items()})
        meta, guides = ci.metadata, ci.guide_assignments
        if ci.embedding is not None:
            id_sets["embedding"] = ci.embedding.index.tolist()
    else:
        for role in ("rna", "protein", "protein_counts", "guide_counts"):
            m = getattr(cfg.inputs, role)
            paths = [cfg.resolve(p) for p in m.paths()]
            if paths and all(p.exists() for p in paths):
                id_sets[role] = _cells_from_matrix_header(paths, m.orientation, m.sep)
        meta = None
        if roles["metadata"]["exists"]:
            meta = pio.read_table(cfg.resolve(cfg.inputs.metadata.file), fmt=cfg.inputs.metadata.format, sep=cfg.inputs.metadata.sep, index_col=cfg.columns.cell_id)
            id_sets["metadata"] = meta.index.tolist()
        guides = None
        if roles["guide_assignments"]["exists"]:
            g = cfg.inputs.guide_assignments
            guides = pio.read_guide_assignments(cfg.resolve(g.file), g.cell_column, g.guides_column, g.list_separator, g.format, g.sep)
            id_sets["guide_assignments"] = guides.index.tolist()
        if roles["embedding"]["exists"]:
            emb = pio.read_table(cfg.resolve(cfg.inputs.embedding.file), fmt=cfg.inputs.embedding.format, sep=cfg.inputs.embedding.sep, index_col=cfg.columns.cell_id)
            id_sets["embedding"] = emb.index.tolist()
    rec["cell_ids"] = id_overlap(id_sets)
    if ci is None and "rna" in id_sets:
        rec["rna_cells_per_file"] = {str(Path(p).name): len(pio.read_header(Path(p), cfg.inputs.rna.sep)) - 1 for p in roles["rna"]["files"]} if cfg.inputs.rna.orientation == "features_x_cells" else None
    # --- metadata / design ----------------------------------------------
    if meta is not None:
        design: Dict[str, Any] = {"n_cells": int(meta.shape[0]), "columns": list(meta.columns)}
        for col in meta.columns:
            vals = meta[col]
            nun = vals.nunique()
            entry: Dict[str, Any] = {"n_unique": int(nun), "n_empty": int((vals == "").sum())}
            if nun <= 30:
                entry["counts"] = {str(k): int(v) for k, v in vals.value_counts().items()}
            else:
                num = pd.to_numeric(vals.replace("", np.nan), errors="coerce")
                if num.notna().mean() > 0.9:
                    entry["numeric_summary"] = {k: float(v) for k, v in num.describe().items()}
                else:
                    entry["examples"] = vals.unique()[:5].tolist()
            design[col] = entry
        mapped = {k: getattr(cfg.columns, k) for k in ("condition", "sample", "donor", "batch", "replicate", "lane", "moi", "guide", "target", "perturbation", "rna_total_counts")}
        design["mapped_columns"] = mapped
        design["missing_mapped_columns"] = [f"{k}={v}" for k, v in mapped.items() if v and v not in meta.columns]
        rec["design"] = design

    # --- perturbations ---------------------------------------------------
    if guides is not None:
        rec["perturbations"] = _guide_stats(guides, cfg, meta)
    elif meta is not None and cfg.columns.guide and cfg.columns.guide in meta.columns:
        ser = meta[cfg.columns.guide].map(lambda v: [v] if v else [])
        rec["perturbations"] = _guide_stats(ser, cfg, meta)

    # --- RNA value state (sampled) --------------------------------------
    if "rna" in id_sets:
        if ci is not None:
            m = ci.rna.subset_cells(ci.rna.cells[:rna_sample_cells].tolist()) if ci.rna.shape[0] > rna_sample_cells else ci.rna
            rs = pio.detect_value_state(m.X).to_dict()
            rs["sampled_cells"], rs["declared_state"], rs["format"] = int(m.shape[0]), ci.rna.state, ci.rna.format
        else:
            paths = [cfg.resolve(p) for p in cfg.inputs.rna.paths()]
            rs, m = _matrix_value_state(paths, cfg.inputs.rna.orientation, cfg.inputs.rna.sep, max_cells=rna_sample_cells, chunk_rows=cfg.compute.chunk_rows)
        feats = [str(f) for f in m.features]
        if feats is not None:
            rs["n_features_total"] = len(feats)
            rs["n_mito_genes"] = int(sum(f.startswith(cfg.rna.mito_prefix) for f in feats))
            rs["n_ribo_genes"] = int(sum(f.startswith(tuple(cfg.rna.ribo_prefixes)) for f in feats))
            rs["duplicate_features"] = int(len(feats) - len(set(feats)))
        if meta is not None and cfg.columns.rna_total_counts and rs.get("log_scale"):
            rs["counts_reconstruction"] = _check_reconstruction(m, cfg, meta, rs["log_scale"])
        rec["rna"] = rs
    # --- protein ---------------------------------------------------------
    try:
        rec["protein"] = _protein_audit(cfg, ci.protein_normalized if ci else None, ci.protein_counts if ci else None)
    except Exception as exc:  # pragma: no cover
        rec["protein"] = {"error": str(exc)}
    # --- guide counts ------------------------------------------------------
    gm = ci.guide_counts if ci is not None else None
    if gm is None and cfg.inputs.guide_counts.paths() and cfg.inputs.guide_counts.format == "dense_csv":
        gi = cfg.inputs.guide_counts
        gm = pio.read_dense_matrix([cfg.resolve(p) for p in gi.paths()], orientation=gi.orientation, sep=gi.sep, chunk_rows=64)
    if gm is not None:
        from .preprocessing.guides import call_guides

        calls = call_guides(gm, gm.cells.tolist(), cfg)
        rec["guide_counts"] = {"n_cells": int(gm.shape[0]), "n_guides": int(gm.shape[1]), "guides": [str(g) for g in gm.features][:200], **{k: v for k, v in calls.info.items() if k not in ("cells_with_guide_counts", "cells_without_guide_counts")}, "assignment_source": cfg.perturbation.assignment.source, "provided_assignments_available": guides is not None}
    return rec


def _check_reconstruction(m: pio.Matrix, cfg: Config, meta: pd.DataFrame, scale: float) -> Dict[str, Any]:
    from .preprocessing.normalize import reconstruct_counts

    tot = pd.to_numeric(meta[cfg.columns.rna_total_counts].reindex(m.cells), errors="coerce").to_numpy()
    counts, info = reconstruct_counts(m.X, tot, scale, tolerance=cfg.rna.reconstruct_tolerance)
    return info


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _kv_table(d: Dict[str, Any], keys: Optional[Sequence[str]] = None) -> str:
    keys = keys or list(d)
    lines = ["| key | value |", "|---|---|"]
    for k in keys:
        if k in d:
            v = d[k]
            if isinstance(v, (dict, list)):
                v = json.dumps(v, default=str)
            lines.append(f"| {k} | {v} |")
    return "\n".join(lines)


def render_markdown(rec: Dict[str, Any]) -> str:
    out: List[str] = [f"# Dataset audit: {rec.get('dataset', '')}", ""]
    inv = rec.get("inventory", rec)
    out.append(f"Input directory: `{rec.get('input_dir', inv.get('data_dir'))}`  ")
    out.append(f"Files: {inv['n_files']}, total {inv['total_size_human']}")
    out.append("")
    out.append("## File inventory")
    out.append("")
    out.append("| file | size | gzip | columns | rows | first col | 2nd-row preview |")
    out.append("|---|---|---|---|---|---|---|")
    for f in inv["files"]:
        rows = f.get("n_data_rows", f.get("n_lines"))
        prev = ", ".join(f.get("row2_values_preview", [])[:3])
        out.append(f"| `{Path(f['path']).relative_to(inv['data_dir']) if f['path'].startswith(inv['data_dir']) else f['name']}` | {f['size_human']} | {f['gzip']} | {f.get('n_columns','')} | {rows if rows is not None else '(not counted)'} | {f.get('first_column','')} | {prev} |")
    out.append("")
    if "cell_ids" in rec:
        out.append("## Cell identifiers")
        out.append("")
        out.append("| set | total | unique | duplicates | examples |")
        out.append("|---|---|---|---|---|")
        for k, v in rec["cell_ids"]["per_set"].items():
            out.append(f"| {k} | {v['total']} | {v['unique']} | {v['duplicates']} | {', '.join(v['examples'])} |")
        out.append("")
        out.append("Pairwise overlaps:")
        out.append("")
        for k, v in rec["cell_ids"]["pairwise"].items():
            out.append(f"- **{k}**: intersection {v['intersection']}; " + "; ".join(f"{kk}={vv}" for kk, vv in v.items() if kk != "intersection" and not kk.endswith("_examples")))
        a = rec["cell_ids"]["all"]
        out.append(f"- **all sets**: intersection {a.get('intersection')} of union {a.get('union')} ({a.get('pct_of_union')}%)")
        out.append("")
    if "design" in rec:
        out.append("## Metadata / experimental design")
        out.append("")
        d = rec["design"]
        out.append(f"Cells: {d['n_cells']}; columns: {', '.join(d['columns'])}")
        out.append("")
        for col in d["columns"]:
            e = d[col]
            desc = json.dumps(e.get("counts", e.get("numeric_summary", e.get("examples"))), ensure_ascii=False, default=str)
            out.append(f"- `{col}`: {e['n_unique']} unique, {e['n_empty']} empty; {desc}")
        out.append("")
        out.append("Mapped columns: " + json.dumps(d["mapped_columns"]))
        if d["missing_mapped_columns"]:
            out.append(f"\n**Missing mapped columns:** {d['missing_mapped_columns']}")
        out.append("")
    if "perturbations" in rec:
        out.append("## Perturbations")
        out.append("")
        out.append(_kv_table(rec["perturbations"]))
        out.append("")
    if "rna" in rec:
        out.append("## RNA matrix (sampled)")
        out.append("")
        out.append(_kv_table(rec["rna"]))
        out.append("")
    if "protein" in rec:
        out.append("## Protein")
        out.append("")
        for key in ("normalized_matrix", "counts_matrix"):
            if key in rec["protein"]:
                p = rec["protein"][key]
                out.append(f"### {key}")
                out.append("")
                out.append(_kv_table(p, [k for k in p if k != "per_feature"]))
                out.append("")
                out.append("| feature | " + " | ".join(next(iter(p["per_feature"].values())).keys()) + " |")
                out.append("|---|" + "---|" * len(next(iter(p["per_feature"].values()))))
                for n, s in p["per_feature"].items():
                    out.append(f"| {n} | " + " | ".join(str(v) for v in s.values()) + " |")
                out.append("")
        if "normalization_check" in rec["protein"]:
            out.append("### Normalization check (provided normalized vs raw counts)")
            out.append("")
            out.append(_kv_table(rec["protein"]["normalization_check"]))
            out.append("")
    return "\n".join(out)
