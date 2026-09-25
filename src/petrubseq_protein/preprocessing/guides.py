"""Guide calling from a guide-count matrix and comparison with provided assignments.

Two rules (``perturbation.assignment.method``):

* ``dominant`` (weili-lab/perturbseq-pipeline rule): a cell is assigned its
  top guide when ``top >= min_umi`` and ``top > dominance_ratio * second``
  (and ``second <= max_second_umi`` when that gate is enabled). A cell whose
  top guide reaches ``min_umi`` but fails the dominance test is ``ambiguous``;
  a cell whose top guide stays below ``min_umi`` is ``unassigned``. (The
  reference labels any non-zero but non-dominant cell ``ambiguous``; here
  ``ambiguous`` requires evidence of at least one real guide.)
* ``threshold``: every guide with ``>= detection_min_umi`` UMIs is assigned,
  giving multi-guide lists for high-MOI screens.

Independently of the rule, the *detected* guide list (``>= detection_min_umi``)
and the dominant call are both recorded per cell, so multi-guide information
is never collapsed into the dominant guide.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import AssignmentConfig, Config
from ..io.readers import Matrix
from .align import classify_target, parse_target

logger = logging.getLogger(__name__)

DIAG_COLUMNS = ["has_guide_counts", "guide_total_counts", "n_guides_detected", "guides_detected", "guide_top", "guide_top_count", "guide_second_count", "guide_dominant_call"]


def top_two(X: sp.csr_matrix, chunk: int = 20_000) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per row: (index of the largest value, largest value, second-largest value)."""
    n, k = X.shape
    top_idx = np.zeros(n, dtype=np.int64)
    top = np.zeros(n, dtype=np.float64)
    second = np.zeros(n, dtype=np.float64)
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        d = X[s:e].toarray().astype(np.float64)
        if k == 1:
            top[s:e] = d[:, 0]
            continue
        part = np.argpartition(d, -2, axis=1)[:, -2:]
        rows = np.arange(d.shape[0])[:, None]
        vals = d[rows, part]
        order = np.argsort(vals, axis=1)
        top_idx[s:e] = part[rows[:, 0], order[:, 1]]
        top[s:e] = vals[rows[:, 0], order[:, 1]]
        second[s:e] = vals[rows[:, 0], order[:, 0]]
    return top_idx, top, second


@dataclass
class GuideCalls:
    #: guide counts aligned to the requested cells (zero rows for cells without counts)
    X: sp.csr_matrix
    guide_ids: np.ndarray
    #: per-cell diagnostics (index = cells), columns ``DIAG_COLUMNS``
    diagnostics: pd.DataFrame
    #: per-cell guide lists from the configured rule (input to harmonize_perturbations)
    lists: pd.Series
    #: cells whose dominant call is ambiguous
    ambiguous: pd.Series
    #: per-guide table
    features: pd.DataFrame
    info: Dict[str, Any]


def call_guides(gm: Matrix, cells: Sequence[str], cfg: Config, guide_target_table: Optional[Dict[str, str]] = None) -> GuideCalls:
    acfg: AssignmentConfig = cfg.perturbation.assignment
    pcfg = cfg.perturbation
    idx = pd.Index(cells)
    pos = pd.Index(gm.cells).get_indexer(idx)
    has = pos >= 0
    n, k = len(idx), gm.shape[1]
    if has.all():
        X = sp.csr_matrix(gm.X[pos], dtype=np.float32)
    else:  # scatter the available rows into a zero matrix via a selection matrix
        rows = np.flatnonzero(has)
        S = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, np.arange(len(rows)))), shape=(n, len(rows)))
        X = sp.csr_matrix(S @ sp.csr_matrix(gm.X[pos[rows]], dtype=np.float32))
    X.sum_duplicates()
    guide_ids = np.asarray(gm.features, dtype=object)
    top_idx, top, second = top_two(X)
    total = np.asarray(X.sum(axis=1)).ravel()
    det_mask = sp.csr_matrix(X >= acfg.detection_min_umi)  # detection_min_umi >= 1, so zeros never pass
    det_mask.eliminate_zeros()
    n_detected = np.asarray(det_mask.sum(axis=1)).ravel().astype(int)
    detected_lists = [list(guide_ids[det_mask.indices[det_mask.indptr[i]:det_mask.indptr[i + 1]]]) for i in range(len(idx))]
    # dominant rule --------------------------------------------------------
    reaches = top >= acfg.min_umi
    dominant = reaches & (top > acfg.dominance_ratio * second)
    if acfg.max_second_umi is not None and acfg.max_second_umi >= 0:
        dominant &= second <= acfg.max_second_umi
    ambiguous = reaches & ~dominant
    call = np.full(len(idx), pcfg.unassigned_label, dtype=object)
    call[ambiguous] = acfg.ambiguous_label
    call[dominant] = guide_ids[top_idx[dominant]]
    top_guide = np.where(top > 0, guide_ids[top_idx], "")
    diag = pd.DataFrame(
        {
            "has_guide_counts": has,
            "guide_total_counts": total.astype(np.float32),
            "n_guides_detected": n_detected,
            "guides_detected": [";".join(l) for l in detected_lists],
            "guide_top": top_guide,
            "guide_top_count": top.astype(np.float32),
            "guide_second_count": second.astype(np.float32),
            "guide_dominant_call": call,
        },
        index=idx,
    )
    if acfg.method == "dominant":
        lists = pd.Series([[g] if d else [] for g, d in zip(call, dominant)], index=idx)
    else:
        lists = pd.Series(detected_lists, index=idx)
    # per-guide table --------------------------------------------------------
    targets = np.asarray([parse_target(g, pcfg.guide_target_regex, guide_target_table) for g in guide_ids], dtype=object)
    classes = np.asarray([classify_target(t, pcfg.control_classes) for t in targets], dtype=object)
    n_dom = pd.Series(call[dominant]).value_counts()
    feats = pd.DataFrame(
        {
            "target": targets,
            "control_class": np.where(classes == "targeting", "", classes),
            "is_control": classes != "targeting",
            "total_umis": np.asarray(X.sum(axis=0)).ravel().astype(np.int64),
            "n_cells_detected": np.asarray(det_mask.sum(axis=0)).ravel().astype(int),
            "n_cells_dominant": [int(n_dom.get(g, 0)) for g in guide_ids],
        },
        index=pd.Index(guide_ids, name="guide"),
    )
    if gm.features_meta is not None:
        for c in gm.features_meta.columns:
            if c not in feats.columns:
                feats[c] = gm.features_meta[c].reindex(guide_ids).to_numpy()
    n_has = int(has.sum())
    info = {
        "method": acfg.method,
        "ambiguous_label": acfg.ambiguous_label,
        "min_umi": acfg.min_umi,
        "dominance_ratio": acfg.dominance_ratio,
        "max_second_umi": acfg.max_second_umi,
        "detection_min_umi": acfg.detection_min_umi,
        "n_guides": int(len(guide_ids)),
        "cells_with_guide_counts": n_has,
        "cells_without_guide_counts": int(len(idx) - n_has),
        "dominant_rule": {"assigned": int(dominant.sum()), "ambiguous": int(ambiguous.sum()), "unassigned": int((~reaches).sum())},
        "threshold_rule": {"cells_with_1_guide": int((n_detected == 1).sum()), "cells_with_2plus_guides": int((n_detected > 1).sum()), "cells_with_0_guides": int((n_detected == 0).sum())},
        "median_guide_umis_per_cell": float(np.median(total[has])) if n_has else float("nan"),
    }
    logger.info("guide calling (%s): %d cells with counts, dominant rule assigned %d / ambiguous %d / unassigned %d; threshold rule: %d single, %d multi, %d none", acfg.method, n_has, *info["dominant_rule"].values(), info["threshold_rule"]["cells_with_1_guide"], info["threshold_rule"]["cells_with_2plus_guides"], info["threshold_rule"]["cells_with_0_guides"])
    return GuideCalls(X, guide_ids, diag, lists, pd.Series(ambiguous, index=idx), feats, info)


def compare_assignments(provided: pd.Series, calls: GuideCalls, cells: Sequence[str]) -> Dict[str, Any]:
    """Agreement between a provided per-cell guide list and the count-derived calls.

    Compares (a) the provided *set* of guides with the detected set (threshold
    rule) and (b) the provided single guide with the dominant call, over cells
    that have guide counts.
    """
    idx = pd.Index(cells)
    prov = provided.reindex(idx).apply(lambda v: v if isinstance(v, list) else [])
    diag = calls.diagnostics
    has = diag["has_guide_counts"].to_numpy()
    det = diag["guides_detected"].map(lambda s: set(s.split(";")) - {""} if s else set())
    prov_sets = prov.map(set)
    same_set = np.asarray([a == b for a, b in zip(prov_sets, det)]) & has
    prov_single = prov.map(len) == 1
    dom = diag["guide_dominant_call"]
    amb = calls.info.get("ambiguous_label", "ambiguous")
    dom_is_guide = dom.isin(calls.guide_ids)
    both_single = prov_single.to_numpy() & dom_is_guide.to_numpy() & has
    agree_single = both_single & (prov.map(lambda l: l[0] if len(l) == 1 else "").to_numpy() == dom.to_numpy())
    prov_has_guides = (prov.map(len) > 0).to_numpy()
    out = {
        "cells_compared": int(has.sum()),
        "set_agreement": {"agree": int(same_set.sum()), "disagree": int((has & ~same_set).sum()), "fraction_agree": float(same_set.sum() / has.sum()) if has.sum() else float("nan")},
        "dominant_vs_provided_single": {"compared": int(both_single.sum()), "agree": int(agree_single.sum()), "disagree": int((both_single & ~agree_single).sum())},
        "provided_assigned_but_counts_ambiguous": int((prov_has_guides & has & (dom == amb).to_numpy()).sum()),
        "provided_assigned_but_counts_unassigned": int((prov_has_guides & has & ~dom_is_guide.to_numpy() & ~(dom == amb).to_numpy()).sum()),
        "provided_unassigned_but_counts_assigned": int((~prov_has_guides & has & dom_is_guide.to_numpy()).sum()),
    }
    return out
