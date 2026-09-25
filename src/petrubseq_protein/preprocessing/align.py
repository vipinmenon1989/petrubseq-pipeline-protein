"""Cell alignment across modalities and perturbation-label harmonization.

Alignment never silently rewrites identifiers: cell IDs are compared exactly
as they appear in each file. Duplicate IDs abort the run (configurable), and
the retained cell set must cover at least ``alignment.min_overlap_fraction``
of the RNA cells.

Perturbation harmonization turns per-cell guide lists into a fixed set of
``obs`` columns that every later stage (and every later dataset) can rely on:

``guides``              all assigned guides, ';'-joined ('' if none)
``n_guides``            number of assigned guides (MOI)
``targets``             unique targets among the guides, ';'-joined
``n_targets``           number of unique targets
``guide``               the guide for single-guide cells, else ''
``target``              target parsed from ``guide`` (regex or explicit table), else ''
``control_class``       control class of the single-guide target ('' for targeting/multi/unassigned)
``perturbation``        target for single-guide cells; 'multi' / 'unassigned' otherwise
``perturbation_class``  single_targeting | single_control | multi_targeting |
                        multi_control | mixed_control_targeting | unassigned
``is_control``          single-guide cell carrying a control guide
``is_targeting``        single-guide cell carrying a gene-targeting guide
``is_single_guide``     n_guides == 1
"""

from __future__ import annotations

import logging
import re
import dataclasses
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..config import Config

logger = logging.getLogger("petrubseq_protein")


class AlignmentError(ValueError):
    pass


# ---------------------------------------------------------------------------
# Cell alignment
# ---------------------------------------------------------------------------


@dataclass
class AlignmentResult:
    cells: List[str]
    counts_before: Dict[str, int]
    duplicates: Dict[str, int]
    only_in: Dict[str, int]
    n_kept: int
    policy: str
    required: List[str]
    optional_missing: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        return {
            "policy": self.policy,
            "required": list(self.required),
            "cells_before": dict(self.counts_before),
            "duplicates": dict(self.duplicates),
            "cells_only_in_one_modality": dict(self.only_in),
            "optional_modality_missing_cells": dict(self.optional_missing),
            "cells_kept": self.n_kept,
        }


def effective_required(required: Sequence[str], id_sets: Dict[str, Sequence[str]]) -> List[str]:
    """Translate the config's ``alignment.required`` to loaded modality names.

    ``protein`` means "a protein representation": it is satisfied by the
    normalized matrix (``protein``) or, when that is absent, by the raw ADT
    counts (``protein_counts``). Everything else is taken literally.
    """
    out: List[str] = []
    for m in required:
        if m == "protein" and "protein" not in id_sets and "protein_counts" in id_sets:
            out.append("protein_counts")
        else:
            out.append(m)
    return out


def align_cells(id_sets: Dict[str, Sequence[str]], cfg: Config, reference: str = "rna", required: Optional[Sequence[str]] = None) -> AlignmentResult:
    """Determine the retained, ordered cell set.

    ``id_sets`` maps modality name -> cell IDs (in file order). The reference
    modality's order is preserved for the output. ``required`` overrides
    ``cfg.alignment.required`` (see :func:`effective_required`).
    """
    acfg = cfg.alignment
    if required is not None:
        acfg = dataclasses.replace(acfg, required=list(required))
    counts_before = {k: len(v) for k, v in id_sets.items()}
    duplicates: Dict[str, int] = {}
    uniq: Dict[str, set] = {}
    for k, v in id_sets.items():
        s = pd.Index(v)
        n_dup = int(s.duplicated().sum())
        duplicates[k] = n_dup
        if n_dup and acfg.fail_on_duplicates:
            ex = s[s.duplicated()][:5].tolist()
            raise AlignmentError(f"Modality '{k}' has {n_dup} duplicate cell IDs (e.g. {ex}); set alignment.fail_on_duplicates=false to override.")
        uniq[k] = set(s)
    required = [m for m in acfg.required if m in id_sets]
    missing_required = [m for m in acfg.required if m not in id_sets]
    if missing_required:
        raise AlignmentError(f"Required modalities not loaded: {missing_required}")
    if reference not in id_sets:
        reference = required[0] if required else next(iter(id_sets))
    if acfg.policy == "intersection":
        keep = set.intersection(*(uniq[m] for m in required)) if required else uniq[reference]
        base_order = id_sets[reference]
    else:  # metadata
        if "metadata" not in id_sets:
            raise AlignmentError("alignment.policy=metadata requires a metadata table")
        keep = uniq["metadata"]
        for m in required:
            lacking = keep - uniq[m]
            if lacking:
                raise AlignmentError(f"alignment.policy=metadata: {len(lacking)} metadata cells missing from '{m}'")
        base_order = id_sets["metadata"]
    seen = set()
    ordered = [c for c in base_order if c in keep and not (c in seen or seen.add(c))]
    only_in = {}
    for k in id_sets:
        others = set.union(*(uniq[o] for o in id_sets if o != k)) if len(id_sets) > 1 else set()
        only_in[k] = len(uniq[k] - others)
    optional_missing = {k: len(keep - uniq[k]) for k in id_sets if k not in required}
    n_ref = len(uniq[reference])
    frac = len(ordered) / max(1, n_ref)
    logger.info("alignment: kept %d cells (%.2f%% of %s)", len(ordered), 100 * frac, reference)
    if frac < acfg.min_overlap_fraction:
        raise AlignmentError(f"Only {len(ordered)} of {n_ref} {reference} cells ({100*frac:.1f}%) are present in all required modalities " f"(alignment.min_overlap_fraction={acfg.min_overlap_fraction}).")
    if not ordered:
        raise AlignmentError("No cells shared across required modalities.")
    return AlignmentResult(ordered, counts_before, duplicates, only_in, len(ordered), acfg.policy, required, optional_missing)


# ---------------------------------------------------------------------------
# Perturbation harmonization
# ---------------------------------------------------------------------------


def parse_target(guide: str, regex: str, table: Optional[Dict[str, str]] = None) -> str:
    if table and guide in table:
        return table[guide]
    m = re.match(regex, guide)
    if m and "target" in m.groupdict() and m.group("target"):
        return m.group("target")
    return guide


def classify_target(target: str, control_classes: Dict[str, List[str]]) -> str:
    for cls, patterns in control_classes.items():
        for pat in patterns:
            if re.search(pat, target, flags=re.IGNORECASE):
                return cls
    return "targeting"


def harmonize_perturbations(
    cells: Sequence[str],
    cfg: Config,
    guide_lists: Optional[pd.Series] = None,
    meta: Optional[pd.DataFrame] = None,
    guide_target_table: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """Build the standard perturbation columns for ``cells``.

    Sources, in order of preference: an explicit per-cell guide list
    (``guide_lists``), otherwise the metadata columns named in ``cfg.columns``.
    """
    pcfg = cfg.perturbation
    idx = pd.Index(cells)
    if guide_lists is not None:
        lists = guide_lists.reindex(idx)
        lists = lists.apply(lambda v: v if isinstance(v, list) else [])
        source = "guide_assignments"
    elif meta is not None and cfg.columns.guide and cfg.columns.guide in meta.columns:
        g = meta[cfg.columns.guide].reindex(idx).fillna("")
        lists = g.map(lambda v: [x.strip() for x in str(v).split(";") if x.strip()] if v else [])
        source = f"metadata:{cfg.columns.guide}"
    elif meta is not None and cfg.columns.perturbation and cfg.columns.perturbation in meta.columns:
        p = meta[cfg.columns.perturbation].reindex(idx).fillna("")
        lists = p.map(lambda v: [str(v)] if v else [])
        source = f"metadata:{cfg.columns.perturbation}"
    else:
        lists = pd.Series([[] for _ in idx], index=idx)
        source = "none"
    n = lists.map(len).astype(int)
    meta_target = None
    if meta is not None and cfg.columns.target and cfg.columns.target in meta.columns:
        meta_target = meta[cfg.columns.target].reindex(idx).fillna("").astype(str)

    def _target(g: str, i: int) -> str:
        if meta_target is not None and n.iloc[i] == 1 and meta_target.iloc[i]:
            return meta_target.iloc[i]
        return parse_target(g, pcfg.guide_target_regex, guide_target_table)

    target_lists = [[_target(g, i) for g in l] for i, l in enumerate(lists)]
    uniq_targets = [list(dict.fromkeys(t)) for t in target_lists]
    class_lists = [[classify_target(t, pcfg.control_classes) for t in t_list] for t_list in uniq_targets]
    guide = lists.map(lambda l: l[0] if len(l) == 1 else "")
    target = pd.Series([t[0] if len(l) == 1 else "" for t, l in zip(target_lists, lists)], index=idx)
    control_class = pd.Series([c[0] if (len(l) == 1 and c[0] != "targeting") else "" for c, l in zip(class_lists, lists)], index=idx)

    def _class(l, cls):
        if len(l) == 0:
            return "unassigned"
        n_ctrl = sum(c != "targeting" for c in cls)
        n_tgt = len(cls) - n_ctrl
        if len(l) == 1:
            return "single_control" if n_ctrl else "single_targeting"
        if n_ctrl and n_tgt:
            return "mixed_control_targeting"
        return "multi_control" if n_ctrl else "multi_targeting"

    pclass = pd.Series([_class(l, c) for l, c in zip(lists, class_lists)], index=idx)
    perturbation = target.where(n == 1, np.where(n == 0, pcfg.unassigned_label, pcfg.multi_guide_label))
    df = pd.DataFrame(
        {
            "guides": lists.map(";".join),
            "n_guides": n,
            "targets": [";".join(t) for t in uniq_targets],
            "n_targets": [len(t) for t in uniq_targets],
            "guide": guide,
            "target": target,
            "control_class": control_class,
            "perturbation": perturbation,
            "perturbation_class": pclass,
            "is_control": pclass == "single_control",
            "is_targeting": pclass == "single_targeting",
            "is_single_guide": n == 1,
        },
        index=idx,
    )
    df.attrs["source"] = source
    df.attrs["control_classes"] = sorted({c for cl in class_lists for c in cl if c != "targeting"})
    # Cross-check against a metadata guide column when both exist.
    if guide_lists is not None and meta is not None and cfg.columns.guide and cfg.columns.guide in meta.columns:
        mg = meta[cfg.columns.guide].reindex(idx).fillna("").astype(str)
        both = (mg != "") & (df["guide"] != "")
        mism = int((mg[both] != df.loc[both, "guide"]).sum())
        df.attrs["metadata_guide_mismatches"] = mism
        df.attrs["metadata_guide_compared"] = int(both.sum())
        if mism:
            logger.warning("%d single-guide cells disagree between guide assignments and metadata '%s'", mism, cfg.columns.guide)
    return df


def build_obs(cells: Sequence[str], cfg: Config, meta: Optional[pd.DataFrame], pert: pd.DataFrame) -> pd.DataFrame:
    """Assemble ``obs``: standard design columns + perturbation columns + kept extras."""
    idx = pd.Index(cells, name="cell_id")
    obs = pd.DataFrame(index=idx)
    std = {
        "condition": cfg.columns.condition,
        "sample": cfg.columns.sample,
        "donor": cfg.columns.donor,
        "batch": cfg.columns.batch,
        "replicate": cfg.columns.replicate,
        "lane": cfg.columns.lane,
    }
    if meta is not None:
        m = meta.reindex(idx)
        for std_name, col in std.items():
            if col and col in m.columns:
                obs[std_name] = m[col].astype(str).values
        if cfg.columns.moi and cfg.columns.moi in m.columns:
            obs["moi_provided"] = pd.to_numeric(m[cfg.columns.moi], errors="coerce").values
        if cfg.columns.rna_total_counts and cfg.columns.rna_total_counts in m.columns:
            obs["rna_total_counts_provided"] = pd.to_numeric(m[cfg.columns.rna_total_counts], errors="coerce").values
        for col in cfg.columns.keep:
            if col in m.columns and col not in obs.columns:
                obs[col] = m[col].values
    for col in pert.columns:
        obs[col] = pert[col].reindex(idx).values
    for col in ("condition", "sample", "donor", "batch", "replicate", "lane", "perturbation", "perturbation_class", "control_class", "target", "guide"):
        if col in obs.columns:
            obs[col] = obs[col].astype("category")
    return obs
