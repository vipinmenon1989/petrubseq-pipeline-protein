"""Permissive prefilter, strict configured filtering and the filtering audit.

Three things are kept strictly apart:

* **prefilter** (``qc.prefilter``): drop obvious empty droplets and
  never-detected genes *before* QC metrics exist. Steps ``prefilter_*``.
* **strict filtering** (``qc.filter``): configured per-modality cell
  thresholds and the perturbation cell mode, applied after the before-filter
  figures. Steps ``rna_*``, ``protein_*``, ``perturbation_*``.
* **coverage flags** (``perturbation.min_cells_per_*``): never remove anything.

Every step is recorded with consistent before/after/removed arithmetic in
``tables/qc_filtering_steps.csv`` and in ``uns['petrubseq_protein']['qc']['filtering_steps']``.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import Config

logger = logging.getLogger("petrubseq_protein")


@dataclass
class FilterStep:
    step: str
    category: str  # input | prefilter | rna | protein | perturbation | final
    threshold: str
    cells_before: int
    cells_after: int
    cells_removed: int
    genes_before: int
    genes_after: int
    genes_removed: int


class FilterAudit:
    """Ordered record of every filtering step of a run."""

    def __init__(self) -> None:
        self.steps: List[FilterStep] = []

    def record(self, step: str, category: str, threshold: str, before: Tuple[int, int], after: Tuple[int, int]) -> FilterStep:
        rec = FilterStep(step, category, threshold, before[0], after[0], before[0] - after[0], before[1], after[1], before[1] - after[1])
        self.steps.append(rec)
        return rec

    def table(self) -> pd.DataFrame:
        return pd.DataFrame([asdict(s) for s in self.steps])

    def to_records(self) -> List[Dict[str, Any]]:
        return [asdict(s) for s in self.steps]

    def check(self) -> None:
        """Internal consistency: each step starts where the previous ended."""
        for a, b in zip(self.steps, self.steps[1:]):
            if (a.cells_after, a.genes_after) != (b.cells_before, b.genes_before):
                raise AssertionError(f"filter audit inconsistent between {a.step} and {b.step}")
        for s in self.steps:
            assert s.cells_removed == s.cells_before - s.cells_after and s.genes_removed == s.genes_before - s.genes_after


def _n_genes_per_cell(adata: ad.AnnData) -> np.ndarray:
    X = adata.X
    if sp.issparse(X):
        return np.diff(sp.csr_matrix(X).indptr)
    return (np.asarray(X) > 0).sum(axis=1)


def _n_cells_per_gene(adata: ad.AnnData) -> np.ndarray:
    X = adata.X
    if sp.issparse(X):
        return np.diff(sp.csc_matrix(X).indptr)
    return (np.asarray(X) > 0).sum(axis=0)


def prefilter(adata: ad.AnnData, cfg: Config, audit: FilterAudit) -> ad.AnnData:
    """Permissive prefilter on detected genes per cell and cells per gene."""
    pcfg = cfg.qc.prefilter
    shape = (adata.n_obs, adata.n_vars)
    if not pcfg.enabled:
        audit.record("prefilter", "prefilter", "disabled", shape, shape)
        logger.info("prefilter disabled")
        return adata
    if pcfg.min_genes_per_cell > 0:
        keep = _n_genes_per_cell(adata) >= pcfg.min_genes_per_cell
        if not keep.all():
            adata = adata[keep].copy()
        audit.record("prefilter_min_genes", "prefilter", f"genes detected >= {pcfg.min_genes_per_cell}", shape, (adata.n_obs, adata.n_vars))
        shape = (adata.n_obs, adata.n_vars)
    if pcfg.min_cells_per_gene > 0:
        keep = _n_cells_per_gene(adata) >= pcfg.min_cells_per_gene
        if not keep.all():
            adata = adata[:, keep].copy()
        audit.record("prefilter_min_cells_per_gene", "prefilter", f"cells with gene >= {pcfg.min_cells_per_gene}", shape, (adata.n_obs, adata.n_vars))
    if adata.n_obs == 0:
        raise ValueError("prefilter removed every cell; lower qc.prefilter.min_genes_per_cell or check the input")
    logger.info("prefilter: %d cells x %d genes retained", adata.n_obs, adata.n_vars)
    return adata


def strict_filter(adata: ad.AnnData, cfg: Config, audit: FilterAudit) -> Tuple[ad.AnnData, pd.Series]:
    """Apply ``qc.filter`` step by step.

    Returns the filtered object and a Series ``removed_by`` over the *input*
    cells (``''`` for retained cells, else the first step that removed them).
    """
    fcfg = cfg.qc.filter
    obs = adata.obs
    removed_by = pd.Series("", index=adata.obs_names, dtype=object)
    alive = np.ones(adata.n_obs, dtype=bool)
    n_genes = adata.n_vars

    def apply(step: str, category: str, threshold: str, fail_mask: np.ndarray) -> None:
        nonlocal alive
        before = int(alive.sum())
        newly = alive & fail_mask
        removed_by[newly] = step
        alive = alive & ~fail_mask
        audit.record(step, category, threshold, (before, n_genes), (int(alive.sum()), n_genes))

    if not fcfg.enabled:
        audit.record("strict_filter", "filter", "disabled (qc.filter.enabled: false)", (adata.n_obs, n_genes), (adata.n_obs, n_genes))
        logger.info("strict filtering disabled: %d cells kept, flags only", adata.n_obs)
        return adata, removed_by

    r = fcfg.rna
    if r.min_genes is not None and "n_genes_by_counts" in obs:
        apply("rna_min_genes", "rna", f"genes detected >= {r.min_genes}", obs["n_genes_by_counts"].to_numpy(float) < r.min_genes)
    if r.min_counts is not None and "total_counts" in obs:
        apply("rna_min_counts", "rna", f"total counts >= {r.min_counts}", obs["total_counts"].to_numpy(float) < r.min_counts)
    if r.max_pct_mt is not None and "pct_counts_mt" in obs:
        # reference: keep pct_counts_mt < max (strict), so cells at exactly the threshold are removed
        apply("rna_max_pct_mt", "rna", f"% mitochondrial < {r.max_pct_mt}", obs["pct_counts_mt"].to_numpy(float) >= r.max_pct_mt)
    if r.max_pct_hb is not None and "pct_counts_hb" in obs:
        apply("rna_max_pct_hb", "rna", f"% haemoglobin < {r.max_pct_hb}", obs["pct_counts_hb"].to_numpy(float) >= r.max_pct_hb)
    p = fcfg.protein
    has_counts = "protein_total_counts" in obs
    if p.min_total_counts is not None and has_counts:
        v = obs["protein_total_counts"].to_numpy(float)
        apply("protein_min_total_counts", "protein", f"total ADT >= {p.min_total_counts}", np.isfinite(v) & (v < p.min_total_counts))
    if p.min_proteins_detected is not None and "protein_n_detected" in obs:
        v = obs["protein_n_detected"].to_numpy(float)
        apply("protein_min_proteins_detected", "protein", f"proteins detected >= {p.min_proteins_detected}", np.isfinite(v) & (v < p.min_proteins_detected))
    if p.max_pct_isotype is not None and "protein_pct_isotype" in obs:
        v = obs["protein_pct_isotype"].to_numpy(float)
        apply("protein_max_pct_isotype", "protein", f"% isotype <= {p.max_pct_isotype}", np.isfinite(v) & (v > p.max_pct_isotype))
    if p.remove_extreme_counts and "protein_extreme_counts" in obs:
        apply("protein_remove_extreme_counts", "protein", f"total ADT <= {cfg.qc.flags.extreme_fold} x q99", obs["protein_extreme_counts"].to_numpy(bool))
    mode = fcfg.perturbation.cells
    if mode == "assigned" and "n_guides" in obs:
        apply("perturbation_assigned_only", "perturbation", "n_guides >= 1", obs["n_guides"].to_numpy(int) < 1)
    elif mode == "single_guide" and "n_guides" in obs:
        apply("perturbation_single_guide", "perturbation", "n_guides == 1", obs["n_guides"].to_numpy(int) != 1)
    if alive.sum() == 0:
        raise ValueError(f"QC filtering removed every cell ({adata.n_obs}); relax qc.filter thresholds or set qc.filter.enabled: false")
    if not alive.all():
        adata = adata[alive].copy()
    # reference: the gene filter is applied again after the cell filters, so genes
    # that fall below min_cells_per_gene once cells are removed leave the matrix
    if r.min_cells_per_gene:
        before = (adata.n_obs, adata.n_vars)
        keep_g = _n_cells_per_gene(adata) >= r.min_cells_per_gene
        if not keep_g.all():
            adata = adata[:, keep_g].copy()
        audit.record("rna_min_cells_per_gene", "rna", f"cells with gene >= {r.min_cells_per_gene}", before, (adata.n_obs, adata.n_vars))
    logger.info("strict filtering: %d -> %d cells, %d -> %d genes", len(alive), int(alive.sum()), n_genes, adata.n_vars)
    return adata, removed_by


def finalize(audit: FilterAudit, adata: ad.AnnData) -> pd.DataFrame:
    shape = (adata.n_obs, adata.n_vars)
    audit.record("final", "final", "-", shape, shape)
    audit.check()
    return audit.table()


def start(audit: FilterAudit, adata: ad.AnnData) -> None:
    shape = (adata.n_obs, adata.n_vars)
    audit.record("input", "input", "-", shape, shape)
