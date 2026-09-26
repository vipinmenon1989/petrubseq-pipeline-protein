"""Perturbation distance vs control, pairwise perturbation distance space, phenotype modules.

Port of the reference ``perturbseq_pipeline/distance.py`` (weili-lab/perturbseq-pipeline,
commit 1c48f9d; audit: docs/reference/PERTURBATION_DISTANCE_AUDIT.md). The reference does
NOT call pertpy at run time: it implements the statistics itself. The energy distance is
the pertpy ``Edistance`` formula (2 E|X-Y| - E|X-X'| - E|Y-Y'| over the full pairwise
Euclidean distance matrices, diagonal zeros included); the permutation test, PCoA and
phenotype modules are the reference's own. Every constant, seed derivation and sampling
rule below reproduces the reference so that the tables agree numerically on the same cells
(docs/reference/PARITY_RESULTS.md).

Three objects, kept distinct:

* **control distance** (``compute_perturbation_distance``): one row per target with
  >= ``min_cells`` single-guide cells: energy distance between the target's cells and the
  control cells in ``obsm[representation]`` (``X_pca``), a seeded label-permutation
  ``DistanceTest`` p-value ``(1 + #perm >= obs) / (1 + B)``, BH-FDR across targets, and the
  optional secondary MMD (Gaussian RBF, median heuristic);
* **pairwise distance** (``compute_distance_space``): the symmetric target x target energy
  distance matrix (zero diagonal), its nearest phenotypic neighbours per target;
* **phenotype space**: classical MDS / PCoA of that matrix (double centring, ``eigh``,
  strictly positive eigenvalues only) and **phenotype modules** = average-linkage
  hierarchical clustering of the matrix cut into ``n_modules`` groups (reference default:
  ``max(2, min(9, K // 4))`` for K >= 8 targets, else K).

Vocabulary: reference ``target_gene`` -> ``target``; reference control ``ntc`` = the
non-targeting class -> ``single_control`` cells of ``control_classes``; reference
``targeting`` -> ``single_targeting``.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import cdist, squareform

from ..config import Config
from ._common import CLASS_TARGETING, bh_fdr
from .perturbation_strength import CONTROL_NTC, CONTROL_OTHER, control_masks

logger = logging.getLogger(__name__)

DISTANCE_COLUMNS = ["target", "n_cells", "n_control", "energy_distance", "pvalue", "fdr", "significant"]


@dataclass
class DistanceResults:
    """Target-level perturbation distance vs control (reference ``DistanceResults``)."""

    table: pd.DataFrame                      # target, n_cells, n_control, energy_distance, pvalue, [mmd_distance], fdr, significant
    skipped: pd.DataFrame = field(default_factory=pd.DataFrame)
    primary_metric: str = "edistance"
    secondary_metric: Optional[str] = None
    control_used: str = CONTROL_NTC
    n_control_cells: int = 0
    representation: str = "X_pca"
    note: str = ""
    info: Dict[str, Any] = field(default_factory=dict)

    @property
    def empty(self) -> bool:
        return self.table.empty

    @property
    def significant_hits(self) -> pd.DataFrame:
        if self.table.empty or "significant" not in self.table.columns:
            return self.table.iloc[0:0]
        return self.table[self.table["significant"].astype(bool)]


@dataclass
class DistanceSpaceResults:
    """Pairwise distance matrix, PCoA coordinates, neighbours, phenotype modules (reference ``DistanceSpaceResults``)."""

    distance_matrix: pd.DataFrame            # targets x targets, symmetric, zero diagonal
    coordinates: pd.DataFrame                # target, PCoA1..PCoAp
    neighbors: pd.DataFrame                  # target, neighbor, distance, rank
    phenotype_modules: pd.DataFrame          # target, phenotype_module
    skipped: pd.DataFrame = field(default_factory=pd.DataFrame)
    eigenvalues: np.ndarray = field(default_factory=lambda: np.zeros(0))
    n_components: int = 0
    metric: str = "edistance"
    linkage_method: str = "average"
    note: str = ""
    info: Dict[str, Any] = field(default_factory=dict)

    @property
    def empty(self) -> bool:
        return self.distance_matrix.empty

    @property
    def module_labels(self) -> List[str]:
        if self.phenotype_modules.empty:
            return []
        return sorted(self.phenotype_modules["phenotype_module"].unique(), key=lambda s: int(str(s)[2:]) if str(s)[2:].isdigit() else 0)


# ---------------------------------------------------------------------------------------- statistics


def energy_distance(X: np.ndarray, Y: np.ndarray) -> float:
    """Empirical energy distance ``2 E|X-Y| - E|X-X'| - E|Y-Y'|`` (reference ``compute_energy_distance``).

    Means run over the full pairwise Euclidean distance matrices, diagonal zeros included
    (pertpy ``Edistance`` convention); clipped at 0 against rounding.
    """
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    n, m = X.shape[0], Y.shape[0]
    if n == 0 or m == 0:
        return float("nan")
    d_xy = cdist(X, Y, metric="euclidean")
    d_xx = cdist(X, X, metric="euclidean")
    d_yy = cdist(Y, Y, metric="euclidean")
    edist = 2.0 * float(np.mean(d_xy)) - float(np.mean(d_xx)) - float(np.mean(d_yy))
    return float(max(edist, 0.0))


def energy_distance_from_cdist(D: np.ndarray, idx_x: np.ndarray, idx_y: np.ndarray) -> float:
    """Energy distance of the split (idx_x, idx_y) of a pooled pairwise-distance matrix (reference)."""
    n, m = len(idx_x), len(idx_y)
    if n == 0 or m == 0:
        return float("nan")
    s_x = float(np.sum(D[np.ix_(idx_x, idx_x)]))
    s_y = float(np.sum(D[np.ix_(idx_y, idx_y)]))
    s_total = float(np.sum(D))
    between = (s_total - s_x - s_y) / (n * m)
    return float(max(between - s_x / (n * n) - s_y / (m * m), 0.0))


def mmd_rbf(X: np.ndarray, Y: np.ndarray, gamma: Optional[float] = None) -> float:
    """Squared MMD with a Gaussian RBF kernel; ``gamma = 1 / median(positive squared distances)`` when None (reference ``compute_mmd``)."""
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    n, m = X.shape[0], Y.shape[0]
    if n == 0 or m == 0:
        return float("nan")
    d2_xx = cdist(X, X, metric="sqeuclidean")
    d2_yy = cdist(Y, Y, metric="sqeuclidean")
    d2_xy = cdist(X, Y, metric="sqeuclidean")
    if gamma is None:
        combined = np.concatenate([d2_xx.ravel(), d2_yy.ravel(), d2_xy.ravel()])
        pos = combined[combined > 0]
        med = float(np.median(pos)) if pos.size > 0 else 1.0
        gamma = 1.0 / max(med, 1e-6)
    mmd2 = float(np.mean(np.exp(-gamma * d2_xx)) + np.mean(np.exp(-gamma * d2_yy)) - 2.0 * np.mean(np.exp(-gamma * d2_xy)))
    return float(max(mmd2, 0.0))


def derive_seed(base_seed: int, identifier: Any) -> int:
    """Reference ``compute.derive_seed``: sha256 of ``"<seed>_<identifier>"``, first 8 hex digits, mod 2**31 - 1."""
    digest = hashlib.sha256(f"{base_seed}_{identifier}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % (2**31 - 1)


def sample_cell_indices(indices: np.ndarray, max_cells: int, rng: np.random.Generator, strata: Optional[np.ndarray] = None) -> np.ndarray:
    """Reference ``_sample_cell_indices``: at most ``max_cells`` indices, proportional per stratum when given, sorted."""
    indices = np.asarray(indices, dtype=np.int64)
    n = len(indices)
    if n <= max_cells:
        return indices
    if strata is not None:
        strata_values = strata[indices]
        unique_strata, strata_counts = np.unique(strata_values, return_counts=True)
        selected = []
        for s_val, s_count in zip(unique_strata, strata_counts):
            s_idx = indices[strata_values == s_val]
            s_quota = max(1, int(np.round(max_cells * (s_count / n))))
            s_quota = min(s_quota, len(s_idx))
            if s_quota > 0:
                selected.append(s_idx if s_quota == len(s_idx) else rng.choice(s_idx, size=s_quota, replace=False))
        if selected:
            combined = np.concatenate(selected)
            if len(combined) > max_cells:
                combined = rng.choice(combined, size=max_cells, replace=False)
            elif len(combined) < max_cells:
                unpicked = np.setdiff1d(indices, combined)
                topup_n = min(max_cells - len(combined), len(unpicked))
                if topup_n > 0:
                    combined = np.concatenate([combined, rng.choice(unpicked, size=topup_n, replace=False)])
            combined.sort()
            return combined.astype(np.int64)
    selected = rng.choice(indices, size=max_cells, replace=False)
    selected.sort()
    return selected.astype(np.int64)


def distance_test_permutation(X: np.ndarray, Y: np.ndarray, n_permutations: int = 1000, seed: int = 123) -> Tuple[float, float]:
    """Reference ``distance_test_permutation``: observed energy distance and the finite-permutation p-value.

    The pooled N x N distance matrix is computed once; every permutation re-uses its row sums
    (``S_Y = S_total + S_X - 2 R_X``). ``p = (1 + #(perm >= obs - 1e-12)) / (1 + B)``; NaN when B = 0.
    """
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    n, m = X.shape[0], Y.shape[0]
    if n == 0 or m == 0:
        return float("nan"), float("nan")
    Z = np.vstack([X, Y])
    N = n + m
    D = cdist(Z, Z, metric="euclidean")
    row_sums = np.sum(D, axis=1)
    s_total = float(np.sum(row_sums))
    idx_x_obs = np.arange(n, dtype=np.int64)
    r_x_obs = float(np.sum(row_sums[idx_x_obs]))
    s_x_obs = float(np.sum(D[np.ix_(idx_x_obs, idx_x_obs)]))
    s_y_obs = s_total + s_x_obs - 2.0 * r_x_obs
    obs_stat = float(max(2.0 * (r_x_obs - s_x_obs) / (n * m) - s_x_obs / (n * n) - s_y_obs / (m * m), 0.0))
    if n_permutations <= 0:
        return obs_stat, float("nan")
    rng = np.random.default_rng(seed)
    all_indices = np.arange(N, dtype=np.int64)
    count = 0
    tol = obs_stat - 1e-12
    inv_nm, inv_nn, inv_mm = 2.0 / (n * m), 1.0 / (n * n), 1.0 / (m * m)
    use_x = n <= m
    for _ in range(n_permutations):
        perm = rng.permutation(all_indices)
        if use_x:
            px = perm[:n]
            r_k = float(np.sum(row_sums[px]))
            s_x = float(np.sum(D[np.ix_(px, px)]))
            s_y = s_total + s_x - 2.0 * r_k
            stat = (r_k - s_x) * inv_nm - s_x * inv_nn - s_y * inv_mm
        else:
            py = perm[n:]
            r_k = float(np.sum(row_sums[py]))
            s_y = float(np.sum(D[np.ix_(py, py)]))
            s_x = s_total + s_y - 2.0 * r_k
            stat = (r_k - s_y) * inv_nm - s_x * inv_nn - s_y * inv_mm
        if stat >= tol:
            count += 1
    return obs_stat, float((1.0 + count) / (1.0 + n_permutations))


def compute_pcoa(dist_matrix: np.ndarray, n_components: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """Classical MDS / PCoA (reference ``compute_pcoa_coordinates``): ``B = -1/2 H D^2 H``, ``eigh``,
    eigenvalues sorted descending, strictly positive ones (> 1e-10) kept, ``coords = V sqrt(lambda)``.
    Eigenvector signs are arbitrary (compare sign-invariantly)."""
    D = np.asarray(dist_matrix, dtype=np.float64)
    K = D.shape[0]
    if K < 2:
        return np.zeros((K, 0)), np.zeros(0)
    D_sq = D**2
    B = -0.5 * (D_sq - np.mean(D_sq, axis=1, keepdims=True) - np.mean(D_sq, axis=0, keepdims=True) + np.mean(D_sq))
    evals, evecs = np.linalg.eigh(B)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    n_pos = int((evals > 1e-10).sum())
    if n_pos == 0:
        logger.warning("PCoA: no positive eigenvalues in the distance matrix")
        return np.zeros((K, 0)), np.zeros(0)
    p = min(n_components, n_pos)
    return evecs[:, :p] * np.sqrt(evals[:p]), evals[:p]


# ---------------------------------------------------------------------------------------- helpers


def _embedding(adata: ad.AnnData, rep: str) -> Optional[np.ndarray]:
    """``obsm[rep]`` as float32 (reference ``get_embedding`` dtype; the statistics cast to float64)."""
    if rep not in adata.obsm:
        return None
    return np.asarray(adata.obsm[rep], dtype=np.float32)


def _targets(adata: ad.AnnData) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    obs = adata.obs
    targets_col = obs["target"].astype(str).to_numpy()
    targeting = obs["perturbation_class"].astype(str).to_numpy() == CLASS_TARGETING
    return targets_col, targeting, sorted(set(targets_col[targeting]))


def _strata(adata: ad.AnnData, stratify_by: Optional[str]) -> Optional[np.ndarray]:
    if stratify_by and stratify_by in adata.obs.columns:
        return adata.obs[stratify_by].astype(str).to_numpy()
    return None


def _target_worker(payload) -> Tuple[Optional[dict], Optional[dict]]:
    (i, target, pert_indices, emb, Y_ctrl, strata, min_cells, max_cells, n_perm, seed, secondary) = payload
    n_pert = int(pert_indices.size)
    if n_pert < min_cells:
        return None, {"target": target, "n_cells": n_pert, "reason": f"fewer than {min_cells} cells ({n_pert})"}
    target_seed = derive_seed(seed, f"{i}_{target}")
    sampled = sample_cell_indices(pert_indices, max_cells, np.random.default_rng(target_seed), strata=strata)
    X_target = emb[sampled]
    edist, pval = distance_test_permutation(X_target, Y_ctrl, n_permutations=n_perm, seed=target_seed)
    row = {"target": target, "n_cells": n_pert, "n_control": int(len(Y_ctrl)), "energy_distance": edist, "pvalue": pval}
    if secondary == "mmd":
        row["mmd_distance"] = mmd_rbf(X_target, Y_ctrl)
    return row, None


def _run(func, tasks, n_jobs: int):
    """Order-preserving execution; joblib (loky) when ``n_jobs > 1``. Results do not depend on n_jobs."""
    if n_jobs <= 1 or len(tasks) <= 1:
        return [func(t) for t in tasks]
    try:
        import joblib

        return joblib.Parallel(n_jobs=n_jobs, backend="loky", return_as="list")(joblib.delayed(func)(t) for t in tasks)
    except Exception as exc:  # pragma: no cover - joblib is installed with scikit-learn
        logger.warning("parallel execution unavailable (%s); running serially", exc)
        return [func(t) for t in tasks]


# ---------------------------------------------------------------------------------------- control distance


def compute_perturbation_distance(adata: ad.AnnData, cfg: Config) -> DistanceResults:
    """Energy distance of every target vs the control cells with a seeded permutation DistanceTest (reference stage 10)."""
    pe = cfg.analysis.perturbation_effects
    dc = pe.distance
    rep = dc.representation
    cols = DISTANCE_COLUMNS + (["mmd_distance"] if dc.secondary_metric == "mmd" else [])
    base_info = {"representation": rep, "primary_metric": dc.primary_metric, "secondary_metric": dc.secondary_metric, "min_cells": dc.min_cells, "max_cells_per_target": dc.max_cells_per_target,
                 "max_control_cells": dc.max_control_cells, "n_permutations": dc.n_permutations, "random_seed": dc.random_seed, "fdr_threshold": dc.fdr_threshold, "stratify_by": dc.stratify_by,
                 "p_value": "(1 + #(perm >= obs)) / (1 + n_permutations); BH-FDR across tested targets", "reference": "perturbseq_pipeline.distance.compute_perturbation_distance"}
    emb = _embedding(adata, rep)
    if emb is None:
        note = f"representation obsm[{rep!r}] missing"
        logger.warning("perturbation distance: %s", note)
        return DistanceResults(pd.DataFrame(columns=cols), pd.DataFrame(), dc.primary_metric, dc.secondary_metric, representation=rep, note=note, info={**base_info, "status": note})
    targets_col, targeting, all_targets = _targets(adata)
    base = control_masks(adata, pe.control_classes)
    available = {k: int(v.sum()) for k, v in base.items()}
    primary = pe.strength.primary_control
    # reference control choice: primary control if it has >= min_cells cells, else ntc, else other, else any non-empty
    ctrl_choice = None
    for cand, need in ((primary, dc.min_cells), (CONTROL_NTC, dc.min_cells), (CONTROL_OTHER, dc.min_cells), (primary, 1), (CONTROL_NTC, 1), (CONTROL_OTHER, 1)):
        if cand in base and available.get(cand, 0) >= need:
            ctrl_choice = cand
            break
    if ctrl_choice is None:
        note = "no control cells available"
        logger.warning("perturbation distance: %s", note)
        return DistanceResults(pd.DataFrame(columns=cols), pd.DataFrame(), dc.primary_metric, dc.secondary_metric, representation=rep, note=note, info={**base_info, "status": note})
    strata = _strata(adata, dc.stratify_by)
    ctrl_all = np.flatnonzero(base[ctrl_choice])
    ctrl_sampled = sample_cell_indices(ctrl_all, dc.max_control_cells, np.random.default_rng(dc.random_seed), strata=strata)
    Y_ctrl = emb[ctrl_sampled]
    logger.info("perturbation distance: representation %s (%d dims), control %s (%d cells, %d used), %d targets, %d permutations", rep, emb.shape[1], ctrl_choice, len(ctrl_all), len(ctrl_sampled), len(all_targets), dc.n_permutations)
    tasks = [(i, t, np.flatnonzero((targets_col == t) & targeting), emb, Y_ctrl, strata, dc.min_cells, dc.max_cells_per_target, dc.n_permutations, dc.random_seed, dc.secondary_metric) for i, t in enumerate(all_targets)]
    results = _run(_target_worker, tasks, cfg.compute.n_jobs)
    rows = [r for r, s in results if r is not None]
    skipped = pd.DataFrame([s for r, s in results if s is not None], columns=["target", "n_cells", "reason"])
    info = {**base_info, "control_used": ctrl_choice, "n_control_cells_available": int(len(ctrl_all)), "n_control_cells_used": int(len(ctrl_sampled)), "n_targets_skipped": int(len(skipped))}
    if not rows:
        note = f"no target has >= {dc.min_cells} cells"
        logger.warning("perturbation distance: %s", note)
        return DistanceResults(pd.DataFrame(columns=cols), skipped, dc.primary_metric, dc.secondary_metric, ctrl_choice, len(ctrl_sampled), rep, note, {**info, "status": note, "n_targets_tested": 0})
    table = pd.DataFrame(rows)
    table["fdr"] = bh_fdr(table["pvalue"].to_numpy())
    table["significant"] = table["fdr"] < dc.fdr_threshold
    table = table.sort_values("energy_distance", ascending=False).reset_index(drop=True)
    n_sig = int(table["significant"].sum())
    info.update({"status": "ok", "n_targets_tested": int(len(table)), "n_significant": n_sig, "top_target": str(table.iloc[0]["target"]), "top_distance": float(table.iloc[0]["energy_distance"])})
    logger.info("perturbation distance: %d targets tested, %d significant at FDR < %g, %d skipped", len(table), n_sig, dc.fdr_threshold, len(skipped))
    return DistanceResults(table, skipped, dc.primary_metric, dc.secondary_metric, ctrl_choice, int(len(ctrl_sampled)), rep, "", info)


# ---------------------------------------------------------------------------------------- distance space


def _pair_worker(task):
    i, j, X_i, X_j, metric = task
    return i, j, (mmd_rbf(X_i, X_j) if metric == "mmd" else energy_distance(X_i, X_j))


def compute_distance_space(adata: ad.AnnData, cfg: Config) -> DistanceSpaceResults:
    """Pairwise target x target distance matrix, PCoA, nearest neighbours and phenotype modules (reference stage 11)."""
    pe = cfg.analysis.perturbation_effects
    ds = pe.distance_space
    rep = ds.representation
    base_info = {"representation": rep, "metric": ds.metric, "n_components": ds.n_components, "nearest_neighbors": ds.nearest_neighbors, "clustering": ds.clustering, "n_modules_requested": ds.n_modules,
                 "cluster_distance_threshold": ds.cluster_distance_threshold, "linkage_method": ds.linkage_method, "min_cells": ds.min_cells, "max_cells_per_target": ds.max_cells_per_target, "random_seed": ds.random_seed,
                 "stratify_by": ds.stratify_by, "pcoa": "double-centred squared distances, eigh, positive eigenvalues, coords = V sqrt(lambda)", "reference": "perturbseq_pipeline.distance.compute_distance_space"}
    empty = pd.DataFrame()

    def _empty(note: str, skipped: pd.DataFrame) -> DistanceSpaceResults:
        logger.warning("distance space: %s", note)
        return DistanceSpaceResults(empty, empty, empty, pd.DataFrame(columns=["target", "phenotype_module"]), skipped, metric=ds.metric, linkage_method=ds.linkage_method, note=note, info={**base_info, "status": note})

    emb = _embedding(adata, rep)
    if emb is None:
        return _empty(f"representation obsm[{rep!r}] missing", empty)
    targets_col, targeting, all_targets = _targets(adata)
    strata = _strata(adata, ds.stratify_by)
    eligible: List[str] = []
    samples: Dict[str, np.ndarray] = {}
    skipped_rows = []
    for i, t in enumerate(all_targets):  # i enumerates ALL targets (skipped ones included), as the reference does
        pert = np.flatnonzero((targets_col == t) & targeting)
        if pert.size < ds.min_cells:
            skipped_rows.append({"target": t, "n_cells": int(pert.size), "reason": f"fewer than {ds.min_cells} cells ({pert.size})"})
            continue
        target_seed = (ds.random_seed + i * 43) % (2**31 - 1)
        samples[t] = emb[sample_cell_indices(pert, ds.max_cells_per_target, np.random.default_rng(target_seed), strata=strata)]
        eligible.append(t)
    skipped = pd.DataFrame(skipped_rows, columns=["target", "n_cells", "reason"])
    K = len(eligible)
    if K < 2:
        return _empty(f"fewer than 2 targets with >= {ds.min_cells} cells (found {K})", skipped)
    logger.info("distance space: pairwise %s over %d targets (%d pairs)", ds.metric, K, K * (K - 1) // 2)
    tasks = [(i, j, samples[eligible[i]], samples[eligible[j]], ds.metric) for i in range(K) for j in range(i + 1, K)]
    M = np.zeros((K, K), dtype=np.float64)
    for i, j, d in _run(_pair_worker, tasks, cfg.compute.n_jobs):
        M[i, j] = M[j, i] = d
    dist_df = pd.DataFrame(M, index=eligible, columns=eligible)
    coords_arr, evals = compute_pcoa(M, n_components=ds.n_components)
    p = coords_arr.shape[1]
    coords = pd.DataFrame(coords_arr, columns=[f"PCoA{c + 1}" for c in range(p)])
    coords.insert(0, "target", eligible)
    k_nn = min(ds.nearest_neighbors, K - 1)
    nrows = []
    for i, t in enumerate(eligible):
        d = M[i].copy()
        d[i] = np.inf
        order = np.argsort(d)
        for rank in range(1, k_nn + 1):
            j = order[rank - 1]
            nrows.append({"target": t, "neighbor": eligible[j], "distance": float(M[i, j]), "rank": rank})
    neighbors = pd.DataFrame(nrows, columns=["target", "neighbor", "distance", "rank"])
    modules = pd.DataFrame(columns=["target", "phenotype_module"])
    n_mod_used: Optional[int] = None
    rule = "disabled"
    if ds.clustering and K >= 2:
        link = linkage(squareform(M, checks=False), method=ds.linkage_method)
        if ds.n_modules is not None:
            n_mod_used = min(ds.n_modules, K)
            clusters = fcluster(link, t=n_mod_used, criterion="maxclust")
            rule = "n_modules"
        elif ds.cluster_distance_threshold is not None:
            clusters = fcluster(link, t=ds.cluster_distance_threshold, criterion="distance")
            rule = "cluster_distance_threshold"
        else:
            n_mod_used = max(2, min(9, K // 4 if K >= 8 else K))
            clusters = fcluster(link, t=n_mod_used, criterion="maxclust")
            rule = "reference default max(2, min(9, K // 4))"
        modules = pd.DataFrame({"target": eligible, "phenotype_module": [f"PM{c}" for c in clusters]})
    n_modules = int(modules["phenotype_module"].nunique()) if not modules.empty else 0
    info = {**base_info, "status": "ok", "n_targets": K, "n_pairs": K * (K - 1) // 2, "n_pcoa_components": p, "pcoa_eigenvalues": [float(v) for v in evals], "n_neighbors_used": k_nn,
            "n_modules": n_modules, "module_rule": rule, "n_modules_cut": n_mod_used, "n_targets_skipped": int(len(skipped))}
    logger.info("distance space: %d x %d matrix, %d PCoA coordinates, %d neighbours per target, %d phenotype modules", K, K, p, k_nn, n_modules)
    return DistanceSpaceResults(dist_df, coords, neighbors, modules, skipped, evals, p, ds.metric, ds.linkage_method, "", info)
