#!/usr/bin/env python
"""Cross-check Stage E against the weili-lab/perturbseq-pipeline reference (E17).

Run inside the reference environment (it provides ``pertps``), with both code
bases on ``sys.path``:

    conda activate perturbseq-pipeline
    python scripts/compare_reference_perturbation.py --reference <perturbseq-pipeline>/src

The same deterministic synthetic AnnData is translated to the reference
vocabulary (obs['target_gene'], obs['perturbation_class'] targeting /
non-targeting, layers['lognorm']) and passed to

* ``pertps.PerturbAnalyzer.calculate_ps_score``  vs ``analysis.ps_score.ps_for_target``
* reference ``lochness`` (its own neighbour graph + ``lochness_score``) vs ``analysis.lochness``
* reference ``modules.build_effect_matrix`` + ``cluster_axis`` vs ``analysis.modules``
  (on the same gene panel, since the reference panel needs Leiden clusters).

Prints the maximum absolute differences; exit code 0 when every check agrees
within 1e-6 (partitions: identical up to label names).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def synthetic(seed: int = 0):
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    n, g = 900, 300
    targets = ["TA", "TB", "TC", "TD", "TE", "TF"]
    lab = np.array(["NT"] * 240 + sum([[t] * 110 for t in targets], []))
    X = rng.poisson(1.2, (n, g)).astype(float)
    for i, t in enumerate(targets):
        m = lab == t
        X[np.ix_(m, np.arange(10 + i * 20, 30 + i * 20))] += rng.poisson(3.0 if i < 3 else 0.5, (m.sum(), 20))
        X[m, i] = rng.poisson(0.2, m.sum())  # knock the target gene down
    X = np.log1p(X / X.sum(1, keepdims=True) * 1e4).astype(np.float32)
    genes = targets + [f"G{i}" for i in range(len(targets), g)]
    a = ad.AnnData(X=sp.csr_matrix(X), var=pd.DataFrame(index=genes))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["target"] = lab
    a.obs["perturbation_class"] = np.where(lab == "NT", "single_control", "single_targeting")
    a.obs["control_class"] = np.where(lab == "NT", "non_targeting", "")
    u, s, _ = np.linalg.svd(X - X.mean(0), full_matrices=False)
    a.obsm["X_pca"] = (u[:, :30] * s[:30]).astype(np.float32)
    return a, targets


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, help="perturbseq-pipeline/src")
    a_ = ap.parse_args()
    sys.path.insert(0, a_.reference)
    sys.path.insert(0, str(REPO / "src"))
    from pertps import PerturbAnalyzer
    from perturbseq_pipeline import lochness as rloch, modules as rmod
    from perturbseq_pipeline.config import Config as RConfig
    from petrubseq_protein.analysis import lochness as oloch, modules as omod, ps_score as ops
    from petrubseq_protein.config import Config

    adata, targets = synthetic()
    ok = True
    # reference-vocabulary copy
    ref = adata.copy()
    ref.layers["lognorm"] = ref.X.copy()
    ref.obs["target_gene"] = np.where(ref.obs["target"] == "NT", "NTC", ref.obs["target"])
    ref.obs["perturbation_class"] = np.where(ref.obs["target"] == "NT", "non-targeting", "targeting")
    # ---- PS ---------------------------------------------------------------------
    work = ref.copy()
    work.X = ref.layers["lognorm"]
    work.obs["gene"] = pd.Categorical(np.where(ref.obs["target"] == "NT", "Non-Targeting", ref.obs["target"]))
    analyzer = PerturbAnalyzer(work, neg_ctrl="Non-Targeting", scale_factor=3.0)
    ctrl = (adata.obs["target"] == "NT").to_numpy()
    max_ps = 0.0
    for t in targets:
        r = analyzer.calculate_ps_score(t, top_n=100)
        pm = (adata.obs["target"] == t).to_numpy()
        idx = np.flatnonzero(pm | ctrl)
        o = ops.ps_for_target(adata.X[idx], pm[idx], adata.var_names, 100, 3.0)
        mine = pd.Series(o["scores"], index=adata.obs_names[idx])
        d = float(np.max(np.abs(mine.reindex(r.index).to_numpy() - r.to_numpy())))
        max_ps = max(max_ps, d)
    print(f"PS: max |ours - pertps| over {len(targets)} targets = {max_ps:.3g}")
    ok &= max_ps < 1e-6
    # ---- lochNESS ---------------------------------------------------------------
    rcfg = RConfig.from_dict({"lochness": {"n_neighbors": 60, "n_pcs": 20, "min_cells_per_target": 10}, "run": {"seed": 0}})
    graph = rloch._build_neighbor_graph(ref, rcfg)
    radj, rcounts = rloch._adjacency(graph)
    ocfg = Config.from_dict({"inputs": {"rna": {"file": "x"}}, "compute": {"seed": 0}, "analysis": {"perturbation_effects": {"enabled": True, "lochness": {"n_neighbors": 60, "n_pcs": 20, "max_k_fraction": 1.0, "n_permutations": 0}}}})
    ores = oloch.compute_lochness(adata, ocfg)
    overall = pd.Series(ref.obs["target_gene"]).value_counts(normalize=True)
    max_l = 0.0
    for t in targets:
        rs = rloch.lochness_score(radj, rcounts, (ref.obs["target_gene"] == t).to_numpy().astype(float), overall[t])
        max_l = max(max_l, float(np.nanmax(np.abs(rs - ores.scores[t].to_numpy()))))
    print(f"lochNESS: k=60, max |ours - reference| = {max_l:.3g} (graphs: identical adjacency = {bool((radj != oloch.neighbor_adjacency(adata, 'X_pca', 60, 20, 0)).nnz == 0)})")
    ok &= max_l < 1e-6
    # ---- modules ----------------------------------------------------------------
    mcfg = Config.from_dict({"inputs": {"rna": {"file": "x"}}, "analysis": {"perturbation_effects": {"enabled": True, "modules": {"n_programs": 4, "n_modules": 3, "min_perturbations": 5, "log2fc_pseudocount": 1e-9, "min_pct_cells_expressing": 0}}}})
    omr = omod.compute_modules(adata, mcfg)
    panel = list(omr.effect.columns)
    rmcfg = RConfig.from_dict({"modules": {"control": "ntc", "min_cells_per_perturbation": 20}})
    reff, rctrl, _ = rmod.build_effect_matrix(ref, panel, sorted(targets), rmcfg)
    d_eff = float(np.max(np.abs(reff.loc[omr.effect.index, panel].to_numpy() - omr.effect.to_numpy())))
    rprog, _, _ = rmod.cluster_axis(reff.T, "pearson", "average", 4, 0.7, "P")
    rmodl, _, _ = rmod.cluster_axis(reff, "spearman", "average", 3, 0.7, "M")
    oprog = omr.gene_programs.set_index("gene")["program"]
    omodl = omr.perturbation_modules.set_index("target")["module"]
    same = lambda a, b: pd.crosstab(a, b.reindex(a.index)).gt(0).sum(axis=1).eq(1).all() and pd.crosstab(b.reindex(a.index), a).gt(0).sum(axis=1).eq(1).all()  # noqa: E731
    print(f"modules: effect matrix {omr.effect.shape}, max |ours - reference| = {d_eff:.3g}; programs identical partition = {same(rprog, oprog)}; modules identical partition = {same(rmodl, omodl)}")
    ok &= d_eff < 1e-6 and same(rprog, oprog) and same(rmodl, omodl)
    print("ALL AGREE" if ok else "DISAGREEMENT")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
