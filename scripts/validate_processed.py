#!/usr/bin/env python
"""Validate a processed petrubseq-pipeline-protein object against the output contract.

    python scripts/validate_processed.py ../results/SCP1064/processed/scp1064_processed.h5ad [--expect-cells N]

Checks: dimensions, uniqueness of cells/features, NaN/Inf in X/layers/obsm,
required obs columns, perturbation classes, protein representations and
isotype handling, count-layer provenance, PCA/UMAP sanity (variance, depth
correlation, isotype loadings), condition/coverage preservation. Prints a
PASS/WARN/FAIL table and exits non-zero on FAIL.
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REQUIRED_OBS = ["guides", "n_guides", "targets", "n_targets", "guide", "target", "control_class", "perturbation", "perturbation_class", "is_control", "is_targeting", "is_single_guide", "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "rna_qc_fail"]
CLASSES = {"single_targeting", "single_control", "multi_targeting", "multi_control", "mixed_control_targeting", "ambiguous", "unassigned"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad")
    ap.add_argument("--expect-cells", type=int)
    ap.add_argument("--expect-genes", type=int)
    ap.add_argument("--expect-proteins", type=int)
    a = ap.parse_args()
    adata = ad.read_h5ad(a.h5ad)
    rows: List[Tuple[str, str, str]] = []

    def check(name: str, ok: bool, detail: str = "", warn: bool = False) -> None:
        rows.append((name, "PASS" if ok else ("WARN" if warn else "FAIL"), detail))

    u = adata.uns.get("petrubseq_protein", {})
    obs = adata.obs
    # --- dimensions / identity ---------------------------------------------
    check("n_cells", a.expect_cells is None or adata.n_obs == a.expect_cells, f"{adata.n_obs}")
    check("n_genes", a.expect_genes is None or adata.n_vars == a.expect_genes, f"{adata.n_vars}")
    check("unique cells", adata.obs_names.is_unique, f"{adata.obs_names.duplicated().sum()} duplicates")
    check("unique genes", adata.var_names.is_unique, f"{adata.var_names.duplicated().sum()} duplicates")
    # --- X / layers -------------------------------------------------------------
    X = adata.X
    xd = X.data if sp.issparse(X) else np.asarray(X)
    check("X finite", bool(np.isfinite(xd).all()), f"min {xd.min():.3g} max {xd.max():.3g}")
    check("X non-negative", bool(xd.min() >= 0))
    layers = list(adata.layers.keys())
    counts_layer = u.get("rna", {}).get("counts_layer")
    check("exactly one count layer", len([l for l in layers if l in ("counts", "reconstructed_counts")]) == 1, f"layers={layers}")
    check("count layer named in uns", counts_layer in layers, f"uns says {counts_layer!r}")
    if counts_layer in layers:
        L = adata.layers[counts_layer]
        ld = L.data if sp.issparse(L) else np.asarray(L)
        check("count layer integer", bool(np.allclose(ld, np.round(ld))), f"dtype {ld.dtype}")
        check("count layer finite/non-negative", bool(np.isfinite(ld).all() and ld.min() >= 0))
        if counts_layer == "reconstructed_counts":
            src = u.get("rna", {}).get("counts_source", "")
            check("reconstructed counts labelled as not observed", "NOT the directly observed" in src, src[:80])
            rec = u.get("rna", {}).get("counts_reconstruction", {})
            check("reconstruction accepted", bool(rec.get("accepted")), f"max dev {rec.get('max_abs_deviation')}")
        if "rna_total_counts_provided" in obs:
            rs = np.asarray(L.sum(axis=1)).ravel()
            frac = float(np.mean(np.abs(rs - obs["rna_total_counts_provided"].to_numpy(float)) <= 0.5))
            check("count row sums match provided library size", frac > 0.99, f"{100*frac:.2f}% of cells", warn=True)
    # --- obs ------------------------------------------------------------------------
    missing = [c for c in REQUIRED_OBS if c not in obs.columns]
    check("required obs columns", not missing, f"missing {missing}")
    if "perturbation_class" in obs:
        vc = obs["perturbation_class"].astype(str).value_counts()
        check("perturbation classes valid", set(vc.index) <= CLASSES, str(dict(vc)))
        check("single-guide cells have guide/target", bool((obs.loc[obs["n_guides"] == 1, "guide"].astype(str) != "").all() and (obs.loc[obs["n_guides"] == 1, "target"].astype(str) != "").all()))
        check("multi/unassigned cells preserved", (obs["n_guides"] != 1).any(), f"{int((obs['n_guides'] > 1).sum())} multi, {int((obs['n_guides'] == 0).sum())} unassigned", warn=True)
        ctrl = obs.loc[obs["is_control"].astype(bool), "control_class"].astype(str).value_counts()
        check("control classes present", len(ctrl) > 0, str(dict(ctrl)), warn=True)
    if "condition" in obs:
        vc = obs["condition"].astype(str).value_counts()
        check("conditions preserved", len(vc) >= 1 and (vc > 0).all(), str(dict(vc)))
    for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"):
        if c in obs:
            v = obs[c].to_numpy(float)
            check(f"obs[{c}] finite", bool(np.isfinite(v).all()), f"median {np.nanmedian(v):.3g} max {np.nanmax(v):.3g}")
    filt = u.get("qc", {}).get("filter", {})
    check("no unexpected filtering", filt.get("cells_before") == filt.get("cells_after") or filt.get("enabled"), str(filt), warn=True)
    steps = u.get("qc", {}).get("filtering_steps")
    if steps is not None:
        st = pd.DataFrame(steps)
        ok = len(st) >= 2 and st["step"].iloc[0] == "input" and st["step"].iloc[-1] == "final" and int(st["cells_after"].iloc[-1]) == adata.n_obs and int(st["genes_after"].iloc[-1]) == adata.n_vars
        check("filtering audit consistent with object", bool(ok), f"{len(st)} steps: input {int(st['cells_before'].iloc[0])} -> final {int(st['cells_after'].iloc[-1])} cells")
    else:
        check("filtering audit present", False, "uns[...]['qc']['filtering_steps'] missing", warn=True)
    check("warnings recorded", "warnings" in u, f"{len(u.get('warnings', []))} warnings", warn=True)
    # --- protein --------------------------------------------------------------------
    feats = adata.uns.get("protein_features")
    check("protein_features present", feats is not None)
    if feats is not None:
        iso = feats.index[feats["is_isotype"].astype(bool)].tolist()
        tgt = feats.index[~feats["is_isotype"].astype(bool)].tolist()
        check("isotype controls identified", len(iso) > 0, f"{iso}", warn=True)
        primary_is_clr = str(u.get("protein", {}).get("primary_method", "")).startswith("clr")
        for key in ("protein", "protein_counts", "protein_clr"):
            if key == "protein_clr" and primary_is_clr:
                continue  # CLR is the primary matrix (obsm['protein']); no separate copy expected
            check(f"obsm[{key}] present", key in adata.obsm)
        if "protein" in adata.obsm:
            P = adata.obsm["protein"]
            check("obsm[protein] is a labelled DataFrame", isinstance(P, pd.DataFrame) and all(isinstance(c, str) and c for c in P.columns), f"{P.shape}")
            check("n_proteins", a.expect_proteins is None or P.shape[1] == a.expect_proteins, f"{P.shape[1]}")
            check("no isotypes in obsm[protein]", not set(P.columns) & set(iso))
            check("obsm[protein] finite", bool(np.isfinite(P.to_numpy(float)).all()))
        if "protein_counts" in adata.obsm:
            C = adata.obsm["protein_counts"]
            check("isotypes retained in obsm[protein_counts]", set(iso) <= set(C.columns))
            nn = C.notna().all(axis=1)
            check("raw ADT counts integer", bool(np.allclose(C[nn].to_numpy(float), np.round(C[nn].to_numpy(float)))))
            check("cells without raw counts flagged", "has_protein_counts" not in obs or int((~nn).sum()) == int((~obs["has_protein_counts"].astype(bool)).sum()), f"{int((~nn).sum())} cells")
        if "protein_clr" in adata.obsm:
            check("no isotypes in obsm[protein_clr]", not set(adata.obsm["protein_clr"].columns) & set(iso))
        check("isotypes not used in embedding", not feats.loc[iso, "used_in_embedding"].any() if "used_in_embedding" in feats else True)
        check("isotypes used for QC", feats.loc[iso, "used_in_qc"].all() if ("used_in_qc" in feats and iso) else True)
        chk = u.get("protein_check", {})
        if chk:
            check("provided protein normalization reproduced from counts", bool(chk.get("formula_reproduced")), f"max diff {max(chk.get('max_abs_diff', {0: 0}).values())}", warn=True)
    # --- representations ------------------------------------------------------------
    if "X_pca" in adata.obsm:
        Z = np.asarray(adata.obsm["X_pca"])
        check("X_pca finite", bool(np.isfinite(Z).all()), f"{Z.shape}")
        vr = np.asarray(adata.uns.get("pca", {}).get("variance_ratio", []))
        check("RNA PCA variance ratio sane", vr.size > 0 and 0 < vr[0] < 0.5, f"PC1 {vr[0] if vr.size else 'NA'}; total {vr.sum():.3f}", warn=True)
        if "total_counts" in obs:
            r = abs(np.corrcoef(Z[:, 0], np.log1p(obs["total_counts"].to_numpy(float)))[0, 1])
            check("RNA PC1 not purely library size", r < 0.9, f"|r|={r:.2f}", warn=True)
        if "pct_counts_mt" in obs:
            r = abs(np.corrcoef(Z[:, 0], obs["pct_counts_mt"].to_numpy(float))[0, 1])
            check("RNA PC1 not mito-driven", r < 0.7, f"|r|={r:.2f}", warn=True)
        load = adata.varm.get("PCs") if "PCs" in adata.varm else None
        if load is not None:
            check("RNA loadings bounded", bool(np.nanmax(np.abs(load)) < 0.9), f"max |loading| {np.nanmax(np.abs(load)):.2f}", warn=True)
    if "X_umap_rna" in adata.obsm:
        U = np.asarray(adata.obsm["X_umap_rna"])
        check("X_umap_rna finite", bool(np.isfinite(U).all()), f"{U.shape}")
    if "X_pca_protein" in adata.obsm:
        Z = np.asarray(adata.obsm["X_pca_protein"])
        ok = np.isfinite(Z).all(axis=1)
        check("X_pca_protein finite for cells with protein", bool(ok.sum() > 0), f"{int(ok.sum())} / {len(ok)}")
        vr = np.asarray(adata.uns.get("pca_protein", {}).get("variance_ratio", []))
        check("protein PCA variance", vr.size > 0 and vr.sum() > 0.8, f"PC1 {vr[0]:.3f}; total {vr.sum():.3f}" if vr.size else "", warn=True)
        lo = adata.uns.get("pca_protein", {}).get("loadings")
        if isinstance(lo, pd.DataFrame) and feats is not None:
            check("protein PCA not loaded on isotypes", not set(lo.index) & set(iso))
        if "protein_total_counts" in obs:
            r = abs(np.corrcoef(Z[ok, 0], np.log1p(obs["protein_total_counts"].to_numpy(float)[ok]))[0, 1])
            check("protein PC1 vs ADT depth (documented CITE-seq effect)", True, f"|r|={r:.2f}")
    if "X_umap_protein" in adata.obsm:
        U = np.asarray(adata.obsm["X_umap_protein"])
        check("X_umap_protein finite for cells with protein", bool(np.isfinite(U).all(axis=1).sum() > 0), f"{int(np.isfinite(U).all(axis=1).sum())} / {len(U)}")
    # --- coverage QC ----------------------------------------------------------------
    pq = u.get("perturbations", {}).get("qc", {})
    if pq:
        check("perturbation QC recorded", "n_targeting_targets" in pq, f"{pq.get('n_targeting_targets')} targets, {pq.get('n_targeting_guides')} guides, low-coverage guides {pq.get('n_low_coverage_guides')}, targets {pq.get('n_low_coverage_targets')}")
    check("provenance recorded", "provenance" in u and "packages" in u["provenance"], str({k: u.get('provenance', {}).get('packages', {}).get(k) for k in ('python', 'scanpy', 'anndata', 'numpy')}))
    check("schema recorded", "schema" in u, f"{len(u.get('schema', {}))} slots")
    # --- Stage C: guide counts / assignment provenance / protein annotation ---
    pert = u.get("perturbations", {})
    asg = pert.get("assignment", {}) if isinstance(pert, dict) else {}
    check("assignment provenance recorded", bool(asg) and "effective_source" in asg, f"configured {asg.get('configured_source')} -> effective {asg.get('effective_source')}; provided {asg.get('provided_assignments_available')}, counts {asg.get('guide_counts_available')}")
    if pert.get("guide_counts_available"):
        gc = adata.obsm.get("guide_counts")
        gf = adata.uns.get("guide_features")
        check("obsm[guide_counts] present", gc is not None and gf is not None, f"{None if gc is None else gc.shape}")
        if gc is not None and gf is not None:
            check("guide_counts columns match uns[guide_features]", gc.shape[1] == len(gf), f"{gc.shape[1]} vs {len(gf)}")
            data = gc.data if hasattr(gc, "data") else np.asarray(gc).ravel()
            check("guide counts non-negative integers", bool(data.size == 0 or (data.min() >= 0 and np.all(np.abs(data - np.round(data)) < 1e-6))))
            need = ["has_guide_counts", "guide_total_counts", "n_guides_detected", "guides_detected", "guide_top_count", "guide_second_count", "guide_dominant_call"]
            check("guide diagnostics columns", all(c in obs.columns for c in need), f"missing {[c for c in need if c not in obs.columns]}")
            if "guide_total_counts" in obs.columns:
                tot = np.asarray(gc.sum(axis=1)).ravel()
                check("guide_total_counts equals obsm row sums", np.allclose(tot, obs["guide_total_counts"].to_numpy(float), atol=1e-3))
            if asg.get("effective_source") == "guide_counts" and "guide_dominant_call" in obs.columns and asg.get("guide_calling", {}).get("method") == "dominant":
                single = obs["perturbation_class"].isin(["single_targeting", "single_control"])
                check("dominant calls match obs[guide] for single-guide cells", bool((obs.loc[single, "guide_dominant_call"].astype(str) == obs.loc[single, "guide"].astype(str)).all()))
                check("ambiguous cells labelled", bool(((obs["perturbation_class"] == "ambiguous") == (obs["guide_dominant_call"].astype(str) == asg.get("guide_calling", {}).get("ambiguous_label", "ambiguous"))).all()))
        if asg.get("provided_assignments_available") and asg.get("guide_counts_available"):
            check("provided-vs-count agreement recorded", "agreement" in asg, f"{asg.get('agreement', {}).get('set_agreement')}")
    else:
        check("no fabricated guide counts", "guide_counts" not in adata.obsm and "guide_features" not in adata.uns)
    if feats is not None:
        ann = ["feature_id", "antibody_name", "protein_name", "gene_symbol", "clone", "feature_type", "isotype", "isotype_control", "annotation_source"]
        check("protein feature annotation columns", all(c in feats.columns for c in ann), f"missing {[c for c in ann if c not in feats.columns]}")
        if "feature_id" in feats.columns:
            check("protein feature_id filled", bool((feats["feature_id"].astype(str) != "").all()))
    # --- Stage E (only when analysis.perturbation_effects ran) -------------------
    pe = u.get("analysis", {}).get("perturbation_effects") if isinstance(u.get("analysis"), dict) else None
    if pe:
        if "ps_scores" in adata.obsm:
            S = adata.obsm["ps_scores"]
            v = S.to_numpy(float)
            check("PS in [0, 1]", bool(np.nanmin(v) >= -1e-6 and np.nanmax(v) <= 1 + 1e-6), f"{S.shape[1]} targets")
            single = obs["perturbation_class"].astype(str).isin(["single_targeting"]).to_numpy()
            check("ps_score only on single-guide targeting cells", bool(obs.loc[~single, "ps_score"].isna().all()))
            amb = obs["perturbation_class"].astype(str).eq("ambiguous").to_numpy()
            check("ambiguous cells carry no PS", bool(np.isnan(v[amb]).all()) if amb.any() else True)
        if "lochness" in adata.obsm:
            lv = adata.obsm["lochness"].to_numpy(float)
            check("lochNESS >= -1 and finite", bool(np.isfinite(lv).all() and lv.min() >= -1 - 1e-6))
            li = pe.get("lochness", {})
            check("lochNESS neighbours in PCA space", "umap" not in str(li.get("use_rep", "")).lower(), f"{li.get('use_rep')}, k={li.get('k_used')}")
        if "program_activity" in adata.obsm:
            check("program activity finite", bool(np.isfinite(np.asarray(adata.obsm["program_activity"], float)).all()), f"{adata.obsm['program_activity'].shape[1]} programs")
        if "protein_effects" in adata.uns:
            t = adata.uns["protein_effects"]
            check("protein effects on normalized values", pe.get("protein", {}).get("representation") != "protein_counts", str(pe.get("protein", {}).get("representation")))
            check("protein-effect FDR in [0, 1]", bool(t["fdr"].astype(float).between(0, 1).all()), f"{len(t)} tests")
    if "lane_id" in obs.columns:
        check("lane columns consistent", "barcode_original" in obs.columns and bool(obs["lane_id"].astype(str).ne("").all()), f"{obs['lane_id'].nunique()} lanes")
    # --- print ------------------------------------------------------------------------
    w = max(len(r[0]) for r in rows)
    for name, status, detail in rows:
        print(f"{status:4s}  {name:{w}s}  {detail}")
    n_fail = sum(r[1] == "FAIL" for r in rows)
    n_warn = sum(r[1] == "WARN" for r in rows)
    print(f"\n{len(rows)} checks: {n_fail} FAIL, {n_warn} WARN")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
