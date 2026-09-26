#!/usr/bin/env python
"""Numerical parity: reference perturbseq-pipeline run versus a petrubseq-protein run
on the SAME cells and SAME counts.

    python scripts/reference_parity.py --reference <ref run dir> --current <current run dir> \
        --out <report.md> [--tol 1e-6]

Both run directories must come from the same input matrix. The reference run is
produced by ``perturbseq-pipeline run`` (its ``tables/`` + processed h5ad); the
current run by ``petrubseq-protein run``. The script needs only anndata / numpy /
pandas / scipy / scikit-learn, so it runs in the petrubseq-protein environment.

Compared, in order: retained cells and genes, QC metrics, normalized expression,
PCA (sign-invariant), Leiden labels (ARI / NMI), perturbation strength (effect
sizes, p-values, FDR, hit calls, ranks), PS (per-cell scores, per-target summary,
quadrant classes), lochNESS (per-cell scores, summaries, ranking, top cluster),
modules (effect matrix, perturbation correlation, module / program partitions,
program activity) and cluster enrichment (odds ratios, p-values, FDR, calls).
Exact equality where the quantity is discrete, absolute tolerance where floating
point applies, ARI / NMI for label partitions. Missing outputs are reported as
MISSING rather than silently skipped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ROWS: List[dict] = []


def rec(section: str, item: str, value, status: str, note: str = "") -> None:
    ROWS.append({"section": section, "item": item, "value": value, "status": status, "note": note})


def _read(path: Path, index_col=None) -> Optional[pd.DataFrame]:
    if not path.is_file():
        return None
    return pd.read_csv(path, index_col=index_col)


def _num_diff(section: str, item: str, a: np.ndarray, b: np.ndarray, tol: float, note: str = "") -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() == 0:
        rec(section, item, "n/a", "MISSING", "no finite pairs")
        return np.nan
    d = float(np.max(np.abs(a[ok] - b[ok])))
    nan_mismatch = int((np.isfinite(a) != np.isfinite(b)).sum())
    st = "MATCH" if d <= tol and nan_mismatch == 0 else "DIFF"
    rec(section, item, f"max|diff|={d:.3g} (n={int(ok.sum())}, nan-mismatch={nan_mismatch})", st, note)
    return d


def _exact(section: str, item: str, a, b, note: str = "") -> None:
    a = pd.Series(a).reset_index(drop=True)
    b = pd.Series(b).reset_index(drop=True)
    eq = (a.astype(str) == b.astype(str))
    n_bad = int((~eq).sum())
    rec(section, item, f"{int(eq.sum())}/{len(eq)} equal", "MATCH" if n_bad == 0 else "DIFF", note)


def _partition(section: str, item: str, a: pd.Series, b: pd.Series, note: str = "") -> None:
    common = a.index.intersection(b.index)
    if len(common) == 0:
        rec(section, item, "n/a", "MISSING", "no common items")
        return
    ari = adjusted_rand_score(a.loc[common].astype(str), b.loc[common].astype(str))
    nmi = normalized_mutual_info_score(a.loc[common].astype(str), b.loc[common].astype(str))
    rec(section, item, f"ARI={ari:.4f} NMI={nmi:.4f} (n={len(common)}, k_ref={a.loc[common].nunique()}, k_cur={b.loc[common].nunique()})", "MATCH" if ari > 0.9999 else ("CLOSE" if ari > 0.9 else "DIFF"), note)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--current", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--label", default="")
    a_ = ap.parse_args()
    tol = a_.tol
    R, C = Path(a_.reference), Path(a_.current)
    rt, ct = R / "tables", C / "tables"
    ref_h5 = next(iter(R.glob("*.h5ad")), None)
    cur_h5 = next(iter((C / "processed").glob("*_processed.h5ad")), None) or next(iter((C / "processed").glob("*.h5ad")), None)
    if ref_h5 is None or cur_h5 is None:
        print("processed h5ad missing", ref_h5, cur_h5)
        return 2
    ra = ad.read_h5ad(ref_h5)
    ca = ad.read_h5ad(cur_h5)

    # ---------------------------------------------------------------- QC
    S = "QC"
    rc, cc_ = set(ra.obs_names), set(ca.obs_names)
    rec(S, "retained cells", f"ref={len(rc)} cur={len(cc_)} common={len(rc & cc_)}", "MATCH" if rc == cc_ else "DIFF")
    rg, cg = set(ra.var_names), set(ca.var_names)
    rec(S, "retained genes", f"ref={len(rg)} cur={len(cg)} common={len(rg & cg)}", "MATCH" if rg == cg else "DIFF", "reference re-applies min_cells_per_gene after the cell filters")
    cells = [c for c in ra.obs_names if c in cc_]
    genes = [g for g in ra.var_names if g in cg]
    for col in ("n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"):
        if col in ra.obs and col in ca.obs:
            _num_diff(S, f"obs[{col}]", ra.obs.loc[cells, col], ca.obs.loc[cells, col], 1e-5, "reference sums in float32 (calculate_qc_metrics on X)")
        else:
            rec(S, f"obs[{col}]", "n/a", "MISSING")
    qs_r, qs_c = _read(rt / "qc_steps.csv"), _read(ct / "qc_filtering_steps.csv")
    if qs_r is not None and qs_c is not None:
        rec(S, "filtering steps", f"ref final {qs_r['cells_after'].iloc[-1]} cells / {qs_r['genes_after'].iloc[-1]} genes; cur final {qs_c['cells_after'].iloc[-1]} / {qs_c['genes_after'].iloc[-1]}",
            "MATCH" if (qs_r['cells_after'].iloc[-1], qs_r['genes_after'].iloc[-1]) == (qs_c['cells_after'].iloc[-1], qs_c['genes_after'].iloc[-1]) else "DIFF")
    # classes
    rk = ra.obs["perturbation_class"].astype(str).map({"targeting": "single_targeting", "non-targeting": "single_control"}).fillna(ra.obs["perturbation_class"].astype(str))
    ck = ca.obs["perturbation_class"].astype(str)
    rvc, cvc = rk.loc[cells].value_counts().to_dict(), ck.loc[cells].value_counts().to_dict()
    rec(S, "perturbation classes", f"ref={rvc} cur={cvc}", "MATCH" if rvc == cvc else "DIFF", "reference labels 0<top<min_umi cells ambiguous; current labels them unassigned (both excluded everywhere)")
    tg_r = ra.obs["target_gene"].astype(str).loc[cells]
    tg_c = ca.obs["target"].astype(str).loc[cells]
    targeting = (rk.loc[cells] == "single_targeting").to_numpy()
    _exact(S, "target label of targeting cells", tg_r[targeting].to_numpy(), tg_c[targeting].to_numpy())

    # ---------------------------------------------------------- normalized X
    S = "NORMALIZATION"
    rl = ra.layers["lognorm"] if "lognorm" in ra.layers else ra.X
    ri = ra.obs_names.get_indexer(cells)
    ci = ca.obs_names.get_indexer(cells)
    rgi = ra.var_names.get_indexer(genes)
    cgi = ca.var_names.get_indexer(genes)
    RX = rl[ri][:, rgi]
    CX = ca.X[ci][:, cgi]
    RX = RX.toarray() if sp.issparse(RX) else np.asarray(RX)
    CX = CX.toarray() if sp.issparse(CX) else np.asarray(CX)
    _num_diff(S, "log-normalized expression (common cells x genes)", RX, CX, 1e-5)
    if "counts" in ra.layers and "counts" in ca.layers:
        RC = ra.layers["counts"][ri][:, rgi]
        CC = ca.layers["counts"][ci][:, cgi]
        RC = RC.toarray() if sp.issparse(RC) else np.asarray(RC)
        CC = CC.toarray() if sp.issparse(CC) else np.asarray(CC)
        _num_diff(S, "raw counts", RC, CC, 0)

    # ------------------------------------------------------------------ PCA
    S = "PCA"
    if "X_pca" in ra.obsm and "X_pca" in ca.obsm:
        P1 = np.asarray(ra.obsm["X_pca"])[ri]
        P2 = np.asarray(ca.obsm["X_pca"])[ci]
        k = min(P1.shape[1], P2.shape[1])
        cors = [abs(np.corrcoef(P1[:, j], P2[:, j])[0, 1]) for j in range(k)]
        signs = np.array([np.sign(np.corrcoef(P1[:, j], P2[:, j])[0, 1]) or 1.0 for j in range(k)])
        d = float(np.max(np.abs(P1[:, :k] - P2[:, :k] * signs)))
        rec(S, "PC components (sign-invariant)", f"min|r|={min(cors):.6f} over {k} PCs; max|diff| after sign alignment={d:.3g}", "MATCH" if min(cors) > 0.9999 else ("CLOSE" if min(cors) > 0.99 else "DIFF"))
        if "pca" in ra.uns and "pca" in ca.uns:
            _num_diff(S, "variance ratio", np.asarray(ra.uns["pca"]["variance_ratio"])[:k], np.asarray(ca.uns["pca"]["variance_ratio"])[:k], 1e-6)
    else:
        rec(S, "X_pca", "n/a", "MISSING")
    hv_r = ra.var["highly_variable"].reindex(genes) if "highly_variable" in ra.var else None
    hv_c = ca.var["highly_variable"].reindex(genes) if "highly_variable" in ca.var else None
    if hv_r is not None and hv_c is not None:
        rec(S, "highly variable genes", f"ref={int(hv_r.sum())} cur={int(hv_c.sum())} shared={int((hv_r & hv_c).sum())}", "MATCH" if (hv_r == hv_c).all() else "DIFF")

    # ------------------------------------------------------------- clusters
    S = "CLUSTERING"
    if "leiden" in ra.obs and "leiden" in ca.obs:
        _partition(S, "Leiden labels", ra.obs["leiden"].astype(str).loc[cells], ca.obs["leiden"].astype(str).loc[cells])
        rec(S, "n clusters", f"ref={ra.obs['leiden'].nunique()} cur={ca.obs['leiden'].nunique()}", "MATCH" if ra.obs['leiden'].nunique() == ca.obs['leiden'].nunique() else "DIFF")
    else:
        rec(S, "Leiden labels", "n/a", "MISSING")

    # -------------------------------------------------- perturbation strength
    S = "PERTURBATION_STRENGTH"
    pr = _read(rt / "perturbation_full.csv")
    pc = _read(ct / "perturbation_strength" / "perturbation_full.csv")
    if pr is None or pc is None:
        rec(S, "perturbation_full table", f"ref={'present' if pr is not None else 'absent'} cur={'present' if pc is not None else 'absent'}", "MISSING")
    else:
        pr = pr.set_index("target_gene")
        pc = pc.set_index("target_gene" if "target_gene" in pc.columns else "target")
        common = pr.index.intersection(pc.index)
        rec(S, "targets tested", f"ref={len(pr)} cur={len(pc)} common={len(common)}", "MATCH" if set(pr.index) == set(pc.index) else "DIFF")
        for arm in ("ntc", "other"):
            for col in (f"log2fc_{arm}", f"pct_knockdown_{arm}", f"ks_stat_{arm}", f"ks_pval_{arm}", f"mwu_pval_less_{arm}", f"ks_fdr_{arm}", f"mwu_fdr_{arm}", f"mean_lognorm_perturbed_{arm}", f"mean_lognorm_control_{arm}"):
                if col in pr and col in pc:
                    _num_diff(S, col, pr.loc[common, col], pc.loc[common, col], tol)
                else:
                    rec(S, col, "n/a", "MISSING")
            if f"is_hit_{arm}" in pr and f"is_hit_{arm}" in pc:
                _exact(S, f"is_hit_{arm}", pr.loc[common, f"is_hit_{arm}"].astype(bool), pc.loc[common, f"is_hit_{arm}"].astype(bool))
        if "rank" in pr and "rank" in pc:
            _exact(S, "rank", pr.loc[common, "rank"], pc.loc[common, "rank"])
        _exact(S, "n_perturbed", pr.loc[common, "n_perturbed"], pc.loc[common, "n_perturbed"])
        sk_r, sk_c = _read(rt / "skipped.csv"), _read(ct / "perturbation_strength" / "skipped.csv")
        if sk_r is not None or sk_c is not None:
            rec(S, "skipped targets", f"ref={sorted(sk_r['target_gene']) if sk_r is not None else []} cur={sorted(sk_c.iloc[:, 0]) if sk_c is not None else []}",
                "MATCH" if (sk_r is not None and sk_c is not None and sorted(sk_r['target_gene']) == sorted(sk_c.iloc[:, 0])) else "DIFF")

    # ------------------------------------------------------------------- PS
    S = "PS"
    ps_r = _read(rt / "ps_score.csv")
    ps_c = _read(ct / "perturbation_effects" / "ps_targets.csv")
    if ps_r is None or ps_c is None:
        rec(S, "ps summary", "n/a", "MISSING")
    else:
        ps_r = ps_r.set_index("target_gene")
        ps_c = ps_c.set_index("target_gene" if "target_gene" in ps_c.columns else "target")
        common = ps_r.index.intersection(ps_c.index)
        rec(S, "targets scored", f"ref={len(ps_r)} cur={len(ps_c)} common={len(common)}", "MATCH" if set(ps_r.index) == set(ps_c.index) else "DIFF")
        for col in ("mean_ps", "median_ps", "pct_high_ps", "pct_successful_kd", "pct_escaper", "pct_non_responder", "pct_low_signal", "pct_controls_called_kd", "net_pct_kd", "expression_cut"):
            if col in ps_r and col in ps_c:
                _num_diff(S, f"summary {col}", ps_r.loc[common, col], ps_c.loc[common, col], 1e-6)
            else:
                rec(S, f"summary {col}", "n/a", "MISSING")
        _exact(S, "summary order (top targets)", list(ps_r.index), list(ps_c.index), "reference sorts by pct_successful_kd")
        # per-cell
        if "ps_scores" in ca.obsm:
            M = ca.obsm["ps_scores"]
            M = M if isinstance(M, pd.DataFrame) else pd.DataFrame(np.asarray(M), index=ca.obs_names)
            worst = 0.0
            n_t = 0
            for t in common:
                col = f"ps_{t}"
                if col in ra.obs and t in M.columns:
                    d = _num_diff(S, f"per-cell ps[{t}]", ra.obs.loc[cells, col].astype(float), M.loc[cells, t].astype(float), 1e-6)
                    worst = max(worst, d if np.isfinite(d) else 0)
                    n_t += 1
            ROWS[:] = [r for r in ROWS if not (r["section"] == S and r["item"].startswith("per-cell ps["))]
            rec(S, "per-cell scores (all targets)", f"max|diff|={worst:.3g} over {n_t} targets", "MATCH" if worst <= 1e-6 else "DIFF")
        else:
            rec(S, "per-cell scores", "n/a", "MISSING")
        if "ps_score" in ra.obs and "ps_score" in ca.obs:
            _num_diff(S, "obs ps_score (own target)", ra.obs.loc[cells, "ps_score"], ca.obs.loc[cells, "ps_score"], 1e-6)
        if "ps_quadrant" in ra.obs and "ps_quadrant" in ca.obs:
            norm = lambda s: s.astype(str).str.replace(" ", "_").str.replace("-", "_")  # noqa: E731
            qa, qb = norm(ra.obs.loc[cells, "ps_quadrant"]), norm(ca.obs.loc[cells, "ps_quadrant"])
            m = qa != "not_applicable"
            _exact(S, "ps_quadrant (own target)", qa[m].to_numpy(), qb[m].to_numpy())
        if "X_lda_umap" in ra.obsm:
            rec(S, "LDA embedding present in current", "yes" if "X_lda_umap" in ca.obsm else "no", "MATCH" if "X_lda_umap" in ca.obsm else "MISSING", "UMAP coordinates are not compared numerically (umap-learn is not bit-reproducible across runs); the LDA label and placed-cell set are")
            if "lda_label" in ra.obs and "lda_label" in ca.obs:
                _exact(S, "lda_label", ra.obs.loc[cells, "lda_label"].astype(str).to_numpy(), ca.obs.loc[cells, "lda_label"].astype(str).to_numpy())
        pv_r, pv_c = _read(rt / "ps_vs_perturbation.csv"), _read(ct / "perturbation_effects" / "ps_vs_perturbation.csv")
        rec(S, "ps_vs_perturbation table", f"ref={'present' if pv_r is not None else 'absent'} cur={'present' if pv_c is not None else 'absent'}", "MATCH" if (pv_r is not None and pv_c is not None and len(pv_r) == len(pv_c)) else "MISSING")

    # ------------------------------------------------------------- lochNESS
    S = "LOCHNESS"
    lr = _read(rt / "lochness.csv")
    lc = _read(ct / "perturbation_effects" / "lochness_targets.csv")
    if lr is None or lc is None:
        rec(S, "lochness summary", "n/a", "MISSING")
    else:
        lr = lr.set_index("target_gene")
        lc = lc.set_index("target_gene" if "target_gene" in lc.columns else "target")
        common = lr.index.intersection(lc.index)
        rec(S, "targets scored", f"ref={len(lr)} cur={len(lc)} common={len(common)}", "MATCH" if set(lr.index) == set(lc.index) else "DIFF")
        for col in ("mean_lochness_in_own_cells", "mean_lochness_all_cells", "max_lochness", "overall_fraction_pct"):
            if col in lr and col in lc:
                _num_diff(S, f"summary {col}", lr.loc[common, col], lc.loc[common, col], 1e-6)
        col_c = "pct_cells_enriched" if "pct_cells_enriched" in lc else "pct_all_cells_enriched"
        if "pct_cells_enriched" in lr and col_c in lc:
            _num_diff(S, "summary pct_cells_enriched", lr.loc[common, "pct_cells_enriched"], lc.loc[common, col_c], 1e-6)
        if "top_cluster" in lr and "top_cluster" in lc:
            _exact(S, "top_cluster", lr.loc[common, "top_cluster"].astype(str), lc.loc[common, "top_cluster"].astype(str), "cluster ids compared as labels; see CLUSTERING for label parity")
        else:
            rec(S, "top_cluster", "n/a", "MISSING")
        rho = stats.spearmanr(lr.loc[common, "mean_lochness_in_own_cells"], lc.loc[common, "mean_lochness_in_own_cells"]).statistic if len(common) > 2 else np.nan
        _exact(S, "ranking (order of targets)", list(lr.index), list(lc.index), f"Spearman rho of own-cell means = {rho:.4f}")
        if "lochness" in ca.obsm:
            M = ca.obsm["lochness"]
            M = M if isinstance(M, pd.DataFrame) else pd.DataFrame(np.asarray(M), index=ca.obs_names)
            worst = 0.0
            n_t = 0
            for t in common:
                col = f"lochness_{t}"
                if col in ra.obs and t in M.columns:
                    d = _num_diff(S, f"per-cell lochness[{t}]", ra.obs.loc[cells, col].astype(float), M.loc[cells, t].astype(float), 1e-6)
                    worst = max(worst, d if np.isfinite(d) else 0)
                    n_t += 1
            ROWS[:] = [r for r in ROWS if not (r["section"] == S and r["item"].startswith("per-cell lochness["))]
            rec(S, "per-cell scores (all targets)", f"max|diff|={worst:.3g} over {n_t} targets", "MATCH" if worst <= 1e-6 else "DIFF")
        else:
            rec(S, "per-cell scores", "n/a", "MISSING")
        if "lochness_self" in ra.obs and "lochness_self" in ca.obs:
            _num_diff(S, "obs lochness_self", ra.obs.loc[cells, "lochness_self"], ca.obs.loc[cells, "lochness_self"], 1e-6)
        bc_r, bc_c = _read(rt / "lochness_by_cluster.csv"), _read(ct / "perturbation_effects" / "lochness_by_cluster.csv")
        rec(S, "by_cluster table", f"ref={'present' if bc_r is not None else 'absent'} cur={'present' if bc_c is not None else 'absent'}", "MATCH" if (bc_r is not None and bc_c is not None) else "MISSING")

    # -------------------------------------------------------------- modules
    S = "MODULES"
    er = _read(rt / "effect_matrix.csv")
    ec = _read(ct / "perturbation_effects" / "perturbation_effect_matrix.csv")
    if er is None or ec is None:
        rec(S, "effect matrix", "n/a", "MISSING")
    else:
        er = er.set_index(er.columns[0])
        ec = ec.set_index(ec.columns[0])
        rec(S, "effect matrix shape", f"ref={er.shape} cur={ec.shape}; shared genes={len(er.columns.intersection(ec.columns))} shared perturbations={len(er.index.intersection(ec.index))}",
            "MATCH" if (set(er.columns) == set(ec.columns) and set(er.index) == set(ec.index)) else "DIFF", "gene panel = union of Leiden cluster markers in the reference")
        gi = er.columns.intersection(ec.columns)
        pi = er.index.intersection(ec.index)
        if len(gi) and len(pi):
            _num_diff(S, "log2FC on shared entries", er.loc[pi, gi].to_numpy(), ec.loc[pi, gi].to_numpy(), 1e-5)
            cr = er.loc[pi, gi].T.corr(method="spearman").to_numpy()
            cc2 = ec.loc[pi, gi].T.corr(method="spearman").to_numpy()
            _num_diff(S, "perturbation Spearman correlation (shared panel)", cr, cc2, 1e-5)
        mr, mc = _read(rt / "cofunctional_modules.csv"), _read(ct / "perturbation_effects" / "perturbation_modules.csv")
        if mr is not None and mc is not None:
            _partition(S, "module assignment", mr.set_index("target_gene")["module"], mc.set_index("target_gene" if "target_gene" in mc.columns else "target")["module"])
            _exact(S, "module labels (leaf-order numbering)", mr.set_index("target_gene")["module"].reindex(pi), mc.set_index("target_gene" if "target_gene" in mc.columns else "target")["module"].reindex(pi))
            r_de = mr.set_index("target_gene")["n_de_genes"].reindex(pi)
            c_de = mc.set_index("target_gene" if "target_gene" in mc.columns else "target")["n_de_genes"].reindex(pi)
            _exact(S, "n_de_genes", r_de, c_de)
        gr, gc = _read(rt / "gene_programs.csv"), _read(ct / "perturbation_effects" / "gene_programs.csv")
        if gr is not None and gc is not None:
            _partition(S, "gene program assignment", gr.set_index("gene")["program"], gc.set_index("gene")["program"])
            _exact(S, "program labels (leaf-order numbering)", gr.set_index("gene")["program"].reindex(gi), gc.set_index("gene")["program"].reindex(gi))
        mp_r, mp_c = _read(rt / "module_program_strength.csv"), _read(ct / "perturbation_effects" / "module_program_strength.csv")
        if mp_r is not None and mp_c is not None:
            mp_r = mp_r.set_index(mp_r.columns[0]); mp_c = mp_c.set_index(mp_c.columns[0])
            ii = mp_r.index.intersection(mp_c.index); jj = mp_r.columns.intersection(mp_c.columns)
            if len(ii) and len(jj):
                _num_diff(S, "module x program strength", mp_r.loc[ii, jj].to_numpy(), mp_c.loc[ii, jj].to_numpy(), 1e-5)
        # program activity per cell
        pa_cols = [c for c in ra.obs.columns if c.startswith("program_P") and c.endswith("_score")]
        if pa_cols and "program_activity" in ca.obsm:
            A = ca.obsm["program_activity"]
            A = A if isinstance(A, pd.DataFrame) else pd.DataFrame(np.asarray(A), index=ca.obs_names)
            worst = 0.0
            for col in pa_cols:
                lab = col[len("program_"):-len("_score")]
                if lab in A.columns:
                    d = float(np.nanmax(np.abs(ra.obs.loc[cells, col].astype(float).to_numpy() - A.loc[cells, lab].astype(float).to_numpy())))
                    worst = max(worst, d)
            rec(S, "program activity per cell", f"max|diff|={worst:.3g} over {len(pa_cols)} programs", "MATCH" if worst <= 1e-5 else "DIFF", "score_genes with the same seed; compared by program label")
        else:
            rec(S, "program activity per cell", "n/a", "MISSING")
        pab_r, pab_c = _read(rt / "program_activity_by_cluster.csv"), _read(ct / "perturbation_effects" / "program_activity_by_cluster.csv")
        rec(S, "program activity by cluster table", f"ref={'present' if pab_r is not None else 'absent'} cur={'present' if pab_c is not None else 'absent'}", "MATCH" if (pab_r is not None and pab_c is not None) else "MISSING")
        for name in ("tf_hubs", "tf_edges", "module_connectivity"):
            tr, tc_ = _read(rt / f"{name}.csv"), _read(ct / "perturbation_effects" / f"{name}.csv")
            if tr is None and tc_ is None:
                rec(S, name, "absent in both", "MATCH")
            elif tr is None or tc_ is None:
                rec(S, name, f"ref={'present' if tr is not None else 'absent'} cur={'present' if tc_ is not None else 'absent'}", "MISSING")
            else:
                rec(S, name, f"rows ref={len(tr)} cur={len(tc_)}", "MATCH" if len(tr) == len(tc_) else "DIFF")

    # ------------------------------------------------------------ enrichment
    S = "ENRICHMENT"
    en_r = _read(rt / "enrichment_full.csv")
    en_c = _read(ct / "cell_states" / "perturbation_cluster_enrichment.csv")
    if en_r is None or en_c is None:
        rec(S, "enrichment table", "n/a", "MISSING")
    else:
        tcol = "target_gene" if "target_gene" in en_c.columns else "target"
        if "control" not in en_c.columns:
            rec(S, "control arms", f"ref={sorted(en_r['control'].unique())} cur=single arm", "DIFF")
            en_c = en_c.assign(control="other" if "n_reference_cells" in en_c.columns else "ntc")
        else:
            rec(S, "control arms", f"ref={sorted(en_r['control'].unique())} cur={sorted(en_c['control'].unique())}", "MATCH" if sorted(en_r['control'].unique()) == sorted(en_c['control'].unique()) else "DIFF")
        en_r["cluster"] = en_r["cluster"].astype(str); en_c["cluster"] = en_c["cluster"].astype(str)
        kr = en_r.set_index(["target_gene", "control", "cluster"])
        kc = en_c.set_index([tcol, "control", "cluster"])
        common = kr.index.intersection(kc.index)
        rec(S, "tested pairs", f"ref={len(kr)} cur={len(kc)} common={len(common)}", "MATCH" if len(kr) == len(kc) == len(common) else "DIFF", "cluster labels compared as strings; requires identical Leiden labels")
        if len(common):
            for rcol, ccol in (("log2_odds_ratio", "log2_odds_ratio" if "log2_odds_ratio" in kc else "log2_or_haldane"), ("pval", "pval" if "pval" in kc else "p_value"), ("fdr", "fdr"), ("pct_of_target", "pct_of_target"), ("pct_of_reference", "pct_of_reference")):
                if rcol in kr and ccol in kc:
                    _num_diff(S, rcol, kr.loc[common, rcol], kc.loc[common, ccol], 1e-6)
                else:
                    rec(S, rcol, "n/a", "MISSING")
            _exact(S, "significant", kr.loc[common, "significant"].astype(bool), kc.loc[common, "significant"].astype(bool))
            _exact(S, "direction", kr.loc[common, "direction"], kc.loc[common, "direction"])
            for rcol in ("guides_concordant", "guides_tested"):
                if rcol in kr and rcol in kc:
                    sig = kr.loc[common, "significant"].astype(bool).to_numpy()
                    _exact(S, rcol + " (significant pairs)", kr.loc[common, rcol][sig].fillna(-1), kc.loc[common, rcol][sig].fillna(-1))
        for name, cur in (("enrichment_composition", "enrichment_composition"), ("enrichment_effect_magnitude", "enrichment_effect_magnitude")):
            tr, tc_ = _read(rt / f"{name}.csv"), _read(ct / "cell_states" / f"{cur}.csv")
            if tr is None or tc_ is None:
                rec(S, name, "n/a", "MISSING")
            else:
                tr = tr.set_index(tr.columns[0]); tc_ = tc_.set_index(tc_.columns[0])
                ii = tr.index.intersection(tc_.index)
                jj = [c for c in tr.columns if c in tc_.columns and pd.api.types.is_numeric_dtype(tr[c])]
                if len(ii) and jj:
                    _num_diff(S, name, tr.loc[ii, jj].to_numpy(dtype=float), tc_.loc[ii, jj].to_numpy(dtype=float), 1e-6)

    # ---------------------------------------------------------------- write
    df = pd.DataFrame(ROWS)
    out = Path(a_.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    counts = df["status"].value_counts().to_dict()
    lines = [f"# Reference parity report {a_.label}", "", f"reference: `{R}`", f"current: `{C}`", "", f"status counts: {counts}", "", "| section | item | value | status | note |", "|---|---|---|---|---|"]
    for _, r in df.iterrows():
        lines.append(f"| {r['section']} | {r['item']} | {r['value']} | **{r['status']}** | {r['note']} |")
    out.write_text("\n".join(lines) + "\n")
    df.to_csv(out.with_suffix(".csv"), index=False)
    print("\n".join(lines))
    return 0 if counts.get("DIFF", 0) == 0 and counts.get("MISSING", 0) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
