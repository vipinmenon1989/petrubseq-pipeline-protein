"""Deterministic synthetic Perturb-CITE-seq fixture in SCP-like file formats.

Writes, under ``outdir``:

* ``rna_a.csv.gz``, ``rna_b.csv.gz``  genes x cells, either raw integer counts
  (``rna_state='raw'``) or ln(TPM+1) (``rna_state='lognorm'``), two chunks
* ``protein_counts.csv.gz``            antibodies (incl. 2 isotypes) x cells, raw counts
* ``protein_norm.csv.gz``              targeting antibodies x cells, isotype-ratio normalized
* ``metadata.csv``                     SCP-style metadata (TYPE row): condition, MOI, sgRNA, UMI_count
* ``guides.txt``                       Cell,sgRNAs (comma list)
* ``umap.csv``                         SCP-style cluster file (X, Y)

Ground truth: two conditions, NO_SITE controls, 4 gene targets x 2 guides,
some multi-guide and unassigned cells, and configurable barcode defects.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

N_GENES = 120
TARGETS = ["GENEA", "GENEB", "GENEC", "GENED"]
PROTEINS = ["CD1", "CD2", "CD3", "CD4", "CD5", "CD6"]
ISOTYPES = ["Mouse_IgG1", "Mouse_IgG2a"]
ISO_MAP = {"CD1": "Mouse_IgG1", "CD2": "Mouse_IgG1", "CD3": "Mouse_IgG2a", "CD4": "Mouse_IgG2a", "CD5": "Mouse_IgG1", "CD6": "Mouse_IgG2a"}


def _write_matrix(path: Path, df: pd.DataFrame) -> None:
    with gzip.open(path, "wt") as fh:
        df.to_csv(fh)


def make_dataset(
    outdir: Path,
    n_cells: int = 300,
    seed: int = 0,
    rna_state: str = "lognorm",
    duplicate_cell: bool = False,
    drop_protein_cells: int = 0,
    protein_norm: bool = True,
) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cells = [f"CELL_{i+1}" for i in range(n_cells)]
    genes = [f"GENE{i:03d}" for i in range(N_GENES)]
    genes[: len(TARGETS)] = TARGETS
    genes[-3:] = ["MT-CO1", "MT-ND1", "RPS1"]
    # design
    condition = np.where(np.arange(n_cells) < n_cells // 2, "Control", "Treated")
    guides = [f"{t}_{k}" for t in TARGETS for k in (1, 2)] + ["NO_SITE_1", "NO_SITE_2"]
    lists: List[List[str]] = []
    for i in range(n_cells):
        r = rng.random()
        if r < 0.08:
            lists.append([])
        elif r < 0.2:
            lists.append(sorted(rng.choice(guides, 2, replace=False), key=guides.index))  # feature order, like count-derived lists
        else:
            lists.append([str(rng.choice(guides))])
    # RNA counts: knock down own target
    base = rng.gamma(2.0, 2.0, size=N_GENES)
    counts = rng.poisson(base[None, :] * rng.uniform(0.5, 1.5, size=(n_cells, 1)))
    for i, l in enumerate(lists):
        if len(l) == 1 and not l[0].startswith("NO_SITE"):
            counts[i, genes.index(l[0].rsplit("_", 1)[0])] = 0
    counts[:, -3] += rng.poisson(5, n_cells)  # some mito
    counts = counts.astype(np.int64)
    total = counts.sum(axis=1)
    total[total == 0] = 1
    if rna_state == "lognorm":
        rna = np.log1p(1e6 * counts / total[:, None]).astype(np.float32)
    else:
        rna = counts
    rna_df = pd.DataFrame(rna.T, index=pd.Index(genes, name="GENE"), columns=cells)
    half = n_cells // 2
    _write_matrix(outdir / "rna_a.csv.gz", rna_df.iloc[:, :half])
    _write_matrix(outdir / "rna_b.csv.gz", rna_df.iloc[:, half:])
    # Protein counts
    pc = rng.poisson(rng.uniform(1, 40, size=(1, len(PROTEINS))), size=(n_cells, len(PROTEINS)))
    iso = rng.poisson(1.0, size=(n_cells, len(ISOTYPES)))
    pcols = PROTEINS[:3] + ISOTYPES + PROTEINS[3:]
    pmat = np.concatenate([pc[:, :3], iso, pc[:, 3:]], axis=1)
    prot_cells = cells[drop_protein_cells:] if drop_protein_cells else cells
    pdf = pd.DataFrame(pmat.T, index=pd.Index(pcols, name=""), columns=cells)[prot_cells]
    if duplicate_cell:
        pdf = pd.concat([pdf, pdf.iloc[:, :1]], axis=1)
    _write_matrix(outdir / "protein_counts.csv.gz", pdf)
    norm = np.zeros((n_cells, len(PROTEINS)), dtype=np.float32)
    for j, p in enumerate(PROTEINS):
        c = pmat[:, pcols.index(p)]
        i = pmat[:, pcols.index(ISO_MAP[p])]
        norm[:, j] = np.maximum(0, np.log((c + 1) / (i + 1)))
    ndf = pd.DataFrame(norm.T, index=pd.Index([f"{p} protein" for p in PROTEINS], name="Protein"), columns=cells)
    if protein_norm:
        _write_matrix(outdir / "protein_norm.csv.gz", ndf)
    # metadata (SCP style)
    meta = pd.DataFrame({"NAME": cells, "condition": condition, "MOI": [len(l) for l in lists], "sgRNA": [l[0] if len(l) == 1 else "" for l in lists], "UMI_count": total.astype(float)})
    with open(outdir / "metadata.csv", "w") as fh:
        fh.write("NAME,condition,MOI,sgRNA,UMI_count\nTYPE,group,numeric,group,numeric\n")
        meta.to_csv(fh, index=False, header=False)
    with open(outdir / "guides.txt", "w") as fh:
        fh.write("Cell,sgRNAs\n")
        for c, l in zip(cells, lists):
            fh.write(f'{c},"{",".join(l)}"\n' if len(l) > 1 else f"{c},{l[0] if l else ''}\n")
    with open(outdir / "umap.csv", "w") as fh:
        fh.write("NAME,X,Y\nTYPE,numeric,numeric\n")
        for c in cells:
            fh.write(f"{c},{rng.normal():.4f},{rng.normal():.4f}\n")
    # Guide UMI counts consistent with ``lists``: assigned guides 8-40 UMIs,
    # background 0-2 UMIs (below the detection threshold of 3), so the
    # threshold rule reproduces ``lists`` exactly; two-guide cells have either
    # one dominant guide (ratio > 2) or two comparable ones (-> ambiguous).
    gc = np.zeros((n_cells, len(guides)), dtype=np.int64)
    noise = rng.random(gc.shape)
    gc[noise > 0.85] = 1
    gc[noise > 0.97] = 2
    for i, l in enumerate(lists):
        for g in l:
            gc[i, guides.index(g)] = int(rng.integers(8, 41))
    gdf = pd.DataFrame(gc.T, index=pd.Index(guides, name="guide"), columns=cells)
    _write_matrix(outdir / "guide_counts.csv.gz", gdf)
    return {"dir": outdir, "cells": cells, "genes": genes, "counts": counts, "rna": rna, "rna_state": rna_state, "lists": lists, "condition": condition, "protein_counts": pmat, "protein_cols": pcols, "protein_norm": norm, "guides": guides, "guide_counts": gc, "total": total}


# ---------------------------------------------------------------------------
# Other serializations of the same experiment (Stage C cross-format tests)
# ---------------------------------------------------------------------------

GENE_IDS = None  # filled lazily


def _gene_ids(genes: List[str]) -> List[str]:
    return [f"ENSG{i:05d}" for i in range(len(genes))]


def _feature_table(d: Dict[str, object], modalities: List[str]) -> pd.DataFrame:
    """10x-style features table (id, name, type) for the requested modalities."""
    rows = []
    genes: List[str] = d["genes"]  # type: ignore[assignment]
    if "rna" in modalities:
        rows += [(i, g, "Gene Expression") for i, g in zip(_gene_ids(genes), genes)]
    if "protein" in modalities:
        rows += [(p, p, "Antibody Capture") for p in d["protein_cols"]]  # type: ignore[union-attr]
    if "guide" in modalities:
        rows += [(g, g.rsplit("_", 1)[0], "CRISPR Guide Capture") for g in d["guides"]]  # type: ignore[union-attr]
    return pd.DataFrame(rows, columns=["id", "name", "type"])


def _stack(d: Dict[str, object], modalities: List[str]) -> np.ndarray:
    blocks = []
    if "rna" in modalities:
        blocks.append(np.asarray(d["counts"]))
    if "protein" in modalities:
        blocks.append(np.asarray(d["protein_counts"]))
    if "guide" in modalities:
        blocks.append(np.asarray(d["guide_counts"]))
    return np.concatenate(blocks, axis=1).astype(np.int64)  # cells x features


def write_mtx(d: Dict[str, object], outdir: Path, modalities: List[str], cells: Optional[List[str]] = None, rows: Optional[np.ndarray] = None) -> Path:
    """Write a 10x MTX directory (gzipped) holding ``modalities`` for ``rows`` cells."""
    import scipy.io
    import scipy.sparse as sp

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    M = _stack(d, modalities)
    if rows is not None:
        M = M[rows]
    cells = list(d["cells"]) if cells is None else cells  # type: ignore[arg-type]
    ft = _feature_table(d, modalities)
    with gzip.open(outdir / "barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(cells) + "\n")
    with gzip.open(outdir / "features.tsv.gz", "wt") as fh:
        ft.to_csv(fh, sep="\t", header=False, index=False)
    with gzip.open(outdir / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, sp.csc_matrix(M.T), field="integer")
    return outdir


def write_10x_h5(d: Dict[str, object], path: Path, modalities: List[str]) -> Path:
    """Write a Cell Ranger v3-style feature-barcode HDF5 file."""
    import h5py
    import scipy.sparse as sp

    M = sp.csr_matrix(_stack(d, modalities))  # cells x features == features x cells CSC
    ft = _feature_table(d, modalities)
    path = Path(path)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=M.data.astype(np.int32))
        g.create_dataset("indices", data=M.indices.astype(np.int64))
        g.create_dataset("indptr", data=M.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array([M.shape[1], M.shape[0]], dtype=np.int32))
        g.create_dataset("barcodes", data=np.array(d["cells"], dtype="S"))
        fg = g.create_group("features")
        fg.create_dataset("id", data=ft["id"].to_numpy().astype("S"))
        fg.create_dataset("name", data=ft["name"].to_numpy().astype("S"))
        fg.create_dataset("feature_type", data=ft["type"].to_numpy().astype("S"))
        fg.create_dataset("genome", data=np.array(["synthetic"] * len(ft), dtype="S"))
        fg.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))
    return path


def write_h5ad_combined(d: Dict[str, object], path: Path) -> Path:
    """One AnnData with every feature type in ``var['feature_types']``."""
    import anndata as ad
    import scipy.sparse as sp

    ft = _feature_table(d, ["rna", "protein", "guide"])
    # gene symbols for RNA (as scanpy's read_10x_* would), feature IDs otherwise
    index = np.where(ft["type"].to_numpy() == "Gene Expression", ft["name"].to_numpy(), ft["id"].to_numpy())
    var = pd.DataFrame({"gene_ids": ft["id"].to_numpy(), "feature_types": ft["type"].to_numpy()}, index=pd.Index(index, name="feature"))
    a = ad.AnnData(X=sp.csr_matrix(_stack(d, ["rna", "protein", "guide"]).astype(np.float32)), var=var)
    a.obs_names = list(d["cells"])  # type: ignore[arg-type]
    a.write_h5ad(path)
    return path


def write_h5ad_slots(d: Dict[str, object], path: Path, guide_obs: bool = True) -> Path:
    """AnnData with RNA counts in ``layers['counts']`` (X = the same counts),
    ADT counts in ``obsm['protein_counts']`` (DataFrame), guide counts in
    ``obsm['guide_counts']`` (array + ``uns['guide_counts_features']``) and the
    provided guide lists in ``obs['guides']`` (``;``-separated)."""
    import anndata as ad
    import scipy.sparse as sp

    genes: List[str] = d["genes"]  # type: ignore[assignment]
    cells: List[str] = d["cells"]  # type: ignore[assignment]
    X = sp.csr_matrix(np.asarray(d["counts"]).astype(np.float32))
    a = ad.AnnData(X=X.copy(), var=pd.DataFrame({"gene_ids": _gene_ids(genes)}, index=pd.Index(genes, name="gene")))
    a.obs_names = cells
    a.layers["counts"] = X
    a.obsm["protein_counts"] = pd.DataFrame(np.asarray(d["protein_counts"]).astype(np.float32), index=a.obs_names, columns=list(d["protein_cols"]))  # type: ignore[arg-type]
    a.obsm["guide_counts"] = sp.csr_matrix(np.asarray(d["guide_counts"]).astype(np.float32))
    a.uns["guide_counts_features"] = list(d["guides"])  # type: ignore[arg-type]
    if guide_obs:
        a.obs["guides"] = [";".join(l) for l in d["lists"]]  # type: ignore[union-attr]
    a.write_h5ad(path)
    return path


def write_lanes(d: Dict[str, object], outdir: Path) -> Dict[str, object]:
    """Two combined MTX lanes whose barcodes repeat across lanes, plus a
    per-lane sample sheet and per-cell metadata / guide lists keyed by the
    global ``<barcode>__<lane>`` IDs."""
    outdir = Path(outdir)
    cells: List[str] = d["cells"]  # type: ignore[assignment]
    n = len(cells)
    half = n // 2
    parts = {"lane1": np.arange(0, half), "lane2": np.arange(half, n)}
    global_ids: List[str] = []
    for lid, rows in parts.items():
        bcs = [f"BC{i:04d}-1" for i in range(len(rows))]
        write_mtx(d, outdir / lid, ["rna", "protein", "guide"], cells=bcs, rows=rows)
        global_ids += [f"{b}__{lid}" for b in bcs]
    with open(outdir / "lane_metadata.csv", "w") as fh:
        fh.write("lane_id,condition,library\nlane1,Control,L1\nlane2,Treated,L2\n")
    meta = pd.DataFrame({"NAME": global_ids, "MOI": [len(l) for l in d["lists"]], "UMI_count": np.asarray(d["total"]).astype(float)})  # type: ignore[union-attr]
    meta.to_csv(outdir / "metadata_lanes.csv", index=False)
    with open(outdir / "guides_lanes.txt", "w") as fh:
        fh.write("Cell,sgRNAs\n")
        for c, l in zip(global_ids, d["lists"]):  # type: ignore[arg-type]
            fh.write(f'{c},"{",".join(l)}"\n' if len(l) > 1 else f"{c},{l[0] if l else ''}\n")
    return {"dir": outdir, "global_ids": global_ids, "lanes": list(parts)}


def write_formats(d: Dict[str, object], outdir: Path) -> Dict[str, Path]:
    """Serialize the experiment in every supported layout; returns the paths."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = {
        "mtx_combined": write_mtx(d, outdir / "mtx_combined", ["rna", "protein", "guide"]),
        "mtx_rna": write_mtx(d, outdir / "mtx_separate" / "rna", ["rna"]),
        "mtx_protein": write_mtx(d, outdir / "mtx_separate" / "protein", ["protein"]),
        "mtx_guide": write_mtx(d, outdir / "mtx_separate" / "guide", ["guide"]),
        "h5_combined": write_10x_h5(d, outdir / "combined.h5", ["rna", "protein", "guide"]),
        "h5ad_combined": write_h5ad_combined(d, outdir / "combined.h5ad"),
        "h5ad_slots": write_h5ad_slots(d, outdir / "slots.h5ad"),
    }
    return out


def base_config(d: Dict[str, object], outdir: Path, **overrides) -> Dict:
    cfg = {
        "dataset": {"name": "synthetic", "input_dir": str(d["dir"])},
        "inputs": {
            "rna": {"files": ["rna_a.csv.gz", "rna_b.csv.gz"]},
            "protein": {"file": "protein_norm.csv.gz", "state": "normalized"},
            "protein_counts": {"file": "protein_counts.csv.gz", "required": False},
            "metadata": {"file": "metadata.csv", "format": "scp_metadata"},
            "guide_assignments": {"file": "guides.txt"},
            "embedding": {"file": "umap.csv", "format": "scp_cluster", "key": "X_umap_provided"},
        },
        "columns": {"cell_id": "NAME", "condition": "condition", "guide": "sgRNA", "moi": "MOI", "rna_total_counts": "UMI_count"},
        "perturbation": {"control_classes": {"non_targeting": ["^NO_SITE$"]}, "min_cells_per_guide": 10, "min_cells_per_target": 20},
        "protein": {"strip_feature_suffix": " protein", "n_pcs": 4},
        "rna": {"hvg": {"n_top_genes": 50}, "n_pcs": 10},
        "umap": {"min_features": 5},
        # synthetic cells have ~100 detected genes: keep the public thresholds from removing everything
        "qc": {"prefilter": {"enabled": True, "min_genes_per_cell": 50, "min_cells_per_gene": 0}, "filter": {"enabled": True, "rna": {"min_genes": 60, "max_pct_mt": 100.0}}},
        "compute": {"n_jobs": 1, "chunk_rows": 40},
        "output": {"dir": str(outdir)},
    }
    def merge(dst, src):
        for k, v in src.items():
            if isinstance(v, dict) and isinstance(dst.get(k), dict):
                merge(dst[k], v)
            else:
                dst[k] = v

    merge(cfg, overrides)
    return cfg
