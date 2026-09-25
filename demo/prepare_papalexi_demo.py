#!/usr/bin/env python
"""Build the deterministic Papalexi ECCITE-seq demo subset (Stage D).

Normal users do NOT need this script: the resulting subset is bundled in the
repository (demo/data/papalexi_eccite/) and ``petrubseq-protein run --config
config/demo_papalexi.yaml`` uses it directly. This script is the provenance /
reconstruction workflow (public GEO source -> fetch -> deterministic
preparation -> the same bundled files) for maintainers and for anyone who
wants to verify the bundle: with seed 0 and n_cells 1800 it reproduces
demo/data/papalexi_eccite byte for byte (demo/check_demo_reproducibility.py).

    python demo/prepare_papalexi_demo.py --source ${DATA_ROOT}/ECCITE-seq \\
        --output ${DATA_ROOT}/demo/papalexi_eccite [--n-cells 1800] [--seed 0]

Source = the GEO GSE153056 supplementary files of the pooled ECCITE-seq screen
(GSM4633614-GSM4633618, see docs/datasets/PAPALEXI_ECCITE_AUDIT.md). The source
directory is only read (checked before and after by size + mtime).

Output (`--output`, created; existing files are overwritten):

    filtered_feature_bc_matrix/{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz}
        10x-style combined matrix: Gene Expression + Antibody Capture + CRISPR
        Guide Capture, raw integer UMI counts, unchanged values
    cells.csv           per-cell metadata: lane, lane_group, HTO call, hto_sample,
                        stimulation, guide-derived stratum used for sampling
    antibodies.csv      protein.feature_table (feature_id, antibody_name,
                        protein_name, gene_symbol; clone left empty)
    guide_targets.csv   guide -> target (from the guide ID; perturbation.guide_target_table)
    selected_cells.txt  the selected cell IDs, one per line
    demo_manifest.json  sources (path, size, sha256), seed, algorithm, strata,
                        derived dimensions, versions

Selection (deterministic, `--seed`):

1. HTO singlets only: top HTO >= 10 UMIs and > 3 x the runner-up.
2. Preliminary guide call with the pipeline's dominant rule (top guide >= 3
   UMIs and > 2 x the runner-up) -> target label, else ``ambiguous`` /
   ``unassigned``; used ONLY to define sampling strata (the pipeline recomputes
   the assignment from the guide counts).
3. Strata = (hto_sample, target label). Each stratum gets
   max(min_per_stratum, share of --n-cells proportional to its size), capped
   at its size; the largest strata are trimmed until the total equals
   --n-cells. Cells inside a stratum are drawn without replacement with
   numpy's default_rng(seed) from the sorted cell list.

Nothing is normalized, filtered by QC, or altered in value.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import platform
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

SOURCES = {
    "rna": "GSM4633614_ECCITE_cDNA_counts.tsv.gz",
    "adt": "GSM4633615_ECCITE_ADT_counts.tsv.gz",
    "adt_barcodes": "GSM4633615_ECCITE_ADT_Barcodes.csv.gz",
    "hto": "GSM4633616_ECCITE_HTO_counts.tsv.gz",
    "gdo": "GSM4633618_ECCITE_GDO_counts.tsv.gz",
    "gdo_barcodes": "GSM4633618_ECCITE_GDO_Barcodes.csv.gz",
}
# Gene symbols of the four antibody targets (public knowledge of the antigen,
# not inferred from the data); clones are not given in the deposit -> left empty.
ANTIBODIES = {
    "CD86": {"protein_name": "CD86 (B7-2)", "gene_symbol": "CD86"},
    "PDL1": {"protein_name": "PD-L1 (CD274)", "gene_symbol": "CD274"},
    "PDL2": {"protein_name": "PD-L2 (CD273)", "gene_symbol": "PDCD1LG2"},
    "CD366": {"protein_name": "TIM-3 (CD366)", "gene_symbol": "HAVCR2"},
}
HTO_MIN_UMI, HTO_RATIO = 10, 3.0
GUIDE_MIN_UMI, GUIDE_RATIO = 3, 2.0


def sha256(path: Path, block: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def snapshot(paths: List[Path]) -> Dict[str, Tuple[int, float]]:
    return {str(p): (p.stat().st_size, p.stat().st_mtime) for p in paths}


def read_small(path: Path) -> pd.DataFrame:
    """features x cells dense TSV -> cells x features int64 DataFrame."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.T


def top_two(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    part = np.argsort(M, axis=1)
    top_idx = part[:, -1]
    top = M[np.arange(len(M)), top_idx]
    second = M[np.arange(len(M)), part[:, -2]] if M.shape[1] > 1 else np.zeros(len(M))
    return top_idx, top, second


def allocate(sizes: pd.Series, n_target: int, min_per_stratum: int) -> pd.Series:
    """Proportional quotas with a floor, capped by stratum size, summing to n_target."""
    sizes = sizes.sort_index()
    share = sizes / sizes.sum() * n_target
    quota = np.maximum(np.floor(share).astype(int), min_per_stratum)
    quota = np.minimum(quota, sizes)
    quota = pd.Series(quota, index=sizes.index)
    # trim the largest strata one cell at a time until the total matches
    while quota.sum() > n_target:
        k = quota.idxmax()
        quota[k] -= 1
    # if the floors under-fill, add to the strata with the largest remainder that still have room
    while quota.sum() < n_target:
        room = sizes - quota
        cand = room[room > 0]
        if cand.empty:
            break
        rem = (share - quota)[cand.index]
        quota[rem.idxmax()] += 1
    return quota.astype(int)


def read_rna_subset(path: Path, cells: List[str], chunk_rows: int = 500) -> Tuple[List[str], sp.csr_matrix]:
    """Stream the dense features x cells TSV and keep only ``cells`` (as columns)."""
    with gzip.open(path, "rt") as fh:
        header = [h.strip('"') for h in fh.readline().rstrip("\n").split("\t")]
    all_cells = header[1:]
    pos = pd.Index(all_cells).get_indexer(cells)
    if (pos < 0).any():
        raise SystemExit(f"{(pos < 0).sum()} selected cells missing from the RNA matrix")
    usecols = [0] + [int(p) + 1 for p in pos]
    genes: List[str] = []
    blocks = []
    for chunk in pd.read_csv(path, sep="\t", header=0, usecols=usecols, chunksize=chunk_rows, dtype={0: str}, compression="gzip"):
        genes.extend(chunk.iloc[:, 0].astype(str).str.strip('"').tolist())
        vals = chunk.iloc[:, 1:].to_numpy()
        if not np.issubdtype(vals.dtype, np.integer):
            if not np.all(np.mod(vals, 1) == 0):
                raise SystemExit("RNA matrix contains non-integer values; expected raw counts")
            vals = vals.astype(np.int64)
        blocks.append(sp.csr_matrix(vals.astype(np.int32)))
    M = sp.vstack(blocks, format="csr")  # genes x selected cells (columns in `cells` order)
    return genes, M


def _gz(path: Path, text: bool):
    """gzip writer with a fixed header timestamp, so identical content gives identical bytes."""
    raw = gzip.GzipFile(filename="", mode="wb", fileobj=open(path, "wb"), mtime=0)
    return io.TextIOWrapper(raw, encoding="utf-8", newline="\n") if text else raw


def write_mtx(outdir: Path, cells: List[str], feats: pd.DataFrame, M: sp.spmatrix) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    with _gz(outdir / "barcodes.tsv.gz", True) as fh:
        fh.write("\n".join(cells) + "\n")
    with _gz(outdir / "features.tsv.gz", True) as fh:
        feats.to_csv(fh, sep="\t", header=False, index=False)
    with _gz(outdir / "matrix.mtx.gz", False) as fh:
        scipy.io.mmwrite(fh, sp.csc_matrix(M, dtype=np.int64), field="integer", comment="Papalexi ECCITE-seq demo subset; raw UMI counts")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True, help="directory with the GSE153056 supplementary files")
    ap.add_argument("--output", required=True, help="demo directory to create")
    ap.add_argument("--n-cells", type=int, default=1800)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-per-stratum", type=int, default=8)
    ap.add_argument("--no-source-hash", action="store_true", help="skip sha256 of the source files (faster)")
    a = ap.parse_args()
    src, out = Path(a.source), Path(a.output)
    paths = {k: src / v for k, v in SOURCES.items()}
    missing = [str(p) for p in paths.values() if not p.is_file()]
    if missing:
        raise SystemExit("missing source files:\n  " + "\n  ".join(missing))
    if out.resolve() == src.resolve() or src.resolve() in out.resolve().parents:
        raise SystemExit("--output must not lie inside --source")
    before = snapshot(list(paths.values()))
    t0 = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    print(f"[prepare] source {src} -> {out}; n_cells={a.n_cells} seed={a.seed}")

    # -- small matrices -------------------------------------------------------
    adt, hto, gdo = read_small(paths["adt"]), read_small(paths["hto"]), read_small(paths["gdo"])
    cells_all = adt.index
    for name, df in (("hto", hto), ("gdo", gdo)):
        if not df.index.equals(cells_all):
            raise SystemExit(f"{name} cells differ from the ADT cells")
    for name, df in (("adt", adt), ("hto", hto), ("gdo", gdo)):
        if not all(np.issubdtype(t, np.integer) for t in df.dtypes):
            raise SystemExit(f"{name}: non-integer values")
    lane = cells_all.str.split("_").str[0]
    lane_group = np.where(lane.isin(["l1", "l2", "l3", "l4"]), "lanes_1-4", "lanes_5-8")
    # -- HTO demultiplexing (simple deterministic rule) -----------------------
    hi, ht, hs = top_two(hto.to_numpy())
    singlet = (ht >= HTO_MIN_UMI) & (ht > HTO_RATIO * hs)
    hto_label = hto.columns.to_numpy()[hi]
    hto_class = np.where(singlet, "singlet", "doublet_or_negative")
    # -- preliminary guide call (strata only) ---------------------------------
    gi, gt, gs = top_two(gdo.to_numpy())
    dominant = (gt >= GUIDE_MIN_UMI) & (gt > GUIDE_RATIO * gs)
    guide_top = gdo.columns.to_numpy()[gi]
    target_top = pd.Series(guide_top).str.replace(r"g\d+$", "", regex=True).to_numpy()
    stratum_target = np.where(dominant, target_top, np.where(gt >= GUIDE_MIN_UMI, "ambiguous", "unassigned"))
    meta = pd.DataFrame(
        {
            "cell": cells_all, "lane": lane, "lane_group": lane_group,
            "hto_top_label": hto_label, "hto_top_count": ht.astype(int), "hto_second_count": hs.astype(int), "hto_class": hto_class,
            "hto_sample": np.where(singlet, hto_label, ""),
            "stimulation": np.where(singlet, np.where(pd.Series(hto_label).str.endswith("-tx"), "IFNg", "control"), ""),
            "sampling_stratum_target": stratum_target,
        }
    ).set_index("cell")
    eligible = meta[meta["hto_class"] == "singlet"]
    strata = eligible.groupby(["hto_sample", "sampling_stratum_target"]).size()
    quota = allocate(strata, a.n_cells, a.min_per_stratum)
    rng = np.random.default_rng(a.seed)
    selected: List[str] = []
    strata_rec = []
    for (sample, target), n in quota.sort_index().items():
        pool = sorted(eligible.index[(eligible["hto_sample"] == sample) & (eligible["sampling_stratum_target"] == target)])
        take = sorted(rng.choice(pool, size=int(n), replace=False).tolist()) if n < len(pool) else pool
        selected.extend(take)
        strata_rec.append({"hto_sample": sample, "target": target, "eligible": int(len(pool)), "selected": int(len(take))})
    selected = [c for c in cells_all if c in set(selected)]  # keep source order
    print(f"[prepare] eligible singlets {len(eligible)} of {len(meta)}; strata {len(quota)}; selected {len(selected)}")

    # -- RNA (streamed, only the selected columns) ----------------------------
    genes, rna = read_rna_subset(paths["rna"], selected)
    if len(set(genes)) != len(genes):
        raise SystemExit("duplicate gene symbols in the RNA matrix")
    adt_sel = adt.loc[selected]
    gdo_sel = gdo.loc[selected]
    feats = pd.concat([
        pd.DataFrame({"id": genes, "name": genes, "type": "Gene Expression"}),
        pd.DataFrame({"id": list(adt.columns), "name": list(adt.columns), "type": "Antibody Capture"}),
        pd.DataFrame({"id": list(gdo.columns), "name": [g for g in gdo.columns], "type": "CRISPR Guide Capture"}),
    ], ignore_index=True)
    M = sp.vstack([rna, sp.csr_matrix(adt_sel.to_numpy().T), sp.csr_matrix(gdo_sel.to_numpy().T)], format="csr")
    out.mkdir(parents=True, exist_ok=True)
    write_mtx(out / "filtered_feature_bc_matrix", selected, feats, M)
    meta.loc[selected].to_csv(out / "cells.csv")
    with open(out / "selected_cells.txt", "w") as fh:
        fh.write("\n".join(selected) + "\n")
    pd.DataFrame([{"feature_id": k, "antibody_name": k, "protein_name": v["protein_name"], "gene_symbol": v["gene_symbol"], "clone": ""} for k, v in ANTIBODIES.items() if k in adt.columns]).to_csv(out / "antibodies.csv", index=False)
    gt_table = pd.DataFrame({"guide": list(gdo.columns), "target": pd.Series(list(gdo.columns)).str.replace(r"g\d+$", "", regex=True)})
    gt_table.to_csv(out / "guide_targets.csv", index=False)
    # -- manifest -------------------------------------------------------------
    after = snapshot(list(paths.values()))
    if before != after:
        raise SystemExit("source files changed during preparation")
    manifest = {
        "demo": "papalexi_eccite", "accession": "GSE153056", "experiment": "pooled ECCITE-seq screen (GSM4633614-GSM4633618)",
        "created": t0, "tool": {"script": Path(__file__).name, "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "parameters": {"n_cells": a.n_cells, "seed": a.seed, "min_per_stratum": a.min_per_stratum, "hto_singlet_rule": f"top >= {HTO_MIN_UMI} UMIs and top > {HTO_RATIO} x second", "guide_stratum_rule": f"top >= {GUIDE_MIN_UMI} UMIs and top > {GUIDE_RATIO} x second (dominant), else ambiguous/unassigned", "sampling": "stratified by (hto_sample, target); proportional quotas with floor, largest strata trimmed to n_cells; numpy default_rng(seed).choice without replacement on sorted cells"},
        "sources": {k: {"path": str(p), "size_bytes": p.stat().st_size, **({} if a.no_source_hash else {"sha256": sha256(p)})} for k, p in paths.items()},
        "source_summary": {"n_cells": int(len(meta)), "n_genes": int(len(genes)), "n_antibodies": int(adt.shape[1]), "n_guides": int(gdo.shape[1]), "cells_per_lane": {k: int(v) for k, v in meta["lane"].value_counts().sort_index().items()}, "hto_classes": {k: int(v) for k, v in meta["hto_class"].value_counts().items()}, "hto_samples_singlets": {k: int(v) for k, v in eligible["hto_sample"].value_counts().items()}, "preliminary_guide_call": {"dominant": int(dominant.sum()), "ambiguous": int(((gt >= GUIDE_MIN_UMI) & ~dominant).sum()), "unassigned": int((gt < GUIDE_MIN_UMI).sum())}},
        "strata": strata_rec,
        "selected": {"n_cells": len(selected), "cells_per_lane": {k: int(v) for k, v in meta.loc[selected, "lane"].value_counts().sort_index().items()}, "hto_samples": {k: int(v) for k, v in meta.loc[selected, "hto_sample"].value_counts().items()}, "targets": {k: int(v) for k, v in meta.loc[selected, "sampling_stratum_target"].value_counts().items()}},
        "derived": {"filtered_feature_bc_matrix": {"n_features": int(M.shape[0]), "n_cells": int(M.shape[1]), "gene_expression": int(len(genes)), "antibody_capture": int(adt.shape[1]), "crispr_guide_capture": int(gdo.shape[1]), "nnz": int(M.nnz)}, "cells_csv_columns": list(meta.columns)},
        "notes": ["All HTOs with signal are 'tx' (IFN-gamma stimulated) samples; stimulation is constant. Replicate naming is ambiguous between the two deposited HTO barcode files; hto_sample keeps the label of the count matrix verbatim.", "No isotype-control antibodies exist in the panel.", "Values are unchanged raw UMI counts; no normalization or QC filtering was applied."],
    }
    with open(out / "demo_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1)
    print(f"[prepare] wrote {out}: {M.shape[0]} features x {M.shape[1]} cells (nnz {M.nnz})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
