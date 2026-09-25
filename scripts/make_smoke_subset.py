#!/usr/bin/env python
"""Build a small, non-destructive SCP1064 subset in the *same file formats* as
the download, so the pipeline's real-input parsing is exercised end to end.

Takes the first N cells (columns) of selected RNA chunk files and subsets the
protein, metadata, guide and embedding files to those cells. Original files
are only read.

    python scripts/make_smoke_subset.py --src ../data/SCP1064 --dest ../data/SCP1064/smoke --n-cells 1500
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
from pathlib import Path

RNA_FILES = [  # one chunk per condition block (Control / IFNg / Co-culture)
    "expression/RNA_expression_subset1a.csv.gz",
    "expression/RNA_expression_subset3.csv.gz",
    "expression/RNA_expression_subset4.csv.gz",
]


def cut_columns(src: Path, dest: Path, n_cells: int) -> list[str]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = f"zcat {src} | cut -d, -f1-{n_cells + 1} | gzip -1 > {dest}"
    subprocess.run(cmd, shell=True, check=True)
    with gzip.open(dest, "rt") as fh:
        return fh.readline().rstrip("\n").split(",")[1:]


def subset_columns(src: Path, dest: Path, cells: set[str]) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if src.suffix == ".gz" else open
    with opener(src, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        keep = [0] + [i for i, c in enumerate(header) if c in cells]
        wopen = gzip.open if dest.suffix == ".gz" else open
        with wopen(dest, "wt") as out:
            out.write(",".join(header[i] for i in keep) + "\n")
            for line in fh:
                toks = line.rstrip("\n").split(",")
                out.write(",".join(toks[i] for i in keep) + "\n")


def subset_rows(src: Path, dest: Path, cells: set[str], scp_type_row: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(src) as fh, open(dest, "w") as out:
        out.write(fh.readline())
        if scp_type_row:
            out.write(fh.readline())
        for line in fh:
            if line.split(",", 1)[0] in cells:
                out.write(line)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--dest", required=True)
    ap.add_argument("--n-cells", type=int, default=1500, help="cells per RNA chunk file")
    a = ap.parse_args()
    src, dest = Path(a.src), Path(a.dest)
    cells: list[str] = []
    for rel in RNA_FILES:
        cells += cut_columns(src / rel, dest / rel, a.n_cells)
        print(f"{rel}: {a.n_cells} cells")
    cs = set(cells)
    subset_columns(src / "documentation/Protein_expression.csv.gz", dest / "documentation/Protein_expression.csv.gz", cs)
    subset_columns(src / "other/raw_CITE_expression.csv.gz", dest / "other/raw_CITE_expression.csv.gz", cs)
    subset_rows(src / "metadata/RNA_metadata.csv", dest / "metadata/RNA_metadata.csv", cs, True)
    subset_rows(src / "documentation/all_sgRNA_assignments.txt", dest / "documentation/all_sgRNA_assignments.txt", cs, False)
    subset_rows(src / "cluster/5fd0e449771a5b0db7207711/RNA_UMAP_cluster.csv", dest / "cluster/5fd0e449771a5b0db7207711/RNA_UMAP_cluster.csv", cs, True)
    print(f"smoke subset written to {dest}: {len(cells)} cells")


if __name__ == "__main__":
    main()
