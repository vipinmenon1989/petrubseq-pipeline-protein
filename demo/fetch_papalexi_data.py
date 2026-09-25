#!/usr/bin/env python
"""Fetch the Papalexi ECCITE-seq source files (GEO GSE153056).

Normal users do NOT need this script: the demo subset is bundled in the
repository (demo/data/papalexi_eccite/). It exists to rebuild that subset from
the public accession (``demo/prepare_papalexi_demo.py``), for provenance and
maintainer verification.

    python demo/fetch_papalexi_data.py --dest ${DATA_ROOT}/ECCITE-seq [--all]

Downloads the supplementary count matrices of the pooled ECCITE-seq screen
(GSM4633614-GSM4633618; ``--all`` adds the CITE-seq and arrayed samples) from
the accession-based GEO FTP layout

    https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4633nnn/<GSM>/suppl/<file>

which is stable (it is derived from the accession, not a session URL).
Existing files with the expected size are skipped, so the script is safe to
re-run and never touches a complete local copy. It downloads into
``<name>.part`` and renames on success. Nothing is decompressed or modified.
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4633nnn/{gsm}/suppl/{name}"
# (GSM, file name, size in bytes as deposited)
POOLED = [
    ("GSM4633614", "GSM4633614_ECCITE_cDNA_counts.tsv.gz", 64135801),
    ("GSM4633615", "GSM4633615_ECCITE_ADT_Barcodes.csv.gz", 116),
    ("GSM4633615", "GSM4633615_ECCITE_ADT_counts.tsv.gz", 196344),
    ("GSM4633616", "GSM4633616_ECCITE_HTO_1-4_Barcodes.csv.gz", 164),
    ("GSM4633616", "GSM4633616_ECCITE_HTO_counts.tsv.gz", 170306),
    ("GSM4633617", "GSM4633617_ECCITE_HTO_5-8_Barcodes.csv.gz", 160),
    ("GSM4633618", "GSM4633618_ECCITE_GDO_Barcodes.csv.gz", 1242),
    ("GSM4633618", "GSM4633618_ECCITE_GDO_counts.tsv.gz", 337680),
]
OTHERS = [
    ("GSM4633605", "GSM4633605_CITE_cDNA_counts.tsv.gz", 22773984),
    ("GSM4633606", "GSM4633606_CITE_ADT_Barcodes.csv.gz", 170),
    ("GSM4633606", "GSM4633606_CITE_ADT_counts.tsv.gz", 104383),
    ("GSM4633607", "GSM4633607_CITE_HTO_Barcodes.csv.gz", 123),
    ("GSM4633607", "GSM4633607_CITE_HTO_counts.tsv.gz", 72853),
    ("GSM4633608", "GSM4633608_ECCITE_Arrayed_cDNA_counts.tsv.gz", 39883220),
    ("GSM4633609", "GSM4633609_ECCITE_Arrayed_ADT_Barcodes.csv.gz", 124),
    ("GSM4633609", "GSM4633609_ECCITE_Arrayed_ADT_counts.tsv.gz", 150657),
    ("GSM4633610", "GSM4633610_ECCITE_Arrayed_HTO1_Barcodes.csv.gz", 187),
    ("GSM4633610", "GSM4633610_ECCITE_Arrayed_HTO_counts.tsv.gz", 187144),
    ("GSM4633611", "GSM4633611_ECCITE_Arrayed_HTO2_Barcodes.csv.gz", 194),
    ("GSM4633612", "GSM4633612_ECCITE_Arrayed_GDO_Barcodes.csv.gz", 270),
    ("GSM4633612", "GSM4633612_ECCITE_Arrayed_GO_CITE03_counts.tsv.gz", 170201),
    ("GSM4633613", "GSM4633613_ECCITE_Arrayed_GO_LENTI_counts.tsv.gz", 104861),
]


def fetch(gsm: str, name: str, size: int, dest: Path, dry: bool) -> str:
    target = dest / name
    if target.is_file() and target.stat().st_size == size:
        return "present"
    url = BASE.format(gsm=gsm, name=name)
    if dry:
        return f"would download {url}"
    part = dest / (name + ".part")
    with urllib.request.urlopen(url, timeout=120) as r, open(part, "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    if part.stat().st_size != size:
        part.unlink()
        raise RuntimeError(f"{name}: downloaded {part.stat().st_size} bytes, expected {size}")
    part.rename(target)
    return "downloaded"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", required=True)
    ap.add_argument("--all", action="store_true", help="also fetch the CITE-seq and arrayed ECCITE samples")
    ap.add_argument("--dry-run", action="store_true", help="only report what would be downloaded")
    a = ap.parse_args()
    dest = Path(a.dest)
    dest.mkdir(parents=True, exist_ok=True)
    files = POOLED + (OTHERS if a.all else [])
    for gsm, name, size in files:
        print(f"{name}: {fetch(gsm, name, size, dest, a.dry_run)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
