#!/usr/bin/env python
"""Prove that two demo preparations are identical.

    python demo/check_demo_reproducibility.py DIR_A DIR_B

Compares sha256 of every file except ``demo_manifest.json`` and ``README.md`` (byte-identical
matrices, barcodes, features, metadata) and the manifests after dropping the
keys that legitimately differ (``created`` and the source/output paths).
Exit code 0 = identical.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

VOLATILE = {"created"}


def hashes(d: Path) -> dict:
    return {str(p.relative_to(d)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(d.rglob("*")) if p.is_file() and p.name not in ("demo_manifest.json", "README.md")}


def manifest(d: Path) -> dict:
    m = json.load(open(d / "demo_manifest.json"))
    for k in VOLATILE:
        m.pop(k, None)
    for rec in m.get("sources", {}).values():
        rec.pop("path", None)
    return m


def main() -> int:
    a, b = Path(sys.argv[1]), Path(sys.argv[2])
    ha, hb = hashes(a), hashes(b)
    diff = sorted(k for k in set(ha) | set(hb) if ha.get(k) != hb.get(k))
    same_manifest = manifest(a) == manifest(b)
    print(f"files compared: {len(set(ha) | set(hb))}; differing: {diff}; manifests identical (excluding {sorted(VOLATILE)} and paths): {same_manifest}")
    return 0 if not diff and same_manifest else 1


if __name__ == "__main__":
    sys.exit(main())
