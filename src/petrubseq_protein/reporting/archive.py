"""Shareable ``<run>_results.tar.gz`` (report, figures, tables, logs; matrices excluded)."""

from __future__ import annotations

import fnmatch
import logging
import tarfile
from pathlib import Path
from typing import List, Optional, Sequence

logger = logging.getLogger("petrubseq_protein")


def archive_results(run_dir: Path, run_name: str, exclude: Sequence[str], name: Optional[str] = None) -> Optional[Path]:
    run_dir = Path(run_dir)
    fname = name or f"{run_name}_results.tar.gz"
    if not fname.endswith((".tar.gz", ".tgz")):
        fname += ".tar.gz"
    dest = run_dir / fname
    members: List[Path] = []
    for p in sorted(run_dir.rglob("*")):
        if not p.is_file() or p == dest or p.name.endswith(".partial"):
            continue
        rel = p.relative_to(run_dir)
        if any(fnmatch.fnmatch(str(rel), pat) or fnmatch.fnmatch(p.name, pat) for pat in exclude):
            continue
        members.append(p)
    if not members:
        logger.warning("nothing to archive in %s", run_dir)
        return None
    tmp = dest.with_suffix(dest.suffix + ".partial")
    try:
        with tarfile.open(tmp, "w:gz") as tar:
            for p in members:
                tar.add(p, arcname=str(Path(run_name) / p.relative_to(run_dir)))
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink()
    logger.info("archived %d files to %s (%.1f MB)", len(members), dest.name, dest.stat().st_size / 1e6)
    return dest
