"""Run provenance: versions, host, git state, inputs, timestamps.

Modelled on ``perturbseq_pipeline.provenance`` but reduced to what v0.1 needs.
The record is JSON-serialisable so it can be stored in ``adata.uns`` and
written next to the report.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import platform
import socket
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from .. import __version__

_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy", "scikit-learn", "umap-learn", "h5py", "matplotlib", "pyyaml")


def _git(args: Iterable[str], cwd: Path) -> Optional[str]:
    try:
        out = subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True, timeout=10)
    except Exception:  # pragma: no cover
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def git_info() -> Dict[str, Optional[str]]:
    here = Path(__file__).resolve().parent
    root = _git(["rev-parse", "--show-toplevel"], here)
    if root is None:
        return {"root": None, "branch": None, "commit": None, "dirty": None}
    status = _git(["status", "--porcelain", "--untracked-files=no"], Path(root))
    return {
        "root": root,
        "branch": _git(["rev-parse", "--abbrev-ref", "HEAD"], Path(root)),
        "commit": _git(["rev-parse", "HEAD"], Path(root)),
        "dirty": bool(status) if status is not None else None,
    }


def package_versions() -> Dict[str, str]:
    import importlib.metadata as md

    out = {"python": platform.python_version(), "petrubseq-pipeline-protein": __version__}
    for name in _PACKAGES:
        try:
            out[name] = md.version(name)
        except Exception:
            out[name] = "not installed"
    return out


def slurm_info() -> Dict[str, Optional[str]]:
    keys = ("SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_JOB_NODELIST", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_JOB_PARTITION")
    return {k.lower(): os.environ.get(k) for k in keys}


def file_record(path: Optional[Path]) -> Optional[Dict[str, Any]]:
    if path is None:
        return None
    p = Path(path)
    rec: Dict[str, Any] = {"path": str(p), "exists": p.exists()}
    if p.exists():
        st = p.stat()
        rec["size_bytes"] = st.st_size
        rec["mtime"] = _dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")
    return rec


def collect(config_path: Optional[str], inputs: Dict[str, Any], extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    rec: Dict[str, Any] = {
        "run_timestamp": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "hostname": socket.gethostname(),
        "cwd": os.getcwd(),
        "user": os.environ.get("USER"),
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
        "config_path": config_path,
        "git": git_info(),
        "slurm": slurm_info(),
        "packages": package_versions(),
        "inputs": inputs,
    }
    if extra:
        rec.update(extra)
    return rec


def write_json(rec: Dict[str, Any], path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(rec, fh, indent=2, default=str)
    return path
