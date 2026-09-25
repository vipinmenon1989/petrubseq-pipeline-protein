"""Command-line interface.

    petrubseq-protein audit    --data /path/to/dataset [--config cfg.yaml] [--out audit.json] [--deep]
    petrubseq-protein validate --config config/scp1064.yaml
    petrubseq-protein run      --config config/scp1064.yaml [--outdir DIR] [--name RUN] [--max-cells N]
    petrubseq-protein init-config my_dataset.yaml
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__
from .config import Config, ConfigError

logger = logging.getLogger("petrubseq_protein")


def _setup_console_logging(verbose: bool) -> None:
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")


def cmd_audit(args: argparse.Namespace) -> int:
    from .audit import audit_config, audit_directory, render_markdown

    if args.config:
        cfg = Config.from_yaml(args.config)
        if args.data:
            cfg.dataset.input_dir = args.data
        rec = audit_config(cfg, deep=args.deep, rna_sample_cells=args.sample_cells)
    else:
        if not args.data:
            print("audit: give --data DIR and/or --config FILE", file=sys.stderr)
            return 2
        rec = audit_directory(args.data, deep=args.deep)
        rec = {"dataset": Path(args.data).name, "input_dir": str(Path(args.data).resolve()), "inventory": rec}
    md = render_markdown(rec)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as fh:
            json.dump(rec, fh, indent=2, default=str)
        md_path = out.with_suffix(".md")
        md_path.write_text(md)
        print(f"audit written: {out} and {md_path}")
    else:
        print(md)
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    from .pipeline import validate_inputs

    cfg = Config.from_yaml(args.config)
    rec = validate_inputs(cfg)
    print(f"config OK: {cfg.source_path}")
    for k, v in rec.items():
        print(f"  {k}: {v['path']} ({v.get('size_bytes', 0)/1e6:.1f} MB)")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from .pipeline import run_pipeline

    cfg = Config.from_yaml(args.config)
    if args.outdir:
        cfg.output.dir = args.outdir
    if args.name:
        cfg.run.name = args.name
    if args.max_cells:
        cfg.compute.max_cells = args.max_cells
    if args.n_jobs:
        cfg.compute.n_jobs = args.n_jobs
    if args.no_umap:
        cfg.umap.enabled = False
    try:
        res = run_pipeline(cfg)
    except Exception as exc:  # full traceback goes to logs/run.log; a clean line to the console
        logger.exception("pipeline failed: %s", exc)
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1
    print(f"done: {res.n_cells} cells x {res.n_genes} genes x {res.n_proteins} proteins in {res.runtime_seconds:.0f}s ({len(res.warnings)} warnings)")
    print(f"  report: {res.report_html}")
    if res.h5ad:
        print(f"  h5ad:   {res.h5ad}")
    return 0


def cmd_init_config(args: argparse.Namespace) -> int:
    """Write the documented defaults as a starting config (inputs left to fill in)."""
    import yaml

    data = Config.defaults_dict()
    data["dataset"] = {"name": "my_dataset", "input_dir": "${DATA_ROOT}/my_dataset", "description": ""}
    data["output"]["dir"] = "${RESULTS_ROOT}/my_dataset"
    data["inputs"]["rna"]["files"] = ["<rna_matrix.csv.gz>"]
    path = Path(args.path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("# petrubseq-pipeline-protein configuration (documented defaults; see docs/DEFAULTS.md).\n")
        fh.write("# Fill in dataset.input_dir, inputs.* and columns.* for your dataset, then:\n")
        fh.write(f"#   petrubseq-protein validate --config {path}\n#   petrubseq-protein run --config {path}\n")
        yaml.safe_dump(data, fh, sort_keys=False, allow_unicode=True)
    print(f"wrote default configuration to {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="petrubseq-protein", description="Perturb-CITE-seq / multimodal Perturb-seq pipeline: input harmonization, QC, representations, perturbation assignment and optional perturbation-effect analyses (PS, lochNESS, gene programs, protein effects).")
    p.add_argument("--version", action="version", version=f"petrubseq-pipeline-protein {__version__}")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = p.add_subparsers(dest="command", required=True)

    a = sub.add_parser("audit", help="inventory and inspect a dataset directory (config-free), or role-aware audit with --config")
    a.add_argument("--data", help="dataset directory")
    a.add_argument("--config", help="run config; enables the role-aware audit (cell overlap, design, protein panel)")
    a.add_argument("--out", help="write JSON (and .md) here instead of printing markdown")
    a.add_argument("--deep", action="store_true", help="count rows of every file, even multi-GB ones")
    a.add_argument("--sample-cells", type=int, default=2000, help="cells sampled from the RNA matrix for value-state checks")
    a.set_defaults(func=cmd_audit)

    v = sub.add_parser("validate", help="parse the config and check that all inputs exist")
    v.add_argument("--config", required=True)
    v.set_defaults(func=cmd_validate)

    r = sub.add_parser("run", help="run the processing pipeline")
    r.add_argument("--config", required=True)
    r.add_argument("--outdir", help="override output.dir")
    r.add_argument("--max-cells", type=int, help="smoke test: load only the first N cells per RNA file")
    r.add_argument("--n-jobs", type=int, help="override compute.n_jobs")
    r.add_argument("--no-umap", action="store_true", help="skip UMAP (faster smoke tests)")
    r.add_argument("--name", help="override run.name (report / archive name)")
    r.set_defaults(func=cmd_run)

    i = sub.add_parser("init-config", help="write a config file with the documented defaults")
    i.add_argument("path")
    i.set_defaults(func=cmd_init_config)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _setup_console_logging(args.verbose)
    try:
        return args.func(args)
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
