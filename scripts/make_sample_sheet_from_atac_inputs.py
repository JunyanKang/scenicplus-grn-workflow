#!/usr/bin/env python
"""Dispatch ATAC sample-sheet generation based on inputs/atac_input_params.tsv."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", default=None, help="Default: $PROJECT_DIR/inputs/atac_input_params.tsv")
    return parser.parse_args()


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(
            f"ATAC input parameter table not found: {path}\n"
            "Create it first with set_atac_input_params.py."
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {str(k).strip(): str(v).strip() for k, v in zip(df.iloc[:, 0], df.iloc[:, 1]) if str(k).strip()}


def main() -> None:
    args = parse_args()
    params_path = Path(args.params).expanduser() if args.params else PROJECT / "inputs" / "atac_input_params.tsv"
    if not params_path.is_absolute():
        params_path = (PROJECT / params_path).resolve()
    params = read_params(params_path)
    layout = params.get("layout", "")
    if layout == "cellranger_outs":
        script = SCRIPT_DIR / "make_sample_sheet_from_cellranger_outs_layout.py"
    elif layout == "split_ge_arc":
        script = SCRIPT_DIR / "make_sample_sheet_from_split_ge_arc_layout.py"
    else:
        raise SystemExit(f"Unsupported ATAC input layout in {params_path}: {layout!r}")
    subprocess.run([sys.executable, str(script), "--params", str(params_path)], check=True)


if __name__ == "__main__":
    main()
