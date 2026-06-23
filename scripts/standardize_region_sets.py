#!/usr/bin/env python
"""Normalize all SCENIC+ region-set BED files to UCSC standard chromosomes."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--region-sets-dir", default=None, help="Default: $PROJECT_DIR/inputs/region_sets")
    parser.add_argument("--allowed-chroms", default=None, help="Default: $CHROMS resolved under $PROJECT_DIR")
    parser.add_argument("--tmp-dir", default=None, help="Default: $PROJECT_DIR/tmp/sort")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_project_path(path_value: str | None, default_rel: str, base: Path) -> Path:
    path = Path(path_value).expanduser() if path_value else base / default_rel
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def shell_quote(path: str | Path) -> str:
    return "'" + str(path).replace("'", "'\"'\"'") + "'"


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required command not found on PATH: {name}")
    return path


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    region_sets_dir = resolve_project_path(args.region_sets_dir, "inputs/region_sets", pdir)
    chroms_value = args.allowed_chroms or os.environ.get("CHROMS")
    if not chroms_value:
        raise SystemExit("ERROR: --allowed-chroms is required unless CHROMS is set in project_env.sh.")
    allowed_chroms = resolve_project_path(chroms_value, chroms_value, pdir)
    tmp_dir = resolve_project_path(args.tmp_dir, "tmp/sort", pdir)
    report_path = pdir / "results" / "pycistopic" / "qc" / "region_set_standardization.tsv"

    if not region_sets_dir.is_dir():
        raise FileNotFoundError(region_sets_dir)
    if not allowed_chroms.exists():
        raise FileNotFoundError(allowed_chroms)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    require_tool("sort")

    normalizer = Path(__file__).resolve().parent / "normalize_bedlike_to_ucsc_standard.py"
    if not normalizer.exists():
        raise FileNotFoundError(normalizer)

    beds = sorted(region_sets_dir.glob("**/*.bed"))
    if not beds:
        raise FileNotFoundError(f"No BED files found under {region_sets_dir}")

    records = []
    for bed in beds:
        before = sum(1 for _ in bed.open())
        with tempfile.NamedTemporaryFile("w", delete=False, dir=str(bed.parent), prefix=bed.name, suffix=".tmp") as handle:
            tmp_path = Path(handle.name)
        cmd = (
            f"{shell_quote(sys.executable)} {shell_quote(normalizer)} "
            f"--allowed-chroms {shell_quote(allowed_chroms)} < {shell_quote(bed)} "
            f"| sort -T {shell_quote(tmp_dir)} -k1,1 -k2,2n > {shell_quote(tmp_path)}"
        )
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        after = sum(1 for _ in tmp_path.open())
        tmp_path.replace(bed)
        records.append({"bed": str(bed), "rows_before": before, "rows_after": after})

    report = pd.DataFrame(records)
    report.to_csv(report_path, sep="\t", index=False)
    print(f"WROTE {report_path}")
    print(report.to_string(index=False))


if __name__ == "__main__":
    main()
