#!/usr/bin/env python
"""Standardize ATAC peak and fragment paths to UCSC primary chromosomes."""
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
    parser.add_argument("--sample-sheet", default=None, help="Default: $PROJECT_DIR/inputs/sample_sheet.tsv")
    parser.add_argument("--allowed-chroms", default=None, help="Default: $CHROMS resolved under $PROJECT_DIR")
    parser.add_argument("--out-sample-sheet", default=None, help="Default: overwrite active sample sheet after validation.")
    parser.add_argument("--standard-peaks-dir", default=None, help="Default: $PROJECT_DIR/work/standard_peaks")
    parser.add_argument("--standard-fragments-dir", default=None, help="Default: $PROJECT_DIR/work/standard_fragments")
    parser.add_argument("--tmp-dir", default=None, help="Default: $PROJECT_DIR/tmp/sort")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required command not found on PATH: {name}")
    return path


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_project_path(path_value: str | None, default_rel: str, base: Path) -> Path:
    path = Path(path_value).expanduser() if path_value else base / default_rel
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def run_pipeline(command: str, out_file: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("wb") as handle:
        subprocess.run(command, shell=True, check=True, stdout=handle, executable="/bin/bash")


def shell_quote(path: str | Path) -> str:
    return "'" + str(path).replace("'", "'\"'\"'") + "'"


def write_table_atomic(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(out_path.parent), prefix=out_path.name, suffix=".tmp") as handle:
        tmp_path = Path(handle.name)
        df.to_csv(handle, sep="\t", index=False)
    tmp_path.replace(out_path)


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    sample_sheet = resolve_project_path(args.sample_sheet, "inputs/sample_sheet.tsv", pdir)
    chroms_value = args.allowed_chroms or os.environ.get("CHROMS")
    if not chroms_value:
        raise SystemExit("ERROR: --allowed-chroms is required unless CHROMS is set in project_env.sh.")
    allowed_chroms = resolve_project_path(chroms_value, chroms_value, pdir)
    if not sample_sheet.exists():
        raise FileNotFoundError(sample_sheet)
    if not allowed_chroms.exists():
        raise FileNotFoundError(allowed_chroms)

    for tool in ["gzip", "sort", "bgzip", "tabix"]:
        require_tool(tool)
    normalizer = Path(__file__).resolve().parent / "normalize_bedlike_to_ucsc_standard.py"
    if not normalizer.exists():
        raise FileNotFoundError(normalizer)

    df = pd.read_csv(sample_sheet, sep="\t", dtype=str).fillna("")
    required = {"sample_id", "fragments_tsv_gz", "atac_peaks_bed"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{sample_sheet} missing columns: {missing}")

    peaks_dir = resolve_project_path(args.standard_peaks_dir, "work/standard_peaks", pdir)
    fragments_dir = resolve_project_path(args.standard_fragments_dir, "work/standard_fragments", pdir)
    tmp_dir = resolve_project_path(args.tmp_dir, "tmp/sort", pdir)
    peaks_dir.mkdir(parents=True, exist_ok=True)
    fragments_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    updated_rows = []
    for row in df.itertuples(index=False):
        row_dict = row._asdict()
        sample_id = str(row_dict["sample_id"])
        peaks_in = resolve_project_path(row_dict["atac_peaks_bed"], row_dict["atac_peaks_bed"], pdir)
        fragments_in = resolve_project_path(row_dict["fragments_tsv_gz"], row_dict["fragments_tsv_gz"], pdir)
        if not peaks_in.exists():
            raise FileNotFoundError(f"Missing ATAC peaks for {sample_id}: {peaks_in}")
        if not fragments_in.exists():
            raise FileNotFoundError(f"Missing fragments for {sample_id}: {fragments_in}")
        if not Path(str(fragments_in) + ".tbi").exists():
            raise FileNotFoundError(f"Missing tabix index for {sample_id}: {fragments_in}.tbi")

        peaks_out = peaks_dir / f"{sample_id}.standard_peaks.bed"
        fragments_out = fragments_dir / f"{sample_id}.fragments.ucsc.standard.tsv.gz"
        if args.overwrite or not peaks_out.exists():
            cmd = (
                f"{shell_quote(sys.executable)} {shell_quote(normalizer)} "
                f"--allowed-chroms {shell_quote(allowed_chroms)} < {shell_quote(peaks_in)} "
                f"| sort -T {shell_quote(tmp_dir)} -k1,1 -k2,2n"
            )
            run_pipeline(cmd, peaks_out)
        if args.overwrite or not fragments_out.exists() or not Path(str(fragments_out) + ".tbi").exists():
            cmd = (
                f"gzip -dc {shell_quote(fragments_in)} "
                f"| {shell_quote(sys.executable)} {shell_quote(normalizer)} "
                f"--allowed-chroms {shell_quote(allowed_chroms)} --min-cols 5 "
                f"| sort -T {shell_quote(tmp_dir)} -k1,1 -k2,2n "
                f"| bgzip -c"
            )
            run_pipeline(cmd, fragments_out)
            subprocess.run(["tabix", "-f", "-p", "bed", str(fragments_out)], check=True)
        row_dict["atac_peaks_bed"] = str(peaks_out)
        row_dict["fragments_tsv_gz"] = str(fragments_out)
        updated_rows.append(row_dict)

    out = pd.DataFrame(updated_rows)
    out_path = resolve_project_path(args.out_sample_sheet, "inputs/sample_sheet.tsv", pdir)
    write_table_atomic(out, out_path)
    print(f"WROTE {out_path}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
