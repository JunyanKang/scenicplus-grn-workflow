#!/usr/bin/env python
"""Validate inputs/sample_sheet.tsv for the active ATAC fragment and peak inputs."""
from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = ["sample_id", "condition", "fragments_tsv_gz", "atac_peaks_bed"]
OUTPUT_COLUMNS = ["sample_id", "condition", "organism", "fragments_tsv_gz", "atac_peaks_bed"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-sheet", default=None, help="Default: $PROJECT_DIR/inputs/sample_sheet.tsv")
    parser.add_argument("--out", default=None, help="Default: overwrite the active sample sheet after validation.")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_existing(path_value: str, base: Path) -> Path:
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def write_table_atomic(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(out_path.parent), prefix=out_path.name, suffix=".tmp") as handle:
        tmp_path = Path(handle.name)
        df.to_csv(handle, sep="\t", index=False)
    tmp_path.replace(out_path)


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    sample_sheet = Path(args.sample_sheet).expanduser() if args.sample_sheet else pdir / "inputs" / "sample_sheet.tsv"
    sample_sheet = sample_sheet.resolve()
    out_path = Path(args.out).expanduser().resolve() if args.out else sample_sheet

    if not sample_sheet.exists():
        raise FileNotFoundError(f"sample_sheet.tsv not found: {sample_sheet}")
    df = pd.read_csv(sample_sheet, sep="\t", dtype=str).fillna("")
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError("sample_sheet.tsv missing columns: " + ", ".join(missing_cols))
    organism = os.environ.get("ORGANISM", "").strip()
    if not organism:
        raise SystemExit("ERROR: ORGANISM is not set. Source project_env.sh from Step 0.")
    if "organism" in df.columns:
        observed = sorted(set(x for x in df["organism"].astype(str) if x))
        if observed and observed != [organism]:
            raise ValueError(
                "sample_sheet.tsv organism values must match Step 0 ORGANISM="
                f"{organism!r}; observed: {observed}"
            )
    df = df[REQUIRED_COLUMNS].copy()

    prepared = []
    report = []
    for row in df.itertuples(index=False):
        sample_id = str(row.sample_id)
        fragments = resolve_existing(row.fragments_tsv_gz, pdir)
        peaks = resolve_existing(row.atac_peaks_bed, pdir)
        checks = {
            "fragments_tsv_gz": fragments.is_file(),
            "fragments_tbi": Path(str(fragments) + ".tbi").is_file(),
            "atac_peaks_bed": peaks.is_file(),
        }
        report.append({"sample_id": sample_id, **{k: str(v) for k, v in checks.items()}})
        if not checks["fragments_tsv_gz"]:
            raise FileNotFoundError(f"missing fragments for {sample_id}: {fragments}")
        if not checks["fragments_tbi"]:
            raise FileNotFoundError(f"missing fragment index for {sample_id}: {fragments}.tbi")
        if not checks["atac_peaks_bed"]:
            raise FileNotFoundError(f"missing ATAC peaks for {sample_id}: {peaks}")

        prepared.append(
            {
                "sample_id": sample_id,
                "condition": str(row.condition),
                "organism": organism,
                "fragments_tsv_gz": str(fragments),
                "atac_peaks_bed": str(peaks),
            }
        )

    out = pd.DataFrame(prepared)
    out = out[OUTPUT_COLUMNS]
    write_table_atomic(out, out_path)
    report_path = pdir / "inputs" / "sample_sheet.validation_report.tsv"
    pd.DataFrame(report).to_csv(report_path, sep="\t", index=False)
    print(f"WROTE {out_path}")
    print(f"WROTE {report_path}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
