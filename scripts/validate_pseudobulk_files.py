#!/usr/bin/env python
"""Validate pycisTopic pseudobulk BED.GZ files and tabix indexes."""
from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk-dir", default="work/pycistopic/pseudobulk_bed")
    parser.add_argument("--out", default="results/pycistopic/qc/pseudobulk_file_validation.tsv")
    return parser.parse_args()


def p(path: str | Path) -> Path:
    out = Path(path).expanduser()
    return out if out.is_absolute() else PROJECT / out


def gzip_ok(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        return True, ""
    except Exception as exc:
        return False, repr(exc)


def main() -> None:
    args = parse_args()
    pseudo_dir = p(args.pseudobulk_dir)
    if not pseudo_dir.exists():
        raise FileNotFoundError(pseudo_dir)
    files = sorted(pseudo_dir.glob("*.bed.gz"))
    if not files:
        raise FileNotFoundError(f"No *.bed.gz files found in {pseudo_dir}")
    rows = []
    for bed_gz in files:
        ok, error = gzip_ok(bed_gz)
        index = bed_gz.with_suffix(bed_gz.suffix + ".tbi")
        rows.append(
            {
                "file": str(bed_gz),
                "size_bytes": bed_gz.stat().st_size,
                "gzip_ok": ok,
                "tabix_index": str(index),
                "tabix_index_exists": index.exists(),
                "error": error,
            }
        )
    out = p(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    print(df.to_string(index=False))
    bad = df.loc[(~df["gzip_ok"]) | (~df["tabix_index_exists"])]
    if not bad.empty:
        raise AssertionError(f"Invalid pseudobulk files detected. See {out}")
    print(f"WROTE {out}")
    print("pseudobulk file validation OK")


if __name__ == "__main__":
    main()
