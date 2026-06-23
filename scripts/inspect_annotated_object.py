#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
import os
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-path", default=os.environ.get("ANNOTATED_OBJECT", ""))
    parser.add_argument("--object-format", default="")
    return parser.parse_args()


def detect_format(path: Path, explicit: str) -> str:
    if explicit:
        return explicit.lower().lstrip(".")
    suffix = path.suffix.lower().lstrip(".")
    if suffix in {"rds", "qs", "h5ad"}:
        return suffix
    raise SystemExit(
        f"Cannot infer annotated object format from {path.name}. "
        "Supported formats: rds, qs, h5ad. Pass --object-format if needed."
    )


def main() -> None:
    args = parse_args()
    if not args.object_path:
        raise SystemExit("Provide --object-path or set ANNOTATED_OBJECT.")
    object_path = Path(args.object_path).expanduser().resolve()
    if not object_path.exists():
        raise FileNotFoundError(object_path)
    fmt = detect_format(object_path, args.object_format)
    if fmt in {"rds", "qs"}:
        cmd = [
            "Rscript",
            str(SCRIPT_DIR / "inspect_annotated_object_for_scenicplus.R"),
            "--object-path",
            str(object_path),
            "--object-format",
            fmt,
        ]
    elif fmt == "h5ad":
        cmd = [
            sys.executable,
            str(SCRIPT_DIR / "inspect_annotated_h5ad_for_scenicplus.py"),
            "--object-path",
            str(object_path),
        ]
    else:
        raise SystemExit(f"Unsupported annotated object format: {fmt}")
    print("RUN", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
