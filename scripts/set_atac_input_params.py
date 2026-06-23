#!/usr/bin/env python
"""Write the ATAC-only input parameter table used by sample-sheet builders."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layout", choices=["cellranger_outs", "split_ge_arc"], default=None)
    parser.add_argument("--atac-data-root", default=None)
    parser.add_argument("--out", default=None, help="Default: $PROJECT_DIR/inputs/atac_input_params.tsv")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def write_project_setting(project: Path, key: str, value: str) -> None:
    config = project / "scenicplus_project.env"
    lines = config.read_text().splitlines() if config.exists() else []
    assignment = f"{key}={value!r}"
    out = []
    replaced = False
    for line in lines:
        if line.startswith(f"{key}="):
            out.append(assignment)
            replaced = True
        else:
            out.append(line)
    if not replaced:
        out.append(assignment)
    config.write_text("\n".join(out) + "\n")


def resolve_root(root: str, base: Path) -> Path:
    path = Path(root).expanduser()
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    layout = args.layout or os.environ.get("ATAC_INPUT_LAYOUT", "")
    atac_data_root = args.atac_data_root or os.environ.get("ATAC_DATA_ROOT", "")
    if layout not in {"cellranger_outs", "split_ge_arc"}:
        raise SystemExit(
            "ERROR: ATAC input layout is required. Set ATAC_INPUT_LAYOUT in project_env.sh "
            "or pass --layout cellranger_outs|split_ge_arc."
        )
    if not atac_data_root:
        raise SystemExit(
            "ERROR: ATAC data root is required. Set ATAC_DATA_ROOT in project_env.sh "
            "or pass --atac-data-root /absolute/path/to/atac_only_root."
        )
    atac_root = resolve_root(atac_data_root, pdir)
    if not atac_root.exists():
        raise FileNotFoundError(f"ATAC data root does not exist: {atac_root}")
    out_path = Path(args.out).expanduser() if args.out else pdir / "inputs" / "atac_input_params.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        [
            {"parameter": "layout", "value": layout},
            {"parameter": "atac_data_root", "value": str(atac_root)},
        ]
    )
    df.to_csv(out_path, sep="\t", index=False)
    write_project_setting(pdir, "ATAC_INPUT_PARAMS", str(out_path))
    write_project_setting(pdir, "ATAC_INPUT_LAYOUT", layout)
    write_project_setting(pdir, "ATAC_DATA_ROOT", str(atac_root))
    print(f"WROTE {out_path}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
