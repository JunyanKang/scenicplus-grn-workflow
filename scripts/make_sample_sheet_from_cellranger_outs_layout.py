#!/usr/bin/env python
"""Create inputs/sample_sheet.tsv from standard Cell Ranger sample/outs folders."""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", default=None, help="Default: $PROJECT_DIR/inputs/atac_input_params.tsv")
    parser.add_argument("--data-root", default=None, help="Override atac_data_root from the parameter table.")
    parser.add_argument("--out", default=None, help="Default: $PROJECT_DIR/inputs/sample_sheet.tsv")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(
            f"ATAC input parameter table not found: {path}\n"
            "Create it first with set_atac_input_params.py."
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.shape[1] < 2:
        raise ValueError(f"{path} must contain at least two tab-separated columns: parameter and value")
    key_col, value_col = df.columns[:2]
    return {str(k).strip(): str(v).strip() for k, v in zip(df[key_col], df[value_col]) if str(k).strip()}


def resolve_input_path(raw_path: str, base: Path) -> Path:
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def choose_existing(candidates: list[Path], label: str, sample_id: str) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    checked = "\n".join(f"  - {p}" for p in candidates)
    raise FileNotFoundError(f"No {label} for {sample_id}. Checked:\n{checked}")


def sample_dirs(root: Path) -> list[Path]:
    direct = [p for p in sorted(root.iterdir()) if p.is_dir() and (p / "outs").is_dir()]
    if direct:
        return direct
    if (root / "outs").is_dir():
        return [root]
    return []


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    params_path = Path(args.params).expanduser() if args.params else pdir / "inputs" / "atac_input_params.tsv"
    params = read_params(params_path)

    layout = params.get("layout", "")
    if layout and layout != "cellranger_outs":
        raise SystemExit(f"ERROR: {params_path} has layout={layout!r}; expected layout=cellranger_outs.")

    raw_root_value = args.data_root or params.get("atac_data_root")
    if not raw_root_value:
        raise SystemExit("ERROR: atac_data_root is missing from the ATAC input parameter table.")
    organism = os.environ.get("ORGANISM")
    if not organism:
        raise SystemExit("ERROR: ORGANISM is not set. Source project_env.sh from Step 0.")

    root = resolve_input_path(raw_root_value, pdir)
    if not root.exists():
        raise FileNotFoundError(f"ATAC data root does not exist: {root}")

    records = []
    dirs = sample_dirs(root)
    if not dirs:
        raise SystemExit(f"No sample folders with an outs/ subdirectory found under {root}")

    for sample_dir in dirs:
        sample_id = sample_dir.name if sample_dir.name != "outs" else root.name
        outs = sample_dir / "outs"
        fragments = choose_existing(
            [
                outs / "atac_fragments.tsv.gz",
                outs / "fragments.tsv.gz",
            ],
            "ATAC fragments file",
            sample_id,
        )
        peaks = choose_existing(
            [
                outs / "atac_peaks.bed",
                outs / "peaks.bed",
            ],
            "ATAC peaks BED",
            sample_id,
        )
        if not Path(str(fragments) + ".tbi").exists():
            raise FileNotFoundError(f"Missing fragment index for {sample_id}: {fragments}.tbi")
        records.append(
            {
                "sample_id": sample_id,
                "condition": re.sub(r"_[0-9]+$", "", sample_id),
                "organism": organism,
                "fragments_tsv_gz": str(fragments),
                "atac_peaks_bed": str(peaks),
            }
        )

    out_path = Path(args.out).expanduser() if args.out else pdir / "inputs" / "sample_sheet.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"WROTE {out_path}")
    print(df.to_string(index=False))
    print("\nEdit only the condition column if needed, then run validate_and_prepare_sample_sheet.py.")


if __name__ == "__main__":
    main()
