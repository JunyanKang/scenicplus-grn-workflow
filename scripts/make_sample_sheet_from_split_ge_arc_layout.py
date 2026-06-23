#!/usr/bin/env python
"""Create inputs/sample_sheet.tsv from split ATAC fragments/peak folders."""
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


def one_file(candidates: list[Path], label: str, sample_id: str) -> Path:
    if len(candidates) == 1:
        return candidates[0]
    checked = "\n".join(f"  - {p}" for p in candidates) or "  - none"
    raise FileNotFoundError(f"Expected one {label} file for {sample_id}; found {len(candidates)}:\n{checked}")


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    params_path = Path(args.params).expanduser() if args.params else pdir / "inputs" / "atac_input_params.tsv"
    params = read_params(params_path)

    layout = params.get("layout", "")
    if layout and layout != "split_ge_arc":
        raise SystemExit(f"ERROR: {params_path} has layout={layout!r}; expected layout=split_ge_arc.")

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
    fragments_root = root / "fragments"
    for arc_dir in sorted(fragments_root.glob("*_arc")):
        arc_id = arc_dir.name
        sample_id = arc_id[:-4]
        frag = one_file(sorted((root / "fragments" / arc_id).glob("*fragments.tsv.gz")), "fragments", sample_id)
        peaks = one_file(sorted((root / "bed" / arc_id).glob("*peaks.bed")), "ATAC peaks BED", sample_id)
        if not Path(str(frag) + ".tbi").exists():
            raise FileNotFoundError(f"Missing fragment index for {sample_id}: {frag}.tbi")

        records.append(
            {
                "sample_id": sample_id,
                "condition": re.sub(r"_[0-9]+$", "", sample_id),
                "organism": organism,
                "fragments_tsv_gz": str(frag),
                "atac_peaks_bed": str(peaks),
            }
        )

    if not records:
        raise SystemExit(f"No *_arc folders found under {fragments_root}")

    out_path = Path(args.out).expanduser() if args.out else pdir / "inputs" / "sample_sheet.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"WROTE {out_path}")
    print(df.to_string(index=False))
    print("\nEdit only the condition column if needed, then run validate_and_prepare_sample_sheet.py.")


if __name__ == "__main__":
    main()
