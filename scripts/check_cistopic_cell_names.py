#!/usr/bin/env python
"""Check cisTopic cell names against the active RNA matrix and metadata."""
from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT / "tmp" / "numba"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)
(PROJECT / "tmp" / "numba").mkdir(parents=True, exist_ok=True)

import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-overlap", type=float, default=None, help="Minimum direct RNA/cisTopic cell-name overlap fraction.")
    parser.add_argument("--out", default="results/pycistopic/qc/cistopic_cell_name_check.tsv")
    return parser.parse_args()


def p(path: str | Path) -> Path:
    out = Path(path).expanduser()
    return out if out.is_absolute() else PROJECT / out


def read_threshold(default: float = 0.50) -> float:
    thresholds = p("inputs/preflight_thresholds.tsv")
    if not thresholds.exists():
        return default
    df = pd.read_csv(thresholds, sep="\t")
    if {"metric", "value"}.issubset(df.columns):
        values = df.set_index("metric")["value"].to_dict()
        for key in ["min_overlap_rna", "min_overlap_atac"]:
            if key in values:
                return float(values[key])
    return default


def cistopic_cell_names(obj) -> list[str]:
    if hasattr(obj, "cell_names"):
        return list(map(str, obj.cell_names))
    if hasattr(obj, "cell_data"):
        return list(map(str, obj.cell_data.index))
    raise ValueError("Cannot find cell names in cisTopic object.")


def main() -> None:
    args = parse_args()
    min_overlap = args.min_overlap if args.min_overlap is not None else read_threshold()
    cistopic_path = p("inputs/cistopic_obj.pkl")
    gex_path = p("inputs/gex.h5ad")
    metadata_path = p("inputs/cell_metadata.tsv")
    sample_sheet_path = p("inputs/sample_sheet.tsv")
    for path in [cistopic_path, gex_path, metadata_path, sample_sheet_path]:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    with cistopic_path.open("rb") as handle:
        cistopic = pickle.load(handle)
    atac_names = cistopic_cell_names(cistopic)
    adata = sc.read_h5ad(gex_path, backed="r")
    rna_names = list(map(str, adata.obs_names))
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna("")
    sample_sheet = pd.read_csv(sample_sheet_path, sep="\t", dtype=str).fillna("")

    atac = set(atac_names)
    rna = set(rna_names)
    meta = set(metadata["cell_id"].astype(str)) if "cell_id" in metadata.columns else set()
    sample_ids = set(sample_sheet["sample_id"].astype(str)) if "sample_id" in sample_sheet.columns else set()
    overlap = atac & rna
    metadata_overlap = atac & meta
    suffix_matches = [
        name for name in atac_names
        if sample_ids and any(name.endswith(f"-{sample_id}") for sample_id in sample_ids)
    ]
    rows = [
        ("n_cistopic_cells", len(atac_names)),
        ("n_rna_cells", len(rna_names)),
        ("n_metadata_cells", len(meta)),
        ("n_direct_rna_cistopic_overlap", len(overlap)),
        ("fraction_cistopic_in_rna", len(overlap) / max(len(atac), 1)),
        ("fraction_rna_in_cistopic", len(overlap) / max(len(rna), 1)),
        ("fraction_cistopic_in_metadata", len(metadata_overlap) / max(len(atac), 1)),
        ("fraction_cistopic_with_sample_suffix", len(suffix_matches) / max(len(atac_names), 1)),
        ("example_cistopic_cells", ",".join(atac_names[:10])),
        ("example_rna_cells", ",".join(rna_names[:10])),
    ]
    out = p(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(out, sep="\t", index=False)
    print(pd.DataFrame(rows, columns=["metric", "value"]).to_string(index=False))

    frac_atac = len(overlap) / max(len(atac), 1)
    frac_rna = len(overlap) / max(len(rna), 1)
    if frac_atac < min_overlap or frac_rna < min_overlap:
        raise AssertionError(
            f"RNA/cisTopic direct cell-name overlap is too low: "
            f"fraction_cistopic_in_rna={frac_atac:.4f}, fraction_rna_in_cistopic={frac_rna:.4f}, "
            f"minimum={min_overlap:.4f}. Rebuild the active RNA or cisTopic object with matching cell IDs."
        )
    print(f"WROTE {out}")
    print("cisTopic cell-name check OK")


if __name__ == "__main__":
    main()
