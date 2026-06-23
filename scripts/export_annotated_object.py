#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
INPUTS = PROJECT / "inputs"
SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cell-label-column",
        dest="cell_label_column",
        default=os.environ.get("CELL_LABEL_COLUMN", ""),
        help="Metadata column in the annotated object to export as standardized cell_label.",
    )
    return parser.parse_args()


def write_project_setting(key: str, value: str) -> None:
    config = PROJECT / "scenicplus_project.env"
    lines: list[str] = []
    if config.exists():
        lines = config.read_text().splitlines()
    assignment = f"{key}={value!r}"
    replaced = False
    new_lines = []
    for line in lines:
        if line.startswith(f"{key}="):
            new_lines.append(assignment)
            replaced = True
        else:
            new_lines.append(line)
    if not replaced:
        new_lines.append(assignment)
    config.write_text("\n".join(new_lines) + "\n")


def update_param_table(path: Path, updates: dict[str, str]) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    params = dict(zip(df["parameter"].astype(str), df["value"].astype(str)))
    params.update(updates)
    ordered = list(df["parameter"].astype(str))
    for key in updates:
        if key not in ordered:
            ordered.append(key)
    out = pd.DataFrame({"parameter": ordered, "value": [params[k] for k in ordered]})
    out.to_csv(path, sep="\t", index=False)
    return params


def update_selected_fields(params: dict[str, str]) -> None:
    path = PROJECT / "results" / "annotated_object" / "annotated_object_selected_fields.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if not {"field", "selected_value"}.issubset(df.columns):
        return
    for key in [
        "assay",
        "layer",
        "sample_col",
        "condition_col",
        "cell_label_column",
        "barcode_col",
        "original_cell_id_col",
        "reduction",
    ]:
        if key in params:
            hit = df["field"].astype(str) == key
            if hit.any():
                df.loc[hit, "selected_value"] = params.get(key, "")
    df.to_csv(path, sep="\t", index=False)
    print("UPDATED", path)


def read_active_params(cell_label_column: str) -> tuple[str, Path, dict[str, str]]:
    h5ad_params = INPUTS / "annotated_h5ad_params.tsv"
    seurat_params = INPUTS / "annotated_object_params.tsv"
    if h5ad_params.exists() and seurat_params.exists():
        h5ad_mtime = h5ad_params.stat().st_mtime
        seurat_mtime = seurat_params.stat().st_mtime
        path = h5ad_params if h5ad_mtime >= seurat_mtime else seurat_params
    elif h5ad_params.exists():
        path = h5ad_params
    elif seurat_params.exists():
        path = seurat_params
    else:
        raise FileNotFoundError(
            "No annotated object parameter table found. Run inspect_annotated_object.py first."
        )
    updates: dict[str, str] = {}
    if cell_label_column:
        updates["cell_label_column"] = cell_label_column
    params = update_param_table(path, updates)
    fmt = params.get("object_format", "").lower().lstrip(".")
    if fmt not in {"rds", "qs", "h5ad"}:
        raise ValueError(f"Unsupported or missing object_format in {path}: {fmt!r}")
    return fmt, path, params


def write_label_summary() -> None:
    metadata = INPUTS / "cell_metadata.tsv"
    if not metadata.exists():
        raise FileNotFoundError(metadata)
    df = pd.read_csv(metadata, sep="\t", dtype=str).fillna("")
    required = {"cell_id", "sample_id", "condition", "cell_label"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{metadata} lacks required columns after export: {', '.join(missing)}")
    if not df["cell_label"].astype(bool).all():
        bad = df.loc[~df["cell_label"].astype(bool), "cell_id"].head(10).tolist()
        raise ValueError(f"cell_label contains missing/empty values; examples: {bad}")
    summary = (
        df.groupby(["sample_id", "condition", "cell_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["sample_id", "condition", "cell_label"])
    )
    summary.to_csv(INPUTS / "grn_label_summary.tsv", sep="\t", index=False)
    print("WROTE", INPUTS / "grn_label_summary.tsv")


def refresh_single_cell_metadata_backup() -> None:
    source = INPUTS / "cell_metadata.tsv"
    backup = INPUTS / "cell_metadata.single_cell.tsv"
    if not source.exists():
        raise FileNotFoundError(source)
    shutil.copyfile(source, backup)
    print("WROTE", backup)


def main() -> None:
    args = parse_args()
    fmt, params_path, params = read_active_params(args.cell_label_column)
    write_project_setting("ANNOTATED_OBJECT", params.get("object_path", ""))
    write_project_setting("ANNOTATED_OBJECT_FORMAT", fmt)
    if not params.get("cell_label_column"):
        raise ValueError(f"{params_path} must contain cell_label_column. Rerun inspect_annotated_object.py.")
    write_project_setting("CELL_LABEL_COLUMN", params["cell_label_column"])
    write_project_setting("ANNOTATED_OBJECT_PARAMS", str(params_path))
    if fmt in {"rds", "qs"}:
        commands = [
            ["Rscript", str(SCRIPT_DIR / "export_annotated_seurat_for_scenicplus.R")],
            [sys.executable, str(SCRIPT_DIR / "make_annotated_seurat_gex_h5ad.py")],
        ]
    elif fmt == "h5ad":
        commands = [[sys.executable, str(SCRIPT_DIR / "export_annotated_h5ad_for_scenicplus.py")]]
    else:
        raise SystemExit(f"Unsupported annotated object format: {fmt}")
    for cmd in commands:
        print("RUN", " ".join(cmd))
        subprocess.run(cmd, check=True)
    update_selected_fields(params)
    refresh_single_cell_metadata_backup()
    write_label_summary()
    write_project_setting("ACTIVE_GEX_H5AD", str(INPUTS / "gex.h5ad"))
    write_project_setting("ACTIVE_CELL_METADATA", str(INPUTS / "cell_metadata.tsv"))
    print("WROTE", PROJECT / "scenicplus_project.env")


if __name__ == "__main__":
    main()
