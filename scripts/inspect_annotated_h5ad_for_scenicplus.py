from __future__ import annotations

import argparse
from pathlib import Path
import os

import anndata as ad
import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
INPUTS = PROJECT / "inputs"
RESULTS = PROJECT / "results" / "annotated_object"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-path", default=os.environ.get("ANNOTATED_OBJECT", ""))
    return parser.parse_args()


def pick_col(columns: list[str], patterns: list[str], default_value: str = "") -> str:
    lower = [c.lower() for c in columns]
    best = ("", 0)
    for col, low in zip(columns, lower):
        score = 0
        for i, pat in enumerate(patterns):
            if pat in low:
                score += len(patterns) - i
        if score > best[1]:
            best = (col, score)
    return best[0] if best[1] > 0 else default_value


def main() -> None:
    args = parse_args()
    if not args.object_path:
        raise SystemExit("Provide --object-path /path/to/annotated.h5ad or set ANNOTATED_OBJECT.")
    path = Path(args.object_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    INPUTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(path, backed="r")
    obs_cols = list(adata.obs.columns)
    sample_col = pick_col(obs_cols, ["sample", "orig.ident", "library", "donor", "replicate"], "sample")
    condition_col = pick_col(obs_cols, ["condition", "genotype", "group", "treatment"], "condition")
    requested_cell_label_column = os.environ.get("CELL_LABEL_COLUMN", "")
    if requested_cell_label_column and requested_cell_label_column in obs_cols:
        cell_label_column = requested_cell_label_column
    else:
        cell_label_column = pick_col(obs_cols, ["cell_label", "cell_type", "celltype", "lineage", "annotation", "cluster"], "cell_label")
    barcode_col = pick_col(obs_cols, ["barcode", "barcodes", "cell_barcode"], "")
    reductions = list(adata.obsm.keys())
    reduction_priority = [
        "X_wnn.umap",
        "wnn.umap",
        "X_umap_wnn",
        "X_rna.umap",
        "rna.umap",
        "X_umap",
    ]
    lower_map = {str(x).lower(): x for x in reductions}
    reduction = ""
    for key in reductions:
        low = str(key).lower()
        if (("final" in low or "eyefinal" in low) and "wnn" in low and "umap" in low):
            reduction = key
            break
    for candidate in reduction_priority:
        if reduction:
            break
        if candidate.lower() in lower_map:
            reduction = lower_map[candidate.lower()]
            break
    if not reduction:
        wnn_hits = [x for x in reductions if "wnn" in str(x).lower() and "umap" in str(x).lower()]
        active_wnn_hits = [
            x
            for x in wnn_hits
            if not any(token in str(x).lower() for token in ["whole", "atlas", "global", "eyefinal", "eyeinitial"])
        ]
        rna_hits = [x for x in reductions if "rna" in str(x).lower() and "umap" in str(x).lower()]
        umap_hits = [x for x in reductions if "umap" in str(x).lower()]
        reduction = (active_wnn_hits or wnn_hits or rna_hits or umap_hits or reductions or [""])[0]

    rows = [
        {"section": "object", "name": "n_cells", "value": str(adata.n_obs)},
        {"section": "object", "name": "n_features", "value": str(adata.n_vars)},
    ]
    rows.extend({"section": "obs_columns", "name": c, "value": "|".join(map(str, adata.obs[c].astype(str).unique()[:8]))} for c in obs_cols)
    rows.extend({"section": "obsm", "name": c, "value": c} for c in reductions)
    pd.DataFrame(rows).to_csv(RESULTS / "annotated_h5ad_candidates.tsv", sep="\t", index=False)

    params = pd.DataFrame(
        [
            ("object_path", str(path)),
            ("object_format", "h5ad"),
            ("sample_col", sample_col),
            ("condition_col", condition_col),
            ("cell_label_column", cell_label_column),
            ("barcode_col", barcode_col),
            ("reduction", reduction),
            ("use_raw", "true"),
        ],
        columns=["parameter", "value"],
    )
    params.to_csv(INPUTS / "annotated_h5ad_params.tsv", sep="\t", index=False)
    selected = pd.DataFrame(
        [
            ("object_format", "h5ad", "Input object format detected from suffix or argument."),
            ("sample_col", sample_col, "Biological sample identifier for matching fragments and pseudobulk grouping."),
            ("condition_col", condition_col, "Condition/group identifier used by downstream differential summaries."),
            ("cell_label_column", cell_label_column, "Metadata column exported as the standardized downstream cell_label."),
            ("barcode_col", barcode_col or "<derive from obs_names>", "Raw barcode column matching fragment-file barcode field, if present."),
            ("reduction", reduction, "Embedding used for inspection and downstream metacell choices."),
            ("use_raw", "true", "Use adata.raw as raw counts when present."),
        ],
        columns=["field", "selected_value", "purpose"],
    )
    selected.to_csv(RESULTS / "annotated_object_selected_fields.tsv", sep="\t", index=False)
    meta_preview = pd.DataFrame(
        {
            "column": obs_cols,
            "dtype": [str(adata.obs[c].dtype) for c in obs_cols],
            "n_unique": [int(adata.obs[c].astype(str).nunique()) for c in obs_cols],
            "example_values": [
                " | ".join(map(str, adata.obs[c].astype(str).unique()[:8]))
                for c in obs_cols
            ],
        }
    )
    meta_preview.to_csv(RESULTS / "annotated_object_metadata_preview.tsv", sep="\t", index=False)
    if cell_label_column in adata.obs.columns:
        preview = adata.obs[cell_label_column].astype(str).value_counts().rename_axis("label").reset_index(name="n_cells")
        preview.to_csv(RESULTS / "annotated_h5ad_label_preview.tsv", sep="\t", index=False)
    report = "\n".join(
        [
            "# Annotated Object Inspection Report",
            "",
            "## Object",
            f"- Path: `{path}`",
            "- Format: `h5ad`",
            f"- Cells: {adata.n_obs}",
            f"- Features: {adata.n_vars}",
            f"- Raw counts present: {adata.raw is not None}",
            f"- obsm entries: {', '.join(reductions)}",
            "",
            "## Automatically Selected Fields",
            f"- sample_col: `{sample_col}`",
            f"- condition_col: `{condition_col}`",
            f"- CELL_LABEL_COLUMN: `{cell_label_column}`",
            f"- barcode_col: `{barcode_col or '<derive from obs_names>'}`",
            f"- reduction: `{reduction}`",
            "",
            "## Files Written",
            "- `inputs/annotated_h5ad_params.tsv`",
            "- `results/annotated_object/annotated_object_selected_fields.tsv`",
            "- `results/annotated_object/annotated_object_metadata_preview.tsv`",
            "- `results/annotated_object/annotated_h5ad_label_preview.tsv`",
        ]
    )
    (RESULTS / "annotated_object_inspection_report.md").write_text(report + "\n")
    print("WROTE", INPUTS / "annotated_h5ad_params.tsv")
    print("WROTE", RESULTS / "annotated_h5ad_candidates.tsv")
    print("WROTE", RESULTS / "annotated_object_selected_fields.tsv")
    print("WROTE", RESULTS / "annotated_object_metadata_preview.tsv")
    print("WROTE", RESULTS / "annotated_object_inspection_report.md")
    print()
    print("ANNOTATED OBJECT SUMMARY")
    print("cells:", adata.n_obs)
    print("features:", adata.n_vars)
    print("selected sample_col:", sample_col)
    print("selected condition_col:", condition_col)
    print("selected CELL_LABEL_COLUMN:", cell_label_column)
    print("selected reduction:", reduction)


if __name__ == "__main__":
    main()
