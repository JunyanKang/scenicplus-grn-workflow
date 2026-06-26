#!/usr/bin/env python
"""Review annotated-object inspection outputs before export."""
from __future__ import annotations

import os
import sys
import argparse
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def p(path: str) -> Path:
    return PROJECT / path


def main() -> None:
    argparse.ArgumentParser(
        description=(
            "Review the annotated-object inspection report and selected export "
            "parameters before writing SCENIC+ RNA inputs. Fails if required "
            "fields such as cell label, assay/layer, sample or condition are missing."
        )
    ).parse_args()
    report = p("results/annotated_object/annotated_object_inspection_report.md")
    params_candidates = [
        p("inputs/annotated_object_params.tsv"),
        p("inputs/annotated_h5ad_params.tsv"),
    ]
    params = next((x for x in params_candidates if x.exists() and x.stat().st_size > 0), None)
    out_dir = p("results/annotated_object")
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_md = out_dir / "annotated_object_pre_export_review.md"

    failures: list[str] = []
    if not report.exists() or report.stat().st_size == 0:
        failures.append(f"Missing inspection report: {report}")
    if params is None:
        failures.append("Missing annotated-object parameter table under inputs/.")

    lines: list[str] = ["# Annotated Object Pre-Export Review", ""]
    lines.append(f"PROJECT_DIR: `{PROJECT}`")
    lines.append("")
    if report.exists():
        lines.append(f"Inspection report: `{report}`")
    if params is not None:
        lines.append(f"Parameter table: `{params}`")
    lines.append("")

    if params is not None:
        df = pd.read_csv(params, sep="\t", dtype=str).fillna("")
        key_col, value_col = df.columns[:2]
        values = {str(k): str(v) for k, v in zip(df[key_col], df[value_col])}
        object_format = values.get("object_format", "").lower().lstrip(".")
        required = ["cell_label_column", "reduction", "sample_col", "condition_col"]
        if object_format in {"rds", "qs"}:
            required.extend(["assay", "layer"])
        rows = []
        for key in required:
            value = values.get(key, "")
            rows.append({"parameter": key, "value": value, "status": "ok" if value else "missing"})
            if not value:
                failures.append(f"Parameter `{key}` is empty in {params}")
        review_tsv = out_dir / "annotated_object_pre_export_review.tsv"
        pd.DataFrame(rows).to_csv(review_tsv, sep="\t", index=False)
        lines.extend(["## Required Fields", ""])
        lines.append(pd.DataFrame(rows).to_markdown(index=False))
        lines.extend(["", f"Machine-readable review: `{review_tsv}`", ""])

    if report.exists() and report.stat().st_size > 0:
        lines.extend(["## Inspection Report", "", report.read_text()])

    if failures:
        lines.extend(["", "## Blocking Issues", ""])
        lines.extend([f"- {x}" for x in failures])
    else:
        lines.extend(["", "## Status", "", "OK: inspection report and required export parameters are present."])
    summary_md.write_text("\n".join(lines) + "\n")
    print(summary_md.read_text())
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
