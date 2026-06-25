#!/usr/bin/env python
"""Check whether required pycisTopic outputs are complete for SCENIC+."""
from __future__ import annotations

import glob
import os
import sys
import argparse
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


REQUIRED_PATTERNS = [
    ("cisTopic object", "inputs/cistopic_obj.pkl", "file", 1),
    ("consensus peaks BED", "work/pycistopic/consensus_peaks.bed", "file", 1),
    ("Topics otsu BEDs", "inputs/region_sets/Topics_otsu/*.bed", "glob", 1),
    ("Topics top BEDs", "inputs/region_sets/Topics_top_3k/*.bed", "glob", 1),
    ("DAR BEDs", "inputs/region_sets/DARs_cell_label/*.bed", "glob", 1),
    ("consensus peaks table", "results/pycistopic/consensus_peaks.tsv", "file", 1),
    ("cell QC metrics", "results/pycistopic/qc/cell_qc_metrics.tsv", "file", 1),
    ("cell QC PDF", "results/pycistopic/qc/cistopic_cell_qc.pdf", "file", 1),
    ("consensus peak QC PDF", "results/pycistopic/qc/consensus_peak_qc.pdf", "file", 1),
    ("cells by sample and label", "results/pycistopic/qc/cells_by_sample_and_label.tsv", "file", 1),
    ("parallelism plan", "results/pycistopic/qc/parallelism_plan.tsv", "file", 1),
    ("pseudobulk resume plan", "results/pycistopic/qc/pseudobulk_resume_plan.tsv", "file", 1),
    ("DAR result tables", "results/pycistopic/dar/*.tsv", "glob", 1),
    ("DAR summary PDF", "results/pycistopic/dar/dar_summary.pdf", "file", 1),
    ("topic model metrics table", "results/pycistopic/model_selection/topic_model_metrics.tsv", "file", 1),
    ("topic model metrics PDF", "results/pycistopic/model_selection/topic_model_metrics.pdf", "file", 1),
    ("topic QC metrics table", "results/pycistopic/model_selection/topic_qc_metrics.tsv", "file", 1),
    ("topic QC metrics PDF", "results/pycistopic/model_selection/topic_qc_metrics.pdf", "file", 1),
    ("topic otsu thresholds PDF", "results/pycistopic/model_selection/topic_otsu_thresholds.pdf", "file", 1),
    ("topic region-set summary table", "results/pycistopic/model_selection/topic_region_set_summary.tsv", "file", 1),
    ("topic region-set summary PDF", "results/pycistopic/model_selection/topic_region_set_summary.pdf", "file", 1),
    ("selected model record", "results/pycistopic/model_selection/selected_model.txt", "file", 1),
    ("pycisTopic manifest", "results/pycistopic/pycistopic_manifest.json", "file", 1),
]


def resolve(pattern: str, kind: str) -> list[Path]:
    if kind == "glob":
        return [Path(x) for x in glob.glob(str(PROJECT / pattern))]
    return [PROJECT / pattern]


def nonempty(paths: list[Path]) -> list[Path]:
    return [x for x in paths if x.exists() and x.is_file() and x.stat().st_size > 0]


def main() -> None:
    argparse.ArgumentParser(
        description=(
            "Check that all pycisTopic outputs required by downstream SCENIC+ "
            "exist and are non-empty. Writes a markdown and TSV completion "
            "report under results/pycistopic/qc."
        )
    ).parse_args()
    rows = []
    for label, pattern, kind, min_count in REQUIRED_PATTERNS:
        paths = sorted(resolve(pattern, kind))
        good = nonempty(paths)
        status = "ok" if len(good) >= min_count else "missing"
        rows.append(
            {
                "check": label,
                "pattern": pattern,
                "min_count": min_count,
                "observed_count": len(good),
                "status": status,
                "example_paths": ";".join(str(x) for x in good[:5]),
            }
        )

    out_dir = PROJECT / "results" / "pycistopic" / "qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    out_tsv = out_dir / "pycistopic_completion_check.tsv"
    out_md = out_dir / "pycistopic_completion_check.md"
    df.to_csv(out_tsv, sep="\t", index=False)
    failed = df[df["status"] != "ok"]
    lines = ["# pycisTopic Completion Check", ""]
    lines.append(f"PROJECT_DIR: `{PROJECT}`")
    lines.append("")
    lines.append(df.loc[:, ["check", "observed_count", "min_count", "status"]].to_markdown(index=False))
    lines.append("")
    lines.append(f"Machine-readable report: `{out_tsv}`")
    if failed.empty:
        lines.extend(["", "Status: OK. Required pycisTopic outputs for SCENIC+ are present."])
    else:
        lines.extend(["", "Status: FAILED. Missing or empty required outputs:", ""])
        for _, row in failed.iterrows():
            lines.append(f"- {row['check']}: `{row['pattern']}`")
    out_md.write_text("\n".join(lines) + "\n")
    print(out_md.read_text())
    if not failed.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
