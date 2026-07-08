#!/usr/bin/env python
"""Run standard SCENIC+ downstream statistics and figures from project params."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run SCENIC+ postprocessing from project parameters."
        )
    )
    parser.add_argument(
        "--layer",
        choices=["direct", "extended", "all"],
        default="all",
        help="eRegulon layer to render. The audit task ignores this value.",
    )
    parser.add_argument(
        "--task",
        choices=["audit", "figures", "stats", "all"],
        default="all",
        help="audit=output tier checks; figures=01-08 PDFs/source tables; stats=09+ condition AUC tables/PDFs; all=complete postprocess.",
    )
    parser.add_argument("--params", default=None, help="Default: $PROJECT_DIR/inputs/postprocess_params.tsv")
    return parser.parse_args()


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        subprocess.run([sys.executable, str(SCRIPT_DIR / "setup_workflow_params.py"), "--section", "postprocess"], check=True)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {str(k): str(v) for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])}


def run_logged(cmd: list[str], log_name: str) -> None:
    logs = PROJECT / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    print("RUN", " ".join(cmd))
    with (logs / log_name).open("w") as log:
        proc = subprocess.Popen(cmd, cwd=PROJECT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        code = proc.wait()
    if code != 0:
        raise SystemExit(code)


def require_inputs(params: dict[str, str], task: str, layers: list[str]) -> None:
    if task == "audit":
        return
    needed = [params.get("metadata", "")]
    if task in {"figures", "all"}:
        needed.append(params.get("tf_to_gene", ""))
        for layer in layers:
            needed.extend([params.get(f"{layer}_auc_h5mu", ""), params.get(f"{layer}_eregulons", "")])
    if task in {"stats", "all"}:
        for layer in layers:
            needed.append(params.get(f"{layer}_auc_h5mu", ""))
    missing = []
    for value in needed:
        if not value:
            continue
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = (PROJECT / path).resolve()
        if not path.exists() or path.stat().st_size == 0:
            missing.append(str(path))
    if missing:
        raise SystemExit("ERROR: SCENIC+ postprocess inputs are not ready. Missing: " + ", ".join(dict.fromkeys(missing)))


def main() -> None:
    args = parse_args()
    params_path = Path(args.params).expanduser() if args.params else PROJECT / "inputs" / "postprocess_params.tsv"
    if not params_path.is_absolute():
        params_path = (PROJECT / params_path).resolve()
    params = read_params(params_path)
    layers = ["direct", "extended"] if args.layer == "all" else [args.layer]
    require_inputs(params, args.task, layers)
    if args.task in {"audit", "all"}:
        run_logged(
            [
                sys.executable,
                str(SCRIPT_DIR / "audit_scenicplus_output_tiers.py"),
            ],
            "audit_scenicplus_output_tiers.log",
        )
        if args.task == "audit":
            return
    def run_figures_layer(layer: str) -> None:
        run_logged(
            [
                sys.executable,
                str(SCRIPT_DIR / "extract_scenicplus_plot_data.py"),
                "--auc-h5mu",
                params[f"{layer}_auc_h5mu"],
                "--metadata",
                params["metadata"],
                "--tf-to-gene",
                params["tf_to_gene"],
                "--eregulons",
                params[f"{layer}_eregulons"],
                "--group-col",
                params.get("group_col", "cell_label"),
                "--condition-col",
                params.get("condition_col", "condition"),
                "--cell-col",
                params.get("cell_col", "cell_id"),
                "--umap-x",
                params.get("umap_x", "umap_1"),
                "--umap-y",
                params.get("umap_y", "umap_2"),
                "--top-n",
                params.get("plot_top_n", "30"),
                "--umap-n",
                params.get("plot_umap_n", "12"),
                "--network-top-tfs",
                params.get("network_top_tfs", "8"),
                "--network-targets-per-tf",
                params.get("network_targets_per_tf", "12"),
                "--regulon-sign-filter",
                params.get("regulon_sign_filter", "tf_positive"),
                "--outdir",
                params[f"{layer}_figures_outdir"],
                "--file-suffix",
                layer,
            ],
            f"extract_scenicplus_plot_data_{layer}.log",
        )
        run_logged(
            [
                "Rscript",
                str(SCRIPT_DIR / "plot_scenicplus_publication_outputs.R"),
                "--outdir",
                params[f"{layer}_figures_outdir"],
                "--file-suffix",
                layer,
                "--plot-style-config",
                params.get("plot_style_config", "results/scenicplus_figures/plot_style_parameters.tsv"),
            ],
            f"plot_scenicplus_publication_outputs_R_{layer}.log",
        )

    def run_stats_layer(layer: str) -> None:
        run_logged(
            [
                sys.executable,
                str(SCRIPT_DIR / "test_eregulon_auc_by_condition.py"),
                "--auc-h5mu",
                params[f"{layer}_auc_h5mu"],
                "--metadata",
                params["metadata"],
                "--group-col",
                params.get("condition_col", "condition"),
                "--label-col",
                params.get("group_col", "cell_label"),
                "--sample-col",
                params.get("sample_col", "sample_id"),
                "--cell-col",
                params.get("cell_col", "cell_id"),
                "--reference-condition",
                params.get("reference_condition", "auto"),
                "--comparison-condition",
                params.get("comparison_condition", "auto"),
                "--regulon-sign-filter",
                params.get("regulon_sign_filter", "tf_positive"),
                "--outdir",
                params[f"{layer}_stats_outdir"],
                "--file-suffix",
                layer,
                "--tables-only",
            ],
            f"test_eregulon_auc_by_condition_{layer}.log",
        )
        run_logged(
            [
                "Rscript",
                str(SCRIPT_DIR / "plot_scenicplus_condition_stats.R"),
                "--outdir",
                params[f"{layer}_stats_outdir"],
                "--file-suffix",
                layer,
                "--plot-style-config",
                params.get("plot_style_config", "results/scenicplus_figures/plot_style_parameters.tsv"),
                "--priority-eregulons",
                params.get("priority_eregulons", ""),
            ],
            f"plot_scenicplus_condition_stats_R_{layer}.log",
        )

    if args.task in {"figures", "all"}:
        for layer in layers:
            run_figures_layer(layer)
    if args.task in {"stats", "all"}:
        for layer in layers:
            run_stats_layer(layer)


if __name__ == "__main__":
    main()
