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
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", choices=["direct", "extended", "all"], default="all")
    parser.add_argument("--task", choices=["stats", "figures", "all"], default="all")
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


def main() -> None:
    args = parse_args()
    params_path = Path(args.params).expanduser() if args.params else PROJECT / "inputs" / "postprocess_params.tsv"
    if not params_path.is_absolute():
        params_path = (PROJECT / params_path).resolve()
    params = read_params(params_path)
    layers = ["direct", "extended"] if args.layer == "all" else [args.layer]
    for layer in layers:
        if args.task in {"stats", "all"}:
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
                    "--sample-col",
                    params.get("sample_col", "sample_id"),
                    "--cell-col",
                    params.get("cell_col", "cell_id"),
                    "--outdir",
                    params[f"{layer}_stats_outdir"],
                ],
                f"test_eregulon_auc_by_condition_{layer}.log",
            )
        if args.task in {"figures", "all"}:
            run_logged(
                [
                    sys.executable,
                    str(SCRIPT_DIR / "plot_scenicplus_publication_outputs.py"),
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
                    "--cell-col",
                    params.get("cell_col", "cell_id"),
                    "--top-n",
                    params.get("plot_top_n", "30"),
                    "--umap-n",
                    params.get("plot_umap_n", "12"),
                    "--network-top-tfs",
                    params.get("network_top_tfs", "8"),
                    "--network-targets-per-tf",
                    params.get("network_targets_per_tf", "12"),
                    "--outdir",
                    params[f"{layer}_figures_outdir"],
                ],
                f"plot_scenicplus_publication_outputs_{layer}.log",
            )


if __name__ == "__main__":
    main()
