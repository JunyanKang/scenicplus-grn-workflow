#!/usr/bin/env python
"""Run SCENIC+ Snakemake using project parameter tables."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["dryrun", "run"], default="dryrun")
    parser.add_argument("--params", default=None, help="Default: $PROJECT_DIR/inputs/snakemake_params.tsv")
    parser.add_argument("--cores", type=int, default=None, help="Override inputs/snakemake_params.tsv.")
    return parser.parse_args()


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        print("Parameter table missing; creating defaults with setup_workflow_params.py", file=sys.stderr)
        subprocess.run([sys.executable, str(SCRIPT_DIR / "setup_workflow_params.py"), "--section", "snakemake"], check=True)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {str(k): str(v) for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])}


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def main() -> None:
    args = parse_args()
    params_path = Path(args.params).expanduser() if args.params else PROJECT / "inputs" / "snakemake_params.tsv"
    if not params_path.is_absolute():
        params_path = (PROJECT / params_path).resolve()
    params = read_params(params_path)
    cores = args.cores if args.cores is not None else int(params.get("cores", "1"))
    snake_dir = PROJECT / "work" / "scenicplus" / "Snakemake"
    if not snake_dir.exists():
        raise FileNotFoundError(snake_dir)
    logs = PROJECT / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log_name = "scenicplus_dryrun.log" if args.mode == "dryrun" else f"scenicplus_run_{datetime.now():%Y%m%d_%H%M%S}.log"
    cmd = ["snakemake", "--cores", str(cores)]
    if args.mode == "dryrun":
        cmd.insert(1, "-n")
        cmd[cmd.index("--cores") + 1] = "1"
    if args.mode == "run" and truthy(params.get("rerun_incomplete", "1")):
        cmd.append("--rerun-incomplete")
    if truthy(params.get("printshellcmds", "1")):
        cmd.append("--printshellcmds")
    if params.get("latency_wait"):
        cmd.extend(["--latency-wait", str(params["latency_wait"])])
    print("RUN", " ".join(cmd))
    with (logs / log_name).open("w") as log:
        proc = subprocess.Popen(cmd, cwd=snake_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        code = proc.wait()
    if code != 0:
        raise SystemExit(code)


if __name__ == "__main__":
    main()
