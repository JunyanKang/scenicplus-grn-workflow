#!/usr/bin/env python
"""Initialize SCENIC+ Snakemake and generate the project config."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
SCRIPT_DIR = Path(__file__).resolve().parent
RUNNER = SCRIPT_DIR.parent / "bin" / "run_python_entrypoint.py"


def patch_snakefile_for_stage_threads(path: Path) -> None:
    text = path.read_text()
    original = text
    text = text.replace(
        'threads: config["params_general"]["n_cpu"]\n    shell:\n        """\n            scenicplus grn_inference motif_enrichment_cistarget',
        'threads: config["params_motif_enrichment"].get("ctx_n_cpu", config["params_general"]["n_cpu"])\n    shell:\n        """\n            scenicplus grn_inference motif_enrichment_cistarget',
    )
    text = text.replace(
        'threads: config["params_general"]["n_cpu"]\n        shell:\n            """\n                scenicplus grn_inference motif_enrichment_dem',
        'threads: config["params_motif_enrichment"].get("dem_n_cpu", config["params_general"]["n_cpu"])\n        shell:\n            """\n                scenicplus grn_inference motif_enrichment_dem',
    )
    if text != original:
        path.write_text(text)
        print(f"PATCHED stage-specific motif-enrichment threads in {path}")


def run_script(script_name: str, *args: str) -> None:
    cmd = [sys.executable]
    if RUNNER.exists():
        cmd.extend([str(RUNNER), str(SCRIPT_DIR / script_name)])
    else:
        cmd.append(str(SCRIPT_DIR / script_name))
    cmd.extend(args)
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)


def main() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(
            "Usage: initialize_scenicplus_snakemake.py\n\n"
            "Initialize the SCENIC+ Snakemake workspace in $PROJECT_DIR/work/scenicplus, "
            "prepare organism files, create default SCENIC+ config parameters and "
            "write the project config."
        )
        return
    out_dir = PROJECT / "work" / "scenicplus"
    snake_dir = out_dir / "Snakemake"
    if not snake_dir.exists():
        scenicplus = shutil.which("scenicplus")
        if scenicplus is None:
            raise FileNotFoundError("scenicplus command not found on PATH")
        env = os.environ.copy()
        env.setdefault("PYTHONWARNINGS", "ignore")
        subprocess.run([scenicplus, "init_snakemake", "--out_dir", str(out_dir)], check=True, env=env)
    else:
        print(f"KEPT {snake_dir}")
    patch_snakefile_for_stage_threads(snake_dir / "workflow" / "Snakefile")
    run_script("prepare_scenicplus_organism_files.py")
    run_script("setup_workflow_params.py", "--section", "scenicplus_config")
    run_script("generate_scenicplus_config.py")


if __name__ == "__main__":
    main()
