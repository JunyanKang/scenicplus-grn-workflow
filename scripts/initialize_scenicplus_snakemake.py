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
        subprocess.run([scenicplus, "init_snakemake", "--out_dir", str(out_dir)], check=True)
    else:
        print(f"KEPT {snake_dir}")
    subprocess.run([sys.executable, str(SCRIPT_DIR / "prepare_scenicplus_organism_files.py")], check=True)
    subprocess.run([sys.executable, str(SCRIPT_DIR / "setup_workflow_params.py"), "--section", "scenicplus_config"], check=True)
    subprocess.run([sys.executable, str(SCRIPT_DIR / "generate_scenicplus_config.py")], check=True)


if __name__ == "__main__":
    main()
