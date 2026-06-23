#!/usr/bin/env python
"""Record software versions and key resource checksums for a SCENIC+ project."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="logs")
    return parser.parse_args()


def p(path: str | Path) -> Path:
    out = Path(path).expanduser()
    return out if out.is_absolute() else PROJECT / out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_text(cmd: list[str], out: Path) -> str:
    exe = shutil.which(cmd[0])
    if exe is None:
        out.write_text(f"{cmd[0]} not found on PATH\n")
        return "not_found"
    proc = subprocess.run([exe, *cmd[1:]], text=True, capture_output=True, check=False)
    out.write_text((proc.stdout or "") + (proc.stderr or ""))
    return f"exit_{proc.returncode}"


def main() -> None:
    args = parse_args()
    out_dir = p(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    packages = ["scenicplus", "scanpy", "anndata", "pycisTopic", "pycistarget", "mudata", "sklearn", "yaml"]
    version_rows = []
    for pkg in packages:
        try:
            mod = importlib.import_module(pkg)
            version_rows.append({"package": pkg, "status": "ok", "version": getattr(mod, "__version__", "version_unknown")})
        except Exception as exc:
            version_rows.append({"package": pkg, "status": "not_importable", "version": repr(exc)})
    versions = pd.DataFrame(version_rows)
    versions.to_csv(out_dir / "python_package_versions.tsv", sep="\t", index=False)
    (out_dir / "python_package_versions.txt").write_text(versions.to_string(index=False) + "\n")

    command_rows = [
        {"command": "pip freeze", "output": str(out_dir / "pip_freeze.txt"), "status": run_text([sys.executable, "-m", "pip", "freeze"], out_dir / "pip_freeze.txt")},
        {"command": "conda list", "output": str(out_dir / "conda_list.txt"), "status": run_text(["conda", "list"], out_dir / "conda_list.txt")},
    ]

    resource_candidates = [
        p("inputs/cistarget_db/motif_annotations.tbl"),
        p("resources/motifs/motifs.txt"),
        p("resources/resource_manifest.json"),
        p("resources/resource_status.tsv"),
        p("work/scenicplus/Snakemake/config/config.yaml"),
        p("work/scenicplus/organism_config.yaml"),
        p("results/cistarget_db/custom_cistarget_db_manifest.tsv"),
    ]
    checksum_rows = []
    for path in resource_candidates:
        checksum_rows.append(
            {
                "path": str(path),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
                "sha256": sha256(path) if path.exists() and path.is_file() else "",
            }
        )
    checksums = pd.DataFrame(checksum_rows)
    checksums.to_csv(out_dir / "resource_sha256.tsv", sep="\t", index=False)
    command_rows.append({"command": "resource sha256", "output": str(out_dir / "resource_sha256.tsv"), "status": "ok"})

    pd.DataFrame(command_rows).to_csv(out_dir / "provenance_commands.tsv", sep="\t", index=False)
    print("WROTE", out_dir / "python_package_versions.tsv")
    print("WROTE", out_dir / "pip_freeze.txt")
    print("WROTE", out_dir / "conda_list.txt")
    print("WROTE", out_dir / "resource_sha256.tsv")
    print("WROTE", out_dir / "provenance_commands.tsv")


if __name__ == "__main__":
    main()
