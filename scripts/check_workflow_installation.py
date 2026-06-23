#!/usr/bin/env python
"""Check that the installed SCENIC+ workflow assets are complete."""
from __future__ import annotations

import argparse
import importlib
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path


REQUIRED_SCRIPTS = [
    "prepare_official_resources.py",
    "init_scenicplus_project.py",
    "set_atac_input_params.py",
    "make_sample_sheet_from_atac_inputs.py",
    "make_sample_sheet_from_cellranger_outs_layout.py",
    "make_sample_sheet_from_split_ge_arc_layout.py",
    "validate_and_prepare_sample_sheet.py",
    "setup_workflow_params.py",
    "normalize_bedlike_to_ucsc_standard.py",
    "standardize_atac_inputs.py",
    "inspect_annotated_object.py",
    "export_annotated_object.py",
    "inspect_annotated_object_for_scenicplus.R",
    "export_annotated_seurat_for_scenicplus.R",
    "inspect_annotated_h5ad_for_scenicplus.py",
    "export_annotated_h5ad_for_scenicplus.py",
    "make_annotated_seurat_gex_h5ad.py",
    "prepare_metacell_inputs_from_seurat.R",
    "make_metacell_gex_h5ad.py",
    "reassign_fragments_to_metacells.py",
    "validate_pseudobulk_files.py",
    "run_pycistopic_workflow.py",
    "check_cistopic_cell_names.py",
    "standardize_region_sets.py",
    "build_custom_cistarget_db.py",
    "prepare_scenicplus_organism_files.py",
    "initialize_scenicplus_snakemake.py",
    "generate_scenicplus_config.py",
    "preflight_scenicplus_inputs.py",
    "run_scenicplus_snakemake.py",
    "record_scenicplus_provenance.py",
    "compare_scenicplus_stability.py",
    "test_eregulon_auc_by_condition.py",
    "plot_scenicplus_publication_outputs.py",
    "run_scenicplus_postprocess.py",
]

REQUIRED_MODULES = [
    "autozyme_runtime.py",
    "autozyme_runtime.R",
    "atac_doublet_scrublet.py",
]

REQUIRED_COMMANDS = [
    "python",
    "scenicplus",
    "snakemake",
    "macs2",
    "cbust",
    "bedtools",
    "samtools",
    "gzip",
    "tabix",
    "bgzip",
]

REQUIRED_IMPORTS = [
    "yaml",
    "pandas",
    "scanpy",
    "anndata",
    "mudata",
    "pycisTopic",
    "pycistarget",
    "scenicplus",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--home", default=os.environ.get("SCENICPLUS_HOME"))
    parser.add_argument("--skip-imports", action="store_true")
    parser.add_argument("--skip-commands", action="store_true")
    return parser.parse_args()


def ok(label: str, detail: str = "") -> None:
    print(f"[OK]   {label}{': ' + detail if detail else ''}")


def fail(label: str, detail: str, failures: list[str]) -> None:
    print(f"[MISS] {label}: {detail}")
    failures.append(f"{label}: {detail}")


def run_version(cmd: list[str]) -> str:
    try:
        env = os.environ.copy()
        env.setdefault("PYTHONWARNINGS", "ignore")
        proc = subprocess.run(cmd, check=False, text=True, capture_output=True, timeout=20, env=env)
        text = (proc.stdout or proc.stderr).strip().splitlines()
        return text[0] if text else "found"
    except Exception as exc:
        return f"found, version check failed: {exc}"


def check_mallet_import(mallet_path: str, failures: list[str]) -> None:
    try:
        with tempfile.TemporaryDirectory(prefix="mallet_check_") as tmp:
            tmp_path = Path(tmp)
            text = tmp_path / "toy.txt"
            out = tmp_path / "toy.mallet"
            text.write_text("doc1\talpha beta gamma\ndoc2\tbeta delta\n")
            proc = subprocess.run(
                [
                    mallet_path,
                    "import-file",
                    "--input",
                    str(text),
                    "--output",
                    str(out),
                    "--keep-sequence",
                    "TRUE",
                ],
                text=True,
                capture_output=True,
                timeout=30,
            )
            if proc.returncode != 0 or not out.exists() or out.stat().st_size == 0:
                detail = (proc.stderr or proc.stdout or "no command output").strip().splitlines()
                fail("command mallet import-file", detail[0] if detail else "smoke test failed", failures)
            else:
                ok("command mallet import-file", "toy corpus converted")
    except Exception as exc:
        fail("command mallet import-file", repr(exc), failures)


def main() -> None:
    warnings.filterwarnings("ignore", category=UserWarning)
    args = parse_args()
    failures: list[str] = []

    home = Path(args.home).resolve() if args.home else None
    if home is None or not home.exists():
        fail("SCENICPLUS_HOME", "set SCENICPLUS_HOME=$CONDA_PREFIX/share/scenicplus-grn after activating the environment", failures)
    else:
        ok("SCENICPLUS_HOME", str(home))

    conda_prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("CONDA_ENV_PREFIX")
    if conda_prefix:
        ok("CONDA_PREFIX", conda_prefix)
    else:
        fail("CONDA_PREFIX", "activate the scenicplus-grn conda environment first", failures)

    if home and home.exists():
        for rel in [
            "README.md",
            "README.zh-CN.md",
            "SCENICPLUS_STEP_BY_STEP.md",
            "SCENICPLUS_STEP_BY_STEP.zh-CN.md",
            "VERSION",
            "CHANGELOG.md",
            "RELEASE_NOTES.md",
            "VERSION_LOCK.md",
            "scenicplus_config_template.yaml",
            "locks/environment-linux-64.solved-lock.yml",
        ]:
            path = home / rel
            ok(rel, str(path)) if path.is_file() and path.stat().st_size > 0 else fail(rel, str(path), failures)
        for script in REQUIRED_SCRIPTS:
            path = home / "scripts" / script
            ok(f"scripts/{script}") if path.is_file() and path.stat().st_size > 0 else fail(f"scripts/{script}", str(path), failures)
        for module in REQUIRED_MODULES:
            path = home / "modules" / module
            ok(f"modules/{module}") if path.is_file() and path.stat().st_size > 0 else fail(f"modules/{module}", str(path), failures)

    if not args.skip_commands:
        commands = list(REQUIRED_COMMANDS)
        if os.environ.get("SCENICPLUS_REQUIRE_MALLET", "1").strip().lower() not in {"0", "false", "no", "off"}:
            commands.append("mallet")
        for command in commands:
            path = shutil.which(command)
            if path:
                detail = path
                if command in {"snakemake", "macs2"}:
                    detail += f" ({run_version([command, '--version'])})"
                ok(f"command {command}", detail)
                if command == "mallet":
                    check_mallet_import(path, failures)
            else:
                fail(f"command {command}", "not found in PATH", failures)

    if not args.skip_imports:
        for module in REQUIRED_IMPORTS:
            try:
                imported = importlib.import_module(module)
                ok(f"python import {module}", getattr(imported, "__version__", "version_unknown"))
            except Exception as exc:
                fail(f"python import {module}", repr(exc), failures)

    if failures:
        print("\nINSTALLATION CHECK FAILED")
        print("Fix the missing items above, then rerun:")
        print("  python $SCENICPLUS_HOME/scripts/check_workflow_installation.py")
        sys.exit(1)

    print("\nINSTALLATION CHECK OK")
    print("Next: continue with the project initialization workflow.")


if __name__ == "__main__":
    main()
