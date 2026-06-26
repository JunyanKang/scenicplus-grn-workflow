#!/usr/bin/env python
"""Initialize a SCENIC+ project directory and write reusable shell variables."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

MODULES_DIR = Path(__file__).resolve().parents[1] / "modules"
sys.path.insert(0, str(MODULES_DIR))

from organism_resources import resolve_organism  # noqa: E402


MOTIF2TF_REFERENCES = ["auto", "human", "mouse", "fly", "chicken"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None, help="Project settings file, for example scenicplus_project.env")
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--organism", default=None, help="Ensembl species name or common alias; validated in resource preparation.")
    parser.add_argument("--autozyme", choices=["on", "off"], default=None)
    parser.add_argument("--conda-root", default=None)
    parser.add_argument("--env-name", default=None)
    parser.add_argument("--ensembl-release", default=None)
    parser.add_argument("--max-memory-gb", default=None, help="Workflow RAM budget in GB, or auto.")
    parser.add_argument("--cell-label-column", default=None)
    parser.add_argument("--motif2tf-reference", choices=MOTIF2TF_REFERENCES, default=None)
    parser.add_argument("--motif2tf-table", default=None)
    parser.add_argument("--env-file", default=None, help="Default: <project-dir>/project_env.sh")
    return parser.parse_args()


def read_config(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    config_path = Path(path).expanduser()
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    values: dict[str, str] = {}
    for raw in config_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.strip().strip('"').strip("'")
        values[key.strip()] = value
    return values


def config_get(config: dict[str, str], key: str, default: str | None = None) -> str | None:
    return config.get(key, default)


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def update_config(path: Path, values: dict[str, str]) -> None:
    existing = path.read_text().splitlines() if path.exists() else []
    seen = set()
    out = []
    for line in existing:
        if "=" not in line or line.lstrip().startswith("#"):
            out.append(line)
            continue
        key = line.split("=", 1)[0].strip()
        if key in values:
            out.append(f"{key}={shell_quote(values[key])}")
            seen.add(key)
        else:
            out.append(line)
    for key, value in values.items():
        if key not in seen:
            out.append(f"{key}={shell_quote(value)}")
    path.write_text("\n".join(out) + "\n")


def read_macs_genome_size(chromsizes: Path) -> str:
    if not chromsizes.exists():
        return ""
    total = 0
    with chromsizes.open() as handle:
        header = next(handle, "")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                total += int(fields[2])
    return str(total) if total > 0 else ""


def normalize_memory_budget(value: str | None) -> str:
    value = (value or "auto").strip()
    if value.lower() == "auto":
        return "auto"
    try:
        number = float(value)
    except ValueError as exc:
        raise SystemExit("ERROR: SCENICPLUS_MAX_MEMORY_GB must be auto or a positive number.") from exc
    if number <= 0:
        raise SystemExit("ERROR: SCENICPLUS_MAX_MEMORY_GB must be auto or a positive number.")
    return str(int(number)) if number.is_integer() else str(number)


def main() -> None:
    args = parse_args()
    config = read_config(args.config)
    project_dir_value = args.project_dir or config_get(config, "PROJECT_DIR")
    organism = args.organism or config_get(config, "ORGANISM")
    ensembl_release = args.ensembl_release or config_get(config, "ENSEMBL_RELEASE", "")
    if organism:
        organism = resolve_organism(organism, ensembl_release or None)
    if not project_dir_value:
        raise SystemExit("ERROR: --project-dir or PROJECT_DIR in --config is required.")
    autozyme = args.autozyme or config_get(config, "AUTOZYME", "on")
    if autozyme not in {"on", "off"}:
        raise SystemExit("ERROR: AUTOZYME must be on or off.")
    conda_root = args.conda_root or config_get(config, "CONDA_ROOT")
    env_name = args.env_name or config_get(config, "ENV_NAME", "scenicplus-grn")
    max_memory_gb = normalize_memory_budget(args.max_memory_gb or config_get(config, "SCENICPLUS_MAX_MEMORY_GB", "auto"))
    cell_label_column = args.cell_label_column or config_get(config, "CELL_LABEL_COLUMN", "")
    atac_input_layout = config_get(config, "ATAC_INPUT_LAYOUT", "")
    atac_data_root = config_get(config, "ATAC_DATA_ROOT", "")
    motif2tf_reference = args.motif2tf_reference or config_get(config, "MOTIF2TF_REFERENCE", "")
    motif2tf_table = args.motif2tf_table or config_get(config, "MOTIF2TF_TABLE", "")
    if motif2tf_reference and motif2tf_reference not in MOTIF2TF_REFERENCES:
        raise SystemExit(f"ERROR: MOTIF2TF_REFERENCE must be one of: {', '.join(MOTIF2TF_REFERENCES)}")

    project_dir = Path(project_dir_value).expanduser().resolve()
    for rel in [
        "inputs",
        "inputs/fragments",
        "inputs/region_sets",
        "inputs/cistarget_db",
        "resources",
        "work",
        "results",
        "logs",
        "tmp",
        "tmp/python_cache",
        "tmp/matplotlib",
        "tmp/numba",
    ]:
        (project_dir / rel).mkdir(parents=True, exist_ok=True)

    chroms = f"resources/chromosomes/{organism}.ucsc.standard.chroms.txt" if organism else ""
    genome = f"resources/{organism}/{organism}.ucsc.standard.fa" if organism else ""
    chromsizes = f"resources/{organism}/{organism}.ucsc.standard.chromsizes.tsv" if organism else ""
    genome_annotation = f"resources/{organism}/{organism}.ucsc.standard.genome_annotation.tsv" if organism else ""
    macs_genome_size = read_macs_genome_size(project_dir / chromsizes) if chromsizes else ""
    autozyme_disabled = "0" if autozyme == "on" else "1"
    env_file = Path(args.env_file).expanduser().resolve() if args.env_file else project_dir / "project_env.sh"
    env_file.parent.mkdir(parents=True, exist_ok=True)
    scenicplus_home = ""
    if conda_root:
        scenicplus_home = str(Path(conda_root).expanduser() / "envs" / str(env_name) / "share" / "scenicplus-grn")
    project_config = project_dir / "scenicplus_project.env"
    if args.config:
        source_config = Path(args.config).expanduser().resolve()
        if source_config != project_config.resolve():
            project_config.write_text(source_config.read_text())
        else:
            project_config.touch(exist_ok=True)
    update_values = {}
    if organism:
        update_values["ORGANISM"] = organism
    if ensembl_release:
        update_values["ENSEMBL_RELEASE"] = str(ensembl_release)
    if update_values:
        update_config(project_config, update_values)

    lines = [
        "# Source this file before running project commands:",
        f"#   source {env_file}",
        f"export PROJECT_DIR={shell_quote(str(project_dir))}",
        f"export SCENICPLUS_PROJECT_CONFIG={shell_quote(str(project_config))}",
        f"export ENV_NAME={shell_quote(str(env_name))}",
        f"export SCENICPLUS_MAX_MEMORY_GB={shell_quote(max_memory_gb)}",
        f"export AUTOZYME={shell_quote(autozyme)}",
        'export MPLCONFIGDIR="$PROJECT_DIR/tmp/matplotlib"',
        'export NUMBA_CACHE_DIR="$PROJECT_DIR/tmp/numba"',
        f"export CHROMS={shell_quote(chroms)}",
        f"export GENOME={shell_quote(genome)}",
        f"export CHROMSIZES={shell_quote(chromsizes)}",
        f"export GENOME_ANNOTATION={shell_quote(genome_annotation)}",
        f"export AUTOZYME_DISABLED={shell_quote(autozyme_disabled)}",
    ]
    if organism:
        lines.insert(3, f"export ORGANISM={shell_quote(organism)}")
    if ensembl_release:
        lines.insert(4 if organism else 3, f"export ENSEMBL_RELEASE={shell_quote(str(ensembl_release))}")
    if conda_root:
        env_prefix = str(Path(conda_root).expanduser() / "envs" / str(env_name))
        lines.insert(3, f"export CONDA_ROOT={shell_quote(str(Path(conda_root).expanduser()))}")
        lines.insert(4, f"export CONDA_ENV_PREFIX={shell_quote(env_prefix)}")
        lines.insert(5, 'export CONDA_PREFIX="$CONDA_ENV_PREFIX"')
        lines.insert(6, 'export PATH="$CONDA_ENV_PREFIX/bin:$PATH"')
    if scenicplus_home:
        lines.insert(6 if conda_root else 3, f"export SCENICPLUS_HOME={shell_quote(scenicplus_home)}")
    if macs_genome_size:
        lines.append(f"export MACS_GENOME_SIZE={shell_quote(macs_genome_size)}")
    else:
        lines.append('export MACS_GENOME_SIZE=""')
    if cell_label_column:
        lines.append(f"export CELL_LABEL_COLUMN={shell_quote(cell_label_column)}")
    if motif2tf_table:
        lines.append(f"export MOTIF2TF_TABLE={shell_quote(motif2tf_table)}")
    if motif2tf_reference:
        lines.append(f"export MOTIF2TF_REFERENCE={shell_quote(motif2tf_reference)}")
    if atac_input_layout:
        lines.append(f"export ATAC_INPUT_LAYOUT={shell_quote(atac_input_layout)}")
    if atac_data_root:
        lines.append(f"export ATAC_DATA_ROOT={shell_quote(atac_data_root)}")
    lines.extend(
        [
            'if [ -f "$SCENICPLUS_PROJECT_CONFIG" ]; then',
            "  set -a",
            '  . "$SCENICPLUS_PROJECT_CONFIG"',
            "  set +a",
            "fi",
        ]
    )
    lines.append('cd "$PROJECT_DIR"')
    env_file.write_text("\n".join(lines) + "\n")

    print(f"WROTE {env_file}")
    print(f"PROJECT_DIR={project_dir}")
    print(f"ORGANISM={organism or 'not set; set in Step 1'}")
    if conda_root:
        print(f"CONDA_ROOT={Path(conda_root).expanduser()}")
    print(f"ENV_NAME={env_name}")
    print(f"ENSEMBL_RELEASE={ensembl_release or 'not set; set in Step 1'}")
    print(f"SCENICPLUS_MAX_MEMORY_GB={max_memory_gb}")
    print(f"MOTIF2TF_REFERENCE={motif2tf_reference or 'not set; prepare_official_resources defaults to auto'}")
    if motif2tf_table:
        print(f"MOTIF2TF_TABLE={motif2tf_table}")
    if cell_label_column:
        print(f"CELL_LABEL_COLUMN={cell_label_column}")
    else:
        print("CELL_LABEL_COLUMN will be set after annotated-object export.")
    print(f"AUTOZYME={'enabled' if autozyme == 'on' else 'disabled'}")
    if macs_genome_size:
        print(f"MACS_GENOME_SIZE={macs_genome_size}")
    else:
        print("MACS_GENOME_SIZE not set yet; prepare_official_resources.py will fill it after resource preparation.")


if __name__ == "__main__":
    main()
