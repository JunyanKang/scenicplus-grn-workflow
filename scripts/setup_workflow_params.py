#!/usr/bin/env python
"""Create or update project parameter tables used by the SCENIC+ workflow."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
INPUTS = PROJECT / "inputs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--section",
        action="append",
        default=[],
        help="Section to write. Repeatable. Default: all.",
    )
    parser.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="SECTION.KEY=VALUE",
        help="Override a default value, for example --set pycistopic.n_iter=300.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing values with defaults plus overrides.")
    return parser.parse_args()


def read_project_env() -> dict[str, str]:
    values: dict[str, str] = {}
    for name in ["project_env.sh", "scenicplus_project.env"]:
        path = PROJECT / name
        if not path.exists():
            continue
        for raw in path.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[len("export ") :]
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip().strip("'").strip('"')
            values[key] = value
    values.update(os.environ)
    return values


def write_project_setting(key: str, value: str) -> None:
    config = PROJECT / "scenicplus_project.env"
    lines = config.read_text().splitlines() if config.exists() else []
    assignment = f"{key}={value!r}"
    out: list[str] = []
    replaced = False
    for line in lines:
        if line.startswith(f"{key}="):
            out.append(assignment)
            replaced = True
        else:
            out.append(line)
    if not replaced:
        out.append(assignment)
    config.write_text("\n".join(out) + "\n")


def parse_overrides(items: list[str]) -> dict[str, dict[str, str]]:
    parsed: dict[str, dict[str, str]] = {}
    for item in items:
        if "=" not in item or "." not in item.split("=", 1)[0]:
            raise SystemExit(f"Invalid --set value: {item}. Use SECTION.KEY=VALUE.")
        left, value = item.split("=", 1)
        section, key = left.split(".", 1)
        parsed.setdefault(section.strip(), {})[key.strip()] = value.strip()
    return parsed


def merge_values(path: Path, defaults: dict[str, str], overrides: dict[str, str], force: bool) -> dict[str, str]:
    values = defaults.copy()
    if path.exists() and not force:
        old = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if old.shape[1] >= 2:
            for key, value in zip(old.iloc[:, 0], old.iloc[:, 1]):
                key = str(key).strip()
                if key:
                    values[key] = str(value).strip()
    values.update(overrides)
    return values


def write_param_table(path: Path, defaults: dict[str, str], overrides: dict[str, str], force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values = merge_values(path, defaults, overrides, force)
    ordered = list(defaults)
    ordered.extend(k for k in values if k not in ordered)
    pd.DataFrame({"parameter": ordered, "value": [values[k] for k in ordered]}).to_csv(
        path,
        sep="\t",
        index=False,
    )
    print(f"WROTE {path}")


def write_two_col_table(path: Path, rows: list[dict[str, str]], key_cols: list[str], force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        old = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if not old.empty:
            print(f"KEPT {path}")
            return
    pd.DataFrame(rows, columns=key_cols).to_csv(path, sep="\t", index=False)
    print(f"WROTE {path}")


def write_sample_threshold_table(path: Path, rows: list[dict[str, str]], force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["sample_id", "min_cistopic_fragments", "min_cistopic_accessible_regions"]
    new_samples = {row["sample_id"] for row in rows}
    if path.exists() and not force:
        old = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if not old.empty and set(old.get("sample_id", pd.Series(dtype=str)).astype(str)) == new_samples:
            print(f"KEPT {path}")
            return
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)
    print(f"WROTE {path}")


def n_cpu_default(env: dict[str, str]) -> str:
    for key in ["SCENICPLUS_N_CPU", "N_CPU"]:
        if env.get(key):
            return env[key]
    count = os.cpu_count() or 4
    return str(max(1, min(16, count)))


def sample_ids() -> list[str]:
    path = INPUTS / "sample_sheet.tsv"
    if not path.exists():
        return ["sample_1", "sample_2"]
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "sample_id" not in df.columns:
        return ["sample_1", "sample_2"]
    ids = [x for x in df["sample_id"].astype(str).tolist() if x]
    return ids or ["sample_1", "sample_2"]


def annotated_param_defaults() -> dict[str, str]:
    path = INPUTS / "annotated_object_params.tsv"
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if {"parameter", "value"}.issubset(df.columns):
            return dict(zip(df["parameter"].astype(str), df["value"].astype(str)))
    except Exception:
        return {}
    return {}


def defaults(env: dict[str, str]) -> dict[str, tuple[Path, dict[str, str]]]:
    n_cpu = n_cpu_default(env)
    max_memory = env.get("SCENICPLUS_MAX_MEMORY_GB", "auto")
    chromsizes = env.get("CHROMSIZES") or env.get("CHROMS") or "resources/organism.ucsc.standard.chromsizes.tsv"
    blacklist = env.get("BLACKLIST", "")
    genome_size = env.get("MACS_GENOME_SIZE", "mm")
    annotated_object = env.get("ANNOTATED_OBJECT", "")
    annotated_format = env.get("ANNOTATED_OBJECT_FORMAT", "").lower().lstrip(".")
    annotated_params = annotated_param_defaults()
    metacell_seurat = annotated_object if (annotated_format in {"rds", "qs"} or annotated_object.lower().endswith((".rds", ".qs"))) else ""
    return {
        "metacell": (
            INPUTS / "metacell_params.tsv",
            {
                "seurat_rds": metacell_seurat,
                "wgcna_name": "scenicplus_metacells",
                "assay": annotated_params.get("assay", "RNA") or "RNA",
                "layer": annotated_params.get("layer", "counts") or "counts",
                "reduction": annotated_params.get("reduction", "wnn.umap") or "wnn.umap",
                "dims": "1 2",
                "k": "12",
                "max_shared": "2",
                "min_cells": "1",
                "target_metacells": "1000000",
                "mode": "average",
            },
        ),
        "pycistopic": (
            INPUTS / "pycistopic_params.tsv",
            {
                "analysis_unit": "metacell",
                "chromsizes": chromsizes,
                "blacklist": blacklist,
                "macs_path": "macs2",
                "genome_size": genome_size,
                "n_cpu": "auto",
                "pseudobulk_n_cpu": "auto",
                "macs_n_cpu": "auto",
                "topic_n_cpu": "auto",
                "ray_temp_dir": "tmp/ray",
                "resume_pseudobulk": "1",
                "resume_macs2": "1",
                "peak_half_width": "250",
                "min_frag": "1000",
                "min_cell": "1",
                "is_acc": "1",
                "lda_backend": "cgs",
                "mallet_path": "mallet",
                "mallet_tmp_dir": "tmp/mallet",
                "mallet_memory_gb": "auto",
                "reuse_mallet_corpus": "1",
                "n_iter": "150",
                "random_state": "555",
                "selected_n_topics": "auto",
                "ntop_regions": "3000",
                "dar_adjpval_thr": "0.05",
                "dar_log2fc_thr": "0.5",
                "expected_doublet_rate": "0.1",
                "doublet_n_prin_comps": "30",
                "doublet_min_counts": "2",
                "doublet_min_cells": "3",
                "write_pseudobulk_bigwig": "0",
            },
        ),
        "cistarget": (
            INPUTS / "cistarget_db_params.tsv",
            {
                "n_cpu": "auto",
                "min_free_memory_gb": "8",
                "memory_gb_per_worker": "3",
                "max_workers": n_cpu,
                "max_cpu_load_fraction": "0.80",
                "use_partial": "auto",
                "partial_n_parts": "auto",
                "target_partial_matrix_gb": "3",
                "max_partial_parts": "64",
                "seed": "555",
            },
        ),
        "scenicplus_config": (
            INPUTS / "scenicplus_config_params.tsv",
            {
                "n_cpu": n_cpu,
                "seed": "666",
                "nr_cells_per_metacells": "10",
                "search_space_upstream": "1000 150000",
                "search_space_downstream": "1000 150000",
                "search_space_extend_tss": "10 10",
                "dem_motif_hit_thr": "3.0",
                "ctx_nes_threshold": "3.0",
                "rho_threshold": "0.05",
                "min_target_genes": "10",
            },
        ),
        "preflight": (
            INPUTS / "preflight_thresholds.tsv",
            {
                "min_overlap_rna": "0.50",
                "min_overlap_atac": "0.50",
                "min_cell_metadata_coverage": "0.95",
                "min_motif_tf_gene_overlap": "0.50",
            },
        ),
        "snakemake": (
            INPUTS / "snakemake_params.tsv",
            {
                "cores": n_cpu,
                "rerun_incomplete": "1",
                "printshellcmds": "1",
                "latency_wait": "60",
            },
        ),
        "postprocess": (
            INPUTS / "postprocess_params.tsv",
            {
                "metadata": "inputs/cell_metadata.tsv",
                "group_col": "cell_label",
                "condition_col": "condition",
                "sample_col": "sample_id",
                "cell_col": "cell_id",
                "tf_to_gene": "work/scenicplus/tf_to_gene_adj.tsv",
                "direct_auc_h5mu": "results/scenicplus/AUCell_direct.h5mu",
                "direct_eregulons": "results/scenicplus/eRegulons_direct.tsv",
                "direct_stats_outdir": "results/scenicplus_stats/auc_by_condition_direct",
                "direct_figures_outdir": "results/scenicplus_figures/direct",
                "extended_auc_h5mu": "results/scenicplus/AUCell_extended.h5mu",
                "extended_eregulons": "results/scenicplus/eRegulons_extended.tsv",
                "extended_stats_outdir": "results/scenicplus_stats/auc_by_condition_extended",
                "extended_figures_outdir": "results/scenicplus_figures/extended",
                "plot_top_n": "30",
                "plot_umap_n": "12",
                "network_top_tfs": "8",
                "network_targets_per_tf": "12",
            },
        ),
    }


def main() -> None:
    args = parse_args()
    env = read_project_env()
    all_defaults = defaults(env)
    sections = args.section or list(all_defaults)
    overrides = parse_overrides(args.overrides)
    unknown = sorted(set(sections) - set(all_defaults))
    if unknown:
        raise SystemExit(f"Unknown section(s): {', '.join(unknown)}")

    for section in sections:
        path, section_defaults = all_defaults[section]
        write_param_table(path, section_defaults, overrides.get(section, {}), args.force)
        write_project_setting(f"{section.upper()}_PARAMS", str(path))

    if "pycistopic" in sections:
        rows = [
            {
                "sample_id": sid,
                "min_cistopic_fragments": "1000",
                "min_cistopic_accessible_regions": "1000",
            }
            for sid in sample_ids()
        ]
        write_sample_threshold_table(
            INPUTS / "atac_qc_thresholds.tsv",
            rows,
            args.force,
        )
        topic_overrides = overrides.get("topic_grid", {})
        n_topics = topic_overrides.get("n_topics", "10,20,30,40,50")
        topic_rows = [{"n_topics": x.strip()} for x in n_topics.replace(";", ",").split(",") if x.strip()]
        write_two_col_table(
            INPUTS / "topic_model_grid.tsv",
            topic_rows,
            ["n_topics"],
            args.force or bool(topic_overrides),
        )
        write_project_setting("ATAC_QC_THRESHOLDS", str(INPUTS / "atac_qc_thresholds.tsv"))
        write_project_setting("TOPIC_MODEL_GRID", str(INPUTS / "topic_model_grid.tsv"))


if __name__ == "__main__":
    main()
