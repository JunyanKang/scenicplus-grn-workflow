#!/usr/bin/env python
"""Generate SCENIC+ Snakemake config.yaml from a project path and template."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
SHARE_DIR = SCRIPT_DIR.parent
DEFAULT_TEMPLATE = SHARE_DIR / "config" / "scenicplus_config_template.yaml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", default=None, help="Project root directory. Default: $PROJECT_DIR")
    parser.add_argument("--template", default=str(DEFAULT_TEMPLATE))
    parser.add_argument("--organism-config", default=None, help="Default: <project-dir>/work/scenicplus/organism_config.yaml")
    parser.add_argument("--out", default=None, help="Default: <project-dir>/work/scenicplus/Snakemake/config/config.yaml")
    parser.add_argument("--params", default=None, help="Default: <project-dir>/inputs/scenicplus_config_params.tsv")
    parser.add_argument("--n-cpu", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--nr-cells-per-metacells", type=int, default=None)
    parser.add_argument("--search-space-upstream", default=None)
    parser.add_argument("--search-space-downstream", default=None)
    parser.add_argument("--search-space-extend-tss", default=None)
    parser.add_argument("--fraction-overlap-w-ctx-database", type=float, default=None)
    parser.add_argument("--fraction-overlap-w-dem-database", type=float, default=None)
    parser.add_argument("--dem-balance-number-of-promoters", default=None)
    parser.add_argument("--dem-promoter-space", type=int, default=None)
    parser.add_argument("--ctx-auc-threshold", type=float, default=None)
    parser.add_argument("--dem-motif-hit-thr", type=float, default=None)
    parser.add_argument("--dem-adj-pval-thr", type=float, default=None)
    parser.add_argument("--dem-log2fc-thr", type=float, default=None)
    parser.add_argument("--dem-mean-fg-thr", type=float, default=None)
    parser.add_argument("--dem-max-bg-regions", type=int, default=None)
    parser.add_argument("--dem-n-cpu", type=int, default=None)
    parser.add_argument("--ctx-nes-threshold", type=float, default=None)
    parser.add_argument("--ctx-rank-threshold", type=float, default=None)
    parser.add_argument("--ctx-n-cpu", type=int, default=None)
    parser.add_argument("--motif-similarity-fdr", type=float, default=None)
    parser.add_argument("--orthologous-identity-threshold", type=float, default=None)
    parser.add_argument("--gsea-n-perm", type=int, default=None)
    parser.add_argument("--quantile-thresholds-region-to-gene", default=None)
    parser.add_argument("--top-n-regiontogenes-per-gene", default=None)
    parser.add_argument("--top-n-regiontogenes-per-region", default=None)
    parser.add_argument("--min-regions-per-gene", type=int, default=None)
    parser.add_argument("--rho-threshold", type=float, default=None)
    parser.add_argument("--min-target-genes", type=int, default=None)
    return parser.parse_args()


def replace_project_dir(value: Any, project_dir: str) -> Any:
    if isinstance(value, str):
        return value.replace("{PROJECT_DIR}", project_dir)
    if isinstance(value, list):
        return [replace_project_dir(x, project_dir) for x in value]
    if isinstance(value, dict):
        return {k: replace_project_dir(v, project_dir) for k, v in value.items()}
    return value


def deep_update(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(
            f"SCENIC+ config parameter table not found: {path}\n"
            "Create inputs/scenicplus_config_params.tsv before running this script."
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.shape[1] < 2:
        raise ValueError(f"{path} must contain two tab-separated columns: parameter and value")
    key_col, value_col = df.columns[:2]
    params = {str(k).strip(): str(v).strip() for k, v in zip(df[key_col], df[value_col]) if str(k).strip()}
    required = [
        "n_cpu",
        "seed",
        "nr_cells_per_metacells",
        "search_space_upstream",
        "search_space_downstream",
        "search_space_extend_tss",
        "dem_motif_hit_thr",
        "dem_n_cpu",
        "ctx_nes_threshold",
        "ctx_n_cpu",
        "rho_threshold",
        "min_target_genes",
    ]
    missing = [x for x in required if x not in params or params[x] == ""]
    if missing:
        raise ValueError(f"{path} missing required parameters: {', '.join(missing)}")
    return params


def get_value(cli_value: Any, params: dict[str, str], key: str, cast):
    if cli_value is not None:
        return cli_value
    return cast(params[key])


def str_to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value}")


def set_optional_param(
    config: dict[str, Any],
    section: str,
    key: str,
    cli_value: Any,
    params: dict[str, str],
    cast,
) -> None:
    if cli_value is not None:
        config[section][key] = cli_value
        return
    if key in params and params[key] != "":
        config[section][key] = cast(params[key])


def main() -> None:
    args = parse_args()
    project_dir = Path(args.project_dir or os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
    template = Path(args.template)
    organism_config_path = Path(args.organism_config) if args.organism_config else project_dir / "work" / "scenicplus" / "organism_config.yaml"
    out = Path(args.out) if args.out else project_dir / "work" / "scenicplus" / "Snakemake" / "config" / "config.yaml"
    params_path = Path(args.params).expanduser() if args.params else project_dir / "inputs" / "scenicplus_config_params.tsv"
    if not params_path.is_absolute():
        params_path = (project_dir / params_path).resolve()
    params = read_params(params_path)

    if not template.exists():
        raise FileNotFoundError(template)
    if not organism_config_path.exists():
        raise FileNotFoundError(organism_config_path)

    config = yaml.safe_load(template.read_text())
    config = replace_project_dir(config, str(project_dir))
    organism_config = yaml.safe_load(organism_config_path.read_text())
    if not isinstance(organism_config, dict):
        raise ValueError(f"Invalid organism config: {organism_config_path}")
    config = deep_update(config, organism_config)

    config["params_general"]["n_cpu"] = get_value(args.n_cpu, params, "n_cpu", int)
    config["params_general"]["seed"] = get_value(args.seed, params, "seed", int)
    config["params_data_preparation"]["nr_cells_per_metacells"] = get_value(args.nr_cells_per_metacells, params, "nr_cells_per_metacells", int)
    config["params_data_preparation"]["search_space_upstream"] = get_value(args.search_space_upstream, params, "search_space_upstream", str)
    config["params_data_preparation"]["search_space_downstream"] = get_value(args.search_space_downstream, params, "search_space_downstream", str)
    config["params_data_preparation"]["search_space_extend_tss"] = get_value(args.search_space_extend_tss, params, "search_space_extend_tss", str)
    set_optional_param(config, "params_motif_enrichment", "fraction_overlap_w_ctx_database", args.fraction_overlap_w_ctx_database, params, float)
    set_optional_param(config, "params_motif_enrichment", "fraction_overlap_w_dem_database", args.fraction_overlap_w_dem_database, params, float)
    set_optional_param(config, "params_motif_enrichment", "dem_balance_number_of_promoters", args.dem_balance_number_of_promoters, params, str_to_bool)
    set_optional_param(config, "params_motif_enrichment", "dem_promoter_space", args.dem_promoter_space, params, int)
    set_optional_param(config, "params_motif_enrichment", "ctx_auc_threshold", args.ctx_auc_threshold, params, float)
    set_optional_param(config, "params_motif_enrichment", "ctx_rank_threshold", args.ctx_rank_threshold, params, float)
    set_optional_param(config, "params_motif_enrichment", "dem_max_bg_regions", args.dem_max_bg_regions, params, int)
    set_optional_param(config, "params_motif_enrichment", "dem_adj_pval_thr", args.dem_adj_pval_thr, params, float)
    set_optional_param(config, "params_motif_enrichment", "dem_log2fc_thr", args.dem_log2fc_thr, params, float)
    set_optional_param(config, "params_motif_enrichment", "dem_mean_fg_thr", args.dem_mean_fg_thr, params, float)
    config["params_motif_enrichment"]["dem_motif_hit_thr"] = get_value(args.dem_motif_hit_thr, params, "dem_motif_hit_thr", float)
    config["params_motif_enrichment"]["dem_n_cpu"] = get_value(args.dem_n_cpu, params, "dem_n_cpu", int)
    config["params_motif_enrichment"]["ctx_nes_threshold"] = get_value(args.ctx_nes_threshold, params, "ctx_nes_threshold", float)
    config["params_motif_enrichment"]["ctx_n_cpu"] = get_value(args.ctx_n_cpu, params, "ctx_n_cpu", int)
    set_optional_param(config, "params_motif_enrichment", "motif_similarity_fdr", args.motif_similarity_fdr, params, float)
    set_optional_param(config, "params_motif_enrichment", "orthologous_identity_threshold", args.orthologous_identity_threshold, params, float)
    set_optional_param(config, "params_inference", "gsea_n_perm", args.gsea_n_perm, params, int)
    set_optional_param(config, "params_inference", "quantile_thresholds_region_to_gene", args.quantile_thresholds_region_to_gene, params, str)
    set_optional_param(config, "params_inference", "top_n_regionTogenes_per_gene", args.top_n_regiontogenes_per_gene, params, str)
    set_optional_param(config, "params_inference", "top_n_regionTogenes_per_region", args.top_n_regiontogenes_per_region, params, str)
    set_optional_param(config, "params_inference", "min_regions_per_gene", args.min_regions_per_gene, params, int)
    config["params_inference"]["rho_threshold"] = get_value(args.rho_threshold, params, "rho_threshold", float)
    config["params_inference"]["min_target_genes"] = get_value(args.min_target_genes, params, "min_target_genes", int)

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(config, sort_keys=False))
    print(f"WROTE {out}")
    print(f"PROJECT_DIR={project_dir}")
    print(f"TEMPLATE={template}")
    print(f"ORGANISM_CONFIG={organism_config_path}")
    print(f"PARAMS={params_path}")


if __name__ == "__main__":
    main()
