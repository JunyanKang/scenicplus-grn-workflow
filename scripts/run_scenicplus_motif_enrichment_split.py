#!/usr/bin/env python
"""Run SCENIC+ motif enrichment in region-set chunks.

This keeps each cisTarget/DEM process bounded in memory while preserving the
official SCENIC+ prepare_menr interface, which accepts multiple motif
enrichment HDF5 files.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import psutil
import yaml


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["ctx", "dem", "both", "prepare-menr", "status"], default="both")
    parser.add_argument("--config", default=None, help="Default: work/scenicplus/Snakemake/config/config.yaml")
    parser.add_argument("--split-root", default="inputs/region_sets_split")
    parser.add_argument("--out-dir", default="work/scenicplus/motif_enrichment_split")
    parser.add_argument("--log-dir", default="logs/motif_enrichment_split")
    parser.add_argument("--max-parallel-chunks", default="auto", help="auto or a positive integer.")
    parser.add_argument("--force", action="store_true", help="Rerun chunks even when output HDF5 exists.")
    parser.add_argument("--skip-prepare-menr", action="store_true")
    parser.add_argument(
        "--dem-relaxed-adjpval-thr",
        type=float,
        default=0.05,
        help="Relaxed DEM adj-p-value threshold when a formal DEM chunk is empty.",
    )
    parser.add_argument(
        "--dem-relaxed-log2fc-thr",
        type=float,
        default=0.0,
        help="Relaxed DEM log2FC threshold when a formal DEM chunk is empty.",
    )
    parser.add_argument(
        "--dem-relaxed-mean-fg-thr",
        type=float,
        default=0.0,
        help="Relaxed DEM mean foreground threshold when a formal DEM chunk is empty.",
    )
    parser.add_argument(
        "--allow-empty-dem",
        action="store_true",
        help="Deprecated: empty DEM chunks are now auto-diagnosed and continue by default.",
    )
    return parser.parse_args()


def p(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else PROJECT / path


def read_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"Invalid SCENIC+ config: {path}")
    return data


def file_size_gb(path: Path) -> float:
    try:
        return path.stat().st_size / 1_000_000_000
    except OSError:
        return 0.0


def resolve_positive_int(value: Any, default: int = 1) -> int:
    try:
        return max(1, int(float(str(value))))
    except Exception:
        return default


def as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(x) for x in value if str(x).strip()]
    text = str(value).strip()
    if not text:
        return []
    return [x for x in text.replace(",", " ").split() if x]


def resolve_parallel_chunks(raw: str, config: dict[str, Any]) -> tuple[int, dict[str, str]]:
    if raw.strip().lower() not in {"", "auto"}:
        return resolve_positive_int(raw), {"mode": "user", "max_parallel_chunks": str(resolve_positive_int(raw))}

    input_data = config["input_data"]
    ctx_db = p(input_data["ctx_db_fname"])
    dem_db = p(input_data["dem_db_fname"])
    vm = psutil.virtual_memory()
    available_gb = vm.available / 1_000_000_000
    total_gb = vm.total / 1_000_000_000

    # pycistarget can retain large ranking/scores tables across multiple region
    # set families in one Python process. Use a conservative per-process memory
    # estimate and reserve headroom for the OS and interactive analysis.
    ctx_worker_gb = max(22.0, file_size_gb(ctx_db) * 1.7 + 5.0)
    dem_worker_gb = max(14.0, file_size_gb(dem_db) * 1.5 + 4.0)
    worker_gb = max(ctx_worker_gb, dem_worker_gb)
    usable_gb = max(1.0, min(available_gb * 0.70, total_gb * 0.55))
    cpu_cap = max(1, min(os.cpu_count() or 1, resolve_positive_int(config["params_general"].get("n_cpu", 1))))
    resolved = max(1, min(cpu_cap, int(usable_gb // worker_gb) or 1))
    report = {
        "mode": "auto_memory_aware",
        "max_parallel_chunks": str(resolved),
        "cpu_cap": str(cpu_cap),
        "total_memory_gb": f"{total_gb:.2f}",
        "available_memory_gb": f"{available_gb:.2f}",
        "usable_memory_gb": f"{usable_gb:.2f}",
        "ctx_db_gb": f"{file_size_gb(ctx_db):.2f}",
        "dem_db_gb": f"{file_size_gb(dem_db):.2f}",
        "estimated_worker_memory_gb": f"{worker_gb:.2f}",
    }
    return resolved, report


def discover_or_create_splits(config: dict[str, Any], split_root: Path) -> list[tuple[str, Path]]:
    region_root = p(config["input_data"]["region_set_folder"])
    if not region_root.exists():
        raise FileNotFoundError(region_root)
    split_root.mkdir(parents=True, exist_ok=True)
    chunks: list[tuple[str, Path]] = []
    children = sorted([x for x in region_root.iterdir() if x.is_dir() and list(x.glob("*.bed"))])
    if not children:
        raise FileNotFoundError(f"No region-set family directories with BED files found in {region_root}")
    for idx, child in enumerate(children, start=1):
        short_name = "DARs" if child.name.startswith("DARs") else child.name
        label = f"{idx:02d}_{short_name}"
        chunk_dir = split_root / label / child.name
        chunk_dir.mkdir(parents=True, exist_ok=True)
        for bed in sorted(child.glob("*.bed")):
            dest = chunk_dir / bed.name
            if not dest.exists() or dest.stat().st_size != bed.stat().st_size:
                shutil.copy2(bed, dest)
        chunks.append((label, split_root / label))
    return chunks


def command_common(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "temp_dir": str(p(config["params_general"]["temp_dir"])),
        "species": config["params_motif_enrichment"]["species"],
        "motif_annotations": str(p(config["input_data"]["path_to_motif_annotations"])),
        "annotation_version": config["params_motif_enrichment"]["annotation_version"],
        "motif_similarity_fdr": str(config["params_motif_enrichment"]["motif_similarity_fdr"]),
        "orthologous_identity_threshold": str(config["params_motif_enrichment"]["orthologous_identity_threshold"]),
        "annotations_to_use": as_list(config["params_motif_enrichment"]["annotations_to_use"]),
        "seed": str(config["params_general"]["seed"]),
    }


def ctx_command(config: dict[str, Any], folder: Path, out: Path, html: Path) -> list[str]:
    common = command_common(config)
    params = config["params_motif_enrichment"]
    n_cpu = resolve_positive_int(params.get("ctx_n_cpu", config["params_general"].get("n_cpu", 1)))
    return [
        "scenicplus",
        "grn_inference",
        "motif_enrichment_cistarget",
        "--region_set_folder",
        str(folder),
        "--cistarget_db_fname",
        str(p(config["input_data"]["ctx_db_fname"])),
        "--output_fname_cistarget_result",
        str(out),
        "--temp_dir",
        common["temp_dir"],
        "--species",
        common["species"],
        "--fr_overlap_w_ctx_db",
        str(params["fraction_overlap_w_ctx_database"]),
        "--auc_threshold",
        str(params["ctx_auc_threshold"]),
        "--nes_threshold",
        str(params["ctx_nes_threshold"]),
        "--rank_threshold",
        str(params["ctx_rank_threshold"]),
        "--path_to_motif_annotations",
        common["motif_annotations"],
        "--annotation_version",
        common["annotation_version"],
        "--motif_similarity_fdr",
        common["motif_similarity_fdr"],
        "--orthologous_identity_threshold",
        common["orthologous_identity_threshold"],
        "--annotations_to_use",
        *common["annotations_to_use"],
        "--write_html",
        "--output_fname_cistarget_html",
        str(html),
        "--n_cpu",
        str(n_cpu),
    ]


def method_n_cpu(method: str, config: dict[str, Any]) -> int:
    params = config["params_motif_enrichment"]
    key = "ctx_n_cpu" if method == "ctx" else "dem_n_cpu"
    return resolve_positive_int(params.get(key, config["params_general"].get("n_cpu", 1)))


def dem_command(
    config: dict[str, Any],
    folder: Path,
    out: Path,
    html: Path,
    dem_adj_pval_thr: float | None = None,
    dem_log2fc_thr: float | None = None,
    dem_mean_fg_thr: float | None = None,
) -> list[str]:
    common = command_common(config)
    params = config["params_motif_enrichment"]
    n_cpu = resolve_positive_int(params.get("dem_n_cpu", config["params_general"].get("n_cpu", 1)))
    adj_pval_thr = params["dem_adj_pval_thr"] if dem_adj_pval_thr is None else dem_adj_pval_thr
    log2fc_thr = params["dem_log2fc_thr"] if dem_log2fc_thr is None else dem_log2fc_thr
    mean_fg_thr = params["dem_mean_fg_thr"] if dem_mean_fg_thr is None else dem_mean_fg_thr
    cmd = [
        "scenicplus",
        "grn_inference",
        "motif_enrichment_dem",
        "--region_set_folder",
        str(folder),
        "--dem_db_fname",
        str(p(config["input_data"]["dem_db_fname"])),
        "--output_fname_dem_result",
        str(out),
        "--temp_dir",
        common["temp_dir"],
        "--species",
        common["species"],
        "--fraction_overlap_w_dem_database",
        str(params["fraction_overlap_w_dem_database"]),
        "--max_bg_regions",
        str(params["dem_max_bg_regions"]),
        "--adjpval_thr",
        str(adj_pval_thr),
        "--log2fc_thr",
        str(log2fc_thr),
        "--mean_fg_thr",
        str(mean_fg_thr),
        "--motif_hit_thr",
        str(params["dem_motif_hit_thr"]),
        "--path_to_motif_annotations",
        common["motif_annotations"],
        "--annotation_version",
        common["annotation_version"],
        "--motif_similarity_fdr",
        common["motif_similarity_fdr"],
        "--orthologous_identity_threshold",
        common["orthologous_identity_threshold"],
        "--annotations_to_use",
        *common["annotations_to_use"],
        "--write_html",
        "--output_fname_dem_html",
        str(html),
        "--seed",
        common["seed"],
        "--n_cpu",
        str(n_cpu),
    ]
    if bool(params.get("dem_balance_number_of_promoters", False)):
        cmd.append("--balance_number_of_promoters")
        cmd.extend(["--genome_annotation", str(p(config["output_data"]["genome_annotation"]))])
        cmd.extend(["--promoter_space", str(params["dem_promoter_space"])])
    return cmd


def _pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    normalized = {str(c).lower().replace(" ", "").replace("_", "").replace("-", ""): c for c in df.columns}
    for candidate in candidates:
        key = candidate.lower().replace(" ", "").replace("_", "").replace("-", "")
        if key in normalized:
            return normalized[key]
    return None


def _summarize_dem_output(hdf5_path: Path, official_adj_thr: float, official_log2fc_thr: float) -> tuple[int, int]:
    total_rows = 0
    pass_rows = 0
    with pd.HDFStore(hdf5_path, mode="r") as store:
        for key in store.keys():
            if not key.endswith("/motif_enrichment"):
                continue
            df = store[key]
            total_rows += len(df)
            adj_col = _pick_col(df, ["adjusted_pval", "adjpval", "q", "qval", "padj", "fdr", "adjustedpval"])
            log2fc_col = _pick_col(df, ["log2fc", "log2_fold_change", "log2fc"])
            if adj_col and log2fc_col:
                adj = pd.to_numeric(df[adj_col], errors="coerce")
                fc = pd.to_numeric(df[log2fc_col], errors="coerce")
                pass_rows += int(((adj <= official_adj_thr) & (fc >= official_log2fc_thr)).sum())
    return total_rows, pass_rows


def run_dem_diagnostic(
    label: str,
    folder: Path,
    config: dict[str, Any],
    out_dir: Path,
    log_dir: Path,
    relaxed_adjpval_thr: float,
    relaxed_log2fc_thr: float,
    relaxed_mean_fg_thr: float,
) -> None:
    diag_out = out_dir / f"{label}_relaxed_dem.hdf5"
    diag_html = out_dir / f"{label}_relaxed_dem.html"
    diag_log = log_dir / f"{label}_relaxed_dem_{datetime.now():%Y%m%d_%H%M%S}.log"
    diag_report_dir = PROJECT / "results" / "scenicplus_diagnostics"
    diag_report_dir.mkdir(parents=True, exist_ok=True)
    diag_summary_tsv = diag_report_dir / f"dem_{label}_relaxed_threshold_diagnostic.tsv"
    diag_summary_md = diag_report_dir / f"dem_{label}_relaxed_threshold_diagnostic.md"

    cmd = dem_command(
        config,
        folder,
        diag_out,
        diag_html,
        dem_adj_pval_thr=relaxed_adjpval_thr,
        dem_log2fc_thr=relaxed_log2fc_thr,
        dem_mean_fg_thr=relaxed_mean_fg_thr,
    )
    env = os.environ.copy()
    thread_count = str(method_n_cpu("dem", config))
    for key in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        env[key] = thread_count
    with diag_log.open("w") as handle:
        handle.write("COMMAND\t" + " ".join(cmd) + "\n")
        proc = subprocess.run(cmd, cwd=PROJECT, stdout=handle, stderr=subprocess.STDOUT, text=True, env=env)

    if proc.returncode != 0:
        pd.DataFrame(
            [
                {
                    "status": "diagnostic_failed",
                    "diagnostic_log": str(diag_log),
                    "diagnostic_hdf5": str(diag_out),
                    "message": f"diagnostic command returned code {proc.returncode}",
                }
            ]
        ).to_csv(diag_summary_tsv, sep="\t", index=False)
        diag_summary_md.write_text(
            "\n".join(
                [
                    "# DEM relaxed-threshold diagnostic",
                    "",
                    f"- diagnostic command failed with code {proc.returncode}",
                    f"- log: {diag_log}",
                    f"- relaxed thresholds: adjpval={relaxed_adjpval_thr}, log2fc={relaxed_log2fc_thr}, mean_fg={relaxed_mean_fg_thr}",
                ]
            )
            + "\n"
        )
        return

    if not diag_out.exists() or diag_out.stat().st_size == 0:
        pd.DataFrame(
            [
                {
                    "status": "diagnostic_empty",
                    "diagnostic_log": str(diag_log),
                    "diagnostic_hdf5": str(diag_out),
                    "message": "no non-empty relaxed-threshold HDF5",
                }
            ]
        ).to_csv(diag_summary_tsv, sep="\t", index=False)
        diag_summary_md.write_text(
            "\n".join(
                [
                    "# DEM relaxed-threshold diagnostic",
                    "",
                    "- no non-empty relaxed-threshold HDF5 was produced",
                    f"- log: {diag_log}",
                    f"- relaxed thresholds: adjpval={relaxed_adjpval_thr}, log2fc={relaxed_log2fc_thr}, mean_fg={relaxed_mean_fg_thr}",
                ]
            )
            + "\n"
        )
        return

    params = config["params_motif_enrichment"]
    official_adj_thr = float(params["dem_adj_pval_thr"])
    official_log2fc_thr = float(params["dem_log2fc_thr"])
    total_rows, pass_rows = _summarize_dem_output(diag_out, official_adj_thr, official_log2fc_thr)
    pd.DataFrame(
        [
            {
                "status": "diagnostic_done",
                "diagnostic_log": str(diag_log),
                "diagnostic_hdf5": str(diag_out),
                "total_rows": total_rows,
                "official_threshold_pass_rows": pass_rows,
                "relaxed_adjpval_thr": relaxed_adjpval_thr,
                "relaxed_log2fc_thr": relaxed_log2fc_thr,
                "relaxed_mean_fg_thr": relaxed_mean_fg_thr,
            }
        ]
    ).to_csv(diag_summary_tsv, sep="\t", index=False)
    diag_summary_md.write_text(
        "\n".join(
            [
                "# DEM relaxed-threshold diagnostic",
                "",
                f"- relaxed thresholds: adjpval={relaxed_adjpval_thr}, log2fc={relaxed_log2fc_thr}, mean_fg={relaxed_mean_fg_thr}",
                f"- diagnostic log: {diag_log}",
                f"- diagnostic html: {diag_html}",
                f"- diagnostic hdf5: {diag_out}",
                f"- total motif rows: {total_rows}",
                f"- official-threshold-pass rows (adj <= {official_adj_thr}, log2fc >= {official_log2fc_thr}): {pass_rows}",
            ]
        )
        + "\n"
    )


def run_one(
    label: str,
    method: str,
    folder: Path,
    config: dict[str, Any],
    out_dir: Path,
    log_dir: Path,
    force: bool,
    relaxed_adjpval_thr: float,
    relaxed_log2fc_thr: float,
    relaxed_mean_fg_thr: float,
) -> Path | None:
    out = out_dir / f"{method}_{label}.hdf5"
    html = out_dir / f"{method}_{label}.html"
    empty_marker = out_dir / f"{method}_{label}.empty.tsv"
    log = log_dir / f"{method}_{label}_{datetime.now():%Y%m%d_%H%M%S}.log"
    if out.exists() and out.stat().st_size > 0 and not force:
        print(f"SKIP {method} {label}: {out}")
        return out
    if method == "dem" and empty_marker.exists() and empty_marker.stat().st_size > 0 and not force:
        print(f"SKIP {method} {label}: previously recorded empty result {empty_marker}, diagnostics already available.")
        return None
    cmd = ctx_command(config, folder, out, html) if method == "ctx" else dem_command(config, folder, out, html)
    log_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"RUN {method} {label}: log={log}")
    env = os.environ.copy()
    thread_count = str(method_n_cpu(method, config))
    for key in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        env[key] = thread_count
    with log.open("w") as handle:
        handle.write("COMMAND\t" + " ".join(cmd) + "\n")
        handle.write(
            "THREAD_ENV\t"
            + "\t".join(f"{key}={env.get(key, '')}" for key in [
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            ])
            + "\n"
        )
        proc = subprocess.run(cmd, cwd=PROJECT, stdout=handle, stderr=subprocess.STDOUT, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"{method} {label} failed; inspect {log}")
    if not out.exists() or out.stat().st_size == 0:
        if method == "dem":
            pd.DataFrame(
                [
                    {
                        "method": method,
                        "label": label,
                        "region_set_folder": str(folder),
                        "log": str(log),
                        "status": "empty",
                        "reason": "DEM run completed but produced no HDF5 under formal thresholds.",
                    }
                ]
            ).to_csv(empty_marker, sep="\t", index=False)
            print(f"WARN {method} {label}: no non-empty HDF5; recorded {empty_marker}")
            run_dem_diagnostic(
                label=label,
                folder=folder,
                config=config,
                out_dir=out_dir,
                log_dir=log_dir,
                relaxed_adjpval_thr=relaxed_adjpval_thr,
                relaxed_log2fc_thr=relaxed_log2fc_thr,
                relaxed_mean_fg_thr=relaxed_mean_fg_thr,
            )
            return None
        raise RuntimeError(f"{method} {label} completed without non-empty output: {out}")
    print(f"DONE {method} {label}: {out}")
    return out


def run_chunks(
    method: str,
    chunks: list[tuple[str, Path]],
    config: dict[str, Any],
    out_dir: Path,
    log_dir: Path,
    max_parallel: int,
    force: bool,
    relaxed_adjpval_thr: float,
    relaxed_log2fc_thr: float,
    relaxed_mean_fg_thr: float,
) -> list[Path]:
    outputs: list[Path] = []
    if max_parallel <= 1:
        for label, folder in chunks:
            out = run_one(
                label,
                method,
                folder,
                config,
                out_dir,
                log_dir,
                force,
                relaxed_adjpval_thr,
                relaxed_log2fc_thr,
                relaxed_mean_fg_thr,
            )
            if out is not None:
                outputs.append(out)
        return outputs
    with cf.ThreadPoolExecutor(max_workers=max_parallel) as pool:
        futures = [
            pool.submit(
                run_one,
                label,
                method,
                folder,
                config,
                out_dir,
                log_dir,
                force,
                relaxed_adjpval_thr,
                relaxed_log2fc_thr,
                relaxed_mean_fg_thr,
            )
            for label, folder in chunks
        ]
        for future in cf.as_completed(futures):
            out = future.result()
            if out is not None:
                outputs.append(out)
    return sorted(outputs)


def ensure_acc_gex(config: dict[str, Any]) -> Path:
    out = p(config["output_data"]["combined_GEX_ACC_mudata"])
    if out.exists() and out.stat().st_size > 0:
        return out
    bc_transform = str(config["params_data_preparation"]["bc_transform_func"]).strip()
    if (bc_transform.startswith('"') and bc_transform.endswith('"')) or (
        bc_transform.startswith("'") and bc_transform.endswith("'")
    ):
        bc_transform = bc_transform[1:-1]
    cmd = [
        "scenicplus",
        "prepare_data",
        "prepare_GEX_ACC",
        "--cisTopic_obj_fname",
        str(p(config["input_data"]["cisTopic_obj_fname"])),
        "--GEX_anndata_fname",
        str(p(config["input_data"]["GEX_anndata_fname"])),
        "--out_file",
        str(out),
        "--bc_transform_func",
        bc_transform,
    ]
    if not bool(config["params_data_preparation"].get("is_multiome", True)):
        cmd.append("--is_not_multiome")
    if config["params_data_preparation"].get("key_to_group_by"):
        cmd.extend(["--key_to_group_by", str(config["params_data_preparation"]["key_to_group_by"])])
    print("RUN", " ".join(cmd))
    subprocess.run(cmd, cwd=PROJECT, check=True)
    return out


def prepare_menr(config: dict[str, Any], motif_results: list[Path]) -> None:
    if not motif_results:
        raise ValueError("No motif enrichment results were provided to prepare_menr.")
    multiome = ensure_acc_gex(config)
    cmd = [
        "scenicplus",
        "prepare_data",
        "prepare_menr",
        "--paths_to_motif_enrichment_results",
        *[str(x) for x in motif_results],
        "--multiome_mudata_fname",
        str(multiome),
        "--out_file_tf_names",
        str(p(config["output_data"]["tf_names"])),
        "--out_file_direct_annotation",
        str(p(config["output_data"]["cistromes_direct"])),
        "--out_file_extended_annotation",
        str(p(config["output_data"]["cistromes_extended"])),
        "--direct_annotation",
        str(config["params_data_preparation"]["direct_annotation"]),
        "--extended_annotation",
        str(config["params_data_preparation"]["extended_annotation"]),
    ]
    print("RUN", " ".join(cmd))
    subprocess.run(cmd, cwd=PROJECT, check=True)


def valid_hdf5(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(8) == b"\x89HDF\r\n\x1a\n"
    except OSError:
        return False


def nonempty_text(path: Path) -> bool:
    try:
        return path.exists() and path.stat().st_size > 0 and bool(path.read_text(errors="ignore").strip())
    except OSError:
        return False


def write_status(config: dict[str, Any], chunks: list[tuple[str, Path]], out_dir: Path) -> None:
    diag_dir = PROJECT / "results" / "scenicplus_diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for label, _folder in chunks:
        for method in ["ctx", "dem"]:
            hdf5 = out_dir / f"{method}_{label}.hdf5"
            empty = out_dir / f"{method}_{label}.empty.tsv"
            diag = diag_dir / f"dem_{label}_relaxed_threshold_diagnostic.tsv"
            if valid_hdf5(hdf5):
                status = "ok"
                message = "valid HDF5"
            elif method == "dem" and empty.exists() and empty.stat().st_size > 0:
                status = "empty_with_diagnostic" if diag.exists() and diag.stat().st_size > 0 else "empty_missing_diagnostic"
                message = "formal DEM empty; diagnostic recorded" if status == "empty_with_diagnostic" else "formal DEM empty but diagnostic missing"
            else:
                status = "missing_or_invalid"
                message = "missing, empty, or invalid HDF5"
            rows.append({"method": method, "label": label, "path": str(hdf5), "status": status, "message": message})

    output_data = config["output_data"]
    for key in ["cistromes_direct", "cistromes_extended"]:
        path = p(output_data[key])
        rows.append({
            "method": "prepare_menr",
            "label": key,
            "path": str(path),
            "status": "ok" if valid_hdf5(path) else "missing_or_invalid",
            "message": "valid HDF5" if valid_hdf5(path) else "missing, empty, or invalid HDF5",
        })
    tf_path = p(output_data["tf_names"])
    rows.append({
        "method": "prepare_menr",
        "label": "tf_names",
        "path": str(tf_path),
        "status": "ok" if nonempty_text(tf_path) else "missing_or_invalid",
        "message": "non-empty text" if nonempty_text(tf_path) else "missing or empty text",
    })

    status_tsv = diag_dir / "motif_enrichment_split_status.tsv"
    status_md = diag_dir / "motif_enrichment_split_status.md"
    df = pd.DataFrame(rows)
    df.to_csv(status_tsv, sep="\t", index=False)

    bad = df[df["status"].isin(["missing_or_invalid", "empty_missing_diagnostic"])]
    dem_rows = df[df["method"] == "dem"]
    all_dem_empty = len(dem_rows) > 0 and not any(dem_rows["status"] == "ok")
    ready = bad.empty and not all_dem_empty
    status_md.write_text(
        "\n".join([
            "# Motif enrichment split status",
            "",
            f"- ready_for_step9_and_step10: {'yes' if ready else 'no'}",
            f"- checked_chunks: {len(chunks)}",
            f"- status_table: {status_tsv}",
            f"- invalid_or_missing_records: {len(bad)}",
            f"- all_dem_chunks_empty: {'yes' if all_dem_empty else 'no'}",
        ])
        + "\n"
    )
    print(status_md.read_text().strip())
    if not ready:
        raise RuntimeError(f"Motif enrichment split is not ready; inspect {status_tsv}")


def main() -> None:
    args = parse_args()
    config_path = p(args.config) if args.config else PROJECT / "work" / "scenicplus" / "Snakemake" / "config" / "config.yaml"
    config = read_config(config_path)
    split_root = p(args.split_root)
    out_dir = p(args.out_dir)
    log_dir = p(args.log_dir)
    chunks = discover_or_create_splits(config, split_root)
    max_parallel, plan = resolve_parallel_chunks(args.max_parallel_chunks, config)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([plan]).to_csv(out_dir / "motif_enrichment_split_resource_plan.tsv", sep="\t", index=False)
    pd.DataFrame([{"label": label, "region_set_folder": str(folder)} for label, folder in chunks]).to_csv(
        out_dir / "motif_enrichment_split_chunks.tsv",
        sep="\t",
        index=False,
    )
    if args.mode == "status":
        write_status(config, chunks, out_dir)
        return
    print(f"Motif enrichment chunks: {len(chunks)}; max_parallel_chunks={max_parallel}")

    motif_results: list[Path] = []
    if args.mode in {"ctx", "both"}:
        motif_results.extend(
            run_chunks(
                "ctx",
                chunks,
                config,
                out_dir,
                log_dir,
                max_parallel,
                args.force,
                args.dem_relaxed_adjpval_thr,
                args.dem_relaxed_log2fc_thr,
                args.dem_relaxed_mean_fg_thr,
            )
        )
    if args.mode in {"dem", "both"}:
        motif_results.extend(
            run_chunks(
                "dem",
                chunks,
                config,
                out_dir,
                log_dir,
                max_parallel,
                args.force,
                args.dem_relaxed_adjpval_thr,
                args.dem_relaxed_log2fc_thr,
                args.dem_relaxed_mean_fg_thr,
            )
        )
    if args.mode == "prepare-menr":
        motif_results = sorted(out_dir.glob("ctx_*.hdf5")) + sorted(out_dir.glob("dem_*.hdf5"))
    if args.mode in {"dem", "both", "prepare-menr"}:
        dem_results = sorted(out_dir.glob("dem_*.hdf5"))
        if not dem_results:
            dem_empty = [out_dir / f"dem_{label}.empty.tsv" for label, _ in chunks]
            if any(x.exists() and x.stat().st_size > 0 for x in dem_empty):
                raise RuntimeError(
                    "All DEM chunks are empty under formal thresholds. "
                    f"Relaxed-threshold diagnostic summaries are in {PROJECT / 'results' / 'scenicplus_diagnostics'} "
                    "and should be reviewed before continuing."
                )
    if args.mode in {"both", "prepare-menr"} and not args.skip_prepare_menr:
        available = sorted(out_dir.glob("ctx_*.hdf5")) + sorted(out_dir.glob("dem_*.hdf5"))
        if not available:
            raise RuntimeError(f"prepare_menr needs at least one motif HDF5 file; found none in {out_dir}")
        prepare_menr(config, available)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
