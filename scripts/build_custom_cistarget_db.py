#!/usr/bin/env python
"""Build custom cisTarget motif databases from the pycisTopic consensus peaks."""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", default=None, help="Default: $PROJECT_DIR/inputs/cistarget_db_params.tsv")
    parser.add_argument("--consensus-bed", default=None, help="Default: $PROJECT_DIR/work/pycistopic/consensus_peaks.bed")
    parser.add_argument("--allowed-chroms", default=None, help="Default: $CHROMS resolved under $PROJECT_DIR")
    parser.add_argument("--genome", default=None, help="Default: $GENOME resolved under $PROJECT_DIR")
    parser.add_argument("--motif-dir", default=None, help="Default: $PROJECT_DIR/resources/motifs/v10nr_clust_public/singletons")
    parser.add_argument("--motifs", default=None, help="Default: $PROJECT_DIR/resources/motifs/motifs.txt")
    parser.add_argument("--out-prefix", default=None, help="Default: $PROJECT_DIR/inputs/cistarget_db/custom")
    parser.add_argument("--work-dir", default=None, help="Default: $PROJECT_DIR/work/cistarget_db")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_project_path(path_value: str | None, default_rel: str, base: Path) -> Path:
    path = Path(path_value).expanduser() if path_value else base / default_rel
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(
            f"cisTarget DB parameter table not found: {path}\n"
            "Create inputs/cistarget_db_params.tsv before running this script."
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.shape[1] < 2:
        raise ValueError(f"{path} must contain two tab-separated columns: parameter and value")
    key_col, value_col = df.columns[:2]
    params = {str(k).strip(): str(v).strip() for k, v in zip(df[key_col], df[value_col]) if str(k).strip()}
    required = ["n_cpu", "seed"]
    missing = [x for x in required if x not in params or params[x] == ""]
    if missing:
        raise ValueError(f"{path} missing required parameters: {', '.join(missing)}")
    return params


def detect_memory_gb() -> tuple[float, float]:
    try:
        import psutil

        mem = psutil.virtual_memory()
        return float(mem.total) / (1024 ** 3), float(mem.available) / (1024 ** 3)
    except Exception:
        pass
    try:
        if hasattr(os, "sysconf"):
            total = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / (1024 ** 3)
            return float(total), float(total) * 0.5
    except Exception:
        pass
    return 16.0, 8.0


def detect_load_fraction() -> float:
    cpu_count = os.cpu_count() or 1
    try:
        import psutil

        return max(0.0, min(1.0, psutil.cpu_percent(interval=1.0) / 100.0))
    except Exception:
        pass
    try:
        load1 = os.getloadavg()[0]
        return max(0.0, min(1.0, float(load1) / float(cpu_count)))
    except Exception:
        return 0.0


def as_float(params: dict[str, str], key: str, default: float) -> float:
    value = str(params.get(key, "")).strip()
    if value == "":
        return default
    if value.lower() == "auto":
        return default
    return float(value)


def as_int(params: dict[str, str], key: str, default: int) -> int:
    value = str(params.get(key, "")).strip()
    if value == "" or value.lower() == "auto":
        return default
    return int(float(value))


def as_bool_auto(params: dict[str, str], key: str, default: str = "auto") -> str:
    value = str(params.get(key, default)).strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return "true"
    if value in {"0", "false", "no", "off"}:
        return "false"
    return "auto"


def resolve_n_cpu(params: dict[str, str]) -> tuple[int, dict[str, str]]:
    raw = str(params.get("n_cpu", "auto")).strip().lower()
    cpu_count = os.cpu_count() or 1
    total_gb, available_gb = detect_memory_gb()
    load_fraction = detect_load_fraction()
    max_memory_raw = str(params.get("max_memory_gb", "auto")).strip().lower()
    if max_memory_raw == "auto" or max_memory_raw == "":
        memory_budget_gb = max(1.0, min(total_gb * 0.65, available_gb))
        memory_source = "auto_min_available_or_65_percent_total"
    else:
        memory_budget_gb = max(1.0, min(float(max_memory_raw), available_gb))
        memory_source = "user_configured_capped_by_available"
    min_free_gb = as_float(params, "min_free_memory_gb", max(4.0, total_gb * 0.10))
    # Cluster-Buster workers are CPU-bound and usually modest in memory; the
    # full cisTarget score/ranking matrices remain the dominant memory cost.
    memory_per_worker_gb = as_float(params, "memory_gb_per_worker", 1.0)
    max_workers_raw = str(params.get("max_workers", "")).strip()
    max_workers = int(max_workers_raw) if max_workers_raw and max_workers_raw.lower() != "auto" else cpu_count
    max_cpu_load_fraction = as_float(params, "max_cpu_load_fraction", 0.80)
    if raw and raw != "auto":
        resolved = max(1, min(int(raw), cpu_count, max_workers))
        mode = "fixed_user_n_cpu"
    else:
        idle_fraction = max(0.05, 1.0 - load_fraction)
        cpu_cap = max(1, int(cpu_count * idle_fraction))
        if load_fraction > max_cpu_load_fraction:
            cpu_cap = 1
        mem_for_workers = max(0.0, memory_budget_gb - min_free_gb)
        memory_cap = max(1, int(mem_for_workers // memory_per_worker_gb)) if memory_per_worker_gb > 0 else max_workers
        resolved = max(1, min(cpu_cap, memory_cap, max_workers, cpu_count))
        mode = "auto_load_memory_aware"
    report = {
        "mode": mode,
        "raw_n_cpu": raw or "auto",
        "resolved_n_cpu": str(resolved),
        "cpu_count": str(cpu_count),
        "load_fraction": f"{load_fraction:.3f}",
        "max_cpu_load_fraction": f"{max_cpu_load_fraction:.3f}",
        "total_memory_gb": f"{total_gb:.2f}",
        "available_memory_gb": f"{available_gb:.2f}",
        "memory_budget_gb": f"{memory_budget_gb:.2f}",
        "memory_source": memory_source,
        "min_free_memory_gb": f"{min_free_gb:.2f}",
        "memory_gb_per_worker": f"{memory_per_worker_gb:.2f}",
        "max_workers": str(max_workers),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    return resolved, report


def resolve_partial_plan(
    params: dict[str, str],
    resource_report: dict[str, str],
    n_regions: int,
    n_motifs: int,
) -> tuple[bool, int, dict[str, str]]:
    matrix_gb = n_regions * n_motifs * 4.0 / (1024 ** 3)
    memory_budget_gb = float(resource_report["memory_budget_gb"])
    use_partial_raw = as_bool_auto(params, "use_partial", "auto")
    target_part_gb = as_float(params, "target_partial_matrix_gb", 3.0)
    max_parts = int(as_float(params, "max_partial_parts", 64.0))
    requested_parts = str(params.get("partial_n_parts", "auto")).strip().lower()

    if requested_parts and requested_parts != "auto":
        n_parts = max(1, int(requested_parts))
    else:
        n_parts = max(1, int(math.ceil(matrix_gb / max(0.5, target_part_gb))))
        n_parts = min(max_parts, n_parts)

    if use_partial_raw == "true":
        use_partial = True
    elif use_partial_raw == "false":
        use_partial = False
    else:
        use_partial = matrix_gb > max(8.0, memory_budget_gb * 0.50)

    if not use_partial:
        n_parts = 1

    report = {
        "full_scores_matrix_gb_float32": f"{matrix_gb:.2f}",
        "use_partial": "1" if use_partial else "0",
        "partial_n_parts": str(n_parts),
        "target_partial_matrix_gb": f"{target_part_gb:.2f}",
        "max_partial_parts": str(max_parts),
    }
    return use_partial, n_parts, report


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required command not found on PATH: {name}")
    return path


def to_ucsc(chrom: str) -> str:
    chrom = chrom.split()[0]
    if chrom in {"MT", "M", "Mt", "mitochondrion_genome"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else "chr" + chrom


def write_named_consensus(consensus_bed: Path, allowed_chroms: Path, out_bed: Path) -> int:
    allowed = set(allowed_chroms.read_text().splitlines())
    rows: list[tuple[str, int, int]] = []
    with consensus_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = to_ucsc(fields[0])
            if chrom not in allowed:
                continue
            start, end = int(fields[1]), int(fields[2])
            if end <= start:
                continue
            rows.append((chrom, start, end))
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with out_bed.open("w") as out:
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    if not rows:
        raise ValueError(f"No valid consensus regions written from {consensus_bed}")
    return len(rows)


def check_motifs(motif_dir: Path, motifs_txt: Path) -> int:
    motifs = [x.strip() for x in motifs_txt.read_text().splitlines() if x.strip()]
    missing = [m for m in motifs if not (motif_dir / f"{m}.cb").exists()]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} motif .cb files; first missing: {missing[:10]}")
    return len(motifs)


def child_process_summary(pid: int) -> str:
    try:
        import psutil

        parent = psutil.Process(pid)
        children = parent.children(recursive=True)
        running = [p for p in children if p.is_running()]
        names: dict[str, int] = {}
        rss = 0
        for proc in running:
            try:
                names[proc.name()] = names.get(proc.name(), 0) + 1
                rss += int(proc.memory_info().rss)
            except Exception:
                continue
        name_text = ",".join(f"{name}:{count}" for name, count in sorted(names.items())) or "none"
        return f"children={len(running)} child_names={name_text} child_rss_gb={rss / (1024 ** 3):.2f}"
    except Exception:
        return "children=unknown"


def watched_path_summary(paths: list[Path]) -> str:
    if not paths:
        return ""
    parts = []
    for path in paths:
        if path.exists():
            parts.append(f"{path.name}=present:{path.stat().st_size}")
        else:
            parts.append(f"{path.name}=missing")
    return " watch=" + ",".join(parts)


def run_checked(
    cmd: list[str],
    heartbeat_seconds: int = 0,
    heartbeat_label: str = "",
    watch_paths: list[Path] | None = None,
) -> None:
    print("RUN " + " ".join(map(str, cmd)), flush=True)
    if heartbeat_seconds <= 0:
        subprocess.run(cmd, check=True)
        return
    started = time.monotonic()
    last_heartbeat = started
    proc = subprocess.Popen(cmd)
    while True:
        return_code = proc.poll()
        if return_code is not None:
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, cmd)
            return
        now = time.monotonic()
        if now - last_heartbeat >= heartbeat_seconds:
            elapsed = (now - started) / 60.0
            label = heartbeat_label or Path(str(cmd[0])).name
            print(
                "HEARTBEAT "
                f"{label} elapsed_min={elapsed:.1f} pid={proc.pid} "
                f"{child_process_summary(proc.pid)}"
                f"{watched_path_summary(watch_paths or [])}",
                flush=True,
            )
            last_heartbeat = now
        time.sleep(min(10, max(1, heartbeat_seconds // 6)))


def final_outputs(out_prefix: Path) -> tuple[Path, Path, Path]:
    rankings = Path(str(out_prefix) + ".regions_vs_motifs.rankings.feather")
    scores = Path(str(out_prefix) + ".regions_vs_motifs.scores.feather")
    motif_scores = Path(str(out_prefix) + ".motifs_vs_regions.scores.feather")
    return rankings, scores, motif_scores


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    params_path = resolve_project_path(args.params, "inputs/cistarget_db_params.tsv", pdir)
    params = read_params(params_path)
    n_cpu, resource_report = resolve_n_cpu(params)
    consensus_bed = resolve_project_path(args.consensus_bed, "work/pycistopic/consensus_peaks.bed", pdir)
    chroms_value = args.allowed_chroms or os.environ.get("CHROMS")
    genome_value = args.genome or os.environ.get("GENOME")
    if not chroms_value:
        raise SystemExit("ERROR: --allowed-chroms is required unless CHROMS is set in project_env.sh.")
    if not genome_value:
        raise SystemExit("ERROR: --genome is required unless GENOME is set in project_env.sh.")
    allowed_chroms = resolve_project_path(chroms_value, chroms_value, pdir)
    genome = resolve_project_path(genome_value, genome_value, pdir)
    motif_dir = resolve_project_path(args.motif_dir, "resources/motifs/v10nr_clust_public/singletons", pdir)
    motifs_txt = resolve_project_path(args.motifs, "resources/motifs/motifs.txt", pdir)
    out_prefix = resolve_project_path(args.out_prefix, "inputs/cistarget_db/custom", pdir)
    work_dir = resolve_project_path(args.work_dir, "work/cistarget_db", pdir)
    named_bed = work_dir / "consensus_regions.named.bed"
    fasta = work_dir / "consensus_regions.fa"
    manifest = pdir / "results" / "cistarget_db" / "custom_cistarget_db_manifest.tsv"

    for path in [consensus_bed, allowed_chroms, genome, motif_dir, motifs_txt]:
        if not path.exists():
            raise FileNotFoundError(path)
    bedtools = require_tool("bedtools")
    cbust = require_tool("cbust")
    conda_prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("CONDA_ENV_PREFIX")
    if not conda_prefix:
        raise SystemExit("ERROR: CONDA_PREFIX or CONDA_ENV_PREFIX must be set.")
    create_db = Path(conda_prefix) / "opt" / "create_cisTarget_databases" / "create_cistarget_motif_databases.py"
    if not create_db.exists():
        raise FileNotFoundError(create_db)

    n_regions = write_named_consensus(consensus_bed, allowed_chroms, named_bed)
    n_motifs = check_motifs(motif_dir, motifs_txt)
    use_partial, n_parts, partial_report = resolve_partial_plan(params, resource_report, n_regions, n_motifs)
    heartbeat_seconds = max(0, as_int(params, "heartbeat_seconds", 600))
    resource_report["heartbeat_seconds"] = str(heartbeat_seconds)
    resource_report.update(partial_report)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    run_checked([bedtools, "getfasta", "-fi", str(genome), "-bed", str(named_bed), "-fo", str(fasta)])
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    resource_path = pdir / "results" / "cistarget_db" / "custom_cistarget_resource_plan.tsv"
    resource_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"parameter": key, "value": value} for key, value in resource_report.items()]
    ).to_csv(resource_path, sep="\t", index=False)
    print(f"WROTE {resource_path}")
    print(f"custom cisTarget resolved n_cpu={n_cpu} ({resource_report['mode']})")
    print(
        "custom cisTarget partial mode="
        f"{'on' if use_partial else 'off'}"
        f" parts={n_parts}"
        f" full_scores_matrix_gb={partial_report['full_scores_matrix_gb_float32']}"
    )

    python_bin = str(Path(conda_prefix) / "bin" / "python")
    common_cmd = [
        python_bin,
        str(create_db),
        "-f",
        str(fasta),
        "-M",
        str(motif_dir),
        "-m",
        str(motifs_txt),
        "-c",
        cbust,
        "-t",
        str(n_cpu),
        "-s",
        str(int(params["seed"])),
    ]
    rankings, scores, motif_scores = final_outputs(out_prefix)
    if rankings.exists() and rankings.stat().st_size > 0 and scores.exists() and scores.stat().st_size > 0:
        print(f"Final cisTarget outputs already exist; skipping build: {rankings}, {scores}")
    elif use_partial:
        partial_dir = work_dir / "partial_scores"
        partial_dir.mkdir(parents=True, exist_ok=True)
        partial_prefix = partial_dir / out_prefix.name
        for part in range(1, n_parts + 1):
            partial_file = Path(
                f"{partial_prefix}.part_{part:04d}_of_{n_parts:04d}.motifs_vs_regions.scores.feather"
            )
            if partial_file.exists() and partial_file.stat().st_size > 0:
                print(f"SKIP existing partial {part}/{n_parts}: {partial_file}")
                continue
            run_checked(
                common_cmd + ["-o", str(partial_prefix), "--partial", str(part), str(n_parts)],
                heartbeat_seconds=heartbeat_seconds,
                heartbeat_label=f"custom_cistarget_partial_{part}_of_{n_parts}",
                watch_paths=[partial_file],
            )

        combine_script = Path(conda_prefix) / "opt" / "create_cisTarget_databases" / "combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py"
        convert_script = Path(conda_prefix) / "opt" / "create_cisTarget_databases" / "convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py"
        for script in [combine_script, convert_script]:
            if not script.exists():
                raise FileNotFoundError(script)
        if not scores.exists() or scores.stat().st_size == 0 or not motif_scores.exists() or motif_scores.stat().st_size == 0:
            run_checked(
                [python_bin, str(combine_script), "-i", str(partial_dir), "-o", str(out_prefix.parent)],
                heartbeat_seconds=heartbeat_seconds,
                heartbeat_label="custom_cistarget_combine_partial_scores",
                watch_paths=[motif_scores, scores],
            )
        if not rankings.exists() or rankings.stat().st_size == 0:
            if not motif_scores.exists() or motif_scores.stat().st_size == 0:
                raise FileNotFoundError(motif_scores)
            run_checked(
                [python_bin, str(convert_script), "-i", str(motif_scores), "-s", str(int(params["seed"]))],
                heartbeat_seconds=heartbeat_seconds,
                heartbeat_label="custom_cistarget_scores_to_rankings",
                watch_paths=[rankings],
            )
    else:
        run_checked(
            common_cmd + ["-o", str(out_prefix)],
            heartbeat_seconds=heartbeat_seconds,
            heartbeat_label="custom_cistarget_full_database",
            watch_paths=[motif_scores, rankings, scores],
        )

    for path in [rankings, scores]:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"item": "named_consensus_bed", "path": str(named_bed), "n": n_regions},
            {"item": "consensus_fasta", "path": str(fasta), "n": n_regions},
            {"item": "motifs_txt", "path": str(motifs_txt), "n": n_motifs},
            {"item": "motif_dir", "path": str(motif_dir), "n": n_motifs},
            {"item": "resolved_n_cpu", "path": resource_report["mode"], "n": n_cpu},
            {"item": "partial_mode", "path": "on" if use_partial else "off", "n": n_parts},
            {"item": "rankings_feather", "path": str(rankings), "n": rankings.stat().st_size},
            {"item": "scores_feather", "path": str(scores), "n": scores.stat().st_size},
        ]
    ).to_csv(manifest, sep="\t", index=False)
    print(f"WROTE {manifest}")
    print(f"WROTE {rankings}")
    print(f"WROTE {scores}")


if __name__ == "__main__":
    main()
