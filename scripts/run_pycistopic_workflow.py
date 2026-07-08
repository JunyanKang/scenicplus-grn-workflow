from __future__ import annotations

from pathlib import Path
from datetime import datetime
import csv
import hashlib
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import traceback
from typing import Mapping, Sequence

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT / "tmp" / "numba"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)
(PROJECT / "tmp" / "numba").mkdir(parents=True, exist_ok=True)
INPUTS = PROJECT / "inputs"
WORK = PROJECT / "work" / "pycistopic"
RESULTS = PROJECT / "results" / "pycistopic"
LOGS = PROJECT / "logs"

if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
    print(
        "Usage: run_pycistopic_workflow.py\n\n"
        "Run pseudobulk peak calling, consensus peak construction, cisTopic object "
        "creation, topic modeling, topic/DAR region set export "
        "and required pycisTopic summary figures for the active project. Requires "
        "$PROJECT_DIR and the parameter tables created by setup_workflow_params.py."
    )
    raise SystemExit(0)

required_start_files = [
    INPUTS / "pycistopic_params.tsv",
    INPUTS / "sample_sheet.tsv",
    INPUTS / "cell_metadata.tsv",
    INPUTS / "topic_model_grid.tsv",
    INPUTS / "atac_qc_thresholds.tsv",
]
missing_start = [str(path) for path in required_start_files if not path.exists() or path.stat().st_size == 0]
if missing_start:
    raise SystemExit(
        "ERROR: pycisTopic inputs are not ready. Run setup/ATAC preparation first. Missing: "
        + ", ".join(missing_start)
    )

with (INPUTS / "pycistopic_params.tsv").open() as handle:
    start_params = {row["parameter"]: row["value"] for row in csv.DictReader(handle, delimiter="\t") if row.get("parameter")}
start_chromsizes = start_params.get("chromsizes", "")
if not start_chromsizes:
    raise SystemExit("ERROR: inputs/pycistopic_params.tsv must define chromsizes.")
start_chromsizes_path = Path(start_chromsizes).expanduser()
if not start_chromsizes_path.is_absolute():
    start_chromsizes_path = PROJECT / start_chromsizes_path
if not start_chromsizes_path.exists() or start_chromsizes_path.stat().st_size == 0:
    raise SystemExit(f"ERROR: chromsizes not found: {start_chromsizes_path}")

import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
try:
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover - optional UI dependency
    class tqdm:
        def __init__(self, total=None, desc=None, unit=None):
            self.total = total
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc, tb):
            return False
        def update(self, n=1):
            return None

_env_prefix = os.environ.get("CONDA_ENV_PREFIX") or os.environ.get("CONDA_PREFIX")
_script_modules = Path(__file__).resolve().parent.parent / "modules"
_default_modules = (
    Path(os.environ["SCENICPLUS_HOME"]) / "modules"
    if os.environ.get("SCENICPLUS_HOME")
    else _script_modules
)
INSTALLER_MODULES = Path(os.environ.get("SCENICPLUS_INSTALLER_MODULES", str(_default_modules)))
if INSTALLER_MODULES.exists():
    sys.path.insert(0, str(INSTALLER_MODULES))
try:
    from autozyme_runtime import activate_autozyme_python
    activate_autozyme_python(("scanpy",))
except Exception as exc:
    print(f"[autozyme] activation helper unavailable: {exc}", file=sys.stderr)


def patch_polars_read_csv_sep_alias() -> None:
    """Keep pycisTopic compatible with current polars read_csv keyword names."""
    try:
        import inspect
        import polars as pl
    except Exception:
        return
    signature = inspect.signature(pl.read_csv)
    if "sep" in signature.parameters or "separator" not in signature.parameters:
        return
    original = pl.read_csv

    def read_csv_compat(*args, **kwargs):
        if "sep" in kwargs and "separator" not in kwargs:
            kwargs["separator"] = kwargs.pop("sep")
        return original(*args, **kwargs)

    pl.read_csv = read_csv_compat


patch_polars_read_csv_sep_alias()

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.diff_features import find_diff_features, impute_accessibility
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.lda_models import evaluate_models, run_cgs_models, run_cgs_models_mallet
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.topic_binarization import binarize_topics


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data: str) -> int:
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def mkdirs() -> None:
    for path in [
        WORK,
        WORK / "pseudobulk_bed",
        WORK / "pseudobulk_bigwig",
        WORK / "macs2",
        WORK / "models",
        INPUTS / "region_sets" / "Topics_otsu",
        INPUTS / "region_sets" / "Topics_top_3k",
        INPUTS / "region_sets" / "DARs_cell_label",
        RESULTS / "model_selection",
        RESULTS / "qc",
        RESULTS / "dar",
        LOGS,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if not {"parameter", "value"}.issubset(df.columns):
        raise ValueError(f"{path} must contain columns: parameter, value")
    return dict(zip(df["parameter"].astype(str), df["value"].astype(str)))


def as_int(params: Mapping[str, str], key: str, default: int) -> int:
    value = params.get(key, "")
    return default if value == "" else int(value)


def as_float(params: Mapping[str, str], key: str, default: float) -> float:
    value = params.get(key, "")
    return default if value == "" else float(value)


def detect_physical_memory_gb() -> float:
    try:
        import psutil

        return float(psutil.virtual_memory().total) / (1024 ** 3)
    except Exception:
        pass
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return float(pages * page_size) / (1024 ** 3)
    except Exception:
        return 16.0


def memory_budget_gb(params: Mapping[str, str]) -> tuple[float, str, float]:
    raw = (
        params.get("max_memory_gb", "").strip()
        or os.environ.get("SCENICPLUS_MAX_MEMORY_GB", "").strip()
        or "auto"
    )
    physical_gb = detect_physical_memory_gb()
    if raw.lower() == "auto":
        return max(1.0, physical_gb * 0.65), "auto_detected_65_percent_physical_ram", physical_gb
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError("max_memory_gb must be auto or a positive number.") from exc
    if value <= 0:
        raise ValueError("max_memory_gb must be auto or a positive number.")
    return value, "user_configured", physical_gb


def resolve_n_cpu(params: Mapping[str, str], memory_gb: float) -> int:
    raw = params.get("n_cpu", "").strip()
    cpu_count = os.cpu_count() or 1
    memory_limited = max(1, int(memory_gb // 12))
    requested = max(1, int(raw)) if raw and raw.lower() != "auto" else cpu_count
    return max(1, min(requested, cpu_count, memory_limited, 8))


def resolve_worker_count(raw: str | None, stage: str, n_cpu: int, memory_gb: float) -> int:
    value = (raw or "").strip()
    if value and value.lower() != "auto":
        return max(1, int(value))
    if stage == "pseudobulk":
        if memory_gb < 96:
            return 1
        if memory_gb < 192:
            return min(n_cpu, 2)
        return min(n_cpu, 4)
    if stage == "macs2":
        return max(1, min(n_cpu, int(memory_gb // 16), 8))
    return n_cpu


def resolve_topic_n_cpu(params: Mapping[str, str], n_topics: Sequence[int], memory_gb: float) -> int:
    raw = params.get("topic_n_cpu", "").strip()
    cpu_count = os.cpu_count() or 1
    if raw and raw.lower() != "auto":
        return max(1, min(int(raw), cpu_count))
    # Topic modeling is usually CPU-bound once the corpus is in memory. Keep the
    # default below the global CPU cap but avoid forcing MALLET to one thread on
    # workstation-scale memory budgets.
    topic_cap_by_memory = max(1, int(memory_gb // 8))
    return max(1, min(cpu_count, max(1, len(n_topics)), topic_cap_by_memory, 8))


def resolve_mallet_memory(params: Mapping[str, str], max_memory_gb: float, physical_memory_gb: float) -> str:
    raw = params.get("mallet_memory_gb", "").strip()
    if raw and raw.lower() != "auto":
        value = float(raw)
        if value <= 0:
            raise ValueError("mallet_memory_gb must be auto or a positive number.")
    else:
        value = min(max_memory_gb * 0.80, physical_memory_gb * 0.75)
        value = max(4.0, value)
    value = min(value, max_memory_gb)
    if value <= 0:
        raise ValueError("Resolved MALLET memory must be positive.")
    return f"{int(value)}g" if float(value).is_integer() else f"{value:.1f}g"


def write_parallelism_plan(
    *,
    physical_memory_gb: float,
    max_memory_gb: float,
    memory_source: str,
    n_cpu: int,
    pseudobulk_n_cpu: int,
    macs_n_cpu: int,
    ray_temp_dir: Path,
) -> None:
    plan = pd.DataFrame(
        [
            ("physical_memory_gb", f"{physical_memory_gb:.2f}"),
            ("max_memory_gb", f"{max_memory_gb:.2f}"),
            ("memory_source", memory_source),
            ("n_cpu", str(n_cpu)),
            ("pseudobulk_n_cpu", str(pseudobulk_n_cpu)),
            ("macs_n_cpu", str(macs_n_cpu)),
            ("ray_temp_dir", str(ray_temp_dir)),
        ],
        columns=["parameter", "value"],
    )
    plan.to_csv(RESULTS / "qc" / "parallelism_plan.tsv", sep="\t", index=False)


def is_valid_gzip(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    gzip_bin = shutil.which("gzip")
    if gzip_bin:
        return subprocess.run(
            [gzip_bin, "-t", str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode == 0
    try:
        import gzip

        with gzip.open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                if not chunk:
                    break
        return True
    except Exception:
        return False


def expected_pseudobulk_bed_paths(labels: Sequence[str]) -> dict[str, str]:
    return {
        sanitize_label(label): str(WORK / "pseudobulk_bed" / f"{sanitize_label(label)}.bed.gz")
        for label in labels
    }


def default_short_ray_temp_dir() -> Path:
    base = Path(os.environ.get("SCENICPLUS_RAY_TMPDIR") or os.environ.get("RAY_TMPDIR") or "/tmp")
    digest = hashlib.sha1(str(PROJECT).encode("utf-8")).hexdigest()[:10]
    uid = os.getuid() if hasattr(os, "getuid") else "user"
    return base / f"spgrn_ray_{uid}_{digest}"


def resolve_ray_temp_dir(value: str) -> Path:
    requested = (value or "").strip()
    if requested.lower() in {"", "auto"}:
        path = default_short_ray_temp_dir()
    else:
        path = Path(requested).expanduser()
        if not path.is_absolute():
            path = PROJECT / path
    path = path.resolve()
    projected_socket = path / "session_YYYY-MM-DD_HH-MM-SS_000000_000000" / "sockets" / "plasma_store"
    if len(str(projected_socket)) > 100:
        fallback = default_short_ray_temp_dir().resolve()
        print(
            "ray_temp_dir path is too long for Linux AF_UNIX sockets; "
            f"using {fallback} instead of {path}",
            file=sys.stderr,
        )
        path = fallback
    path.mkdir(parents=True, exist_ok=True)
    os.environ["RAY_TMPDIR"] = str(path)
    return path


def expected_macs_narrow_peak_paths(labels: Sequence[str]) -> dict[str, Path]:
    return {
        sanitize_label(label): WORK / "macs2" / f"{sanitize_label(label)}_peaks.narrowPeak"
        for label in labels
    }


def read_narrow_peak_as_pyranges(path: Path) -> pr.PyRanges:
    narrow_peak = pd.read_csv(path, sep="\t", header=None)
    if narrow_peak.shape[1] < 10:
        raise ValueError(f"narrowPeak file has fewer than 10 columns: {path}")
    narrow_peak = narrow_peak.iloc[:, :10].copy()
    narrow_peak.columns = [
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Score",
        "Strand",
        "FC_summit",
        "-log10_pval",
        "-log10_qval",
        "Summit",
    ]
    return pr.PyRanges(narrow_peak)


def require_file(path: str | Path, label: str) -> Path:
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = PROJECT / p
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


def sanitize_label(label: str) -> str:
    label = str(label).strip()
    label = re.sub(r"[^A-Za-z0-9_.]+", "_", label)
    return label.strip("_") or "unlabeled"


def raw_barcode_from_cell_id(cell_id: str, sample_id: str) -> str:
    suffix = f"-{sample_id}"
    if not str(cell_id).endswith(suffix):
        raise ValueError(
            f"cell_id '{cell_id}' does not end with expected sample suffix '{suffix}'. "
            "Use barcode-sample_id before running pycisTopic."
        )
    return str(cell_id)[: -len(suffix)]


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, str]]:
    params = read_params(INPUTS / "pycistopic_params.tsv")
    sample_sheet = pd.read_csv(INPUTS / "sample_sheet.tsv", sep="\t", dtype=str).fillna("")
    cell_meta = pd.read_csv(INPUTS / "cell_metadata.tsv", sep="\t", dtype=str).fillna("")
    topic_grid = pd.read_csv(INPUTS / "topic_model_grid.tsv", sep="\t")
    atac_qc = pd.read_csv(INPUTS / "atac_qc_thresholds.tsv", sep="\t")

    required_sample_cols = {"sample_id", "fragments_tsv_gz"}
    missing_sample = sorted(required_sample_cols - set(sample_sheet.columns))
    if missing_sample:
        raise ValueError(f"inputs/sample_sheet.tsv missing columns: {missing_sample}")

    required_cell_cols = {"cell_id", "sample_id", "cell_label"}
    missing_cell = sorted(required_cell_cols - set(cell_meta.columns))
    if missing_cell:
        raise ValueError(f"inputs/cell_metadata.tsv missing columns: {missing_cell}")

    if "n_topics" not in topic_grid.columns:
        raise ValueError("inputs/topic_model_grid.tsv must contain column: n_topics")
    topic_grid["n_topics"] = topic_grid["n_topics"].astype(int)

    sample_sheet["sample_id"] = sample_sheet["sample_id"].astype(str)
    cell_meta["sample_id"] = cell_meta["sample_id"].astype(str)
    cell_meta["cell_id"] = cell_meta["cell_id"].astype(str)
    cell_meta["cell_label"] = cell_meta["cell_label"].map(sanitize_label)

    known_samples = set(sample_sheet["sample_id"])
    unknown = sorted(set(cell_meta["sample_id"]) - known_samples)
    if unknown:
        raise ValueError(f"cell_metadata.tsv contains samples absent from sample_sheet.tsv: {unknown}")

    cell_meta["barcode"] = [
        raw_barcode_from_cell_id(cell_id, sample_id)
        for cell_id, sample_id in zip(cell_meta["cell_id"], cell_meta["sample_id"])
    ]
    if cell_meta["cell_id"].duplicated().any():
        dup = cell_meta.loc[cell_meta["cell_id"].duplicated(), "cell_id"].head().tolist()
        raise ValueError(f"Duplicate cell_id values in inputs/cell_metadata.tsv: {dup}")

    for row in sample_sheet.itertuples(index=False):
        frag = require_file(row.fragments_tsv_gz, f"fragments_tsv_gz for sample {row.sample_id}")
        require_file(str(frag) + ".tbi", f"tabix index for sample {row.sample_id}")

    chromsizes = params.get("chromsizes", "")
    if chromsizes == "":
        raise ValueError("inputs/pycistopic_params.tsv must define chromsizes")
    require_file(chromsizes, "chromsizes")

    blacklist = params.get("blacklist", "")
    if blacklist:
        require_file(blacklist, "blacklist")

    if not {"sample_id", "min_cistopic_fragments", "min_cistopic_accessible_regions"}.issubset(atac_qc.columns):
        raise ValueError(
            "inputs/atac_qc_thresholds.tsv must contain sample_id, "
            "min_cistopic_fragments, min_cistopic_accessible_regions"
        )
    atac_qc["sample_id"] = atac_qc["sample_id"].astype(str)
    missing_qc = sorted(known_samples - set(atac_qc["sample_id"]))
    if missing_qc:
        raise ValueError(f"atac_qc_thresholds.tsv lacks thresholds for samples: {missing_qc}")

    return sample_sheet, cell_meta, topic_grid, atac_qc, params


def fragment_paths_by_sample(sample_sheet: pd.DataFrame) -> dict[str, str]:
    return dict(zip(sample_sheet["sample_id"].astype(str), sample_sheet["fragments_tsv_gz"].astype(str)))


def write_cell_accounting(cell_meta: pd.DataFrame) -> None:
    summary = (
        cell_meta.groupby(["sample_id", "cell_label"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["sample_id", "cell_label"])
    )
    summary.to_csv(RESULTS / "qc" / "cells_by_sample_and_label.tsv", sep="\t", index=False)


def require_nonempty_pdf(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"required PDF was not created: {path}")


def export_cell_qc_pdf(qc: pd.DataFrame) -> None:
    pdf_path = RESULTS / "qc" / "cistopic_cell_qc.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.ravel()
        axes[0].hist(np.log10(qc["cisTopic_nr_frag"].astype(float) + 1), bins=80, color="#4C78A8", alpha=0.85)
        axes[0].set_xlabel("log10 fragments + 1")
        axes[0].set_ylabel("cells")
        axes[0].set_title("Fragment depth")
        axes[1].hist(np.log10(qc["cisTopic_nr_acc"].astype(float) + 1), bins=80, color="#72B7B2", alpha=0.85)
        axes[1].set_xlabel("log10 accessible regions + 1")
        axes[1].set_ylabel("cells")
        axes[1].set_title("Accessible regions")
        axes[2].scatter(
            np.log10(qc["cisTopic_nr_frag"].astype(float) + 1),
            np.log10(qc["cisTopic_nr_acc"].astype(float) + 1),
            s=3,
            alpha=0.25,
            color="#333333",
            linewidths=0,
        )
        axes[2].set_xlabel("log10 fragments + 1")
        axes[2].set_ylabel("log10 accessible regions + 1")
        axes[2].set_title("Cell QC phase space")
        pass_rate = qc.groupby("sample_id")["pass_atac_qc"].mean().sort_index()
        axes[3].bar(pass_rate.index.astype(str), pass_rate.values, color="#59A14F", alpha=0.9)
        axes[3].set_ylim(0, 1)
        axes[3].set_ylabel("fraction passing ATAC QC")
        axes[3].set_title("Per-sample pass rate")
        axes[3].tick_params(axis="x", rotation=45)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    require_nonempty_pdf(pdf_path)


def export_consensus_peak_qc_pdf(consensus_df: pd.DataFrame) -> None:
    pdf_path = RESULTS / "qc" / "consensus_peak_qc.pdf"
    peak_width = consensus_df["End"].astype(int) - consensus_df["Start"].astype(int)
    chrom_counts = consensus_df["Chromosome"].astype(str).value_counts().sort_index()
    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].hist(peak_width, bins=80, color="#9C755F", alpha=0.85)
        axes[0].set_xlabel("peak width (bp)")
        axes[0].set_ylabel("peaks")
        axes[0].set_title("Consensus peak widths")
        axes[1].bar(chrom_counts.index, chrom_counts.values, color="#4E79A7", alpha=0.9)
        axes[1].set_ylabel("consensus peaks")
        axes[1].set_title("Peaks per chromosome")
        axes[1].tick_params(axis="x", rotation=90, labelsize=7)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    require_nonempty_pdf(pdf_path)


def export_topic_qc_pdf(cistopic_obj) -> None:
    from pycisTopic.topic_qc import compute_topic_metrics

    metrics = compute_topic_metrics(cistopic_obj, return_metrics=True)
    metrics.to_csv(RESULTS / "model_selection" / "topic_qc_metrics.tsv", sep="\t")
    numeric = metrics.select_dtypes(include=[np.number])
    pdf_path = RESULTS / "model_selection" / "topic_qc_metrics.pdf"
    with PdfPages(pdf_path) as pdf:
        if numeric.shape[1] == 0:
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.text(0.5, 0.5, "No numeric topic QC metrics returned", ha="center", va="center")
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
        else:
            for col in numeric.columns:
                fig, ax = plt.subplots(figsize=(8, 4))
                ax.bar(metrics.index.astype(str), numeric[col].values, color="#4E79A7", alpha=0.9)
                ax.set_ylabel(col)
                ax.set_xlabel("topic")
                ax.set_title(f"Topic QC: {col}")
                ax.tick_params(axis="x", rotation=90, labelsize=7)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
    require_nonempty_pdf(pdf_path)


def export_dar_summary_pdf(summary: pd.DataFrame) -> None:
    pdf_path = RESULTS / "dar" / "dar_summary.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(summary)), 4))
        ax.bar(summary["contrast"].astype(str), summary["n_dar"].astype(int), color="#F28E2B", alpha=0.9)
        ax.set_ylabel("DAR count")
        ax.set_title("DARs by cell-label contrast")
        ax.tick_params(axis="x", rotation=90, labelsize=7)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    require_nonempty_pdf(pdf_path)


def run_pseudobulk_and_peak_calling(
    sample_sheet: pd.DataFrame,
    cell_meta: pd.DataFrame,
    params: Mapping[str, str],
) -> Path:
    max_memory_gb, memory_source, physical_memory_gb = memory_budget_gb(params)
    n_cpu = resolve_n_cpu(params, max_memory_gb)
    pseudobulk_n_cpu = resolve_worker_count(
        params.get("pseudobulk_n_cpu", "auto"),
        "pseudobulk",
        n_cpu,
        max_memory_gb,
    )
    macs_n_cpu = resolve_worker_count(
        params.get("macs_n_cpu", "auto"),
        "macs2",
        n_cpu,
        max_memory_gb,
    )
    chromsizes = params["chromsizes"]
    fragments = fragment_paths_by_sample(sample_sheet)
    ray_temp_dir = resolve_ray_temp_dir(params.get("ray_temp_dir", "auto"))
    write_parallelism_plan(
        physical_memory_gb=physical_memory_gb,
        max_memory_gb=max_memory_gb,
        memory_source=memory_source,
        n_cpu=n_cpu,
        pseudobulk_n_cpu=pseudobulk_n_cpu,
        macs_n_cpu=macs_n_cpu,
        ray_temp_dir=ray_temp_dir,
    )
    write_bigwig = params.get("write_pseudobulk_bigwig", "0").strip().lower() in {"1", "true", "yes", "on"}
    resume_pseudobulk = params.get("resume_pseudobulk", "1").strip().lower() in {"1", "true", "yes", "on"}

    pseudobulk_meta = cell_meta.set_index("cell_id", drop=False)[
        ["sample_id", "cell_label", "barcode"]
    ].copy()
    label_order = pd.Series(pseudobulk_meta["cell_label"].astype(str).unique()).tolist()
    expected_bed_paths = expected_pseudobulk_bed_paths(label_order)
    existing_ok = {
        label: is_valid_gzip(Path(path))
        for label, path in expected_bed_paths.items()
    }
    labels_to_build = [
        original
        for original in label_order
        if not (resume_pseudobulk and existing_ok.get(sanitize_label(original), False))
    ]
    pseudobulk_status = pd.DataFrame(
        [
            {
                "label": label,
                "safe_label": sanitize_label(label),
                "bed_path": expected_bed_paths[sanitize_label(label)],
                "status": "reuse_existing" if resume_pseudobulk and existing_ok.get(sanitize_label(label), False) else "build",
            }
            for label in label_order
        ]
    )
    pseudobulk_status.to_csv(RESULTS / "qc" / "pseudobulk_resume_plan.tsv", sep="\t", index=False)
    ray_kwargs = {"_temp_dir": str(ray_temp_dir)}
    bw_paths = {}
    if labels_to_build:
        for label in labels_to_build:
            candidate = Path(expected_bed_paths[sanitize_label(label)])
            if candidate.exists() and not is_valid_gzip(candidate):
                candidate.unlink()
        input_for_export = pseudobulk_meta.loc[
            pseudobulk_meta["cell_label"].astype(str).isin(labels_to_build)
        ].copy()
        bw_paths, built_bed_paths = export_pseudobulk(
            input_data=input_for_export,
            variable="cell_label",
            chromsizes=chromsizes,
            bed_path=str(WORK / "pseudobulk_bed"),
            bigwig_path=str(WORK / "pseudobulk_bigwig") if write_bigwig else None,
            path_to_fragments=fragments,
            sample_id_col="sample_id",
            n_cpu=pseudobulk_n_cpu,
            normalize_bigwig=True,
            remove_duplicates=True,
            split_pattern="-",
            **(ray_kwargs if pseudobulk_n_cpu > 1 else {}),
        )
        for label, path in built_bed_paths.items():
            expected_bed_paths[sanitize_label(label)] = path
    bed_paths = expected_bed_paths
    (RESULTS / "pseudobulk_bed_paths.json").write_text(json.dumps(bed_paths, indent=2))
    (RESULTS / "pseudobulk_bigwig_paths.json").write_text(json.dumps(bw_paths, indent=2))
    missing_bed = [
        f"{group}: {path}"
        for group, path in bed_paths.items()
        if not Path(path).exists() or Path(path).stat().st_size == 0
    ]
    if missing_bed:
        raise FileNotFoundError(
            "Pseudobulk BED files were not created or are empty: " + "; ".join(missing_bed)
        )

    macs_path = params.get("macs_path", "macs2") or "macs2"
    if shutil.which(macs_path) is None:
        raise FileNotFoundError(
            f"MACS2 executable not found: {macs_path}. Activate the SCENIC+ environment "
            "or set macs_path in inputs/pycistopic_params.tsv."
        )

    resume_macs2 = params.get("resume_macs2", "1").strip().lower() in {"1", "true", "yes", "on"}
    expected_narrow = expected_macs_narrow_peak_paths(list(bed_paths.keys()))
    narrow_peaks = {}
    bed_paths_to_call = {}
    macs_status_rows = []
    for label, bed_path in bed_paths.items():
        narrow_path = expected_narrow[label]
        if resume_macs2 and narrow_path.exists() and narrow_path.stat().st_size > 0:
            narrow_peaks[label] = read_narrow_peak_as_pyranges(narrow_path)
            status = "reuse_existing"
        else:
            bed_paths_to_call[label] = bed_path
            status = "call_macs2"
        macs_status_rows.append(
            {"label": label, "bed_path": bed_path, "narrow_peak": str(narrow_path), "status": status}
        )
    pd.DataFrame(macs_status_rows).to_csv(
        RESULTS / "qc" / "macs2_resume_plan.tsv",
        sep="\t",
        index=False,
    )
    if bed_paths_to_call:
        called_peaks = peak_calling(
            macs_path=macs_path,
            bed_paths=bed_paths_to_call,
            outdir=str(WORK / "macs2"),
            genome_size=params.get("genome_size", "mm") or "mm",
            n_cpu=macs_n_cpu,
            input_format="BEDPE",
            shift=73,
            ext_size=146,
            keep_dup="all",
            q_value=0.05,
            nolambda=True,
            **(ray_kwargs if macs_n_cpu > 1 else {}),
        )
        narrow_peaks.update(called_peaks)
    consensus = get_consensus_peaks(
        narrow_peaks,
        peak_half_width=as_int(params, "peak_half_width", 250),
        chromsizes=chromsizes,
        path_to_blacklist=params.get("blacklist", "") or None,
    )
    consensus_bed = WORK / "consensus_peaks.bed"
    consensus_df = consensus.df.loc[:, ["Chromosome", "Start", "End"]].drop_duplicates()
    consensus_df = consensus_df.sort_values(["Chromosome", "Start", "End"])
    consensus_df.to_csv(consensus_bed, sep="\t", header=False, index=False)
    consensus_df.to_csv(RESULTS / "consensus_peaks.tsv", sep="\t", index=False)
    export_consensus_peak_qc_pdf(consensus_df)
    return consensus_bed


def build_cistopic_object(
    sample_sheet: pd.DataFrame,
    cell_meta: pd.DataFrame,
    atac_qc: pd.DataFrame,
    params: Mapping[str, str],
    consensus_bed: Path,
):
    max_memory_gb, _, _ = memory_budget_gb(params)
    n_cpu = resolve_n_cpu(params, max_memory_gb)
    min_frag = as_int(params, "min_frag", 1000)
    min_cell = as_int(params, "min_cell", 1)
    is_acc = as_int(params, "is_acc", 1)
    fragments = fragment_paths_by_sample(sample_sheet)

    objects = []
    for sample_id in sample_sheet["sample_id"].astype(str):
        sample_cells = cell_meta.loc[cell_meta["sample_id"] == sample_id].copy()
        if sample_cells.empty:
            raise ValueError(f"No cells in cell_metadata.tsv for sample {sample_id}")
        obj = create_cistopic_object_from_fragments(
            path_to_fragments=fragments[sample_id],
            path_to_regions=str(consensus_bed),
            path_to_blacklist=params.get("blacklist", "") or None,
            valid_bc=sample_cells["barcode"].tolist(),
            n_cpu=n_cpu,
            min_frag=min_frag,
            min_cell=min_cell,
            is_acc=is_acc,
            project=sample_id,
            split_pattern="-",
            use_polars=True,
        )
        objects.append(obj)

    cistopic_obj = objects[0] if len(objects) == 1 else merge(objects, project="cisTopic_merge", split_pattern="-")
    cistopic_obj.add_cell_data(cell_meta.set_index("cell_id", drop=False), split_pattern="-")

    qc = cistopic_obj.cell_data.copy()
    if "cell_id" in qc.columns:
        qc = qc.drop(columns=["cell_id"])
    qc.index.name = "cell_id"
    qc = qc.reset_index()
    threshold_map = atac_qc.set_index("sample_id").to_dict(orient="index")
    qc["pass_cistopic_fragment_threshold"] = [
        float(row.cisTopic_nr_frag) >= float(threshold_map[str(row.sample_id)]["min_cistopic_fragments"])
        for row in qc.itertuples(index=False)
    ]
    qc["pass_cistopic_accessible_region_threshold"] = [
        float(row.cisTopic_nr_acc) >= float(threshold_map[str(row.sample_id)]["min_cistopic_accessible_regions"])
        for row in qc.itertuples(index=False)
    ]
    qc["pass_atac_qc"] = qc["pass_cistopic_fragment_threshold"] & qc["pass_cistopic_accessible_region_threshold"]
    qc.to_csv(RESULTS / "qc" / "cell_qc_metrics.tsv", sep="\t", index=False)
    export_cell_qc_pdf(qc)

    with open(WORK / "cistopic_obj_pre_models.pkl", "wb") as handle:
        pickle.dump(cistopic_obj, handle)
    return cistopic_obj


def write_model_metrics(models: Sequence[object]) -> None:
    rows = []
    for model in models:
        n_topic = getattr(model, "n_topic", getattr(model, "n_topics", None))
        metrics = getattr(model, "metrics", {})
        if isinstance(metrics, Mapping):
            for metric, value in metrics.items():
                rows.append({"n_topics": n_topic, "metric": metric, "value": value})
        else:
            df = pd.DataFrame(metrics).copy()
            df["n_topics"] = n_topic
            rows.extend(df.to_dict(orient="records"))
    pd.DataFrame(rows).to_csv(
        RESULTS / "model_selection" / "topic_model_metrics.tsv",
        sep="\t",
        index=False,
    )


def run_topic_models(cistopic_obj, topic_grid: pd.DataFrame, params: Mapping[str, str]):
    n_topics = topic_grid["n_topics"].astype(int).tolist()
    max_memory_gb, _, physical_memory_gb = memory_budget_gb(params)
    n_cpu = resolve_n_cpu(params, max_memory_gb)
    topic_n_cpu = resolve_topic_n_cpu(params, n_topics, max_memory_gb)
    backend = params.get("lda_backend", "cgs").strip().lower()
    if backend == "cgs":
        models = run_cgs_models(
            cistopic_obj,
            n_topics=n_topics,
            n_cpu=topic_n_cpu,
            n_iter=as_int(params, "n_iter", 150),
            random_state=as_int(params, "random_state", 555),
            save_path=str(WORK / "models"),
        )
    elif backend == "mallet":
        mallet_value = params.get("mallet_path", "mallet") or "mallet"
        mallet_path = shutil.which(mallet_value) or (mallet_value if Path(mallet_value).exists() else None)
        if mallet_path is None:
            raise FileNotFoundError(
                "lda_backend=mallet was requested, but the MALLET executable was not found. "
                "Install MALLET in the environment or set pycistopic.mallet_path to the executable path. "
                "Use pycistopic.lda_backend=cgs for the built-in backend."
            )
        mallet_tmp = Path(params.get("mallet_tmp_dir", "tmp/mallet") or "tmp/mallet")
        if not mallet_tmp.is_absolute():
            mallet_tmp = PROJECT / mallet_tmp
        mallet_tmp.mkdir(parents=True, exist_ok=True)
        reuse_corpus = params.get("reuse_mallet_corpus", "1").strip().lower() in {"1", "true", "yes", "on"}
        mallet_memory = resolve_mallet_memory(params, max_memory_gb, physical_memory_gb)
        old_mallet_memory = os.environ.get("MALLET_MEMORY")
        os.environ["MALLET_MEMORY"] = mallet_memory
        (RESULTS / "model_selection" / "mallet_runtime.tsv").write_text(
            "parameter\tvalue\n"
            f"mallet_path\t{mallet_path}\n"
            f"mallet_memory\t{mallet_memory}\n"
            f"mallet_tmp_dir\t{mallet_tmp}\n"
            f"reuse_mallet_corpus\t{int(reuse_corpus)}\n"
            f"topic_n_cpu\t{topic_n_cpu}\n"
        )
        try:
            models = run_cgs_models_mallet(
                path_to_mallet_binary=str(mallet_path),
                cistopic_obj=cistopic_obj,
                n_topics=n_topics,
                n_cpu=topic_n_cpu,
                n_iter=as_int(params, "n_iter", 150),
                random_state=as_int(params, "random_state", 555),
                tmp_path=str(mallet_tmp),
                save_path=str(WORK / "models"),
                reuse_corpus=reuse_corpus,
            )
        finally:
            if old_mallet_memory is None:
                os.environ.pop("MALLET_MEMORY", None)
            else:
                os.environ["MALLET_MEMORY"] = old_mallet_memory
    else:
        raise ValueError("lda_backend must be 'cgs' or 'mallet'.")
    selected_value = params.get("selected_n_topics", "auto").strip().lower()
    selected_n_topics = None if selected_value in {"", "auto"} else int(selected_value)
    selected_model = evaluate_models(
        models,
        select_model=selected_n_topics,
        return_model=True,
        plot=True,
        save=str(RESULTS / "model_selection" / "topic_model_metrics.pdf"),
    )
    cistopic_obj.add_LDA_model(selected_model)
    write_model_metrics(models)
    export_topic_qc_pdf(cistopic_obj)
    selected_n = getattr(selected_model, "n_topic", getattr(selected_model, "n_topics", "unknown"))
    (RESULTS / "model_selection" / "selected_model.txt").write_text(f"selected_n_topics\t{selected_n}\n")
    return cistopic_obj


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, coords = str(region).split(":", 1)
    start, end = coords.replace("-", ":").split(":")[:2]
    return chrom, int(start), int(end)


def write_region_set(region_names: Sequence[str], out_file: Path) -> None:
    rows = [parse_region(region) for region in region_names]
    if rows:
        df = pd.DataFrame(rows).drop_duplicates().sort_values([0, 1, 2])
    else:
        df = pd.DataFrame(columns=[0, 1, 2])
    df.to_csv(out_file, sep="\t", header=False, index=False)


def export_topic_region_sets(cistopic_obj, params: Mapping[str, str]) -> None:
    otsu_sets = binarize_topics(
        cistopic_obj,
        target="region",
        method="otsu",
        plot=True,
        save=str(RESULTS / "model_selection" / "topic_otsu_thresholds.pdf"),
    )
    top_sets = binarize_topics(
        cistopic_obj,
        target="region",
        method="ntop",
        ntop=as_int(params, "ntop_regions", 3000),
        plot=False,
    )
    summary = []
    for topic, regions in otsu_sets.items():
        out_file = INPUTS / "region_sets" / "Topics_otsu" / f"{topic}.bed"
        write_region_set(list(regions.index), out_file)
        summary.append({"region_set_type": "Topics_otsu", "name": str(topic), "n_regions": int(len(regions.index)), "bed": str(out_file)})
    for topic, regions in top_sets.items():
        out_file = INPUTS / "region_sets" / "Topics_top_3k" / f"{topic}.bed"
        write_region_set(list(regions.index), out_file)
        summary.append({"region_set_type": "Topics_top_3k", "name": str(topic), "n_regions": int(len(regions.index)), "bed": str(out_file)})
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULTS / "model_selection" / "topic_region_set_summary.tsv", sep="\t", index=False)
    pdf_path = RESULTS / "model_selection" / "topic_region_set_summary.pdf"
    with PdfPages(pdf_path) as pdf:
        for region_set_type, df in summary_df.groupby("region_set_type", sort=True):
            fig, ax = plt.subplots(figsize=(max(6, 0.25 * len(df)), 4))
            ax.bar(df["name"].astype(str), df["n_regions"].astype(int), color="#76B7B2", alpha=0.9)
            ax.set_ylabel("regions")
            ax.set_title(region_set_type)
            ax.tick_params(axis="x", rotation=90, labelsize=7)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    require_nonempty_pdf(pdf_path)


def export_dar_region_sets(cistopic_obj, params: Mapping[str, str]) -> None:
    levels = sorted(cistopic_obj.cell_data["cell_label"].dropna().astype(str).unique())
    if len(levels) < 2:
        raise ValueError("DAR calling requires at least two cell_label groups")
    contrasts = [[[label], [x for x in levels if x != label]] for label in levels]
    pycistopic_keys = ["_".join(x[0]) + "_VS_" + "_".join(x[1]) for x in contrasts]
    output_names = [f"{label}_VS_rest" for label in levels]

    imputed = impute_accessibility(cistopic_obj)
    dar = find_diff_features(
        cistopic_obj,
        imputed,
        variable="cell_label",
        contrasts=contrasts,
        adjpval_thr=as_float(params, "dar_adjpval_thr", 0.05),
        log2fc_thr=as_float(params, "dar_log2fc_thr", 0.5),
        n_cpu=resolve_n_cpu(params, memory_budget_gb(params)[0]),
        split_pattern="-",
    )
    if not isinstance(dar, Mapping):
        raise TypeError("find_diff_features did not return a contrast-to-table mapping")

    summary = []
    for key, out_name in zip(pycistopic_keys, output_names):
        if key not in dar:
            raise KeyError(f"DAR result missing expected pycisTopic contrast: {key}")
        table = dar[key]
        out_tsv = RESULTS / "dar" / f"{out_name}.tsv"
        out_bed = INPUTS / "region_sets" / "DARs_cell_label" / f"{out_name}.bed"
        table.to_csv(out_tsv, sep="\t")
        write_region_set(list(table.index.astype(str)), out_bed)
        summary.append({
            "contrast": out_name,
            "pycistopic_contrast": key,
            "n_dar": int(table.shape[0]),
            "bed": str(out_bed),
            "table": str(out_tsv),
        })
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULTS / "dar" / "dar_summary.tsv", sep="\t", index=False)
    export_dar_summary_pdf(summary_df)


def write_manifest(
    sample_sheet: pd.DataFrame,
    cell_meta: pd.DataFrame,
    topic_grid: pd.DataFrame,
    params: Mapping[str, str],
    consensus_bed: Path,
) -> None:
    region_set_files = sorted(str(p) for p in (INPUTS / "region_sets").glob("**/*.bed"))
    manifest = {
        "sample_sheet": "inputs/sample_sheet.tsv",
        "cell_metadata": "inputs/cell_metadata.tsv",
        "params": dict(params),
        "samples": sample_sheet["sample_id"].astype(str).tolist(),
        "n_cells_input": int(cell_meta.shape[0]),
        "topic_grid": topic_grid["n_topics"].astype(int).tolist(),
        "consensus_peaks": str(consensus_bed),
        "region_sets": region_set_files,
        "required_outputs": [
            "inputs/cistopic_obj.pkl",
            "work/pycistopic/consensus_peaks.bed",
            "inputs/region_sets/Topics_otsu/*.bed",
            "inputs/region_sets/Topics_top_3k/*.bed",
            "inputs/region_sets/DARs_cell_label/*.bed",
            "results/pycistopic/qc/cell_qc_metrics.tsv",
            "results/pycistopic/qc/cistopic_cell_qc.pdf",
            "results/pycistopic/qc/consensus_peak_qc.pdf",
            "results/pycistopic/dar/dar_summary.tsv",
            "results/pycistopic/dar/dar_summary.pdf",
            "results/pycistopic/model_selection/topic_model_metrics.tsv",
            "results/pycistopic/model_selection/topic_model_metrics.pdf",
            "results/pycistopic/model_selection/topic_qc_metrics.tsv",
            "results/pycistopic/model_selection/topic_qc_metrics.pdf",
            "results/pycistopic/model_selection/topic_otsu_thresholds.pdf",
            "results/pycistopic/model_selection/topic_region_set_summary.tsv",
            "results/pycistopic/model_selection/topic_region_set_summary.pdf",
            "results/pycistopic/model_selection/selected_model.txt",
        ],
    }
    (RESULTS / "pycistopic_manifest.json").write_text(json.dumps(manifest, indent=2))


def final_checks() -> None:
    required_files = [
        INPUTS / "cistopic_obj.pkl",
        WORK / "consensus_peaks.bed",
        RESULTS / "qc" / "cell_qc_metrics.tsv",
        RESULTS / "qc" / "cistopic_cell_qc.pdf",
        RESULTS / "qc" / "consensus_peak_qc.pdf",
        RESULTS / "dar" / "dar_summary.tsv",
        RESULTS / "dar" / "dar_summary.pdf",
        RESULTS / "model_selection" / "topic_model_metrics.tsv",
        RESULTS / "model_selection" / "topic_model_metrics.pdf",
        RESULTS / "model_selection" / "topic_qc_metrics.tsv",
        RESULTS / "model_selection" / "topic_qc_metrics.pdf",
        RESULTS / "model_selection" / "topic_otsu_thresholds.pdf",
        RESULTS / "model_selection" / "topic_region_set_summary.tsv",
        RESULTS / "model_selection" / "topic_region_set_summary.pdf",
        RESULTS / "model_selection" / "selected_model.txt",
        RESULTS / "pycistopic_manifest.json",
    ]
    missing = [str(p) for p in required_files if not p.exists() or p.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing required pycisTopic outputs: " + ", ".join(missing))
    required_dirs = [
        INPUTS / "region_sets" / "Topics_otsu",
        INPUTS / "region_sets" / "Topics_top_3k",
        INPUTS / "region_sets" / "DARs_cell_label",
    ]
    for d in required_dirs:
        if not list(d.glob("*.bed")):
            raise FileNotFoundError(f"No BED files written in {d}")


def run_workflow() -> None:
    mkdirs()
    with tqdm(total=8, desc="pycisTopic workflow", unit="stage") as progress:
        sample_sheet, cell_meta, topic_grid, atac_qc, params = load_inputs()
        write_cell_accounting(cell_meta)
        progress.update()
        consensus_bed = run_pseudobulk_and_peak_calling(sample_sheet, cell_meta, params)
        progress.update()
        cistopic_obj = build_cistopic_object(sample_sheet, cell_meta, atac_qc, params, consensus_bed)
        progress.update()
        cistopic_obj = run_topic_models(cistopic_obj, topic_grid, params)
        progress.update()
        export_topic_region_sets(cistopic_obj, params)
        progress.update()
        export_dar_region_sets(cistopic_obj, params)
        progress.update()
        with open(INPUTS / "cistopic_obj.pkl", "wb") as handle:
            pickle.dump(cistopic_obj, handle)
        write_manifest(sample_sheet, cell_meta, topic_grid, params, consensus_bed)
        progress.update()
        final_checks()
        progress.update()
    print("pycisTopic workflow completed")
    print(f"cells_in_metadata\t{cell_meta.shape[0]}")
    print(f"consensus_peaks\t{sum(1 for _ in open(consensus_bed))}")


def main() -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    log_path = LOGS / f"run_pycistopic_workflow_{datetime.now():%Y%m%d_%H%M%S}.log"
    with log_path.open("w") as log_handle:
        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = Tee(old_stdout, log_handle)
        sys.stderr = Tee(old_stderr, log_handle)
        try:
            print(f"pycisTopic workflow log\t{log_path}")
            run_workflow()
        except Exception:
            traceback.print_exc()
            raise SystemExit(1)
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


if __name__ == "__main__":
    main()
