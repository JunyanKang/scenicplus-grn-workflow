"""Reproducible ATAC doublet calls from a cisTopic peak-by-cell matrix."""
from __future__ import annotations

from pathlib import Path
import os
from typing import Mapping, Sequence

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scrublet as scr


def _as_float(params: Mapping[str, str], key: str, default: float) -> float:
    value = params.get(key, "")
    return default if value == "" else float(value)


def _as_int(params: Mapping[str, str], key: str, default: int) -> int:
    value = params.get(key, "")
    return default if value == "" else int(value)


def run_scrublet_atac_doublets(
    count_matrix,
    cell_names: Sequence[str],
    cell_metadata: pd.DataFrame,
    params: Mapping[str, str],
    out_tsv: str | Path,
    out_pdf: str | Path,
) -> pd.DataFrame:
    """Run Scrublet on cells x regions ATAC counts and write calls plus PDF.

    Parameters
    ----------
    count_matrix
        Sparse or dense cells x regions count matrix.
    cell_names
        Cell IDs in the same order as rows of count_matrix.
    cell_metadata
        Metadata indexed by cell ID; sample_id/cell_label are copied when present.
    params
        String-valued parameter mapping. Uses expected_doublet_rate,
        doublet_n_prin_comps, doublet_min_counts, doublet_min_cells, random_state.
    out_tsv, out_pdf
        Required output files.
    """
    out_tsv = Path(out_tsv)
    out_pdf = Path(out_pdf)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    scrub = scr.Scrublet(
        count_matrix,
        expected_doublet_rate=_as_float(params, "expected_doublet_rate", 0.1),
        random_state=_as_int(params, "random_state", 555),
    )
    scores, predicted = scrub.scrub_doublets(
        min_counts=_as_int(params, "doublet_min_counts", 2),
        min_cells=_as_int(params, "doublet_min_cells", 3),
        n_prin_comps=_as_int(params, "doublet_n_prin_comps", 30),
        verbose=True,
    )

    doublets = pd.DataFrame({
        "cell_id": list(cell_names),
        "doublet_score": scores,
        "predicted_doublet": predicted.astype(bool),
    }).set_index("cell_id")
    metadata_cols = [c for c in ["sample_id", "cell_label"] if c in cell_metadata.columns]
    if metadata_cols:
        doublets = doublets.join(cell_metadata.loc[:, metadata_cols], how="left")
    doublets.reset_index().to_csv(out_tsv, sep="\t", index=False)

    with PdfPages(out_pdf) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].hist(scores, bins=80, color="#E15759", alpha=0.75, label="observed cells")
        sim_scores = getattr(scrub, "doublet_scores_sim_", None)
        if sim_scores is not None:
            axes[0].hist(sim_scores, bins=80, color="#BAB0AC", alpha=0.55, label="simulated doublets")
        threshold = getattr(scrub, "threshold_", None)
        if threshold is not None:
            axes[0].axvline(threshold, color="black", lw=1, ls="--", label=f"threshold={threshold:.3f}")
        axes[0].set_xlabel("Scrublet doublet score")
        axes[0].set_ylabel("cells")
        axes[0].set_title("ATAC doublet score distribution")
        axes[0].legend(frameon=False, fontsize=8)

        if "sample_id" in doublets.columns:
            by_sample = doublets.groupby("sample_id")["predicted_doublet"].mean().sort_index()
        else:
            by_sample = pd.Series(dtype=float)
        axes[1].bar(by_sample.index.astype(str), by_sample.values, color="#B07AA1", alpha=0.9)
        y_max = max(0.05, float(by_sample.max()) * 1.2 if len(by_sample) else 0.05)
        axes[1].set_ylim(0, y_max)
        axes[1].set_ylabel("predicted doublet fraction")
        axes[1].set_title("Per-sample doublet rate")
        axes[1].tick_params(axis="x", rotation=45)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    if not out_pdf.exists() or out_pdf.stat().st_size == 0:
        raise FileNotFoundError(f"doublet PDF was not created: {out_pdf}")
    return doublets
