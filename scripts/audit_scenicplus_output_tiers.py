#!/usr/bin/env python
"""Audit SCENIC+/pycisTopic outputs by scientific evidence tier."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(PROJECT / "tmp" / "xdg_cache"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT / "tmp" / "numba_cache"))
for _path in ["MPLCONFIGDIR", "XDG_CACHE_HOME", "NUMBA_CACHE_DIR"]:
    Path(os.environ[_path]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", default=str(PROJECT))
    parser.add_argument("--outdir", default="results/scenicplus_output_tiers")
    return parser.parse_args()


def resolve(path: str, project: Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else project / p


def exists_nonempty(path: Path) -> bool:
    if any(ch in str(path) for ch in ["*", "?", "["]):
        return any(p.is_file() and p.stat().st_size > 0 for p in path.parent.glob(path.name))
    return path.is_file() and path.stat().st_size > 0


def dir_has_files(path: Path, pattern: str = "*") -> bool:
    return path.is_dir() and any(p.is_file() and p.stat().st_size > 0 for p in path.glob(pattern))


def tier_rows(project: Path) -> list[dict[str, str]]:
    checks: list[dict[str, object]] = [
        {
            "tier": "0_input_qc",
            "scientific_question": "Are the active RNA/ATAC cells, metadata and ATAC QC reliable enough for GRN inference?",
            "required": [
                "inputs/gex.h5ad",
                "inputs/cell_metadata.tsv",
                "inputs/sample_sheet.tsv",
                "results/pycistopic/qc/cistopic_cell_qc.pdf",
            ],
            "recommended_outputs": "cell/fragment QC PDF and sample-by-label cell accounting; doublet filtering belongs to the upstream single-cell QC object",
        },
        {
            "tier": "1_chromatin_topics",
            "scientific_question": "What are the major chromatin-accessibility states in the active cells?",
            "required": [
                "inputs/cistopic_obj.pkl",
                "results/pycistopic/model_selection/topic_model_metrics.pdf",
                "results/pycistopic/model_selection/topic_qc_metrics.pdf",
                "results/pycistopic/model_selection/selected_model.txt",
            ],
            "recommended_outputs": "topic model selection, topic QC, topic-by-cell-state heatmap, topic UMAP overlays",
        },
        {
            "tier": "2_region_sets_and_DARs",
            "scientific_question": "Which accessible-region programs and differential regions define cell states or perturbations?",
            "required": [
                "work/pycistopic/consensus_peaks.bed",
                "inputs/region_sets/Topics_otsu",
                "inputs/region_sets/Topics_top_3k",
                "inputs/region_sets/DARs_cell_label",
                "results/pycistopic/dar/dar_summary.tsv",
                "results/pycistopic/dar/dar_summary.pdf",
            ],
            "recommended_outputs": "topic region BEDs, DAR BEDs/tables, DAR summary PDF",
        },
        {
            "tier": "3_motif_cistarget",
            "scientific_question": "Which TF motifs are enriched in the project-specific region universe?",
            "required": [
                "inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather",
                "inputs/cistarget_db/custom.regions_vs_motifs.scores.feather",
                "inputs/cistarget_db/motif_annotations.tbl",
            ],
            "recommended_outputs": "custom cisTarget rankings/scores database and motif2TF annotation manifest",
        },
        {
            "tier": "4_eregulon_activity",
            "scientific_question": "Which TF-region-gene eRegulons are active in each cell state?",
            "required": [
                "results/scenicplus/AUCell_direct.h5mu",
                "results/scenicplus/eRegulons_direct.tsv",
                "work/scenicplus/tf_to_gene_adj.tsv",
                "work/scenicplus/region_to_gene_adj.tsv",
                "results/scenicplus_figures/01_eregulon_auc_heatmap_direct.pdf",
                "results/scenicplus_figures/01_eregulon_auc_heatmap_condition_direct.pdf",
                "results/scenicplus_figures/03_eregulon_dot_heatmap_direct.pdf",
                "results/scenicplus_figures/03_eregulon_dot_heatmap_condition_direct.pdf",
            ],
            "recommended_outputs": "AUC heatmap, RSS/specificity, dot heatmap, selected eRegulon UMAP overlays",
        },
        {
            "tier": "5_condition_effects",
            "scientific_question": "Which expressed eRegulons change between experimental conditions at sample level?",
            "required": [
                "results/scenicplus_figures/sample_mean_auc_direct.tsv",
                "results/scenicplus_figures/condition_eregulon_auc_statistics_direct.tsv",
                "results/scenicplus_figures/by_cell_label_condition_eregulon_auc_statistics_direct.tsv",
                "results/scenicplus_figures/09_condition_sample_counts_direct.pdf",
                "results/scenicplus_figures/12_condition_by_cell_label_effect_heatmap_direct.pdf",
            ],
            "recommended_outputs": "overall and cell-label-stratified sample-level AUC statistics with one-panel vector PDFs",
        },
        {
            "tier": "6_mechanistic_views",
            "scientific_question": "Which candidate TF-target mechanisms deserve locus-level or network-level presentation?",
            "required": [
                "results/scenicplus_figures/08_tf_target_network_direct.pdf",
                "results/scenicplus_figures/source_tf_target_network_edges_direct.tsv",
            ],
            "recommended_outputs": "focused TF-target network, locus coverage/arcs when pseudobulk bigWigs and interactions are available",
        },
    ]
    rows: list[dict[str, str]] = []
    for check in checks:
        missing: list[str] = []
        present: list[str] = []
        for item in check["required"]:  # type: ignore[index]
            path = resolve(str(item), project)
            ok = dir_has_files(path, "*.bed") or exists_nonempty(path)
            (present if ok else missing).append(str(item))
        rows.append(
            {
                "tier": str(check["tier"]),
                "scientific_question": str(check["scientific_question"]),
                "status": "ready" if not missing else "missing_inputs",
                "n_required": str(len(check["required"])),  # type: ignore[arg-type]
                "n_present": str(len(present)),
                "missing": ";".join(missing),
                "recommended_outputs": str(check["recommended_outputs"]),
            }
        )
    return rows


def write_pdf(df: pd.DataFrame, out_pdf: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.6,
        }
    )
    with PdfPages(out_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(8.2, max(3.5, 0.55 * len(df))))
        ax.axis("off")
        display = df[["tier", "status", "n_present", "n_required"]].copy()
        display["missing_summary"] = df["missing"].map(lambda x: "none" if not x else f"{len(x.split(';'))} missing")
        table = ax.table(
            cellText=display.values,
            colLabels=display.columns,
            cellLoc="left",
            loc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1, 1.35)
        for (row, col), cell in table.get_celld().items():
            cell.set_linewidth(0.25)
            if row == 0:
                cell.set_facecolor("#EAEFF2")
                cell.set_text_props(weight="bold")
            elif display.iloc[row - 1]["status"] == "ready":
                cell.set_facecolor("#EEF5EF")
            else:
                cell.set_facecolor("#F7EFEA")
        ax.set_title("SCENIC+ output tier audit", fontsize=9, pad=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    project = Path(args.project_dir).expanduser().resolve()
    outdir = resolve(args.outdir, project)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(tier_rows(project))
    tsv = outdir / "scenicplus_output_tier_audit.tsv"
    pdf = outdir / "scenicplus_output_tier_audit.pdf"
    df.to_csv(tsv, sep="\t", index=False)
    write_pdf(df, pdf)
    print(f"WROTE {tsv}")
    print(f"WROTE {pdf}")
    if (df["status"] != "ready").any():
        print("Some scientific output tiers are not ready; inspect the missing column in the TSV.")


if __name__ == "__main__":
    main()
