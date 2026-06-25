#!/usr/bin/env python
"""SCENIC+ result figures from AUCell and eRegulon tables."""
from __future__ import annotations

import argparse
import ast
import os
import re
import textwrap
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import anndata as ad
import mudata as md
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.spatial.distance import jensenshannon


GROUP_COLORS = {}
FALLBACK_COLORS = ["#2F4858", "#B57457", "#667A9A", "#9A8067", "#7D9477", "#8E5A63"]
DIV_CMAP = sns.diverging_palette(230, 20, s=70, l=55, center="light", as_cmap=True)
SEQ_CMAP = sns.light_palette("#2F4858", as_cmap=True)
OVERLAP_CMAP = sns.light_palette("#8E5A63", as_cmap=True)
UMAP_CMAP = sns.light_palette("#B57457", as_cmap=True)

STYLE_DEFAULTS = {
    "font_family": "Arial,Helvetica,DejaVu Sans,sans-serif",
    "font_size": "7",
    "axes_label_size": "6.5",
    "axes_title_size": "7",
    "tick_label_size": "5.2",
    "legend_font_size": "5.5",
    "figure_legend_font_size": "5.0",
    "axes_line_width": "0.7",
    "heatmap_cbar_shrink": "0.52",
    "heatmap_left_margin": "0.34",
    "heatmap_right_margin": "0.88",
    "heatmap_bottom_margin": "0.20",
    "heatmap_top_margin": "0.88",
    "dot_min_size": "8",
    "dot_size_scale": "95",
    "umap_point_size": "2.0",
    "umap_point_alpha": "0.88",
    "embedding_point_size": "3.0",
    "embedding_point_alpha": "0.75",
    "max_panels_per_row": "3",
    "condition_background_alpha": "0.18",
    "condition_background_point_size": "1.2",
    "network_tf_node_size": "150",
    "network_target_node_size": "18",
    "network_edge_alpha": "0.24",
    "color_fallback_1": "#2F4858",
    "color_fallback_2": "#B57457",
    "color_fallback_3": "#667A9A",
    "color_fallback_4": "#9A8067",
    "color_fallback_5": "#7D9477",
    "color_fallback_6": "#8E5A63",
}
STYLE = STYLE_DEFAULTS.copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--auc-h5mu", default="results/scenicplus/AUCell_direct.h5mu")
    parser.add_argument("--metadata", default="inputs/cell_metadata.tsv")
    parser.add_argument("--tf-to-gene", default="work/scenicplus/tf_to_gene_adj.tsv")
    parser.add_argument("--region-to-gene", default="work/scenicplus/region_to_gene_adj.tsv")
    parser.add_argument("--eregulons", default="results/scenicplus/eRegulons_direct.tsv")
    parser.add_argument("--group-col", default="cell_label")
    parser.add_argument("--condition-col", default="condition")
    parser.add_argument("--cell-col", default="cell_id")
    parser.add_argument("--umap-x", default="umap_1")
    parser.add_argument("--umap-y", default="umap_2")
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--umap-n", type=int, default=12)
    parser.add_argument("--network-top-tfs", type=int, default=8)
    parser.add_argument("--network-targets-per-tf", type=int, default=12)
    parser.add_argument("--outdir", default="results/scenicplus_figures")
    parser.add_argument("--file-suffix", default="")
    parser.add_argument("--plot-style-config", default="results/scenicplus_figures/plot_style_parameters.tsv")
    return parser.parse_args()


def project_dir() -> Path:
    return PROJECT


def resolve_project_path(path_value: str, base: Path) -> str:
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = (base / path).resolve()
    return str(path)


def resolve_args_paths(args: argparse.Namespace) -> argparse.Namespace:
    base = project_dir()
    for attr in ["auc_h5mu", "metadata", "tf_to_gene", "region_to_gene", "eregulons", "outdir", "plot_style_config"]:
        setattr(args, attr, resolve_project_path(getattr(args, attr), base))
    return args


def style_float(name: str, default: float) -> float:
    try:
        return float(STYLE.get(name, default))
    except Exception:
        return default


def style_int(name: str, default: int) -> int:
    try:
        return int(float(STYLE.get(name, default)))
    except Exception:
        return default


def load_or_write_style_config(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if {"parameter", "value"}.issubset(df.columns):
            STYLE.update({str(row.parameter): str(row.value) for row in df.itertuples(index=False)})
        return
    rows = [
        {
            "parameter": key,
            "value": value,
            "scope": "all_figures",
            "description": "Edit value and rerun postprocess figures/stats to redraw with updated styling.",
        }
        for key, value in STYLE_DEFAULTS.items()
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def group_color(label: object, idx: int = 0) -> str:
    key = f"color_{label}"
    if key in STYLE:
        return STYLE[key]
    if str(label) in GROUP_COLORS:
        return GROUP_COLORS[str(label)]
    fallback = [STYLE.get(f"color_fallback_{i}", FALLBACK_COLORS[(i - 1) % len(FALLBACK_COLORS)]) for i in range(1, 7)]
    return fallback[idx % len(fallback)]


def set_publication_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [x.strip() for x in STYLE.get("font_family", STYLE_DEFAULTS["font_family"]).split(",")],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": style_float("axes_line_width", 0.7),
            "axes.labelsize": style_float("axes_label_size", 6.5),
            "axes.titlesize": style_float("axes_title_size", 7),
            "xtick.labelsize": style_float("tick_label_size", 5.2),
            "ytick.labelsize": style_float("tick_label_size", 5.2),
            "legend.fontsize": style_float("legend_font_size", 5.5),
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.transparent": False,
            "legend.frameon": False,
        }
    )
    sns.set_theme(style="white", rc={"axes.linewidth": style_float("axes_line_width", 0.7)})


def require_pdf(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"required PDF was not created: {path}")


def clean_regulon_name(name: object) -> str:
    return re.sub(r"\s+", "_", str(name))


def tf_from_regulon(name: object) -> str:
    text = clean_regulon_name(name)
    text = re.sub(r"_(direct|extended).*$", "", text)
    text = re.sub(r"_[+-]/[+-].*$", "", text)
    return re.split(r"[_()]", text)[0]


def regulon_display_label(name: object) -> str:
    text = clean_regulon_name(name)
    tf = tf_from_regulon(text)
    genes = re.search(r"\((\d+)g\)", text)
    return f"{tf} ({genes.group(1)} targets)" if genes else tf


def display_index(values: pd.Index | list[str]) -> list[str]:
    return [regulon_display_label(x) for x in values]


def short_text(value: object, width: int = 18) -> str:
    return textwrap.shorten(str(value), width=width, placeholder="...")


def save_pdf(fig: plt.Figure, pdf: Path) -> None:
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    require_pdf(pdf)


def named(stem: str, suffix: str, ext: str = ".pdf") -> str:
    return f"{stem}_{suffix}{ext}" if suffix else f"{stem}{ext}"


def add_figure_legend(fig: plt.Figure, text: str) -> None:
    wrapped = "\n".join(textwrap.wrap(text, width=135))
    fig.text(0.02, 0.015, wrapped, ha="left", va="bottom", fontsize=style_float("figure_legend_font_size", 5.0))


def ordered_unique(values: pd.Series) -> list[str]:
    present = [str(x) for x in pd.unique(values.astype(str)) if str(x)]
    return present


def add_joint_group(meta: pd.DataFrame, group_col: str, condition_col: str) -> tuple[pd.Series, list[str]]:
    if not condition_col or condition_col not in meta.columns:
        return meta[group_col].astype(str), ordered_unique(meta[group_col])
    group_order = ordered_unique(meta[group_col])
    condition_order = ordered_unique(meta[condition_col])
    joint = meta[group_col].astype(str) + "_" + meta[condition_col].astype(str)
    order = [f"{group}_{condition}" for group in group_order for condition in condition_order if (joint == f"{group}_{condition}").any()]
    return joint, order


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.astype(float)
    mean = values.mean(axis=1)
    sd = values.std(axis=1).replace(0, np.nan)
    return values.sub(mean, axis=0).div(sd, axis=0).fillna(0)


def cluster_order(values: pd.DataFrame, axis: int) -> list:
    matrix = values.astype(float).replace([np.inf, -np.inf], np.nan).fillna(0)
    if matrix.shape[axis] < 3:
        return list(matrix.index if axis == 0 else matrix.columns)
    data = matrix.values if axis == 0 else matrix.values.T
    if np.all(np.nanstd(data, axis=1) == 0):
        return list(matrix.index if axis == 0 else matrix.columns)
    return list((matrix.index if axis == 0 else matrix.columns)[leaves_list(linkage(pdist(data), method="average"))])


def cluster_frame(values: pd.DataFrame, rows: bool = True, cols: bool = True) -> pd.DataFrame:
    ordered = values.copy()
    if rows:
        ordered = ordered.loc[cluster_order(ordered, axis=0)]
    if cols:
        ordered = ordered.loc[:, cluster_order(ordered, axis=1)]
    return ordered


def choose_column(columns: pd.Index, candidates: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return str(lower[candidate.lower()])
    return None


def find_auc_mod(mdata: md.MuData, meta_cell_ids: set[str]) -> pd.DataFrame:
    candidates = []
    for name, adata in mdata.mod.items():
        obs_overlap = len(set(map(str, adata.obs_names)) & meta_cell_ids)
        candidates.append((obs_overlap, adata.n_vars, name, adata, False))
        var_overlap = len(set(map(str, adata.var_names)) & meta_cell_ids)
        candidates.append((var_overlap, adata.n_obs, name, adata, True))
    candidates.sort(reverse=True, key=lambda x: (x[0], x[1]))
    overlap, _, name, adata, transpose = candidates[0]
    if overlap == 0:
        raise ValueError("No AUCell modality has cell IDs overlapping metadata.")
    print("selected_auc_mod", name, "transpose", transpose, "overlap", overlap)
    matrix = adata.to_df()
    return matrix.T if transpose else matrix


def read_auc_and_metadata(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    required = [args.cell_col, args.group_col]
    if args.condition_col:
        required.append(args.condition_col)
    for col in required:
        if col not in meta.columns:
            raise ValueError(f"metadata missing required column: {col}")
    meta = meta.drop_duplicates(args.cell_col).set_index(args.cell_col)
    mdata = md.read_h5mu(args.auc_h5mu)
    auc = find_auc_mod(mdata, set(meta.index.astype(str)))
    auc.index = auc.index.astype(str)
    auc.columns = [clean_regulon_name(x) for x in auc.columns]
    overlap = auc.index.intersection(meta.index)
    if len(overlap) == 0:
        raise ValueError("AUCell matrix and metadata have no overlapping cells.")
    auc = auc.loc[overlap].apply(pd.to_numeric, errors="coerce").fillna(0)
    meta = meta.loc[overlap].copy()
    meta[args.group_col] = meta[args.group_col].astype(str)
    if args.condition_col:
        meta[args.condition_col] = meta[args.condition_col].astype(str)
    return auc, meta


def calculate_rss(auc: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    rss = pd.DataFrame(index=sorted(groups.astype(str).unique()), columns=auc.columns, dtype=float)
    values = auc.clip(lower=0).astype(float)
    for regulon in values.columns:
        p_reg = values[regulon].values + 1e-12
        p_reg = p_reg / p_reg.sum()
        for group in rss.index:
            p_group = (groups.values == group).astype(float) + 1e-12
            p_group = p_group / p_group.sum()
            rss.loc[group, regulon] = 1.0 - float(jensenshannon(p_reg, p_group, base=2))
    return rss


def select_top_regulons(mean_by_group: pd.DataFrame, rss: pd.DataFrame, top_n: int) -> list[str]:
    variable = mean_by_group.var(axis=0).sort_values(ascending=False)
    specific = rss.max(axis=0).sort_values(ascending=False)
    ranked = pd.concat([variable.rank(ascending=False), specific.rank(ascending=False)], axis=1)
    ranked.columns = ["variance_rank", "rss_rank"]
    ranked["combined_rank"] = ranked.mean(axis=1)
    return ranked.sort_values("combined_rank").head(top_n).index.tolist()


def heatmap_pdf(
    data: pd.DataFrame,
    pdf: Path,
    cmap,
    center: float | None,
    cbar_label: str,
    title: str,
    figsize: tuple[float, float] = (4.3, 5.7),
    x_rotation: int = 0,
    legend: str = "",
    cluster_cols: bool = True,
) -> None:
    plot_data = cluster_frame(data, rows=True, cols=cluster_cols)
    plot_data.index = display_index(plot_data.index)
    plot_data.columns = [short_text(x, 16) for x in plot_data.columns]
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        plot_data,
        cmap=cmap,
        center=center,
        ax=ax,
        linewidths=0,
        cbar_kws={"label": cbar_label, "shrink": style_float("heatmap_cbar_shrink", 0.52), "pad": 0.03},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, pad=4)
    ax.tick_params(length=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=x_rotation, ha="right" if x_rotation else "center")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if legend:
        add_figure_legend(fig, legend)
    fig.subplots_adjust(
        left=style_float("heatmap_left_margin", 0.34),
        right=style_float("heatmap_right_margin", 0.88),
        bottom=style_float("heatmap_bottom_margin", 0.20),
        top=style_float("heatmap_top_margin", 0.88),
    )
    save_pdf(fig, pdf)


def plot_auc_heatmap(mean_by_group: pd.DataFrame, top: list[str], outdir: Path, suffix: str, condition: bool = False) -> None:
    data = zscore_rows(mean_by_group[top].T)
    data.to_csv(outdir / named("source_auc_heatmap_zscore", suffix, ".tsv"), sep="\t")
    heatmap_pdf(
        data,
        outdir / named("01_eregulon_auc_heatmap", suffix),
        cmap=DIV_CMAP,
        center=0,
        cbar_label="AUC z-score",
        title="eRegulon activity" if not condition else "eRegulon activity by state and condition",
        figsize=(4.4, 5.8),
        legend=(
            "Each row is an eRegulon and each column is a cell-state/condition stratum. Values are row-scaled AUCell means; "
            "columns are ordered by cell label then condition so condition shifts can be compared within the same biological state."
            if condition
            else "Each row is an eRegulon and each column is a cell label. Values are row-scaled AUCell means, emphasizing state-biased activity patterns."
        ),
        cluster_cols=not condition,
    )


def plot_rss_heatmap(rss: pd.DataFrame, top: list[str], outdir: Path, suffix: str, condition: bool = False) -> None:
    data = rss[top].T
    data.to_csv(outdir / named("source_rss_specificity", suffix, ".tsv"), sep="\t")
    heatmap_pdf(
        data,
        outdir / named("02_eregulon_specificity_heatmap", suffix),
        cmap=SEQ_CMAP,
        center=None,
        cbar_label="specificity",
        title="eRegulon specificity" if not condition else "eRegulon specificity by state and condition",
        figsize=(4.4, 5.8),
        legend=(
            "Regulon specificity score is calculated against cell-state/condition strata. Columns are ordered by cell label then condition; rows are clustered by specificity profile."
            if condition
            else "Regulon specificity score summarizes how selectively an eRegulon is active in each cell label. Rows and columns are clustered by profile similarity."
        ),
        cluster_cols=not condition,
    )


def plot_dot_heatmap(
    auc: pd.DataFrame,
    meta: pd.DataFrame,
    mean_by_group: pd.DataFrame,
    top: list[str],
    group_values: pd.Series,
    outdir: Path,
    suffix: str,
    condition: bool = False,
) -> None:
    active = auc[top].gt(auc[top].quantile(0.75), axis=1)
    frac = active.groupby(group_values).mean()
    z = zscore_rows(mean_by_group[top].T).T
    z_plot = cluster_frame(z[top], rows=not condition, cols=True)
    groups = list(z_plot.index)
    top = list(z_plot.columns)
    rows = []
    for group in groups:
        for regulon in top:
            rows.append(
                {
                    "group": group,
                    "eregulon": regulon,
                    "z_auc": float(z.loc[group, regulon]),
                    "active_fraction": float(frac.loc[group, regulon]) if group in frac.index else 0.0,
                }
            )
    dot = pd.DataFrame(rows)
    dot.to_csv(outdir / named("source_dot_heatmap", suffix, ".tsv"), sep="\t", index=False)
    x_map = {value: idx for idx, value in enumerate(top)}
    y_map = {value: idx for idx, value in enumerate(groups)}
    vmax = max(1.5, float(np.nanmax(np.abs(dot["z_auc"]))))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    sc = ax.scatter(
        [x_map[x] for x in dot["eregulon"]],
        [y_map[y] for y in dot["group"]],
        c=dot["z_auc"],
        s=style_float("dot_min_size", 8) + style_float("dot_size_scale", 95) * dot["active_fraction"],
        cmap=DIV_CMAP,
        vmin=-vmax,
        vmax=vmax,
        edgecolors="none",
        alpha=0.92,
    )
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(display_index(top), rotation=90, ha="center")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("eRegulon activity dot heatmap" if not condition else "eRegulon activity by state and condition", pad=4)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(sc, ax=ax, label="AUC z-score", shrink=0.60, pad=0.015)
    add_figure_legend(
        fig,
        (
            "Each dot summarizes one eRegulon in one cell-state/condition stratum. Dot color shows row-scaled mean AUCell activity and dot size shows the fraction of cells above the eRegulon's 75th percentile."
            if condition
            else "Each dot summarizes one eRegulon in one cell label. Dot color shows row-scaled mean AUCell activity and dot size shows the fraction of active cells."
        ),
    )
    fig.subplots_adjust(left=0.10, right=0.91, bottom=0.56, top=0.84)
    pdf = outdir / named("03_eregulon_dot_heatmap", suffix)
    save_pdf(fig, pdf)


def plot_umap(auc: pd.DataFrame, meta: pd.DataFrame, top: list[str], args: argparse.Namespace, outdir: Path) -> None:
    pdf = outdir / named("04_eregulon_auc_umap", args.file_suffix)
    if args.umap_x not in meta.columns or args.umap_y not in meta.columns:
        pairs = []
        for col in meta.columns:
            if str(col).endswith("_1"):
                mate = f"{str(col)[:-2]}_2"
                if mate in meta.columns and "umap" in str(col).lower():
                    pairs.append((str(col), mate))
        if pairs:
            args.umap_x, args.umap_y = pairs[0]
    if args.umap_x not in meta.columns or args.umap_y not in meta.columns:
        with PdfPages(pdf) as pp:
            fig, ax = plt.subplots(figsize=(3.5, 2.2))
            ax.text(0.5, 0.5, "UMAP columns not found", ha="center", va="center")
            ax.axis("off")
            add_figure_legend(fig, "UMAP coordinates were not present in the metadata, so eRegulon activity could not be projected onto the annotated embedding.")
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        require_pdf(pdf)
        return
    coords = meta[[args.umap_x, args.umap_y]].apply(pd.to_numeric, errors="coerce")
    keep = coords.notna().all(axis=1)
    with PdfPages(pdf) as pp:
        for start in range(0, min(len(top), args.umap_n), 6):
            subset = top[start : start + 6]
            fig, axes = plt.subplots(2, 3, figsize=(8.2, 4.6))
            axes = axes.ravel()
            for ax, regulon in zip(axes, subset):
                values = auc.loc[coords.index, regulon].astype(float)
                q1, q99 = np.nanpercentile(values, [1, 99])
                sc = ax.scatter(
                    coords.loc[keep, args.umap_x],
                    coords.loc[keep, args.umap_y],
                    c=values.loc[keep],
                    s=style_float("umap_point_size", 2.0),
                    cmap=UMAP_CMAP,
                    vmin=q1,
                    vmax=q99,
                    linewidths=0,
                    alpha=style_float("umap_point_alpha", 0.88),
                )
                ax.set_title(regulon_display_label(regulon), fontsize=6.5, pad=2)
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.015)
            for ax in axes[len(subset) :]:
                ax.axis("off")
            add_figure_legend(
                fig,
                "Each panel projects single-cell eRegulon AUCell activity onto the annotated embedding. Color limits are clipped to the 1st and 99th percentiles for each eRegulon to keep sparse high-activity cells visible without changing source values.",
            )
            fig.subplots_adjust(left=0.04, right=0.955, bottom=0.12, top=0.91, wspace=0.32, hspace=0.28)
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    require_pdf(pdf)


def plot_umap_by_condition(auc: pd.DataFrame, meta: pd.DataFrame, top: list[str], args: argparse.Namespace, outdir: Path) -> None:
    if not args.condition_col or args.condition_col not in meta.columns:
        return
    pdf = outdir / named("04_eregulon_auc_umap", f"condition_{args.file_suffix}")
    if args.umap_x not in meta.columns or args.umap_y not in meta.columns:
        pairs = []
        for col in meta.columns:
            if str(col).endswith("_1"):
                mate = f"{str(col)[:-2]}_2"
                if mate in meta.columns and "umap" in str(col).lower():
                    pairs.append((str(col), mate))
        if pairs:
            args.umap_x, args.umap_y = pairs[0]
    if args.umap_x not in meta.columns or args.umap_y not in meta.columns:
        return
    coords = meta[[args.umap_x, args.umap_y]].apply(pd.to_numeric, errors="coerce")
    keep = coords.notna().all(axis=1)
    conditions = ordered_unique(meta[args.condition_col])
    ncols = min(style_int("max_panels_per_row", 3), max(1, len(conditions)))
    nrows = int(np.ceil(len(conditions) / ncols))
    with PdfPages(pdf) as pp:
        for regulon in top[: args.umap_n]:
            values = auc.loc[coords.index, regulon].astype(float)
            q1, q99 = np.nanpercentile(values, [1, 99])
            fig, axes = plt.subplots(nrows, ncols, figsize=(2.65 * ncols, 2.45 * nrows), squeeze=False)
            axes_flat = axes.ravel()
            for ax, condition in zip(axes_flat, conditions):
                mask = keep & (meta[args.condition_col].astype(str) == condition)
                sc_handle = ax.scatter(
                    coords.loc[mask, args.umap_x],
                    coords.loc[mask, args.umap_y],
                    c=values.loc[mask],
                    s=style_float("umap_point_size", 2.0),
                    cmap=UMAP_CMAP,
                    vmin=q1,
                    vmax=q99,
                    linewidths=0,
                    alpha=style_float("umap_point_alpha", 0.88),
                )
                ax.set_title(f"{regulon_display_label(regulon)} | {condition}", fontsize=6.5, pad=2)
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                fig.colorbar(sc_handle, ax=ax, fraction=0.035, pad=0.015)
            for ax in axes_flat[len(conditions) :]:
                ax.axis("off")
            add_figure_legend(
                fig,
                "Condition-resolved panels display the same eRegulon on the same embedding, split by condition. Color limits are shared across condition panels for the displayed eRegulon.",
            )
            fig.subplots_adjust(left=0.04, right=0.96, bottom=0.13, top=0.88, wspace=0.30, hspace=0.34)
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    require_pdf(pdf)


def plot_activity_embedding(auc: pd.DataFrame, meta: pd.DataFrame, group_col: str, condition_col: str, outdir: Path, suffix: str) -> None:
    pdf = outdir / named("05_eregulon_activity_embedding", suffix)
    adata = ad.AnnData(auc.astype(float))
    adata.obs[group_col] = meta.loc[adata.obs_names, group_col].astype(str).values
    if condition_col and condition_col in meta.columns:
        adata.obs[condition_col] = meta.loc[adata.obs_names, condition_col].astype(str).values
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=min(15, max(2, adata.n_obs - 1)))
    sc.tl.umap(adata, random_state=7)
    coords = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["eregulon_umap_1", "eregulon_umap_2"])
    coords[group_col] = adata.obs[group_col].values
    if condition_col and condition_col in adata.obs:
        coords[condition_col] = adata.obs[condition_col].values
    coords.to_csv(outdir / named("source_eregulon_activity_embedding", suffix, ".tsv"), sep="\t")
    groups = ordered_unique(coords[group_col])
    palette = {g: group_color(g, i) for i, g in enumerate(groups)}
    conditions = ordered_unique(coords[condition_col]) if condition_col and condition_col in coords.columns else []
    panels = ["All"] + conditions
    ncols = min(style_int("max_panels_per_row", 3), len(panels))
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.95 * ncols, 2.75 * nrows), squeeze=False)
    axes_flat = axes.ravel()
    for ax, panel in zip(axes_flat, panels):
        if panel == "All":
            visible = coords.index
            title = "All cells"
        else:
            visible = coords.index[coords[condition_col].astype(str) == panel]
            title = str(panel)
            ax.scatter(
                coords["eregulon_umap_1"],
                coords["eregulon_umap_2"],
                s=style_float("condition_background_point_size", 1.2),
                color="#D0D0D0",
                alpha=style_float("condition_background_alpha", 0.18),
                linewidths=0,
            )
        for i, group in enumerate(groups):
            sub = coords.loc[coords.index.intersection(visible)]
            sub = sub.loc[sub[group_col] == group]
            if sub.empty:
                continue
            ax.scatter(
                sub["eregulon_umap_1"],
                sub["eregulon_umap_2"],
                s=style_float("embedding_point_size", 3.0),
                color=palette[group],
                alpha=style_float("embedding_point_alpha", 0.75),
                linewidths=0,
                label=group if panel == "All" else None,
            )
        ax.set_xlabel("eRegulon UMAP 1")
        ax.set_ylabel("eRegulon UMAP 2")
        ax.set_title(title, pad=4)
        for spine in ax.spines.values():
            spine.set_visible(False)
    for ax in axes_flat[len(panels) :]:
        ax.axis("off")
    handles, labels = axes_flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, frameon=False, loc="center right", bbox_to_anchor=(0.995, 0.52), markerscale=2)
    add_figure_legend(
        fig,
        "The embedding is computed from the eRegulon AUCell matrix. The first panel shows all cells by cell label; condition panels reuse the same coordinates and display each condition separately with the same color mapping.",
    )
    fig.subplots_adjust(left=0.08, right=0.86 if handles else 0.97, bottom=0.14, top=0.89, wspace=0.34, hspace=0.36)
    save_pdf(fig, pdf)


def parse_distance(value: object) -> float:
    text = str(value)
    nums = re.findall(r"-?\d+\.?\d*", text)
    if not nums:
        return np.nan
    return float(nums[0])


def plot_region_gene_link_structure(args: argparse.Namespace, outdir: Path) -> None:
    pdf = outdir / named("06_region_gene_link_structure", args.file_suffix)
    path = Path(args.region_to_gene)
    with PdfPages(pdf) as pp:
        if not path.exists() or path.stat().st_size == 0:
            fig, ax = plt.subplots(figsize=(3.5, 2.2))
            ax.text(0.5, 0.5, "region-to-gene table not found", ha="center", va="center")
            ax.axis("off")
            add_figure_legend(fig, "No region-to-gene adjacency table was available, so link-structure diagnostics could not be drawn.")
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            pd.DataFrame().to_csv(outdir / named("source_region_gene_link_structure", args.file_suffix, ".tsv"), sep="\t", index=False)
        else:
            df = pd.read_csv(path, sep="\t")
            gene_col = choose_column(df.columns, ["target", "gene", "Gene", "target_gene"])
            region_col = choose_column(df.columns, ["region", "Region"])
            importance_col = choose_column(df.columns, ["importance", "importance_x_abs_rho", "importance_x_rho"])
            distance_col = choose_column(df.columns, ["Distance", "distance"])
            if not gene_col or not region_col:
                raise ValueError("region-to-gene table must contain gene and region columns.")
            counts = df.groupby(gene_col)[region_col].nunique().rename("n_regions").reset_index()
            counts.to_csv(outdir / named("source_regions_per_gene", args.file_suffix, ".tsv"), sep="\t", index=False)
            summary = {
                "n_links": len(df),
                "n_genes": counts[gene_col].nunique(),
                "median_regions_per_gene": float(counts["n_regions"].median()),
                "fraction_genes_1_10_regions": float(counts["n_regions"].between(1, 10).mean()),
            }
            pd.DataFrame([summary]).to_csv(outdir / named("source_region_gene_link_structure", args.file_suffix, ".tsv"), sep="\t", index=False)
            fig, axes = plt.subplots(1, 2, figsize=(6.2, 2.6))
            axes[0].hist(counts["n_regions"], bins=np.arange(1, min(50, counts["n_regions"].max() + 2)), color="#667A9A", alpha=0.9)
            axes[0].set_xlabel("regions per gene")
            axes[0].set_ylabel("genes")
            axes[0].set_title("Enhancer load per gene", pad=4)
            if distance_col and importance_col:
                sample = df[[distance_col, importance_col]].copy()
                sample["abs_distance"] = sample[distance_col].map(parse_distance).abs()
                sample["importance"] = pd.to_numeric(sample[importance_col], errors="coerce")
                sample = sample.dropna().sample(min(len(sample.dropna()), 200000), random_state=7)
                axes[1].hexbin(np.log10(sample["abs_distance"] + 1), sample["importance"], gridsize=45, cmap=SEQ_CMAP, mincnt=1)
                axes[1].set_xlabel("log10 distance to gene + 1")
                axes[1].set_ylabel("region-gene importance")
            else:
                axes[1].axis("off")
                axes[1].text(0.5, 0.5, "distance or importance column not found", ha="center", va="center")
            axes[1].set_title("Distance vs link score", pad=4)
            for ax in axes:
                for spine in ax.spines.values():
                    spine.set_visible(False)
            add_figure_legend(
                fig,
                "Region-to-gene links are model-level SCENIC+ associations rather than condition-specific subsets. The left panel summarizes enhancer load per gene; the right panel relates genomic distance to link importance where both fields are available.",
            )
            fig.subplots_adjust(left=0.10, right=0.97, bottom=0.28, top=0.82, wspace=0.34)
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    require_pdf(pdf)


def plot_target_region_overlap(args: argparse.Namespace, top: list[str], outdir: Path) -> None:
    pdf = outdir / named("07_target_region_overlap_heatmap", args.file_suffix)
    path = Path(args.eregulons)
    if not path.exists() or path.stat().st_size == 0:
        pd.DataFrame().to_csv(outdir / named("source_target_region_overlap", args.file_suffix, ".tsv"), sep="\t")
        heatmap_pdf(pd.DataFrame([[0]], index=["missing"], columns=["missing"]), pdf, SEQ_CMAP, None, "Jaccard", "Target-region overlap")
        return
    df = pd.read_csv(path, sep="\t")
    name_col = choose_column(df.columns, ["eRegulon_name", "eregulon", "Gene_signature_name", "Region_signature_name"])
    region_col = choose_column(df.columns, ["Region", "region"])
    if not name_col or not region_col:
        raise ValueError("eRegulon table must contain eRegulon and Region columns.")
    selected = [x for x in top if x in set(df[name_col].astype(str))]
    if len(selected) < 2:
        selected = list(pd.unique(df[name_col].astype(str)))[: min(20, df[name_col].nunique())]
    selected = selected[: min(25, len(selected))]
    region_sets = {name: set(df.loc[df[name_col].astype(str) == name, region_col].astype(str)) for name in selected}
    mat = pd.DataFrame(index=selected, columns=selected, dtype=float)
    for a in selected:
        for b in selected:
            union = len(region_sets[a] | region_sets[b])
            mat.loc[a, b] = len(region_sets[a] & region_sets[b]) / union if union else 0.0
    mat.to_csv(outdir / named("source_target_region_overlap", args.file_suffix, ".tsv"), sep="\t")
    heatmap_pdf(
        mat,
        pdf,
        OVERLAP_CMAP,
        None,
        "Jaccard",
        "Target-region overlap",
        figsize=(5.8, 5.8),
        x_rotation=90,
        legend="Each cell is the Jaccard overlap between two eRegulons' linked regulatory-region sets. This is a model-level structural view; condition effects are evaluated with AUCell activity statistics.",
    )


def parse_list_like(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value)
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, (list, tuple, set)):
            return [str(x) for x in parsed]
    except Exception:
        pass
    return [x.strip() for x in re.split(r"[;,|]", text) if x.strip()]


def read_network_edges(args: argparse.Namespace, top_tfs: list[str]) -> pd.DataFrame:
    tf_path = Path(args.tf_to_gene)
    edges = []
    if tf_path.exists() and tf_path.stat().st_size > 0:
        df = pd.read_csv(tf_path, sep="\t")
        tf_col = choose_column(df.columns, ["TF", "tf", "source", "transcription_factor"])
        gene_col = choose_column(df.columns, ["target", "gene", "Gene", "target_gene"])
        weight_col = choose_column(df.columns, ["importance", "weight", "rho", "score", "adj_importance"])
        if tf_col and gene_col:
            df = df.loc[df[tf_col].astype(str).isin(top_tfs)].copy()
            df["_weight"] = pd.to_numeric(df[weight_col], errors="coerce").fillna(0).abs() if weight_col else 1.0
            df = df.sort_values([tf_col, "_weight"], ascending=[True, False])
            for tf, sub in df.groupby(tf_col):
                for _, row in sub.head(args.network_targets_per_tf).iterrows():
                    edges.append((str(row[tf_col]), str(row[gene_col]), float(row["_weight"])))
    if not edges:
        ereg_path = Path(args.eregulons)
        if ereg_path.exists() and ereg_path.stat().st_size > 0:
            df = pd.read_csv(ereg_path, sep="\t")
            tf_col = choose_column(df.columns, ["TF", "tf", "transcription_factor"])
            gene_col = choose_column(df.columns, ["target_genes", "Target_genes", "genes", "Genes", "gene"])
            name_col = choose_column(df.columns, ["eRegulon_name", "eregulon", "Gene_signame", "Region_signame"])
            for _, row in df.iterrows():
                tf = str(row[tf_col]) if tf_col else tf_from_regulon(row[name_col]) if name_col else ""
                if tf not in top_tfs:
                    continue
                genes = parse_list_like(row[gene_col]) if gene_col else []
                for gene in genes[: args.network_targets_per_tf]:
                    edges.append((tf, gene, 1.0))
    return pd.DataFrame(edges, columns=["tf", "target_gene", "weight"]).drop_duplicates()


def plot_network(args: argparse.Namespace, top: list[str], outdir: Path) -> None:
    top_tfs = []
    for regulon in top:
        tf = tf_from_regulon(regulon)
        if tf and tf not in top_tfs:
            top_tfs.append(tf)
        if len(top_tfs) >= args.network_top_tfs:
            break
    edge_df = read_network_edges(args, top_tfs)
    edge_df.to_csv(outdir / named("source_tf_target_network_edges", args.file_suffix, ".tsv"), sep="\t", index=False)
    pdf = outdir / named("08_tf_target_network", args.file_suffix)
    with PdfPages(pdf) as pp:
        if edge_df.empty:
            fig, ax = plt.subplots(figsize=(3.5, 2.2))
            ax.text(0.5, 0.5, "No TF-target edges found", ha="center", va="center")
            ax.axis("off")
            add_figure_legend(fig, "No TF-target edges were available for the selected eRegulons.")
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        else:
            graph = nx.DiGraph()
            for row in edge_df.itertuples(index=False):
                graph.add_edge(row.tf, row.target_gene, weight=row.weight)
            pos = nx.spring_layout(graph, seed=7, k=0.75)
            tf_nodes = sorted(edge_df["tf"].unique())
            target_nodes = [node for node in graph.nodes if node not in tf_nodes]
            fig, ax = plt.subplots(figsize=(6.8, 5.2))
            widths = [0.35 + min(1.8, graph[u][v].get("weight", 1.0)) for u, v in graph.edges]
            nx.draw_networkx_edges(graph, pos, ax=ax, alpha=style_float("network_edge_alpha", 0.24), arrows=True, width=widths, edge_color="#777777")
            nx.draw_networkx_nodes(graph, pos, nodelist=target_nodes, node_size=style_float("network_target_node_size", 18), node_color="#9AA6B2", ax=ax, linewidths=0)
            nx.draw_networkx_nodes(graph, pos, nodelist=tf_nodes, node_size=style_float("network_tf_node_size", 150), node_color="#B57457", ax=ax, linewidths=0.3, edgecolors="#6F3A2A")
            labels = {node: node for node in tf_nodes}
            labels.update({node: short_text(node, 12) for node in target_nodes if graph.degree[node] > 1})
            nx.draw_networkx_labels(graph, pos, labels=labels, font_size=5.5, ax=ax)
            ax.set_title("TF-target network", pad=4)
            ax.axis("off")
            add_figure_legend(
                fig,
                "The network displays the highest-ranked TF-target links among selected eRegulons. TFs are highlighted as larger nodes; target genes are smaller nodes; edge width follows the available link weight.",
            )
            fig.subplots_adjust(left=0.02, right=0.98, bottom=0.10, top=0.91)
            pp.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    require_pdf(pdf)


def main() -> None:
    args = resolve_args_paths(parse_args())
    set_publication_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    auc, meta = read_auc_and_metadata(args)
    mean_by_group = auc.groupby(meta[args.group_col]).mean()
    mean_by_group.to_csv(outdir / named("source_auc_mean_by_group", args.file_suffix, ".tsv"), sep="\t")
    rss = calculate_rss(auc, meta[args.group_col])
    top = select_top_regulons(mean_by_group, rss, args.top_n)
    pd.DataFrame({"eregulon": top, "display_label": display_index(top)}).to_csv(
        outdir / named("source_selected_top_eregulons", args.file_suffix, ".tsv"), sep="\t", index=False
    )
    plot_auc_heatmap(mean_by_group, top, outdir, args.file_suffix)
    plot_rss_heatmap(rss, top, outdir, args.file_suffix)
    plot_dot_heatmap(auc, meta, mean_by_group, top, args.group_col, outdir, args.file_suffix)
    plot_activity_embedding(auc, meta, args.group_col, outdir, args.file_suffix)
    plot_umap(auc, meta, top, args, outdir)
    plot_region_gene_link_structure(args, outdir)
    plot_target_region_overlap(args, top, outdir)
    plot_network(args, top, outdir)
    print("SCENIC+ result figures written to", outdir)


if __name__ == "__main__":
    main()
