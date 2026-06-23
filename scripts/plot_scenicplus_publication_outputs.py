#!/usr/bin/env python
"""SCENIC+ result figures from AUCell and eRegulon tables."""
from __future__ import annotations

import argparse
import ast
import os
import re
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import mudata as md
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import jensenshannon


GROUP_COLORS = [
    "#2F4858",
    "#B57457",
    "#667A9A",
    "#9A8067",
    "#7D9477",
    "#8E5A63",
    "#4F6D7A",
    "#C49A6C",
    "#5F6F52",
    "#8C6A5D",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--auc-h5mu", default="results/scenicplus/AUCell_direct.h5mu")
    parser.add_argument("--metadata", default="inputs/cell_metadata.tsv")
    parser.add_argument("--tf-to-gene", default="work/scenicplus/tf_to_gene_adj.tsv")
    parser.add_argument("--eregulons", default="results/scenicplus/eRegulons_direct.tsv")
    parser.add_argument("--group-col", default="cell_label")
    parser.add_argument("--cell-col", default="cell_id")
    parser.add_argument("--umap-x", default="umap_1")
    parser.add_argument("--umap-y", default="umap_2")
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--umap-n", type=int, default=12)
    parser.add_argument("--network-top-tfs", type=int, default=8)
    parser.add_argument("--network-targets-per-tf", type=int, default=12)
    parser.add_argument("--outdir", default="results/scenicplus_figures/direct")
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
    for attr in ["auc_h5mu", "metadata", "tf_to_gene", "eregulons", "outdir"]:
        setattr(args, attr, resolve_project_path(getattr(args, attr), base))
    return args


def set_publication_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.6,
            "axes.labelsize": 6.5,
            "axes.titlesize": 7,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
            "figure.dpi": 300,
            "savefig.dpi": 300,
        }
    )
    sns.set_theme(style="white", rc={"axes.linewidth": 0.6})


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


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.astype(float)
    mean = values.mean(axis=1)
    sd = values.std(axis=1).replace(0, np.nan)
    return values.sub(mean, axis=0).div(sd, axis=0).fillna(0)


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
    for col in [args.cell_col, args.group_col]:
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


def heatmap_pdf(data: pd.DataFrame, pdf: Path, cmap: str, center: float | None, cbar_label: str, title: str) -> None:
    height = max(3.2, 0.16 * data.shape[0])
    width = max(3.6, 0.22 * data.shape[1])
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        data,
        cmap=cmap,
        center=center,
        ax=ax,
        linewidths=0,
        cbar_kws={"label": cbar_label, "shrink": 0.55},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, pad=4)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout(pad=0.6)
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    require_pdf(pdf)


def plot_auc_heatmap(mean_by_group: pd.DataFrame, top: list[str], outdir: Path) -> None:
    data = zscore_rows(mean_by_group[top].T)
    data.to_csv(outdir / "source_auc_heatmap_zscore.tsv", sep="\t")
    heatmap_pdf(
        data,
        outdir / "eregulon_auc_heatmap.pdf",
        cmap="vlag",
        center=0,
        cbar_label="AUC z-score",
        title="eRegulon activity",
    )


def plot_rss_heatmap(rss: pd.DataFrame, top: list[str], outdir: Path) -> None:
    data = rss[top].T
    data.to_csv(outdir / "source_rss_specificity.tsv", sep="\t")
    heatmap_pdf(
        data,
        outdir / "eregulon_specificity_heatmap.pdf",
        cmap="mako",
        center=None,
        cbar_label="specificity",
        title="eRegulon specificity",
    )


def plot_dot_heatmap(auc: pd.DataFrame, meta: pd.DataFrame, mean_by_group: pd.DataFrame, top: list[str], group_col: str, outdir: Path) -> None:
    active = auc[top].gt(auc[top].quantile(0.75), axis=1)
    frac = active.groupby(meta[group_col]).mean()
    z = zscore_rows(mean_by_group[top].T).T
    rows = []
    for group in z.index:
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
    dot.to_csv(outdir / "source_dot_heatmap.tsv", sep="\t", index=False)
    groups = list(z.index)
    x_map = {value: idx for idx, value in enumerate(top)}
    y_map = {value: idx for idx, value in enumerate(groups)}
    vmax = max(1.5, float(np.nanmax(np.abs(dot["z_auc"]))))
    fig, ax = plt.subplots(figsize=(max(4.2, 0.19 * len(top)), max(2.8, 0.28 * len(groups))))
    sc = ax.scatter(
        [x_map[x] for x in dot["eregulon"]],
        [y_map[y] for y in dot["group"]],
        c=dot["z_auc"],
        s=8 + 95 * dot["active_fraction"],
        cmap="vlag",
        vmin=-vmax,
        vmax=vmax,
        edgecolors="none",
        alpha=0.92,
    )
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top, rotation=90, ha="center")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("eRegulon activity dot heatmap", pad=4)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(sc, ax=ax, label="AUC z-score", shrink=0.55)
    fig.tight_layout(pad=0.6)
    pdf = outdir / "eregulon_dot_heatmap.pdf"
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    require_pdf(pdf)


def plot_umap(auc: pd.DataFrame, meta: pd.DataFrame, top: list[str], args: argparse.Namespace, outdir: Path) -> None:
    pdf = outdir / "eregulon_auc_umap.pdf"
    if args.umap_x not in meta.columns or args.umap_y not in meta.columns:
        with PdfPages(pdf) as pp:
            fig, ax = plt.subplots(figsize=(3.5, 2.2))
            ax.text(0.5, 0.5, "UMAP columns not found", ha="center", va="center")
            ax.axis("off")
            pp.savefig(fig)
            plt.close(fig)
        require_pdf(pdf)
        return
    coords = meta[[args.umap_x, args.umap_y]].apply(pd.to_numeric, errors="coerce")
    keep = coords.notna().all(axis=1)
    with PdfPages(pdf) as pp:
        for start in range(0, min(len(top), args.umap_n), 6):
            subset = top[start : start + 6]
            fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.6))
            axes = axes.ravel()
            for ax, regulon in zip(axes, subset):
                values = auc.loc[coords.index, regulon].astype(float)
                q1, q99 = np.nanpercentile(values, [1, 99])
                sc = ax.scatter(
                    coords.loc[keep, args.umap_x],
                    coords.loc[keep, args.umap_y],
                    c=values.loc[keep],
                    s=2.0,
                    cmap="viridis",
                    vmin=q1,
                    vmax=q99,
                    linewidths=0,
                    alpha=0.88,
                )
                ax.set_title(regulon, fontsize=6.5, pad=2)
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
            for ax in axes[len(subset) :]:
                ax.axis("off")
            fig.tight_layout(pad=0.4)
            pp.savefig(fig)
            plt.close(fig)
    require_pdf(pdf)


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
    edge_df.to_csv(outdir / "source_tf_target_network_edges.tsv", sep="\t", index=False)
    pdf = outdir / "tf_target_network.pdf"
    with PdfPages(pdf) as pp:
        if edge_df.empty:
            fig, ax = plt.subplots(figsize=(3.5, 2.2))
            ax.text(0.5, 0.5, "No TF-target edges found", ha="center", va="center")
            ax.axis("off")
            pp.savefig(fig)
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
            nx.draw_networkx_edges(graph, pos, ax=ax, alpha=0.24, arrows=True, width=widths, edge_color="#777777")
            nx.draw_networkx_nodes(graph, pos, nodelist=target_nodes, node_size=18, node_color="#9AA6B2", ax=ax, linewidths=0)
            nx.draw_networkx_nodes(graph, pos, nodelist=tf_nodes, node_size=150, node_color="#B57457", ax=ax, linewidths=0.3, edgecolors="#6F3A2A")
            labels = {node: node for node in tf_nodes}
            labels.update({node: node for node in target_nodes if graph.degree[node] > 1})
            nx.draw_networkx_labels(graph, pos, labels=labels, font_size=5.5, ax=ax)
            ax.set_title("TF-target network", pad=4)
            ax.axis("off")
            fig.tight_layout(pad=0.4)
            pp.savefig(fig)
            plt.close(fig)
    require_pdf(pdf)


def main() -> None:
    args = resolve_args_paths(parse_args())
    set_publication_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    auc, meta = read_auc_and_metadata(args)
    mean_by_group = auc.groupby(meta[args.group_col]).mean()
    mean_by_group.to_csv(outdir / "source_auc_mean_by_group.tsv", sep="\t")
    rss = calculate_rss(auc, meta[args.group_col])
    top = select_top_regulons(mean_by_group, rss, args.top_n)
    pd.Series(top, name="eregulon").to_csv(outdir / "source_selected_top_eregulons.tsv", sep="\t", index=False)
    plot_auc_heatmap(mean_by_group, top, outdir)
    plot_rss_heatmap(rss, top, outdir)
    plot_dot_heatmap(auc, meta, mean_by_group, top, args.group_col, outdir)
    plot_umap(auc, meta, top, args, outdir)
    plot_network(args, top, outdir)
    print("SCENIC+ result figures written to", outdir)


if __name__ == "__main__":
    main()
