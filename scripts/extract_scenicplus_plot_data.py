#!/usr/bin/env python
"""Extract SCENIC+ plotting source tables without drawing figures."""
from __future__ import annotations

import argparse
import ast
import os
import re
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()

import anndata as ad
import mudata as md
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import jensenshannon, pdist


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract source tables used by the R SCENIC+ figure renderer.")
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
    parser.add_argument(
        "--regulon-sign-filter",
        choices=["tf_positive", "all"],
        default="tf_positive",
        help="Default keeps only SCENIC+ +/+ eRegulons for source tables and figures. Use all for exploratory signed regulon output.",
    )
    return parser.parse_args()


def resolve_project_path(path_value: str) -> Path:
    path = Path(path_value).expanduser()
    return path if path.is_absolute() else (PROJECT / path).resolve()


def named(stem: str, suffix: str, ext: str = ".tsv") -> str:
    return f"{stem}_{suffix}{ext}" if suffix else f"{stem}{ext}"


def clean_regulon_name(name: object) -> str:
    return re.sub(r"\s+", "_", str(name))


def tf_from_regulon(name: object) -> str:
    text = clean_regulon_name(name)
    text = re.sub(r"_(direct|extended).*$", "", text)
    text = re.sub(r"_[+-]/[+-].*$", "", text)
    return re.split(r"[_()]", text)[0]


def regulon_sign(name: object) -> str | None:
    match = re.search(r"([+-]/[+-])", clean_regulon_name(name))
    return match.group(1) if match else None


def filter_regulon_columns(auc: pd.DataFrame, mode: str) -> pd.DataFrame:
    if mode == "all":
        return auc
    signs = pd.Series({col: regulon_sign(col) for col in auc.columns}, dtype=object)
    signed = signs.notna()
    if not signed.any():
        print("No SCENIC+ eRegulon sign labels detected; keeping all regulons.")
        return auc
    keep = signs.index[signs == "+/+"].tolist()
    if not keep:
        raise ValueError("No +/+ eRegulons found after applying --regulon-sign-filter=tf_positive.")
    print(f"Keeping {len(keep)} +/+ eRegulons; filtered {int(signed.sum()) - len(keep)} non +/+ signed eRegulons.")
    return auc.loc[:, keep]


def filter_eregulon_rows(df: pd.DataFrame, name_col: str, mode: str) -> pd.DataFrame:
    if mode == "all":
        return df
    signs = df[name_col].map(regulon_sign)
    if signs.notna().any():
        return df.loc[signs == "+/+"].copy()
    return df


def regulon_display_label(name: object) -> str:
    text = clean_regulon_name(name)
    tf = tf_from_regulon(text)
    genes = re.search(r"\((\d+)g\)", text)
    signs = re.search(r"([+-]/[+-])", text)
    sign_label = f" {signs.group(1)}" if signs and signs.group(1) != "+/+" else ""
    return f"{tf}{sign_label} ({genes.group(1)} targets)" if genes else f"{tf}{sign_label}"


def ordered_unique(values: pd.Series) -> list[str]:
    present = [str(x) for x in pd.unique(values.astype(str)) if str(x)]
    return present


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
    meta = pd.read_csv(resolve_project_path(args.metadata), sep="\t", dtype=str).fillna("")
    required = [args.cell_col, args.group_col]
    if args.condition_col:
        required.append(args.condition_col)
    for col in required:
        if col not in meta.columns:
            raise ValueError(f"metadata missing required column: {col}")
    meta = meta.drop_duplicates(args.cell_col).set_index(args.cell_col)
    mdata = md.read_h5mu(resolve_project_path(args.auc_h5mu))
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


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.astype(float)
    mean = values.mean(axis=1)
    sd = values.std(axis=1).replace(0, np.nan)
    return values.sub(mean, axis=0).div(sd, axis=0).fillna(0)


def calculate_rss(auc: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    rss = pd.DataFrame(index=ordered_unique(groups), columns=auc.columns, dtype=float)
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


def make_joint_groups(meta: pd.DataFrame, group_col: str, condition_col: str) -> tuple[pd.Series, list[str]]:
    group_order = ordered_unique(meta[group_col])
    condition_order = ordered_unique(meta[condition_col])
    joint = meta[group_col].astype(str) + "_" + meta[condition_col].astype(str)
    order = [f"{group}_{condition}" for group in group_order for condition in condition_order if (joint == f"{group}_{condition}").any()]
    return joint, order


def write_matrix_long(matrix: pd.DataFrame, out: Path, value_name: str) -> None:
    data = matrix.copy()
    data.insert(0, "eregulon", data.index.astype(str))
    long = data.melt(id_vars="eregulon", var_name="group", value_name=value_name)
    long.insert(1, "display_label", long["eregulon"].map(regulon_display_label))
    long.to_csv(out, sep="\t", index=False)


def write_dot_table(auc: pd.DataFrame, groups: pd.Series, mean_by_group: pd.DataFrame, top: list[str], out: Path, preserve_group_order: list[str] | None = None) -> None:
    active = auc[top].gt(auc[top].quantile(0.75), axis=1)
    frac = active.groupby(groups).mean()
    z = zscore_rows(mean_by_group[top].T).T
    z_plot = cluster_frame(z[top], rows=preserve_group_order is None, cols=True)
    group_order = preserve_group_order if preserve_group_order is not None else list(z_plot.index)
    regulon_order = list(z_plot.columns)
    rows = []
    for group in group_order:
        if group not in z.index:
            continue
        for regulon in regulon_order:
            rows.append(
                {
                    "group": group,
                    "eregulon": regulon,
                    "display_label": regulon_display_label(regulon),
                    "z_auc": float(z.loc[group, regulon]),
                    "active_fraction": float(frac.loc[group, regulon]) if group in frac.index else 0.0,
                    "group_order": group_order.index(group) + 1,
                    "eregulon_order": regulon_order.index(regulon) + 1,
                }
            )
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


def parse_distance(value: object) -> float:
    nums = re.findall(r"-?\d+\.?\d*", str(value))
    return float(nums[0]) if nums else np.nan


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
    tf_path = resolve_project_path(args.tf_to_gene)
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
        ereg_path = resolve_project_path(args.eregulons)
        if ereg_path.exists() and ereg_path.stat().st_size > 0:
            df = pd.read_csv(ereg_path, sep="\t")
            tf_col = choose_column(df.columns, ["TF", "tf", "transcription_factor"])
            gene_col = choose_column(df.columns, ["target_genes", "Target_genes", "genes", "Genes", "gene"])
            name_col = choose_column(df.columns, ["eRegulon_name", "eregulon", "Gene_signame", "Region_signame"])
            if name_col:
                df = filter_eregulon_rows(df, name_col, args.regulon_sign_filter)
            for _, row in df.iterrows():
                tf = str(row[tf_col]) if tf_col else tf_from_regulon(row[name_col]) if name_col else ""
                if tf not in top_tfs:
                    continue
                genes = parse_list_like(row[gene_col]) if gene_col else []
                for gene in genes[: args.network_targets_per_tf]:
                    edges.append((tf, gene, 1.0))
    return pd.DataFrame(edges, columns=["tf", "target_gene", "weight"]).drop_duplicates()


def main() -> None:
    args = parse_args()
    outdir = resolve_project_path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    auc, meta = read_auc_and_metadata(args)
    auc = filter_regulon_columns(auc, args.regulon_sign_filter)

    group_order = ordered_unique(meta[args.group_col])
    mean_by_group = auc.groupby(meta[args.group_col]).mean().reindex(group_order)
    rss = calculate_rss(auc, meta[args.group_col])
    top = select_top_regulons(mean_by_group, rss, args.top_n)
    pd.DataFrame({"eregulon": top, "display_label": [regulon_display_label(x) for x in top]}).to_csv(
        outdir / named("source_selected_top_eregulons", args.file_suffix), sep="\t", index=False
    )

    auc_z = cluster_frame(zscore_rows(mean_by_group[top].T), rows=True, cols=True)
    write_matrix_long(auc_z, outdir / named("source_auc_heatmap_zscore", args.file_suffix), "z_auc")
    write_matrix_long(cluster_frame(rss[top].T, rows=True, cols=True), outdir / named("source_rss_specificity", args.file_suffix), "specificity")
    write_dot_table(auc, meta[args.group_col], mean_by_group, top, outdir / named("source_dot_heatmap", args.file_suffix))

    if args.condition_col and args.condition_col in meta.columns:
        joint, joint_order = make_joint_groups(meta, args.group_col, args.condition_col)
        mean_by_joint = auc.groupby(joint).mean().reindex(joint_order)
        rss_joint = calculate_rss(auc, joint).reindex(joint_order)
        auc_joint_z = zscore_rows(mean_by_joint[top].T)
        auc_joint_z = auc_joint_z.loc[cluster_order(auc_joint_z, axis=0), joint_order]
        write_matrix_long(auc_joint_z, outdir / named("source_auc_heatmap_zscore_condition", args.file_suffix), "z_auc")
        rss_joint_plot = rss_joint[top].T
        rss_joint_plot = rss_joint_plot.loc[cluster_order(rss_joint_plot, axis=0), joint_order]
        write_matrix_long(rss_joint_plot, outdir / named("source_rss_specificity_condition", args.file_suffix), "specificity")
        write_dot_table(auc, joint, mean_by_joint, top, outdir / named("source_dot_heatmap_condition", args.file_suffix), joint_order)

    coords_cols = [args.umap_x, args.umap_y]
    if args.umap_x not in meta.columns or args.umap_y not in meta.columns:
        pairs = []
        for col in meta.columns:
            if str(col).endswith("_1"):
                mate = f"{str(col)[:-2]}_2"
                if mate in meta.columns and "umap" in str(col).lower():
                    pairs.append((str(col), mate))
        if pairs:
            coords_cols = list(pairs[0])
    if all(col in meta.columns for col in coords_cols):
        coords = meta[coords_cols].apply(pd.to_numeric, errors="coerce")
        keep = coords.notna().all(axis=1)
        umap_rows = []
        for regulon in top[: args.umap_n]:
            tmp = pd.DataFrame(
                {
                    "cell_id": coords.index[keep],
                    "umap_1": coords.loc[keep, coords_cols[0]].values,
                    "umap_2": coords.loc[keep, coords_cols[1]].values,
                    "cell_label": meta.loc[coords.index[keep], args.group_col].values,
                    "condition": meta.loc[coords.index[keep], args.condition_col].values if args.condition_col in meta.columns else "",
                    "eregulon": regulon,
                    "display_label": regulon_display_label(regulon),
                    "auc": auc.loc[coords.index[keep], regulon].astype(float).values,
                }
            )
            umap_rows.append(tmp)
        pd.concat(umap_rows, ignore_index=True).to_csv(outdir / named("source_eregulon_auc_umap", args.file_suffix), sep="\t", index=False)

    adata = ad.AnnData(auc.astype(float))
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=min(15, max(2, adata.n_obs - 1)))
    sc.tl.umap(adata, random_state=7)
    emb = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["eregulon_umap_1", "eregulon_umap_2"])
    emb.insert(0, "cell_id", emb.index)
    emb["cell_label"] = meta.loc[emb.index, args.group_col].astype(str).values
    emb["condition"] = meta.loc[emb.index, args.condition_col].astype(str).values if args.condition_col in meta.columns else ""
    emb.to_csv(outdir / named("source_eregulon_activity_embedding", args.file_suffix), sep="\t", index=False)

    rg_path = resolve_project_path(args.region_to_gene)
    if rg_path.exists() and rg_path.stat().st_size > 0:
        df = pd.read_csv(rg_path, sep="\t")
        gene_col = choose_column(df.columns, ["target", "gene", "Gene", "target_gene"])
        region_col = choose_column(df.columns, ["region", "Region"])
        importance_col = choose_column(df.columns, ["importance", "importance_x_abs_rho", "importance_x_rho"])
        distance_col = choose_column(df.columns, ["Distance", "distance"])
        if gene_col and region_col:
            counts = df.groupby(gene_col)[region_col].nunique().rename("n_regions").reset_index().rename(columns={gene_col: "gene"})
            counts.to_csv(outdir / named("source_regions_per_gene", args.file_suffix), sep="\t", index=False)
            summary = {
                "n_links": len(df),
                "n_genes": counts["gene"].nunique(),
                "median_regions_per_gene": float(counts["n_regions"].median()),
                "fraction_genes_1_10_regions": float(counts["n_regions"].between(1, 10).mean()),
            }
            pd.DataFrame([summary]).to_csv(outdir / named("source_region_gene_link_structure", args.file_suffix), sep="\t", index=False)
            if distance_col and importance_col:
                sample = df[[distance_col, importance_col]].copy()
                sample["abs_distance"] = sample[distance_col].map(parse_distance).abs()
                sample["importance"] = pd.to_numeric(sample[importance_col], errors="coerce")
                sample = sample.dropna().sample(min(len(sample.dropna()), 200000), random_state=7)
                sample[["abs_distance", "importance"]].to_csv(outdir / named("source_region_distance_importance", args.file_suffix), sep="\t", index=False)

    ereg_path = resolve_project_path(args.eregulons)
    if ereg_path.exists() and ereg_path.stat().st_size > 0:
        df = pd.read_csv(ereg_path, sep="\t")
        name_col = choose_column(df.columns, ["eRegulon_name", "eregulon", "Gene_signature_name", "Region_signature_name"])
        region_col = choose_column(df.columns, ["Region", "region"])
        if name_col and region_col:
            df = filter_eregulon_rows(df, name_col, args.regulon_sign_filter)
            selected = [x for x in top if x in set(df[name_col].astype(str))]
            if len(selected) < 2:
                selected = list(pd.unique(df[name_col].astype(str)))[: min(20, df[name_col].nunique())]
            selected = selected[: min(25, len(selected))]
            region_sets = {name: set(df.loc[df[name_col].astype(str) == name, region_col].astype(str)) for name in selected}
            rows = []
            for a in selected:
                for b in selected:
                    union = len(region_sets[a] | region_sets[b])
                    rows.append({"eregulon_a": a, "eregulon_b": b, "display_a": regulon_display_label(a), "display_b": regulon_display_label(b), "jaccard": len(region_sets[a] & region_sets[b]) / union if union else 0.0})
            pd.DataFrame(rows).to_csv(outdir / named("source_target_region_overlap", args.file_suffix), sep="\t", index=False)

    top_tfs = []
    for regulon in top:
        tf = tf_from_regulon(regulon)
        if tf and tf not in top_tfs:
            top_tfs.append(tf)
        if len(top_tfs) >= args.network_top_tfs:
            break
    read_network_edges(args, top_tfs).to_csv(outdir / named("source_tf_target_network_edges", args.file_suffix), sep="\t", index=False)
    print("SCENIC+ plotting source tables written to", outdir)


if __name__ == "__main__":
    main()
