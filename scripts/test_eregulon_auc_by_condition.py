from pathlib import Path
import argparse
import os
import re

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import mudata as md
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "font.size": 7,
        "axes.labelsize": 6.5,
        "axes.titlesize": 7,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "legend.frameon": False,
    }
)

parser = argparse.ArgumentParser()
parser.add_argument("--auc-h5mu", default="results/scenicplus/AUCell_direct.h5mu")
parser.add_argument("--metadata", default="inputs/cell_metadata.tsv")
parser.add_argument("--group-col", default="condition")
parser.add_argument("--label-col", default="")
parser.add_argument("--sample-col", default="sample_id")
parser.add_argument("--cell-col", default="cell_id")
parser.add_argument("--reference-condition", default="auto")
parser.add_argument("--comparison-condition", default="auto")
parser.add_argument("--outdir", default="results/scenicplus_figures/condition")
parser.add_argument("--file-suffix", default="")
parser.add_argument("--tables-only", action="store_true", help="Write condition statistics tables only; PDF plotting is handled by the R renderer.")
args = parser.parse_args()

def resolve_project_path(path_value: str) -> Path:
    path = Path(path_value).expanduser()
    return path if path.is_absolute() else PROJECT / path


outdir = resolve_project_path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    idx = np.where(ok)[0][order]
    q[idx] = adjusted
    return q


def clean_regulon_name(name):
    return re.sub(r"\s+", "_", str(name))


def tf_from_regulon(name):
    text = clean_regulon_name(name)
    text = re.sub(r"_(direct|extended).*$", "", text)
    text = re.sub(r"_[+-]/[+-].*$", "", text)
    return re.split(r"[_()]", text)[0]


def regulon_display_label(name):
    text = clean_regulon_name(name)
    tf = tf_from_regulon(text)
    genes = re.search(r"\((\d+)g\)", text)
    signs = re.search(r"([+-]/[+-])", text)
    sign_label = f" {signs.group(1)}" if signs else ""
    return f"{tf}{sign_label} ({genes.group(1)} targets)" if genes else f"{tf}{sign_label}"


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "label"


def save_pdf(fig, path):
    fig.savefig(path)
    plt.close(fig)
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)


def named(stem, ext=".pdf"):
    return f"{stem}_{args.file_suffix}{ext}" if args.file_suffix else f"{stem}{ext}"


def cluster_order(values, axis):
    matrix = values.astype(float).replace([np.inf, -np.inf], np.nan).fillna(0)
    if matrix.shape[axis] < 3:
        return list(matrix.index if axis == 0 else matrix.columns)
    data = matrix.values if axis == 0 else matrix.values.T
    if np.all(np.nanstd(data, axis=1) == 0):
        return list(matrix.index if axis == 0 else matrix.columns)
    return list((matrix.index if axis == 0 else matrix.columns)[leaves_list(linkage(pdist(data), method="average"))])


def cluster_frame(values, rows=True, cols=True):
    ordered = values.copy()
    if rows:
        ordered = ordered.loc[cluster_order(ordered, axis=0)]
    if cols:
        ordered = ordered.loc[:, cluster_order(ordered, axis=1)]
    return ordered


def resolve_condition_pair(conditions):
    conditions = [str(x) for x in conditions if str(x)]
    if len(conditions) != 2:
        return None, None, sorted(conditions)
    if args.reference_condition != "auto" and args.comparison_condition != "auto":
        if args.reference_condition not in conditions or args.comparison_condition not in conditions:
            raise ValueError("reference/comparison condition not found in sample metadata.")
        return args.reference_condition, args.comparison_condition, [args.reference_condition, args.comparison_condition]
    control_like = [x for x in conditions if re.search(r"(^|[_-])(ctrl|control|wt|wildtype|wild_type)($|[_-])", x, flags=re.I)]
    reference = control_like[0] if control_like else sorted(conditions)[0]
    if args.reference_condition != "auto":
        if args.reference_condition not in conditions:
            raise ValueError("reference condition not found in sample metadata.")
        reference = args.reference_condition
    remaining = [x for x in conditions if x != reference]
    comparison = remaining[0]
    if args.comparison_condition != "auto":
        if args.comparison_condition not in conditions:
            raise ValueError("comparison condition not found in sample metadata.")
        comparison = args.comparison_condition
    return reference, comparison, [reference, comparison]

def find_auc_mod(mdata, meta_cell_ids):
    candidates = []
    for name, ad in mdata.mod.items():
        obs_overlap = len(set(map(str, ad.obs_names)) & meta_cell_ids)
        candidates.append((obs_overlap, ad.n_vars, name, ad, False))
        var_overlap = len(set(map(str, ad.var_names)) & meta_cell_ids)
        candidates.append((var_overlap, ad.n_obs, name, ad, True))
    candidates.sort(reverse=True, key=lambda x: (x[0], x[1]))
    overlap, _, name, ad, transpose = candidates[0]
    if overlap == 0:
        raise ValueError("No AUCell modality has cell IDs overlapping metadata.")
    print("selected_auc_mod", name, "transpose", transpose, "overlap", overlap)
    df = ad.to_df()
    return df.T if transpose else df

meta = pd.read_csv(resolve_project_path(args.metadata), sep="\t", dtype=str).fillna("")
required_cols = [args.cell_col, args.sample_col, args.group_col]
if args.label_col:
    required_cols.append(args.label_col)
for col in required_cols:
    if col not in meta.columns:
        raise ValueError(f"metadata missing required column: {col}")
meta = meta.drop_duplicates(args.cell_col).set_index(args.cell_col)
meta_cell_ids = set(meta.index.astype(str))

mdata = md.read_h5mu(resolve_project_path(args.auc_h5mu))
auc = find_auc_mod(mdata, meta_cell_ids)
auc.index = auc.index.astype(str)
overlap = auc.index.intersection(meta.index)
if len(overlap) == 0:
    raise ValueError("AUCell matrix and metadata have no overlapping cells.")
auc = auc.loc[overlap]
meta = meta.loc[overlap]

def compute_condition_stats(sample_auc):
    reference, comparison, conditions = resolve_condition_pair(sample_auc[args.group_col].unique())
    rows = []
    regulons = [x for x in sample_auc.columns if x not in {args.group_col, args.label_col}]
    for regulon in regulons:
        values_by_condition = [
            sample_auc.loc[sample_auc[args.group_col] == condition, regulon].astype(float).values
            for condition in conditions
        ]
        row = {"eregulon": regulon, "n_conditions": len(conditions)}
        for condition, values in zip(conditions, values_by_condition):
            row[f"mean_auc__{condition}"] = float(np.mean(values)) if len(values) else np.nan
            row[f"n_samples__{condition}"] = int(len(values))
        if len(conditions) == 2:
            ref_values, comp_values = values_by_condition
            row["reference_condition"] = reference
            row["comparison_condition"] = comparison
            row["contrast"] = f"{comparison} - {reference}"
            row["delta_mean_auc"] = row[f"mean_auc__{comparison}"] - row[f"mean_auc__{reference}"]
            if len(ref_values) >= 2 and len(comp_values) >= 2:
                row["p_value"] = float(stats.ttest_ind(comp_values, ref_values, equal_var=False, nan_policy="omit").pvalue)
                row["test"] = "welch_t_test_on_sample_means"
            else:
                row["p_value"] = np.nan
                row["test"] = "not_tested_less_than_two_samples_per_condition"
        else:
            valid = [x for x in values_by_condition if len(x) > 0]
            row["delta_mean_auc"] = np.nan
            if len(conditions) < 2:
                row["contrast"] = "single_condition"
                row["p_value"] = np.nan
                row["test"] = "descriptive_only_single_condition"
            elif len(valid) >= 2 and all(len(x) >= 2 for x in valid):
                row["contrast"] = "multi_condition_omnibus"
                row["p_value"] = float(stats.kruskal(*valid).pvalue)
                row["test"] = "kruskal_wallis_on_sample_means"
            else:
                row["contrast"] = "multi_condition_omnibus"
                row["p_value"] = np.nan
                row["test"] = "not_tested_insufficient_samples"
        rows.append(row)
    res = pd.DataFrame(rows)
    res.insert(1, "display_label", res["eregulon"].map(regulon_display_label))
    res["fdr"] = bh_fdr(res["p_value"].values)
    sort_cols = ["fdr", "p_value"] if res["p_value"].notna().any() else ["eregulon"]
    return res.sort_values(sort_cols, na_position="last")


def plot_sample_counts(sample_auc, path, title):
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    counts = sample_auc[args.group_col].value_counts().sort_index()
    fallback = ["#2F4858", "#B57457", "#667A9A", "#9A8067", "#7D9477", "#8E5A63"]
    colors = {name: fallback[i % len(fallback)] for i, name in enumerate(counts.index.astype(str))}
    ax.bar(counts.index.astype(str), counts.values, color=[colors[x] for x in counts.index.astype(str)])
    ax.set_ylabel("samples")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=0)
    fig.subplots_adjust(left=0.18, right=0.95, bottom=0.20, top=0.82)
    save_pdf(fig, path)


def plot_volcano(res, path, title):
    fig, ax = plt.subplots(figsize=(4.2, 3.1))
    if "delta_mean_auc" in res and res["p_value"].notna().any():
        x = res["delta_mean_auc"].astype(float)
        y = -np.log10(res["p_value"].astype(float).clip(lower=np.nextafter(0, 1)))
        sig = res["fdr"].astype(float) < 0.05
        ax.scatter(x[~sig], y[~sig], s=9, alpha=0.50, color="#8F8F8F", linewidths=0)
        ax.scatter(x[sig], y[sig], s=12, alpha=0.88, color="#B57457", linewidths=0)
        for row in res.loc[sig].head(3).itertuples(index=False):
            ax.text(row.delta_mean_auc, -np.log10(max(row.p_value, np.nextafter(0, 1))), row.display_label, fontsize=5.0)
        contrast = res["contrast"].dropna().iloc[0] if "contrast" in res and res["contrast"].notna().any() else "condition contrast"
        ax.set_xlabel(f"delta mean AUC ({contrast})")
        ax.set_ylabel("-log10 P")
        ax.set_title(title)
    else:
        ax.text(0.5, 0.5, "Insufficient sample replicates for P values", ha="center", va="center")
        ax.axis("off")
    fig.subplots_adjust(left=0.17, right=0.96, bottom=0.18, top=0.86)
    save_pdf(fig, path)


def plot_sample_heatmap(sample_auc, res, path, title):
    top = res.head(min(25, res.shape[0]))["eregulon"].tolist()
    if not top:
        return
    mat = sample_auc.set_index(args.group_col, append=True)[top].T
    mat = cluster_frame(mat)
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    im = ax.imshow(mat.values, aspect="auto", cmap="Blues")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([regulon_display_label(x) for x in top], fontsize=5.3)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([f"{s}|{g}" for s, g in mat.columns], rotation=90, fontsize=5.8)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="AUC", shrink=0.62, pad=0.02)
    fig.subplots_adjust(left=0.34, right=0.90, bottom=0.30, top=0.90)
    save_pdf(fig, path)


pdf_index = 9


def numbered_pdf(stem):
    global pdf_index
    path = outdir / named(f"{pdf_index:02d}_{stem}")
    pdf_index += 1
    return path


sample_groups = meta[[args.sample_col, args.group_col]].drop_duplicates().set_index(args.sample_col)[args.group_col]
sample_auc = auc.groupby(meta[args.sample_col]).mean()
sample_auc.insert(0, args.group_col, sample_groups.loc[sample_auc.index].astype(str).values)
sample_auc.to_csv(outdir / named("sample_mean_auc", ".tsv"), sep="\t")
res = compute_condition_stats(sample_auc)
res.to_csv(outdir / named("condition_eregulon_auc_statistics", ".tsv"), sep="\t", index=False)
if not args.tables_only:
    plot_sample_counts(sample_auc, numbered_pdf("condition_sample_counts"), "Samples per condition")
    plot_volcano(res, numbered_pdf("condition_overall_eregulon_auc_volcano"), "Condition-level eRegulon AUC")
    plot_sample_heatmap(sample_auc, res, numbered_pdf("condition_overall_top_sample_mean_auc_heatmap"), "Top eRegulons: sample mean AUC")

if args.label_col:
    label_sample_tables = []
    label_stats = []
    label_plots = []
    for label in sorted(x for x in meta[args.label_col].unique() if x):
        label_cells = meta.index[meta[args.label_col] == label]
        label_sample_auc = auc.loc[label_cells].groupby(meta.loc[label_cells, args.sample_col]).mean()
        if label_sample_auc.empty:
            continue
        label_groups = meta.loc[label_cells, [args.sample_col, args.group_col]].drop_duplicates().set_index(args.sample_col)[args.group_col]
        label_sample_auc.insert(0, args.group_col, label_groups.loc[label_sample_auc.index].astype(str).values)
        label_sample_auc.insert(0, args.label_col, label)
        label_sample_tables.append(label_sample_auc)
        if label_sample_auc[args.group_col].nunique() < 2:
            continue
        label_res = compute_condition_stats(label_sample_auc)
        label_res.insert(0, args.label_col, label)
        label_stats.append(label_res)
        label_plots.append((label, label_res))
    if label_sample_tables:
        by_label_sample = pd.concat(label_sample_tables)
        by_label_sample.to_csv(outdir / named("by_cell_label_sample_mean_auc", ".tsv"), sep="\t")
    if label_stats:
        by_label = pd.concat(label_stats, ignore_index=True)
        by_label.to_csv(outdir / named("by_cell_label_condition_eregulon_auc_statistics", ".tsv"), sep="\t", index=False)
        top_regs = (
            by_label.assign(abs_delta=by_label["delta_mean_auc"].abs())
            .sort_values("abs_delta", ascending=False)
            .drop_duplicates("eregulon")
            .head(30)["eregulon"]
            .tolist()
        )
        mat = by_label.pivot_table(index="eregulon", columns=args.label_col, values="delta_mean_auc", aggfunc="mean").reindex(top_regs)
        mat = cluster_frame(mat)
        if not args.tables_only:
            fig, ax = plt.subplots(figsize=(4.9, 5.4))
            vmax = np.nanmax(np.abs(mat.values)) if np.isfinite(mat.values).any() else 1.0
            im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            ax.set_yticks(range(mat.shape[0]))
            ax.set_yticklabels([regulon_display_label(x) for x in mat.index], fontsize=5.1)
            ax.set_xticks(range(mat.shape[1]))
            ax.set_xticklabels(mat.columns.astype(str), rotation=45, ha="right", fontsize=5.8)
            ax.set_title("Condition effect by cell label")
            fig.colorbar(im, ax=ax, label="delta mean AUC", shrink=0.65, pad=0.02)
            fig.subplots_adjust(left=0.37, right=0.88, bottom=0.18, top=0.90)
            save_pdf(fig, numbered_pdf("condition_by_cell_label_effect_heatmap"))
            for label, label_res in label_plots:
                plot_volcano(label_res, numbered_pdf(f"condition_by_cell_label_eregulon_auc_volcano__{safe_name(label)}"), f"{label}: condition eRegulon AUC")

print(res.head(20).to_string(index=False))
