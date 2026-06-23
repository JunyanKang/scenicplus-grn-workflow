from pathlib import Path
import argparse
import os

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import mudata as md
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument("--auc-h5mu", default="results/scenicplus/AUCell_direct.h5mu")
parser.add_argument("--metadata", default="inputs/cell_metadata.tsv")
parser.add_argument("--group-col", default="condition")
parser.add_argument("--sample-col", default="sample_id")
parser.add_argument("--cell-col", default="cell_id")
parser.add_argument("--outdir", default="results/scenicplus_stats/auc_by_condition")
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
for col in [args.cell_col, args.sample_col, args.group_col]:
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

sample_groups = (
    meta[[args.sample_col, args.group_col]]
    .drop_duplicates()
    .set_index(args.sample_col)[args.group_col]
)
sample_auc = auc.groupby(meta[args.sample_col]).mean()
sample_auc.insert(0, args.group_col, sample_groups.loc[sample_auc.index].astype(str).values)
sample_auc.to_csv(outdir / "sample_mean_auc.tsv", sep="\t")

conditions = sorted(sample_auc[args.group_col].unique())
rows = []
for regulon in sample_auc.columns.drop(args.group_col):
    values_by_condition = [
        sample_auc.loc[sample_auc[args.group_col] == condition, regulon].astype(float).values
        for condition in conditions
    ]
    row = {"eregulon": regulon, "n_conditions": len(conditions)}
    for condition, values in zip(conditions, values_by_condition):
        row[f"mean_auc__{condition}"] = float(np.mean(values)) if len(values) else np.nan
        row[f"n_samples__{condition}"] = int(len(values))
    if len(conditions) == 2:
        a, b = values_by_condition
        row["delta_mean_auc"] = row[f"mean_auc__{conditions[1]}"] - row[f"mean_auc__{conditions[0]}"]
        if len(a) >= 2 and len(b) >= 2:
            row["p_value"] = float(stats.ttest_ind(b, a, equal_var=False, nan_policy="omit").pvalue)
            row["test"] = "welch_t_test_on_sample_means"
        else:
            row["p_value"] = np.nan
            row["test"] = "not_tested_less_than_two_samples_per_condition"
    else:
        valid = [x for x in values_by_condition if len(x) > 0]
        row["delta_mean_auc"] = np.nan
        if len(valid) >= 2 and all(len(x) >= 2 for x in valid):
            row["p_value"] = float(stats.kruskal(*valid).pvalue)
            row["test"] = "kruskal_wallis_on_sample_means"
        else:
            row["p_value"] = np.nan
            row["test"] = "not_tested_insufficient_samples"
    rows.append(row)

res = pd.DataFrame(rows)
res["fdr"] = bh_fdr(res["p_value"].values)
sort_cols = ["fdr", "p_value"] if res["p_value"].notna().any() else ["eregulon"]
res = res.sort_values(sort_cols, na_position="last")
res.to_csv(outdir / "condition_eregulon_auc_statistics.tsv", sep="\t", index=False)

pdf = outdir / "condition_eregulon_auc_statistics.pdf"
with PdfPages(pdf) as pp:
    fig, ax = plt.subplots(figsize=(7, 4))
    counts = sample_auc[args.group_col].value_counts().sort_index()
    ax.bar(counts.index.astype(str), counts.values, color="#4E79A7")
    ax.set_ylabel("samples")
    ax.set_title("Samples per condition")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    if "delta_mean_auc" in res and res["p_value"].notna().any():
        x = res["delta_mean_auc"].astype(float)
        y = -np.log10(res["p_value"].astype(float).clip(lower=np.nextafter(0, 1)))
        ax.scatter(x, y, s=10, alpha=0.65, color="#333333", linewidths=0)
        ax.set_xlabel("delta mean AUC")
        ax.set_ylabel("-log10 P")
        ax.set_title("Condition-level eRegulon AUC")
    else:
        ax.text(0.5, 0.5, "Insufficient sample replicates for P values", ha="center", va="center")
        ax.axis("off")
    fig.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

    top = res.head(min(25, res.shape[0]))["eregulon"].tolist()
    if top:
        mat = sample_auc.set_index(args.group_col, append=True)[top].T
        fig, ax = plt.subplots(figsize=(max(7, 0.25 * mat.shape[1]), max(4, 0.18 * len(top))))
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top, fontsize=6)
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels([f"{s}|{g}" for s, g in mat.columns], rotation=90, fontsize=6)
        ax.set_title("Top eRegulons: sample mean AUC")
        fig.colorbar(im, ax=ax, label="AUC")
        fig.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

if not pdf.exists() or pdf.stat().st_size == 0:
    raise FileNotFoundError(pdf)
print(res.head(20).to_string(index=False))
