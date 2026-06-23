from pathlib import Path
import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--run-a", required=True)
parser.add_argument("--run-b", required=True)
parser.add_argument("--out", default="results/scenicplus_stability/stability_summary.tsv")
args = parser.parse_args()

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_project_path(path_value: str) -> Path:
    path = Path(path_value).expanduser()
    return path if path.is_absolute() else PROJECT / path


def edge_set(path):
    if not Path(path).exists():
        return set()
    df = pd.read_csv(path, sep="\t")
    lower = {c.lower(): c for c in df.columns}
    if {"tf", "gene"}.issubset(lower):
        cols = [lower["tf"], lower["gene"]]
    elif {"region", "gene"}.issubset(lower):
        cols = [lower["region"], lower["gene"]]
    elif {"target", "region"}.issubset(lower):
        cols = [lower["region"], lower["target"]]
    else:
        cols = list(df.columns[:min(3, df.shape[1])])
    return set(map(tuple, df[cols].astype(str).itertuples(index=False, name=None)))

rows = []
for name in ["eRegulons_direct.tsv", "eRegulons_extended.tsv", "region_to_gene_adj.tsv"]:
    a = edge_set(resolve_project_path(args.run_a) / name)
    b = edge_set(resolve_project_path(args.run_b) / name)
    union = len(a | b)
    rows.append({
        "file": name,
        "run_a_edges": len(a),
        "run_b_edges": len(b),
        "shared_edges": len(a & b),
        "jaccard": len(a & b) / union if union else 0,
    })
summary = pd.DataFrame(rows)
out = resolve_project_path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(out, sep="\t", index=False)
print(summary.to_string(index=False))
