from __future__ import annotations

from pathlib import Path
import os
import sys

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT / "tmp" / "numba"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)
(PROJECT / "tmp" / "numba").mkdir(parents=True, exist_ok=True)

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse


INPUTS = PROJECT / "inputs"
RESULTS = PROJECT / "results" / "annotated_object"


def read_params(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return dict(zip(df["parameter"].astype(str), df["value"].astype(str)))


def main() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(
            "Usage: export_annotated_h5ad_for_scenicplus.py\n\n"
            "Export active RNA counts and metadata from an inspected AnnData h5ad "
            "object to inputs/cell_metadata.tsv and inputs/gex.h5ad. Requires "
            "$PROJECT_DIR/inputs/annotated_h5ad_params.tsv."
        )
        return
    INPUTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    params = read_params(INPUTS / "annotated_h5ad_params.tsv")
    path = Path(params["object_path"]).expanduser().resolve()
    adata = ad.read_h5ad(path)
    sample_col = params.get("sample_col", "sample")
    condition_col = params.get("condition_col", "condition")
    cell_label_column = params.get("cell_label_column", "")
    if not cell_label_column:
        raise ValueError("annotated_h5ad_params.tsv must define cell_label_column. Rerun inspect_annotated_object.py.")
    barcode_col = params.get("barcode_col", "")
    reduction = params.get("reduction", "")
    use_raw = params.get("use_raw", "true").lower() in {"1", "true", "yes", "on"}
    for col in [sample_col, condition_col, cell_label_column]:
        if col not in adata.obs.columns:
            raise ValueError(f"AnnData obs lacks required column: {col}")
    obs_names = pd.Index(adata.obs_names.astype(str))
    sample_id = adata.obs[sample_col].astype(str).to_numpy()
    condition = adata.obs[condition_col].astype(str).to_numpy()
    cell_label = adata.obs[cell_label_column].astype(str).to_numpy()
    if barcode_col and barcode_col in adata.obs.columns:
        barcode = adata.obs[barcode_col].astype(str).to_numpy()
    else:
        barcode = obs_names.to_numpy().copy()
        for sid in sorted(set(sample_id)):
            mask = sample_id == sid
            barcode[mask] = pd.Series(barcode[mask]).str.replace(f"^{sid}_", "", regex=True).to_numpy()
    cell_id = pd.Index([f"{bc}-{sid}" for bc, sid in zip(barcode, sample_id)])
    if cell_id.has_duplicates:
        raise ValueError("Duplicated cell_id values after barcode-sample_id conversion.")
    cell_meta = pd.DataFrame(
        {
            "cell_id": cell_id,
            "original_cell_id": obs_names.to_numpy(),
            "barcode": barcode,
            "sample_id": sample_id,
            "condition": condition,
            "cell_label": cell_label,
            "source_label": cell_label,
            "analysis_unit": "cell",
        }
    )
    if reduction and reduction in adata.obsm:
        emb = adata.obsm[reduction]
        for i in range(emb.shape[1]):
            cell_meta[f"{reduction}_{i + 1}"] = emb[:, i]

    matrix_source = adata.raw if use_raw and adata.raw is not None else adata
    X = matrix_source.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    out = ad.AnnData(
        X=X.copy(),
        obs=cell_meta.set_index("cell_id"),
        var=pd.DataFrame(index=pd.Index(matrix_source.var_names.astype(str), name="gene")),
    )
    out.var_names_make_unique()
    raw_counts = out.copy()
    sc.pp.filter_genes(out, min_cells=3)
    if out.n_vars == 0:
        out = raw_counts.copy()
        sc.pp.filter_genes(out, min_cells=1)
    if out.n_vars == 0:
        raise ValueError("No expressed genes remain in the selected annotated h5ad cells.")
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    out.raw = raw_counts
    cell_meta.to_csv(INPUTS / "cell_metadata.tsv", sep="\t", index=False)
    out.write_h5ad(INPUTS / "gex.h5ad")
    pd.crosstab(cell_meta["sample_id"], cell_meta["cell_label"]).to_csv(
        RESULTS / "cells_by_sample_and_label.tsv",
        sep="\t",
    )
    print("WROTE", INPUTS / "cell_metadata.tsv")
    print("WROTE", INPUTS / "gex.h5ad")


if __name__ == "__main__":
    main()
