from __future__ import annotations

from pathlib import Path
import os

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT / "tmp" / "numba"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)
(PROJECT / "tmp" / "numba").mkdir(parents=True, exist_ok=True)

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import io


INPUTS = PROJECT / "inputs"
WORK = PROJECT / "work" / "annotated_seurat"


def main() -> None:
    matrix_path = WORK / "rna_counts.genes_by_cells.mtx.gz"
    genes_path = WORK / "genes.tsv"
    cells_path = WORK / "cells.tsv"
    metadata_path = INPUTS / "cell_metadata.tsv"
    for path in [matrix_path, genes_path, cells_path, metadata_path]:
        if not path.exists():
            raise FileNotFoundError(path)
    matrix = io.mmread(matrix_path).tocsr().T
    genes = pd.read_csv(genes_path, sep="\t", header=None, names=["gene"], dtype=str)
    cells = pd.read_csv(cells_path, sep="\t", header=None, names=["cell_id"], dtype=str)
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna("")
    if matrix.shape != (cells.shape[0], genes.shape[0]):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match cells x genes "
            f"{cells.shape[0]} x {genes.shape[0]}"
        )
    if set(cells["cell_id"]) != set(metadata["cell_id"]):
        raise ValueError("RNA cell IDs and metadata cell IDs do not match.")
    metadata = metadata.set_index("cell_id").loc[cells["cell_id"]]
    adata = ad.AnnData(
        X=matrix,
        obs=metadata,
        var=pd.DataFrame(index=pd.Index(genes["gene"].astype(str), name="gene")),
    )
    adata.var_names_make_unique()
    raw_counts = adata.copy()
    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_vars == 0:
        adata = raw_counts.copy()
        sc.pp.filter_genes(adata, min_cells=1)
    if adata.n_vars == 0:
        raise ValueError("No expressed genes remain in the selected annotated Seurat cells.")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = raw_counts
    adata.write_h5ad(INPUTS / "gex.h5ad")
    print(adata)
    print("raw exists", adata.raw is not None)
    print("WROTE", INPUTS / "gex.h5ad")


if __name__ == "__main__":
    main()
