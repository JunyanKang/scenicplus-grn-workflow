from pathlib import Path
from datetime import datetime
import os
import sys

PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "tmp" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT / "tmp" / "numba"))
(PROJECT / "tmp" / "matplotlib").mkdir(parents=True, exist_ok=True)
(PROJECT / "tmp" / "numba").mkdir(parents=True, exist_ok=True)

import pickle
import pandas as pd
import scanpy as sc
import yaml

if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
    print(
        "Usage: preflight_scenicplus_inputs.py\n\n"
        "Validate the project SCENIC+ inputs before running Snakemake. Requires "
        "$PROJECT_DIR with inputs/gex.h5ad, inputs/cell_metadata.tsv, "
        "inputs/cistopic_obj.pkl, region sets, cisTarget databases and motif annotations."
    )
    raise SystemExit(0)


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


log_dir = PROJECT / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
log_handle = (log_dir / f"preflight_scenicplus_inputs_{datetime.now():%Y%m%d_%H%M%S}.log").open("w")
sys.stdout = Tee(sys.stdout, log_handle)
sys.stderr = Tee(sys.stderr, log_handle)


def p(path: str) -> Path:
    out = Path(path).expanduser()
    return out if out.is_absolute() else PROJECT / out


threshold_table = pd.read_csv(p("inputs/preflight_thresholds.tsv"), sep="\t", dtype=str).fillna("")
key_col = "metric" if "metric" in threshold_table.columns else "parameter"
if key_col not in threshold_table.columns or "value" not in threshold_table.columns:
    raise ValueError("inputs/preflight_thresholds.tsv must contain parameter/value or metric/value columns")
thresholds = threshold_table.set_index(key_col)["value"].astype(float).to_dict()
required = [
    "inputs/gex.h5ad",
    "inputs/cell_metadata.tsv",
    "inputs/cistopic_obj.pkl",
    "inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather",
    "inputs/cistarget_db/custom.regions_vs_motifs.scores.feather",
    "inputs/cistarget_db/motif_annotations.tbl",
]
for f in required:
    assert p(f).is_file() and p(f).stat().st_size > 0, f"missing or empty: {p(f)}"
assert p("inputs/region_sets").is_dir(), "missing inputs/region_sets"
assert list(p("inputs/region_sets").glob("**/*.bed")), "no BED files in inputs/region_sets"

adata = sc.read_h5ad(p("inputs/gex.h5ad"))
assert adata.raw is not None, "gex.h5ad does not contain adata.raw"
assert adata.obs_names.is_unique, "RNA cell names are not unique"
assert adata.raw.n_vars >= adata.n_vars, "raw gene count should be >= X gene count"

meta = pd.read_csv(p("inputs/cell_metadata.tsv"), sep="\t")
for col in ["cell_id", "sample_id", "condition", "cell_label"]:
    assert col in meta.columns, f"inputs/cell_metadata.tsv missing column: {col}"
assert meta["cell_id"].is_unique, "cell_metadata cell_id values are not unique"
assert meta["cell_label"].notna().all(), "cell_label contains missing values"
meta_coverage = len(set(adata.obs_names) & set(meta["cell_id"])) / max(adata.n_obs, 1)
print("cell_metadata/RNA coverage", meta_coverage)
assert meta_coverage >= thresholds.get("min_cell_metadata_coverage", 0.95), "cell_metadata does not cover RNA cells"

cistopic = pickle.load(open(p("inputs/cistopic_obj.pkl"), "rb"))
if hasattr(cistopic, "cell_names"):
    atac = set(cistopic.cell_names)
elif hasattr(cistopic, "cell_data"):
    atac = set(cistopic.cell_data.index)
else:
    raise ValueError("Cannot find cell names in cisTopic object.")
rna = set(adata.obs_names)
overlap = len(rna & atac)
overlap_rna = overlap / max(len(rna), 1)
overlap_atac = overlap / max(len(atac), 1)
print("RNA", len(rna), "ATAC", len(atac), "overlap", overlap)
print("overlap/RNA", overlap_rna)
print("overlap/ATAC", overlap_atac)
assert overlap_rna >= thresholds.get("min_overlap_rna", 0.50), "RNA/ATAC barcode overlap too low for RNA side"
assert overlap_atac >= thresholds.get("min_overlap_atac", 0.50), "RNA/ATAC barcode overlap too low for ATAC side"

motif_dir = p("resources/motifs/v10nr_clust_public/singletons")
motifs = [x.strip() for x in p("resources/motifs/motifs.txt").read_text().splitlines() if x.strip()]
missing = [m for m in motifs if not (motif_dir / f"{m}.cb").exists()]
assert not missing, f"motif directory/list mismatch: {missing[:20]}"

motif_annot = pd.read_csv(p("inputs/cistarget_db/motif_annotations.tbl"), sep="\t")
if "gene_name" in motif_annot.columns:
    motif_tfs = set(motif_annot["gene_name"].dropna().astype(str))
    raw_genes = set(map(str, adata.raw.var_names))
    tf_overlap = len(motif_tfs & raw_genes) / max(len(motif_tfs), 1)
    print("motif TF/RNA raw gene overlap", tf_overlap)
    assert tf_overlap >= thresholds.get("min_motif_tf_gene_overlap", 0.50), "motif2TF gene names poorly match gex.h5ad.raw.var_names"

chromsizes_candidates = [
    p("work/scenicplus/chromsizes.tsv"),
    p(os.environ.get("CHROMSIZES", "resources/organism.ucsc.standard.chromsizes.tsv")),
]
chromsizes_path = next((x for x in chromsizes_candidates if x.exists()), None)
if chromsizes_path is None:
    raise FileNotFoundError(
        "Cannot locate active UCSC-standard chromsizes. Run prepare_official_resources.py "
        "and initialize_scenicplus_snakemake.py before preflight."
    )
allowed_chroms = set(pd.read_csv(chromsizes_path, sep="\t", header=None, dtype=str).iloc[:, 0])
for bed in p("inputs/region_sets").glob("**/*.bed"):
    with bed.open() as f:
        for i, line in enumerate(f):
            if i >= 1000:
                break
            if not line.strip() or line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            assert chrom in allowed_chroms, f"non-standard chromosome in {bed}: {chrom}"

cfg_path = p("work/scenicplus/Snakemake/config/config.yaml")
if cfg_path.exists():
    cfg = yaml.safe_load(cfg_path.read_text())
    assert cfg["params_data_preparation"]["bc_transform_func"] == '"lambda x: x"'
    assert float(cfg["params_motif_enrichment"]["dem_motif_hit_thr"]) == 3.0
    assert cfg["params_inference"]["tf_to_gene_importance_method"] == "GBM"
    assert cfg["params_inference"]["region_to_gene_importance_method"] == "GBM"
    assert cfg["params_inference"]["region_to_gene_correlation_method"] == "SR"
print("preflight OK")
