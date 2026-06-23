from pathlib import Path
import pandas as pd
from pybiomart import Server

human_motif2tf = "resources/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
orthology_cache = Path("resources/cyno/human_cynomolgus_orthology_biomart.tsv")
out_file = "inputs/cistarget_db/motif_annotations.tbl"
audit_file = "resources/cyno/cyno_motif2tf_orthology_audit.tsv"
required = [
    "external_gene_name",
    "mfascicularis_homolog_associated_gene_name",
    "mfascicularis_homolog_orthology_type",
    "mfascicularis_homolog_perc_id",
    "mfascicularis_homolog_perc_id_r1",
]

if orthology_cache.exists():
    raw_orth = pd.read_csv(orthology_cache, sep="\t")
else:
    server = Server(host="http://www.ensembl.org")
    mart = server["ENSEMBL_MART_ENSEMBL"]
    dataset = mart["hsapiens_gene_ensembl"]
    missing = [x for x in required if x not in dataset.attributes]
    if missing:
        available = [x for x in dataset.attributes if "mfascicularis_homolog" in x]
        raise ValueError("Missing BioMart attributes: " + ", ".join(missing) + "\n" + "\n".join(available[:80]))
    raw_orth = dataset.query(attributes=required)
    raw_orth.to_csv(orthology_cache, sep="\t", index=False)

raw_orth = raw_orth.copy()
raw_orth.columns = ["human_gene", "cyno_gene", "orthology_type", "human_perc_id", "cyno_perc_id"]
orth = raw_orth.dropna(subset=["human_gene", "cyno_gene"])
orth = orth.query("human_gene != '' and cyno_gene != ''")
orth = orth.query("orthology_type == 'ortholog_one2one'").copy()
orth = orth.sort_values(["human_gene", "human_perc_id", "cyno_perc_id"], ascending=[True, False, False])
orth = orth.drop_duplicates("human_gene")

motif = pd.read_csv(human_motif2tf, sep="\t")
if "gene_name" not in motif.columns:
    raise ValueError("Expected column gene_name not found in motif2TF table.")
converted = motif.merge(orth[["human_gene", "cyno_gene"]], left_on="gene_name", right_on="human_gene", how="inner")
converted["gene_name"] = converted["cyno_gene"]
converted = converted[motif.columns].drop_duplicates()
converted.to_csv(out_file, sep="\t", index=False)

summary = pd.DataFrame([
    {"metric": "human_motif_rows", "value": motif.shape[0]},
    {"metric": "human_unique_tfs", "value": motif["gene_name"].nunique()},
    {"metric": "raw_orthology_rows", "value": raw_orth.shape[0]},
    {"metric": "one_to_one_ortholog_genes", "value": orth.shape[0]},
    {"metric": "cynomolgus_motif_rows", "value": converted.shape[0]},
    {"metric": "cynomolgus_unique_tfs", "value": converted["gene_name"].nunique()},
])
summary.to_csv(audit_file, sep="\t", index=False)
print(summary.to_string(index=False))
print("orthology_cache", orthology_cache)
