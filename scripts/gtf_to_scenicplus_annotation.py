import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--gtf", required=True)
parser.add_argument("--allowed-chroms", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()
allowed = set(open(args.allowed_chroms).read().splitlines())

def attr_value(attrs, keys):
    for key in keys:
        m = re.search(rf'{key} "([^"]+)"', attrs)
        if m:
            return m.group(1)
    return None

rows = []
with open(args.gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = fields
        if chrom not in allowed or feature != "transcript":
            continue
        biotype = attr_value(attrs, ["transcript_biotype", "gene_biotype", "transcript_type", "gene_type"])
        if biotype != "protein_coding":
            continue
        gene = attr_value(attrs, ["gene_name", "gene", "Name", "gene_id"])
        if not gene:
            continue
        start_i, end_i = int(start), int(end)
        tss = start_i if strand == "+" else end_i
        rows.append({
            "Chromosome": chrom,
            "Start": start_i,
            "End": end_i,
            "Strand": strand,
            "Gene": gene,
            "Transcription_Start_Site": tss,
            "Transcript_type": "protein_coding",
        })

df = pd.DataFrame(rows).drop_duplicates()
df.to_csv(args.out, sep="\t", index=False)
print(df.shape)
print(df.head())
