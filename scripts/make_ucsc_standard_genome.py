import argparse
import gzip
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--fasta-gz", required=True)
parser.add_argument("--gtf-gz", required=True)
parser.add_argument("--allowed-chroms", required=True)
parser.add_argument("--out-prefix", required=True)
args = parser.parse_args()

allowed = set(Path(args.allowed_chroms).read_text().splitlines())

def to_ucsc(chrom):
    chrom = chrom.split()[0]
    if chrom in {"MT", "M", "Mt", "mitochondrion_genome"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else "chr" + chrom

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

out_fasta = Path(args.out_prefix + ".ucsc.standard.fa")
out_gtf = Path(args.out_prefix + ".ucsc.standard.gtf")
out_chromsizes = Path(args.out_prefix + ".ucsc.standard.chromsizes.tsv")

with open_text(args.fasta_gz) as fin, out_fasta.open("w") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            raw = line[1:].strip().split()[0]
            ucsc = to_ucsc(raw)
            keep = ucsc in allowed
            if keep:
                fout.write(f">{ucsc}\n")
        elif keep:
            fout.write(line)

chrom = None
length = 0
rows = []
with out_fasta.open() as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if chrom is not None:
                rows.append((chrom, 0, length))
            chrom = line[1:].split()[0]
            length = 0
        else:
            length += len(line)
    if chrom is not None:
        rows.append((chrom, 0, length))
with out_chromsizes.open("w") as out:
    out.write("Chromosome\tStart\tEnd\n")
    for chrom, start, end in rows:
        out.write(f"{chrom}\t{start}\t{end}\n")

with open_text(args.gtf_gz) as fin, out_gtf.open("w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        fields[0] = to_ucsc(fields[0])
        if fields[0] not in allowed:
            continue
        fout.write("\t".join(fields) + "\n")

print(out_fasta)
print(out_gtf)
print(out_chromsizes)
