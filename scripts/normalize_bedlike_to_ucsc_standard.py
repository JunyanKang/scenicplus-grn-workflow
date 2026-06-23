import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--allowed-chroms", required=True)
parser.add_argument("--min-cols", type=int, default=3)
args = parser.parse_args()
allowed = set(open(args.allowed_chroms).read().splitlines())

def to_ucsc(chrom):
    if chrom in {"MT", "M", "Mt", "mitochondrion_genome"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else "chr" + chrom

for line in sys.stdin:
    if line.startswith("#") or line.strip() == "":
        continue
    fields = line.rstrip("\n").split("\t")
    if len(fields) < args.min_cols:
        continue
    fields[0] = to_ucsc(fields[0])
    if fields[0] not in allowed:
        continue
    print("\t".join(fields))
