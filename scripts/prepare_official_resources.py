#!/usr/bin/env python
"""Prepare and validate official genome/motif resources for SCENIC+ projects."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd

MODULES_DIR = Path(__file__).resolve().parents[1] / "modules"
sys.path.insert(0, str(MODULES_DIR))

from organism_resources import (  # noqa: E402
    DEFAULT_REF_BY_SPECIES,
    DIRECT_MOTIF2TF_BY_SPECIES,
    infer_assembly_from_url,
    infer_biomart_orthology_prefix,
    listed_ensembl_species,
    motif2tf_plan,
    pick_ensembl_file,
    resolve_organism,
)

ENSEMBL_RELEASE = "115"
MOTIF_COLLECTION = "v10nr_clust_public"
MOTIF2TF_VERSION = "v10nr_clust"
DIRECT_MOTIF2TF = {
    "human": {
        "url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "biomart_host": "http://www.ensembl.org",
        "mart": "ENSEMBL_MART_ENSEMBL",
        "dataset": "hsapiens_gene_ensembl",
        "orthology_prefix": "hsapiens",
    },
    "mouse": {
        "url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        "file": "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        "biomart_host": "http://www.ensembl.org",
        "mart": "ENSEMBL_MART_ENSEMBL",
        "dataset": "mmusculus_gene_ensembl",
        "orthology_prefix": "mmusculus",
    },
    "fly": {
        "url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl",
        "file": "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl",
        "biomart_host": "http://metazoa.ensembl.org",
        "mart": "ENSEMBL_MART_ENSEMBL",
        "dataset": "dmelanogaster_gene_ensembl",
        "orthology_prefix": "dmelanogaster",
    },
    "chicken": {
        "url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.chicken-m0.001-o0.0.tbl",
        "file": "motifs-v10nr_clust-nr.chicken-m0.001-o0.0.tbl",
        "biomart_host": "http://www.ensembl.org",
        "mart": "ENSEMBL_MART_ENSEMBL",
        "dataset": "ggallus_gene_ensembl",
        "orthology_prefix": "ggallus",
    },
}
SPECIES = {
    "human": {
        "organism": "homo_sapiens",
        "assembly": "GRCh38",
        "prefix": "human",
        "chroms": [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "direct",
    },
    "mouse": {
        "organism": "mus_musculus",
        "assembly": "GRCm39",
        "prefix": "mouse",
        "chroms": [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/mus_musculus/Mus_musculus.GRCm39.{release}.gtf.gz",
        "motif2tf_reference": "mouse",
        "motif2tf_mode": "direct",
    },
    "cyno": {
        "organism": "macaca_fascicularis",
        "assembly": "Macaca_fascicularis_6.0",
        "prefix": "cyno",
        "chroms": [f"chr{i}" for i in range(1, 21)] + ["chrX"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "mfascicularis",
        "orthology_label": "cynomolgus",
    },
    "rhesus": {
        "organism": "macaca_mulatta",
        "assembly": "Mmul_10",
        "prefix": "rhesus",
        "chroms": [f"chr{i}" for i in range(1, 21)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "mmulatta",
        "orthology_label": "rhesus",
    },
    "rat": {
        "organism": "rattus_norvegicus",
        "assembly": "GRCr8",
        "prefix": "rat",
        "chroms": [f"chr{i}" for i in range(1, 21)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.{release}.gtf.gz",
        "motif2tf_reference": "mouse",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "rnorvegicus",
        "orthology_label": "rat",
    },
    "rabbit": {
        "organism": "oryctolagus_cuniculus",
        "assembly": "OryCun2.0",
        "prefix": "rabbit",
        "chroms": [f"chr{i}" for i in range(1, 22)] + ["chrX"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "ocuniculus",
        "orthology_label": "rabbit",
    },
    "pig": {
        "organism": "sus_scrofa",
        "assembly": "Sscrofa11.1",
        "prefix": "pig",
        "chroms": [f"chr{i}" for i in range(1, 19)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "sscrofa",
        "orthology_label": "pig",
    },
    "cow": {
        "organism": "bos_taurus",
        "assembly": "ARS-UCD2.0",
        "prefix": "cow",
        "chroms": [f"chr{i}" for i in range(1, 30)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/bos_taurus/Bos_taurus.ARS-UCD2.0.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "btaurus",
        "orthology_label": "cow",
    },
    "dog": {
        "organism": "canis_lupus_familiaris",
        "assembly": "ROS_Cfam_1.0",
        "prefix": "dog",
        "chroms": [f"chr{i}" for i in range(1, 39)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.{release}.gtf.gz",
        "motif2tf_reference": "human",
        "motif2tf_mode": "mapped",
        "orthology_prefix": "clfamiliaris",
        "orthology_label": "dog",
    },
    "chicken": {
        "organism": "gallus_gallus",
        "assembly": "bGalGal1.mat.broiler.GRCg7b",
        "prefix": "chicken",
        "chroms": [f"chr{i}" for i in range(1, 40)] + ["chrZ", "chrW"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.{release}.gtf.gz",
        "motif2tf_reference": "chicken",
        "motif2tf_mode": "direct",
    },
    "zebrafish": {
        "organism": "danio_rerio",
        "assembly": "GRCz11",
        "prefix": "zebrafish",
        "chroms": [f"chr{i}" for i in range(1, 26)],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/danio_rerio/Danio_rerio.GRCz11.{release}.gtf.gz",
        "motif2tf_reference": None,
        "motif2tf_mode": "restricted",
        "orthology_prefix": "drerio",
        "orthology_label": "zebrafish",
        "restricted_reason": "teleost genome duplication makes simple one-to-one motif2TF transfer incomplete; use a curated species-specific table or explicit exploratory mapping.",
    },
}
MOTIF_COLLECTION_URL = "https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/"
CURATED_BY_ENSEMBL = {cfg["organism"]: cfg for cfg in SPECIES.values()}
RESOLVED_SPECIES: dict[str, dict] = {}


def resolve_species_config(query: str, release: str, motif2tf_reference: str = "auto", target_orthology_prefix: str = "") -> tuple[str, dict]:
    species = resolve_organism(query, release)
    if species in RESOLVED_SPECIES:
        return species, RESOLVED_SPECIES[species]
    curated = CURATED_BY_ENSEMBL.get(species, {})
    fasta_url = curated.get("fasta_url", "").format(release=release) if curated.get("fasta_url") else pick_ensembl_file(release, species, "fasta")
    gtf_url = curated.get("gtf_url", "").format(release=release) if curated.get("gtf_url") else pick_ensembl_file(release, species, "gtf")
    if not fasta_url or not gtf_url:
        raise ValueError(
            f"Could not resolve Ensembl FASTA/GTF for organism={query!r} "
            f"(resolved as {species!r}) in release {release}. "
            "Run spgrn-query-organism-resources --list first."
        )
    motif_status, motif_strategy, motif_ref = motif2tf_plan(species, motif2tf_reference)
    motif_mode = "direct" if motif_strategy == "direct" else ("mapped" if motif_ref else "restricted")
    cfg = {
        "organism": species,
        "assembly": curated.get("assembly") or infer_assembly_from_url(fasta_url, release),
        "prefix": species,
        "chroms": curated.get("chroms", []),
        "fasta_url": fasta_url,
        "gtf_url": gtf_url,
        "motif2tf_reference": motif_ref,
        "motif2tf_mode": motif_mode,
        "orthology_prefix": curated.get("orthology_prefix") or target_orthology_prefix or infer_biomart_orthology_prefix(species),
        "orthology_label": curated.get("orthology_label") or species.replace("_", " "),
        "restricted_reason": curated.get(
            "restricted_reason",
            "no direct Aerts motif2TF table is available for this species; use --ref human|mouse|fly|chicken or MOTIF2TF_TABLE.",
        ),
        "motif2tf_status": motif_status,
    }
    RESOLVED_SPECIES[species] = cfg
    return species, cfg


def cfg_for(species: str) -> dict:
    if species in RESOLVED_SPECIES:
        return RESOLVED_SPECIES[species]
    if species in CURATED_BY_ENSEMBL:
        cfg = dict(CURATED_BY_ENSEMBL[species])
        cfg["prefix"] = species
        RESOLVED_SPECIES[species] = cfg
        return cfg
    if species in SPECIES:
        cfg = dict(SPECIES[species])
        cfg["prefix"] = cfg["organism"]
        RESOLVED_SPECIES[cfg["organism"]] = cfg
        return cfg
    raise KeyError(f"Unresolved organism: {species}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", default=os.environ.get("ORGANISM"), help="Ensembl species name or common alias, for example mus_musculus or mouse.")
    parser.add_argument("--resources-dir", default=None, help="Default: $PROJECT_DIR/resources")
    parser.add_argument("--inputs-dir", default=None, help="Default: $PROJECT_DIR/inputs")
    parser.add_argument("--ensembl-release", default=os.environ.get("ENSEMBL_RELEASE", ENSEMBL_RELEASE))
    parser.add_argument("--mode", choices=["prepare", "check", "status"], default="prepare")
    parser.add_argument("--force", action="store_true", help="Regenerate derived files and re-download direct resources.")
    parser.add_argument("--skip-motifs", action="store_true", help="Skip motif collection and motif2TF preparation.")
    parser.add_argument("--motifs-only", action="store_true", help="Prepare only motif collection and motif2TF files for the selected organism.")
    parser.add_argument(
        "--motif2tf-reference",
        choices=["auto"] + sorted(DIRECT_MOTIF2TF),
        default=os.environ.get("MOTIF2TF_REFERENCE", "auto"),
        help="Reference Aerts motif2TF table for orthology mapping. Direct species ignore this unless a non-matching reference is explicit.",
    )
    parser.add_argument(
        "--ref",
        choices=sorted(DIRECT_MOTIF2TF),
        default=None,
        help="Short alias for --motif2tf-reference human|mouse|fly|chicken.",
    )
    parser.add_argument(
        "--orthology-policy",
        choices=["one2one", "paralog-aware"],
        default=os.environ.get("ORTHOLOGY_POLICY", "one2one"),
        help="Orthology mapping policy for non-direct motif2TF preparation.",
    )
    parser.add_argument("--install", action="store_true", help="Accepted for clarity; single-organism runs always install motif_annotations.tbl.")
    parser.add_argument(
        "--motif2tf-table",
        default=os.environ.get("MOTIF2TF_TABLE", ""),
        help="User-supplied species-specific motif2TF table. If set, it is copied to inputs/cistarget_db/motif_annotations.tbl.",
    )
    parser.add_argument(
        "--generate-motif2tf",
        action="store_true",
        help=(
            "Generate an audited species-specific motif2TF table by matching target gene symbols "
            "against Aerts direct human/mouse/fly/chicken motif2TF evidence."
        ),
    )
    parser.add_argument(
        "--target-orthology-prefix",
        default=os.environ.get("TARGET_ORTHOLOGY_PREFIX", ""),
        help="Advanced: target BioMart orthology prefix for explicit --ref mapping when automatic inference is insufficient.",
    )
    parser.add_argument("--log", default=None, help="Log file path. Default: logs/prepare_official_resources_<timestamp>.log")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_project_path(path_value: str | None, default_rel: str, base: Path) -> Path:
    path = Path(path_value).expanduser() if path_value else base / default_rel
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def setup_logging(log_path: str | None, base: Path) -> Path:
    if log_path is None:
        (base / "logs").mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = base / "logs" / f"prepare_official_resources_{stamp}.log"
    else:
        path = Path(log_path).expanduser()
        if not path.is_absolute():
            path = (base / path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(path), logging.StreamHandler(sys.stdout)],
    )
    return path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except Exception:
        return False


def file_ok(path: Path, min_size: int = 1) -> bool:
    return path.exists() and path.stat().st_size >= min_size


def run(cmd: list[str]) -> None:
    logging.info("RUN %s", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required command not found in PATH: {name}")
    return path


def curl_download(url: str, out: Path, force: bool = False) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0 and not force:
        if out.suffix == ".gz" and not gzip_ok(out):
            logging.warning("Existing gzip failed validation, re-downloading: %s", out)
        else:
            logging.info("SKIP existing download: %s", out)
            return
    require_tool("curl")
    run(["curl", "-L", "--retry", "3", "--retry-delay", "5", "-C", "-", "-o", str(out), url])
    if out.suffix == ".gz" and not gzip_ok(out):
        raise ValueError(f"Downloaded gzip failed validation: {out}")
    if not file_ok(out):
        raise ValueError(f"Downloaded file is empty: {out}")


def write_chrom_allowlist(resources: Path, species: str) -> Path:
    cfg = cfg_for(species)
    chrom_dir = resources / "chromosomes"
    chrom_dir.mkdir(parents=True, exist_ok=True)
    out = chrom_dir / f"{species}.ucsc.standard.chroms.txt"
    chroms = cfg.get("chroms") or []
    if not chroms:
        fasta_gz = resources / species / "genome.ensembl.fa.gz"
        chroms = infer_standard_chroms_from_fasta(fasta_gz, resources / species / f"{species}_chromosome_audit.tsv")
        cfg["chroms"] = chroms
    out.write_text("\n".join(chroms) + "\n")
    return out


def to_ucsc(chrom: str) -> str:
    chrom = chrom.split()[0]
    if chrom in {"MT", "M", "Mt", "mitochondrion_genome"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else "chr" + chrom


def infer_standard_chroms_from_fasta(fasta_gz: Path, audit_path: Path) -> list[str]:
    if not fasta_gz.exists():
        raise FileNotFoundError(f"Cannot infer chromosomes; FASTA is missing: {fasta_gz}")
    rows = []
    keep_chroms: list[str] = []
    bad_tokens = ("_", ".", "GL", "KI", "JH", "NW_", "random", "alt", "hap", "fix", "patch", "scaffold", "unplaced")
    with open_text(fasta_gz) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            raw = line[1:].strip().split()[0]
            ucsc = to_ucsc(raw)
            simple = re.fullmatch(r"(chr)?([0-9]{1,2}|[XYWZ]|M|MT)", raw, flags=re.IGNORECASE)
            blocked = any(token.lower() in raw.lower() for token in bad_tokens)
            keep = bool(simple) and not blocked
            reason = "standard_candidate" if keep else "excluded_non_standard_or_ambiguous"
            rows.append({"raw_chrom": raw, "ucsc_chrom": ucsc, "keep": keep, "reason": reason})
            if keep:
                keep_chroms.append(ucsc)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(audit_path, sep="\t", index=False)
    keep_chroms = sorted(set(keep_chroms), key=lambda x: (not x[3:].isdigit(), int(x[3:]) if x[3:].isdigit() else x))
    if not keep_chroms:
        raise ValueError(f"No standard chromosome candidates inferred from {fasta_gz}; audit written to {audit_path}")
    return keep_chroms


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def convert_genome(species: str, resources: Path, force: bool = False) -> dict[str, Path]:
    cfg = cfg_for(species)
    prefix = cfg["prefix"]
    species_dir = resources / species
    allowed_path = write_chrom_allowlist(resources, species)
    allowed = set(allowed_path.read_text().splitlines())
    fasta_gz = species_dir / "genome.ensembl.fa.gz"
    gtf_gz = species_dir / "genes.ensembl.gtf.gz"
    out_fasta = species_dir / f"{prefix}.ucsc.standard.fa"
    out_gtf = species_dir / f"{prefix}.ucsc.standard.gtf"
    out_chromsizes = species_dir / f"{prefix}.ucsc.standard.chromsizes.tsv"
    if out_fasta.exists() and out_gtf.exists() and out_chromsizes.exists() and not force:
        logging.info("SKIP existing UCSC-converted genome resources for %s", species)
        return {"fasta": out_fasta, "gtf": out_gtf, "chromsizes": out_chromsizes, "allowed": allowed_path}

    with open_text(fasta_gz) as fin, out_fasta.open("w") as fout:
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
    with out_fasta.open() as handle:
        for line in handle:
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

    with open_text(gtf_gz) as fin, out_gtf.open("w") as fout:
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

    if shutil.which("samtools"):
        run(["samtools", "faidx", str(out_fasta)])
    else:
        logging.warning("samtools not found; FASTA .fai index was not created.")
    return {"fasta": out_fasta, "gtf": out_gtf, "chromsizes": out_chromsizes, "allowed": allowed_path}


def attr_value(attrs: str, keys: Iterable[str]) -> str | None:
    for key in keys:
        m = re.search(rf'{key} "([^"]+)"', attrs)
        if m:
            return m.group(1)
    return None


def build_annotation(species: str, resources: Path, force: bool = False) -> Path:
    cfg = cfg_for(species)
    prefix = cfg["prefix"]
    species_dir = resources / species
    gtf = species_dir / f"{prefix}.ucsc.standard.gtf"
    allowed = set((resources / "chromosomes" / f"{species}.ucsc.standard.chroms.txt").read_text().splitlines())
    out_file = species_dir / f"{prefix}.ucsc.standard.genome_annotation.tsv"
    if out_file.exists() and out_file.stat().st_size > 0 and not force:
        logging.info("SKIP existing genome annotation: %s", out_file)
        return out_file
    rows = []
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            if chrom not in allowed or feature != "transcript":
                continue
            biotype = attr_value(attrs, ["transcript_biotype", "gene_biotype", "transcript_type", "gene_type"])
            if biotype != "protein_coding":
                continue
            gene = attr_value(attrs, ["gene_name", "gene", "Name", "gene_id"])
            if not gene:
                continue
            start_i, end_i = int(start), int(end)
            rows.append(
                {
                    "Chromosome": chrom,
                    "Start": start_i,
                    "End": end_i,
                    "Strand": strand,
                    "Gene": gene,
                    "Transcription_Start_Site": start_i if strand == "+" else end_i,
                    "Transcript_type": "protein_coding",
                }
            )
    df = pd.DataFrame(rows).drop_duplicates()
    if df.empty:
        raise ValueError(f"No protein-coding transcript annotations written for {species}")
    df.to_csv(out_file, sep="\t", index=False)
    return out_file


def prepare_motif_collection(resources: Path, force: bool = False) -> Path:
    motifs_root = resources / "motifs"
    singletons = motifs_root / MOTIF_COLLECTION / "singletons"
    motifs_txt = motifs_root / "motifs.txt"
    if singletons.exists() and list(singletons.glob("*.cb")) and motifs_txt.exists() and not force:
        logging.info("SKIP existing motif collection: %s", singletons)
        return singletons
    motifs_root.mkdir(parents=True, exist_ok=True)
    collection_dir = motifs_root / MOTIF_COLLECTION
    zip_path = collection_dir / f"{MOTIF_COLLECTION}.zip"
    curl_download(MOTIF_COLLECTION_URL + f"{MOTIF_COLLECTION}.zip", zip_path, force=force)
    singletons.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        cb_members = [name for name in archive.namelist() if name.endswith(".cb") and "/singletons/" in name]
        if not cb_members:
            raise FileNotFoundError(f"No singleton .cb motif files found in {zip_path}")
        for member in cb_members:
            out_file = singletons / Path(member).name
            if out_file.exists() and out_file.stat().st_size > 0 and not force:
                continue
            with archive.open(member) as src, out_file.open("wb") as dst:
                shutil.copyfileobj(src, dst)
    cb_files = sorted(singletons.glob("*.cb"))
    if not cb_files:
        raise FileNotFoundError(f"No .cb motif files found under {singletons}")
    motifs_txt.write_text("\n".join(p.stem for p in cb_files) + "\n")
    (motifs_root / "motifs.txt.sha256").write_text(f"{sha256(motifs_txt)}  motifs.txt\n")
    return singletons


def resolve_reference(species: str, requested: str) -> str | None:
    cfg = cfg_for(species)
    if requested != "auto":
        return requested
    return cfg.get("motif2tf_reference")


def reference_source(resources: Path, reference: str, force: bool = False) -> Path:
    ref = DIRECT_MOTIF2TF[reference]
    motif2tf_dir = resources / "motif2tf"
    motif2tf_dir.mkdir(parents=True, exist_ok=True)
    source = motif2tf_dir / ref["file"]
    curl_download(ref["url"], source, force=force)
    return source


def install_motif2tf(species_out: Path, inputs: Path) -> None:
    inputs_db = inputs / "cistarget_db"
    inputs_db.mkdir(parents=True, exist_ok=True)
    target = inputs_db / "motif_annotations.tbl"
    for conflict in inputs_db.glob("motif_annotations*_Conflict.tbl"):
        conflict.unlink()
    tmp_target = inputs_db / f".motif_annotations.{os.getpid()}.tmp"
    with species_out.open("rb") as src, tmp_target.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    os.replace(tmp_target, target)


def motif2tf_output_path(resources: Path, species: str, strategy: str, detail: str = "") -> Path:
    safe_detail = re.sub(r"[^A-Za-z0-9_.-]+", "_", detail).strip("_")
    suffix = f".{safe_detail}" if safe_detail else ""
    return resources / "motif2tf" / f"motif_annotations.{species}.{strategy}{suffix}.tbl"


def write_active_motif2tf_pointer(resources: Path, species: str, strategy: str, source: Path, installed_to_inputs: bool) -> None:
    pointer = resources / "motif2tf" / f"motif_annotations.{species}.active.tsv"
    pointer.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "species": species,
                "strategy": strategy,
                "source_table": str(source),
                "installed_as": "inputs/cistarget_db/motif_annotations.tbl" if installed_to_inputs else "",
            }
        ]
    ).to_csv(pointer, sep="\t", index=False)


def use_user_motif2tf_table(species: str, resources: Path, inputs: Path, table: str, install_to_inputs: bool = True) -> Path:
    source = Path(table).expanduser().resolve()
    if not source.exists() or source.stat().st_size == 0:
        raise FileNotFoundError(f"User-supplied motif2TF table is missing or empty: {source}")
    df = pd.read_csv(source, sep="\t", nrows=5)
    if "gene_name" not in df.columns:
        raise ValueError("User-supplied motif2TF table must contain a gene_name column.")
    motif2tf_dir = resources / "motif2tf"
    motif2tf_dir.mkdir(parents=True, exist_ok=True)
    species_out = motif2tf_output_path(resources, species, "user_supplied")
    shutil.copy2(source, species_out)
    audit = resources / species / f"{species}_motif2tf_user_table_audit.tsv"
    audit.parent.mkdir(parents=True, exist_ok=True)
    full = pd.read_csv(source, sep="\t")
    pd.DataFrame(
        [
            {"metric": "strategy", "value": "user_supplied_species_specific_table"},
            {"metric": "source_path", "value": str(source)},
            {"metric": "motif_rows", "value": full.shape[0]},
            {"metric": "unique_tfs", "value": full["gene_name"].nunique()},
        ]
    ).to_csv(audit, sep="\t", index=False)
    if install_to_inputs:
        install_motif2tf(species_out, inputs)
    write_active_motif2tf_pointer(resources, species, "user_supplied", species_out, install_to_inputs)
    return species_out


def generate_motif2tf_from_symbol_evidence(
    species: str,
    resources: Path,
    inputs: Path,
    force: bool = False,
    install_to_inputs: bool = True,
) -> Path:
    cfg = cfg_for(species)
    prefix = cfg["prefix"]
    annotation = resources / species / f"{prefix}.ucsc.standard.genome_annotation.tsv"
    if not annotation.exists() or annotation.stat().st_size == 0:
        raise FileNotFoundError(
            f"Genome annotation is required before generating motif2TF by symbol evidence: {annotation}. "
            "Run spgrn-prepare-official-resources without --motifs-only first."
        )
    motif2tf_dir = resources / "motif2tf"
    motif2tf_dir.mkdir(parents=True, exist_ok=True)
    out_file = motif2tf_output_path(resources, species, "generated", "symbol_evidence")
    species_dir = resources / species
    species_dir.mkdir(parents=True, exist_ok=True)
    audit_file = species_dir / f"{species}_motif2tf_generated_symbol_audit.tsv"
    unmatched_file = species_dir / f"{species}_motif2tf_generated_unmatched_reference_tfs.tsv"
    sources_file = species_dir / f"{species}_motif2tf_generated_sources.tsv"
    if out_file.exists() and audit_file.exists() and not force:
        logging.info("SKIP existing generated motif2TF table: %s", out_file)
        if install_to_inputs:
            install_motif2tf(out_file, inputs)
        return out_file

    ann = pd.read_csv(annotation, sep="\t")
    if "Gene" not in ann.columns:
        raise ValueError(f"Genome annotation lacks Gene column: {annotation}")
    target_genes = set(ann["Gene"].dropna().astype(str))
    if not target_genes:
        raise ValueError(f"No target gene symbols found in genome annotation: {annotation}")

    frames = []
    source_rows = []
    for reference in sorted(DIRECT_MOTIF2TF):
        source = reference_source(resources, reference, force=force)
        table = pd.read_csv(source, sep="\t")
        if "gene_name" not in table.columns:
            raise ValueError(f"Expected gene_name column in motif2TF table: {source}")
        table = table.copy()
        table["_reference"] = reference
        table["_source_path"] = str(source)
        frames.append(table)
        source_rows.append(
            {
                "reference": reference,
                "source_path": str(source),
                "rows": table.shape[0],
                "unique_tfs": table["gene_name"].nunique(),
            }
        )
    all_ref = pd.concat(frames, ignore_index=True)
    all_ref["gene_name"] = all_ref["gene_name"].astype(str)
    generated = all_ref.loc[all_ref["gene_name"].isin(target_genes)].copy()
    unmatched = sorted(set(all_ref["gene_name"]) - target_genes)
    pd.DataFrame({"reference_gene": unmatched}).to_csv(unmatched_file, sep="\t", index=False)
    pd.DataFrame(source_rows).to_csv(sources_file, sep="\t", index=False)

    source_cols = ["_reference", "_source_path"]
    generated_out = generated.drop(columns=[c for c in source_cols if c in generated.columns]).drop_duplicates()
    summary = pd.DataFrame(
        [
            {"metric": "strategy", "value": "generated_from_aerts_direct_symbol_evidence"},
            {"metric": "target_species", "value": species},
            {"metric": "target_protein_coding_genes", "value": len(target_genes)},
            {"metric": "reference_sources", "value": ",".join(sorted(DIRECT_MOTIF2TF))},
            {"metric": "reference_rows_total", "value": all_ref.shape[0]},
            {"metric": "reference_unique_tfs_total", "value": all_ref["gene_name"].nunique()},
            {"metric": "matched_rows", "value": generated_out.shape[0]},
            {"metric": "matched_unique_tfs", "value": generated_out["gene_name"].nunique() if "gene_name" in generated_out.columns else 0},
            {
                "metric": "matched_tf_fraction_of_reference",
                "value": round(generated_out["gene_name"].nunique() / all_ref["gene_name"].nunique(), 4)
                if all_ref["gene_name"].nunique()
                else 0,
            },
            {"metric": "audit_note", "value": "strict gene-symbol match to target genome annotation; not de novo motif discovery"},
        ]
    )
    summary.to_csv(audit_file, sep="\t", index=False)
    if generated_out.empty:
        raise ValueError(
            f"Generated motif2TF table for {species} is empty. "
            f"Audit written to {audit_file}; unmatched TFs written to {unmatched_file}."
        )
    generated_out.to_csv(out_file, sep="\t", index=False)
    if install_to_inputs:
        install_motif2tf(out_file, inputs)
    write_active_motif2tf_pointer(resources, species, "generated_symbol_evidence", out_file, install_to_inputs)
    return out_file


def prepare_motif2tf(
    species: str,
    resources: Path,
    inputs: Path,
    force: bool = False,
    install_to_inputs: bool = True,
    motif2tf_reference: str = "auto",
    motif2tf_table: str = "",
    orthology_policy: str = "one2one",
    generate_motif2tf: bool = False,
) -> Path:
    cfg = cfg_for(species)
    if motif2tf_table:
        return use_user_motif2tf_table(species, resources, inputs, motif2tf_table, install_to_inputs=install_to_inputs)
    if generate_motif2tf:
        return generate_motif2tf_from_symbol_evidence(species, resources, inputs, force=force, install_to_inputs=install_to_inputs)

    reference = resolve_reference(species, motif2tf_reference)
    if not reference:
        reason = cfg.get("restricted_reason", "no default motif2TF reference is defined")
        raise ValueError(
            f"{species} has no safe default motif2TF mapping strategy: {reason} "
            "Provide MOTIF2TF_TABLE=/path/to/motif_annotations.tbl for species-specific analysis, "
            "set --motif2tf-reference human|mouse|fly|chicken for explicit mapping, "
            "or run with --generate-motif2tf for an audited symbol-evidence table."
        )
    source = reference_source(resources, reference, force=force)
    motif2tf_dir = resources / "motif2tf"
    strategy = "direct" if cfg["motif2tf_mode"] == "direct" and reference == cfg.get("motif2tf_reference") else "mapped"
    species_out = motif2tf_output_path(resources, species, strategy, reference)
    if cfg["motif2tf_mode"] == "direct" and reference == cfg.get("motif2tf_reference"):
        shutil.copy2(source, species_out)
        audit_direct_motif2tf(species, resources, reference, source, species_out)
    else:
        convert_reference_motif2tf_to_target(
            species,
            resources,
            reference,
            source,
            species_out,
            force=force,
            orthology_policy=orthology_policy,
        )
    if install_to_inputs:
        install_motif2tf(species_out, inputs)
    write_active_motif2tf_pointer(resources, species, strategy, species_out, install_to_inputs)
    return species_out


def audit_direct_motif2tf(species: str, resources: Path, reference: str, source: Path, out_file: Path) -> None:
    species_dir = resources / species
    species_dir.mkdir(parents=True, exist_ok=True)
    motif = pd.read_csv(out_file, sep="\t")
    pd.DataFrame(
        [
            {"metric": "strategy", "value": "direct_aerts_public_table"},
            {"metric": "reference", "value": reference},
            {"metric": "source_path", "value": str(source)},
            {"metric": "motif_rows", "value": motif.shape[0]},
            {"metric": "unique_tfs", "value": motif["gene_name"].nunique() if "gene_name" in motif.columns else ""},
        ]
    ).to_csv(species_dir / f"{species}_motif2tf_direct_audit.tsv", sep="\t", index=False)


def convert_reference_motif2tf_to_target(
    species: str,
    resources: Path,
    reference: str,
    reference_motif2tf: Path,
    out_file: Path,
    force: bool = False,
    orthology_policy: str = "one2one",
) -> Path:
    cfg = cfg_for(species)
    ref = DIRECT_MOTIF2TF[reference]
    label = cfg["orthology_label"]
    ref_prefix = ref["orthology_prefix"]
    orth_prefix = cfg["orthology_prefix"]
    species_dir = resources / species
    species_dir.mkdir(parents=True, exist_ok=True)
    orthology_cache = species_dir / f"{reference}_{species}_orthology_biomart.tsv"
    audit_file = species_dir / f"{species}_motif2tf_{reference}_orthology_audit.tsv"
    missing_tfs_file = species_dir / f"{species}_motif2tf_{reference}_unmapped_reference_tfs.tsv"
    non_one_to_one_file = species_dir / f"{species}_motif2tf_{reference}_non_one_to_one_orthology.tsv"
    if out_file.exists() and audit_file.exists() and orthology_cache.exists() and not force:
        logging.info("SKIP existing %s motif2TF conversion from %s: %s", label, reference, out_file)
        return out_file
    try:
        from pybiomart import Server
    except Exception as exc:
        raise ImportError(f"pybiomart is required for {label} motif2TF orthology conversion from {reference}.") from exc

    required = [
        "external_gene_name",
        f"{orth_prefix}_homolog_associated_gene_name",
        f"{orth_prefix}_homolog_orthology_type",
        f"{orth_prefix}_homolog_perc_id",
        f"{orth_prefix}_homolog_perc_id_r1",
    ]
    if orthology_cache.exists() and not force:
        raw_orth = pd.read_csv(orthology_cache, sep="\t")
    else:
        server = Server(host=ref["biomart_host"])
        mart = server[ref["mart"]]
        dataset = mart[ref["dataset"]]
        missing = [x for x in required if x not in dataset.attributes]
        if missing:
            available = [x for x in dataset.attributes if f"{orth_prefix}_homolog" in x]
            raise ValueError(
                f"BioMart does not expose a usable {reference}->{species} orthology mapping. "
                "Missing attributes: " + ", ".join(missing) + "\n" + "\n".join(available[:80])
            )
        raw_orth = dataset.query(attributes=required)
        raw_orth.to_csv(orthology_cache, sep="\t", index=False)

    raw_orth = raw_orth.copy()
    raw_orth.columns = ["reference_gene", "target_gene", "orthology_type", "reference_perc_id", "target_perc_id"]
    orth = raw_orth.dropna(subset=["reference_gene", "target_gene"])
    orth = orth.query("reference_gene != '' and target_gene != ''")
    non_one = orth.query("orthology_type != 'ortholog_one2one'").copy()
    if not non_one.empty:
        non_one.to_csv(non_one_to_one_file, sep="\t", index=False)
    if orthology_policy == "one2one":
        orth = orth.query("orthology_type == 'ortholog_one2one'").copy()
        orth = orth.sort_values(["reference_gene", "reference_perc_id", "target_perc_id"], ascending=[True, False, False])
        orth = orth.drop_duplicates("reference_gene")
    else:
        orth = orth.sort_values(["reference_gene", "target_gene", "reference_perc_id", "target_perc_id"], ascending=[True, True, False, False])
        orth = orth.drop_duplicates(["reference_gene", "target_gene"])
    motif = pd.read_csv(reference_motif2tf, sep="\t")
    if "gene_name" not in motif.columns:
        raise ValueError("Expected column gene_name not found in motif2TF table.")
    missing_tfs = sorted(set(motif["gene_name"].dropna().astype(str)) - set(orth["reference_gene"].dropna().astype(str)))
    pd.DataFrame({"reference_gene": missing_tfs}).to_csv(missing_tfs_file, sep="\t", index=False)
    converted = motif.merge(orth[["reference_gene", "target_gene"]], left_on="gene_name", right_on="reference_gene", how="inner")
    converted["gene_name"] = converted["target_gene"]
    converted = converted[motif.columns].drop_duplicates()
    if converted.empty:
        raise ValueError(f"{label.capitalize()} motif2TF conversion from {reference} produced zero rows.")
    converted.to_csv(out_file, sep="\t", index=False)
    ref_unique = motif["gene_name"].nunique()
    target_unique = converted["gene_name"].nunique()
    summary = pd.DataFrame(
        [
            {"metric": "strategy", "value": "reference_to_target_orthology"},
            {"metric": "orthology_policy", "value": orthology_policy},
            {"metric": "reference", "value": reference},
            {"metric": "reference_dataset", "value": ref["dataset"]},
            {"metric": "reference_orthology_prefix", "value": ref_prefix},
            {"metric": "target_orthology_prefix", "value": orth_prefix},
            {"metric": f"{reference}_motif_rows", "value": motif.shape[0]},
            {"metric": f"{reference}_unique_tfs", "value": ref_unique},
            {"metric": "raw_orthology_rows", "value": raw_orth.shape[0]},
            {"metric": "one_to_one_ortholog_genes", "value": orth.shape[0]},
            {"metric": "non_one_to_one_ortholog_rows", "value": non_one.shape[0]},
            {"metric": "unmapped_reference_tfs", "value": len(missing_tfs)},
            {"metric": f"{species}_motif_rows", "value": converted.shape[0]},
            {"metric": f"{species}_unique_tfs", "value": target_unique},
            {"metric": "reference_tf_coverage", "value": round(target_unique / ref_unique, 4) if ref_unique else 0},
        ]
    )
    summary.to_csv(audit_file, sep="\t", index=False)
    return out_file


def status_row(name: str, path: Path, kind: str, required: bool = True, min_size: int = 1) -> dict[str, object]:
    exists = path.exists()
    ok = file_ok(path, min_size=min_size)
    if kind == "gzip" and exists:
        ok = gzip_ok(path)
    return {
        "item": name,
        "path": str(path),
        "kind": kind,
        "required": required,
        "exists": exists,
        "ok": ok,
        "size": path.stat().st_size if exists else 0,
        "sha256": sha256(path) if ok and path.is_file() else "",
    }


def collect_status(
    organism_list: list[str],
    resources: Path,
    inputs: Path,
    skip_motifs: bool,
    require_project_motif_annotation: bool,
    require_genome_resources: bool = True,
) -> pd.DataFrame:
    rows = []
    for species in organism_list:
        prefix = cfg_for(species)["prefix"]
        species_dir = resources / species
        rows.extend(
            [
                status_row(f"{species}.genome_ensembl_fasta_gz", species_dir / "genome.ensembl.fa.gz", "gzip", required=False),
                status_row(f"{species}.genes_ensembl_gtf_gz", species_dir / "genes.ensembl.gtf.gz", "gzip", required=False),
                status_row(f"{species}.allowed_chroms", resources / "chromosomes" / f"{species}.ucsc.standard.chroms.txt", "text", required=require_genome_resources),
                status_row(f"{species}.ucsc_fasta", species_dir / f"{prefix}.ucsc.standard.fa", "text", required=require_genome_resources),
                status_row(f"{species}.ucsc_gtf", species_dir / f"{prefix}.ucsc.standard.gtf", "text", required=require_genome_resources),
                status_row(f"{species}.chromsizes", species_dir / f"{prefix}.ucsc.standard.chromsizes.tsv", "table", required=require_genome_resources),
                status_row(f"{species}.genome_annotation", species_dir / f"{prefix}.ucsc.standard.genome_annotation.tsv", "table", required=require_genome_resources),
            ]
        )
    if not skip_motifs:
        rows.extend(
            [
                status_row("motifs.singletons_dir", resources / "motifs" / MOTIF_COLLECTION / "singletons", "directory"),
                status_row("motifs.motifs_txt", resources / "motifs" / "motifs.txt", "text"),
            ]
        )
        for species in organism_list:
            rows.append(status_row(f"{species}.motif_annotations_active_pointer", resources / "motif2tf" / f"motif_annotations.{species}.active.tsv", "table"))
        rows.append(status_row("project.motif_annotations", inputs / "cistarget_db" / "motif_annotations.tbl", "table", required=require_project_motif_annotation))
    return pd.DataFrame(rows)


def validate_status(df: pd.DataFrame) -> None:
    bad = df.loc[df["required"].astype(bool) & ~df["ok"].astype(bool)]
    if not bad.empty:
        raise FileNotFoundError("Missing or incomplete resources:\n" + bad[["item", "path", "exists", "ok", "size"]].to_string(index=False))


def write_manifest(
    organism_list: list[str],
    resources: Path,
    inputs: Path,
    release: str,
    status: pd.DataFrame,
    log_path: Path,
    motif2tf_reference: str,
    motif2tf_table: str,
) -> Path:
    organisms = {}
    for sp in organism_list:
        cfg = cfg_for(sp)
        organisms[sp] = {k: v for k, v in cfg.items() if k not in {"chroms"}}
        organisms[sp]["chromosomes"] = cfg.get("chroms", [])
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "ensembl_release": release,
        "organisms": organisms,
        "motif_collection": MOTIF_COLLECTION,
        "motif_collection_url": MOTIF_COLLECTION_URL,
        "motif2tf_version": MOTIF2TF_VERSION,
        "direct_motif2tf_references": DIRECT_MOTIF2TF,
        "motif2tf_reference": motif2tf_reference,
        "motif2tf_table": motif2tf_table,
        "resources_dir": str(resources),
        "inputs_dir": str(inputs),
        "log": str(log_path),
        "status": status.to_dict(orient="records"),
    }
    out = resources / "resource_manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    status.to_csv(resources / "resource_status.tsv", sep="\t", index=False)
    return out


def read_macs_genome_size(chromsizes: Path) -> str:
    if not chromsizes.exists():
        return ""
    total = 0
    with chromsizes.open() as handle:
        header = next(handle, "")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                total += int(fields[2])
    return str(total) if total > 0 else ""


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def update_project_env_macs_genome_size(species: str, resources: Path) -> None:
    env_path = project_dir() / "project_env.sh"
    if not env_path.exists():
        logging.info("project_env.sh not found under PROJECT_DIR; skip MACS_GENOME_SIZE update")
        return
    prefix = cfg_for(species)["prefix"]
    chromsizes = resources / species / f"{prefix}.ucsc.standard.chromsizes.tsv"
    genome_size = read_macs_genome_size(chromsizes)
    if not genome_size:
        logging.info("chromsizes not ready; skip MACS_GENOME_SIZE update: %s", chromsizes)
        return
    lines = env_path.read_text().splitlines()
    updated = False
    new_line = f"export MACS_GENOME_SIZE={shell_quote(genome_size)}"
    for i, line in enumerate(lines):
        if line.startswith("export MACS_GENOME_SIZE="):
            lines[i] = new_line
            updated = True
            break
    if not updated:
        insert_at = max(0, len(lines) - 1) if lines and lines[-1] == 'cd "$PROJECT_DIR"' else len(lines)
        lines.insert(insert_at, new_line)
    env_path.write_text("\n".join(lines) + "\n")
    logging.info("UPDATED %s with MACS_GENOME_SIZE=%s", env_path, genome_size)


def prepare_species(
    species: str,
    resources: Path,
    inputs: Path,
    release: str,
    force: bool,
    skip_motifs: bool,
    install_motif_to_inputs: bool,
    motif2tf_reference: str,
    motif2tf_table: str,
    orthology_policy: str,
    generate_motif2tf: bool,
) -> None:
    cfg = cfg_for(species)
    species_dir = resources / species
    species_dir.mkdir(parents=True, exist_ok=True)
    curl_download(cfg["fasta_url"].format(release=release), species_dir / "genome.ensembl.fa.gz", force=force)
    curl_download(cfg["gtf_url"].format(release=release), species_dir / "genes.ensembl.gtf.gz", force=force)
    convert_genome(species, resources, force=force)
    build_annotation(species, resources, force=force)
    if not skip_motifs:
        prepare_motif_collection(resources, force=force)
        prepare_motif2tf(
            species,
            resources,
            inputs,
            force=force,
            install_to_inputs=install_motif_to_inputs,
            motif2tf_reference=motif2tf_reference,
            motif2tf_table=motif2tf_table,
            orthology_policy=orthology_policy,
            generate_motif2tf=generate_motif2tf,
        )


def main() -> None:
    args = parse_args()
    if args.ref:
        args.motif2tf_reference = args.ref
    pdir = project_dir()
    if not args.organism:
        raise SystemExit("ERROR: --organism is required unless ORGANISM is set in the project environment.")
    log_path = setup_logging(args.log, pdir)
    resources = resolve_project_path(args.resources_dir, "resources", pdir)
    inputs = resolve_project_path(args.inputs_dir, "inputs", pdir)
    resources.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(parents=True, exist_ok=True)
    if args.organism == "all":
        organism_list = []
        for cfg in SPECIES.values():
            species, _ = resolve_species_config(cfg["organism"], args.ensembl_release, args.motif2tf_reference, args.target_orthology_prefix)
            organism_list.append(species)
    else:
        species, _ = resolve_species_config(args.organism, args.ensembl_release, args.motif2tf_reference, args.target_orthology_prefix)
        organism_list = [species]
    single_organism = len(organism_list) == 1
    logging.info(
        "mode=%s organism=%s ensembl_release=%s motif2tf_reference=%s orthology_policy=%s motif2tf_table=%s resources=%s",
        args.mode,
        ",".join(organism_list),
        args.ensembl_release,
        args.motif2tf_reference,
        args.orthology_policy,
        args.motif2tf_table,
        resources,
    )
    if args.mode == "prepare":
        for species in organism_list:
            if args.motifs_only:
                prepare_motif_collection(resources, force=args.force)
                prepare_motif2tf(
                    species,
                    resources,
                    inputs,
                    force=args.force,
                    install_to_inputs=single_organism,
                    motif2tf_reference=args.motif2tf_reference,
                    motif2tf_table=args.motif2tf_table,
                    orthology_policy=args.orthology_policy,
                    generate_motif2tf=args.generate_motif2tf,
                )
            else:
                prepare_species(
                    species,
                    resources,
                    inputs,
                    args.ensembl_release,
                    args.force,
                    args.skip_motifs,
                    install_motif_to_inputs=single_organism,
                    motif2tf_reference=args.motif2tf_reference,
                    motif2tf_table=args.motif2tf_table,
                    orthology_policy=args.orthology_policy,
                    generate_motif2tf=args.generate_motif2tf,
                )
    status = collect_status(
        organism_list,
        resources,
        inputs,
        args.skip_motifs,
        require_project_motif_annotation=single_organism and not args.skip_motifs,
        require_genome_resources=not args.motifs_only,
    )
    status_path = resources / "resource_status.tsv"
    status.to_csv(status_path, sep="\t", index=False)
    manifest_path = write_manifest(
        organism_list,
        resources,
        inputs,
        args.ensembl_release,
        status,
        log_path,
        args.motif2tf_reference,
        args.motif2tf_table,
    )
    logging.info("WROTE %s", status_path)
    logging.info("WROTE %s", manifest_path)
    print(status[["item", "required", "ok", "size", "path"]].to_string(index=False))
    if args.mode in {"prepare", "check"}:
        validate_status(status)
        logging.info("resource validation OK")
    if args.mode == "prepare" and single_organism:
        update_project_env_macs_genome_size(organism_list[0], resources)


if __name__ == "__main__":
    main()
