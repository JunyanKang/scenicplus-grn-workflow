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


ENSEMBL_RELEASE = "115"
MOTIF_COLLECTION = "v10nr_clust_public"
MOTIF2TF_VERSION = "v10nr_clust"
SPECIES = {
    "human": {
        "organism": "homo_sapiens",
        "assembly": "GRCh38",
        "prefix": "human",
        "chroms": [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz",
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_mode": "direct",
    },
    "mouse": {
        "organism": "mus_musculus",
        "assembly": "GRCm39",
        "prefix": "mouse",
        "chroms": [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/mus_musculus/Mus_musculus.GRCm39.{release}.gtf.gz",
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        "motif2tf_mode": "direct",
    },
    "cyno": {
        "organism": "macaca_fascicularis",
        "assembly": "Macaca_fascicularis_6.0",
        "prefix": "cyno",
        "chroms": [f"chr{i}" for i in range(1, 21)] + ["chrX"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.{release}.gtf.gz",
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_mode": "human_to_target_orthology",
        "orthology_prefix": "mfascicularis",
        "orthology_label": "cynomolgus",
    },
    "rat": {
        "organism": "rattus_norvegicus",
        "assembly": "GRCr8",
        "prefix": "rat",
        "chroms": [f"chr{i}" for i in range(1, 21)] + ["chrX", "chrY"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.{release}.gtf.gz",
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_mode": "human_to_target_orthology",
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
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_mode": "human_to_target_orthology",
        "orthology_prefix": "ocuniculus",
        "orthology_label": "rabbit",
    },
    "chicken": {
        "organism": "gallus_gallus",
        "assembly": "bGalGal1.mat.broiler.GRCg7b",
        "prefix": "chicken",
        "chroms": [f"chr{i}" for i in range(1, 40)] + ["chrZ", "chrW"],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.{release}.gtf.gz",
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_mode": "human_to_target_orthology",
        "orthology_prefix": "ggallus",
        "orthology_label": "chicken",
    },
    "zebrafish": {
        "organism": "danio_rerio",
        "assembly": "GRCz11",
        "prefix": "zebrafish",
        "chroms": [f"chr{i}" for i in range(1, 26)],
        "fasta_url": "https://ftp.ensembl.org/pub/release-{release}/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-{release}/gtf/danio_rerio/Danio_rerio.GRCz11.{release}.gtf.gz",
        "motif2tf_url": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_file": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "motif2tf_mode": "human_to_target_orthology",
        "orthology_prefix": "drerio",
        "orthology_label": "zebrafish",
    },
}
MOTIF_COLLECTION_URL = "https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", choices=sorted(SPECIES) + ["all"], default=os.environ.get("ORGANISM"))
    parser.add_argument("--resources-dir", default=None, help="Default: $PROJECT_DIR/resources")
    parser.add_argument("--inputs-dir", default=None, help="Default: $PROJECT_DIR/inputs")
    parser.add_argument("--ensembl-release", default=os.environ.get("ENSEMBL_RELEASE", ENSEMBL_RELEASE))
    parser.add_argument("--mode", choices=["prepare", "check", "status"], default="prepare")
    parser.add_argument("--force", action="store_true", help="Regenerate derived files and re-download direct resources.")
    parser.add_argument("--skip-motifs", action="store_true", help="Skip motif collection and motif2TF preparation.")
    parser.add_argument("--motifs-only", action="store_true", help="Prepare only motif collection and motif2TF files for the selected organism.")
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
    cfg = SPECIES[species]
    chrom_dir = resources / "chromosomes"
    chrom_dir.mkdir(parents=True, exist_ok=True)
    out = chrom_dir / f"{species}.ucsc.standard.chroms.txt"
    out.write_text("\n".join(cfg["chroms"]) + "\n")
    return out


def to_ucsc(chrom: str) -> str:
    chrom = chrom.split()[0]
    if chrom in {"MT", "M", "Mt", "mitochondrion_genome"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else "chr" + chrom


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def convert_genome(species: str, resources: Path, force: bool = False) -> dict[str, Path]:
    cfg = SPECIES[species]
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
    cfg = SPECIES[species]
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


def prepare_motif2tf(species: str, resources: Path, inputs: Path, force: bool = False, install_to_inputs: bool = True) -> Path:
    cfg = SPECIES[species]
    motif2tf_dir = resources / "motif2tf"
    motif2tf_dir.mkdir(parents=True, exist_ok=True)
    inputs_db = inputs / "cistarget_db"
    inputs_db.mkdir(parents=True, exist_ok=True)
    source = motif2tf_dir / cfg["motif2tf_file"]
    curl_download(cfg["motif2tf_url"], source, force=force)
    species_out = motif2tf_dir / f"motif_annotations.{species}.tbl"
    if cfg["motif2tf_mode"] == "direct":
        shutil.copy2(source, species_out)
    else:
        convert_human_motif2tf_to_target(species, resources, source, species_out, force=force)
    if install_to_inputs:
        target = inputs_db / "motif_annotations.tbl"
        for conflict in inputs_db.glob("motif_annotations*_Conflict.tbl"):
            conflict.unlink()
        tmp_target = inputs_db / f".motif_annotations.{os.getpid()}.tmp"
        with species_out.open("rb") as src, tmp_target.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        os.replace(tmp_target, target)
    return species_out


def convert_human_motif2tf_to_target(species: str, resources: Path, human_motif2tf: Path, out_file: Path, force: bool = False) -> Path:
    cfg = SPECIES[species]
    label = cfg["orthology_label"]
    orth_prefix = cfg["orthology_prefix"]
    species_dir = resources / species
    species_dir.mkdir(parents=True, exist_ok=True)
    orthology_cache = species_dir / f"human_{species}_orthology_biomart.tsv"
    audit_file = species_dir / f"{species}_motif2tf_orthology_audit.tsv"
    if out_file.exists() and audit_file.exists() and orthology_cache.exists() and not force:
        logging.info("SKIP existing %s motif2TF conversion: %s", label, out_file)
        return out_file
    try:
        from pybiomart import Server
    except Exception as exc:
        raise ImportError(f"pybiomart is required for {label} motif2TF orthology conversion.") from exc

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
        server = Server(host="http://www.ensembl.org")
        mart = server["ENSEMBL_MART_ENSEMBL"]
        dataset = mart["hsapiens_gene_ensembl"]
        missing = [x for x in required if x not in dataset.attributes]
        if missing:
            available = [x for x in dataset.attributes if f"{orth_prefix}_homolog" in x]
            raise ValueError("Missing BioMart attributes: " + ", ".join(missing) + "\n" + "\n".join(available[:80]))
        raw_orth = dataset.query(attributes=required)
        raw_orth.to_csv(orthology_cache, sep="\t", index=False)

    raw_orth = raw_orth.copy()
    raw_orth.columns = ["human_gene", "target_gene", "orthology_type", "human_perc_id", "target_perc_id"]
    orth = raw_orth.dropna(subset=["human_gene", "target_gene"])
    orth = orth.query("human_gene != '' and target_gene != ''")
    orth = orth.query("orthology_type == 'ortholog_one2one'").copy()
    orth = orth.sort_values(["human_gene", "human_perc_id", "target_perc_id"], ascending=[True, False, False])
    orth = orth.drop_duplicates("human_gene")
    motif = pd.read_csv(human_motif2tf, sep="\t")
    if "gene_name" not in motif.columns:
        raise ValueError("Expected column gene_name not found in motif2TF table.")
    converted = motif.merge(orth[["human_gene", "target_gene"]], left_on="gene_name", right_on="human_gene", how="inner")
    converted["gene_name"] = converted["target_gene"]
    converted = converted[motif.columns].drop_duplicates()
    if converted.empty:
        raise ValueError(f"{label.capitalize()} motif2TF conversion produced zero rows.")
    converted.to_csv(out_file, sep="\t", index=False)
    summary = pd.DataFrame(
        [
            {"metric": "human_motif_rows", "value": motif.shape[0]},
            {"metric": "human_unique_tfs", "value": motif["gene_name"].nunique()},
            {"metric": "raw_orthology_rows", "value": raw_orth.shape[0]},
            {"metric": "one_to_one_ortholog_genes", "value": orth.shape[0]},
            {"metric": f"{species}_motif_rows", "value": converted.shape[0]},
            {"metric": f"{species}_unique_tfs", "value": converted["gene_name"].nunique()},
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


def collect_status(organism_list: list[str], resources: Path, inputs: Path, skip_motifs: bool, require_project_motif_annotation: bool) -> pd.DataFrame:
    rows = []
    for species in organism_list:
        prefix = SPECIES[species]["prefix"]
        species_dir = resources / species
        rows.extend(
            [
                status_row(f"{species}.genome_ensembl_fasta_gz", species_dir / "genome.ensembl.fa.gz", "gzip", required=False),
                status_row(f"{species}.genes_ensembl_gtf_gz", species_dir / "genes.ensembl.gtf.gz", "gzip", required=False),
                status_row(f"{species}.allowed_chroms", resources / "chromosomes" / f"{species}.ucsc.standard.chroms.txt", "text"),
                status_row(f"{species}.ucsc_fasta", species_dir / f"{prefix}.ucsc.standard.fa", "text"),
                status_row(f"{species}.ucsc_gtf", species_dir / f"{prefix}.ucsc.standard.gtf", "text"),
                status_row(f"{species}.chromsizes", species_dir / f"{prefix}.ucsc.standard.chromsizes.tsv", "table"),
                status_row(f"{species}.genome_annotation", species_dir / f"{prefix}.ucsc.standard.genome_annotation.tsv", "table"),
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
            rows.append(status_row(f"{species}.motif_annotations", resources / "motif2tf" / f"motif_annotations.{species}.tbl", "table"))
        rows.append(status_row("project.motif_annotations", inputs / "cistarget_db" / "motif_annotations.tbl", "table", required=require_project_motif_annotation))
    return pd.DataFrame(rows)


def validate_status(df: pd.DataFrame) -> None:
    bad = df.loc[df["required"].astype(bool) & ~df["ok"].astype(bool)]
    if not bad.empty:
        raise FileNotFoundError("Missing or incomplete resources:\n" + bad[["item", "path", "exists", "ok", "size"]].to_string(index=False))


def write_manifest(organism_list: list[str], resources: Path, inputs: Path, release: str, status: pd.DataFrame, log_path: Path) -> Path:
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "ensembl_release": release,
        "organisms": {sp: {k: v for k, v in SPECIES[sp].items() if k not in {"chroms"}} | {"chromosomes": SPECIES[sp]["chroms"]} for sp in organism_list},
        "motif_collection": MOTIF_COLLECTION,
        "motif_collection_url": MOTIF_COLLECTION_URL,
        "motif2tf_version": MOTIF2TF_VERSION,
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
    prefix = SPECIES[species]["prefix"]
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


def prepare_species(species: str, resources: Path, inputs: Path, release: str, force: bool, skip_motifs: bool, install_motif_to_inputs: bool) -> None:
    cfg = SPECIES[species]
    species_dir = resources / species
    species_dir.mkdir(parents=True, exist_ok=True)
    curl_download(cfg["fasta_url"].format(release=release), species_dir / "genome.ensembl.fa.gz", force=force)
    curl_download(cfg["gtf_url"].format(release=release), species_dir / "genes.ensembl.gtf.gz", force=force)
    convert_genome(species, resources, force=force)
    build_annotation(species, resources, force=force)
    if not skip_motifs:
        prepare_motif_collection(resources, force=force)
        prepare_motif2tf(species, resources, inputs, force=force, install_to_inputs=install_motif_to_inputs)


def main() -> None:
    args = parse_args()
    pdir = project_dir()
    if not args.organism:
        raise SystemExit("ERROR: --organism is required unless ORGANISM is set in the project environment.")
    log_path = setup_logging(args.log, pdir)
    resources = resolve_project_path(args.resources_dir, "resources", pdir)
    inputs = resolve_project_path(args.inputs_dir, "inputs", pdir)
    resources.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(parents=True, exist_ok=True)
    organism_list = list(SPECIES) if args.organism == "all" else [args.organism]
    single_organism = len(organism_list) == 1
    logging.info("mode=%s organism=%s ensembl_release=%s resources=%s", args.mode, ",".join(organism_list), args.ensembl_release, resources)
    if args.mode == "prepare":
        for species in organism_list:
            if args.motifs_only:
                prepare_motif_collection(resources, force=args.force)
                prepare_motif2tf(species, resources, inputs, force=args.force, install_to_inputs=single_organism)
            else:
                prepare_species(species, resources, inputs, args.ensembl_release, args.force, args.skip_motifs, install_motif_to_inputs=single_organism)
    status = collect_status(organism_list, resources, inputs, args.skip_motifs, require_project_motif_annotation=single_organism and not args.skip_motifs)
    status_path = resources / "resource_status.tsv"
    status.to_csv(status_path, sep="\t", index=False)
    manifest_path = write_manifest(organism_list, resources, inputs, args.ensembl_release, status, log_path)
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
