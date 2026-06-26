#!/usr/bin/env python
"""Query Ensembl organism resource availability for this SCENIC+ workflow."""
from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from pathlib import Path


MODULES_DIR = Path(__file__).resolve().parents[1] / "modules"
sys.path.insert(0, str(MODULES_DIR))

from organism_resources import (  # noqa: E402
    AERTS_MOTIF_COLLECTION,
    DIRECT_MOTIF2TF_BY_SPECIES,
    infer_assembly_from_url,
    listed_ensembl_species,
    motif2tf_plan,
    pick_ensembl_file,
    resolve_organism,
)


def url_exists(url: str) -> bool:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=20):
            return True
    except Exception:
        try:
            with urllib.request.urlopen(url, timeout=20) as handle:
                handle.read(1)
            return True
        except Exception:
            return False


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def print_one(query: str, release: str, check_remote: bool) -> int:
    species_list = listed_ensembl_species(release)
    ensembl_species = resolve_organism(query, release, species_list)
    if ensembl_species not in species_list and ensembl_species not in DIRECT_MOTIF2TF_BY_SPECIES:
        print(f"organism_query\t{query}")
        print(f"ensembl_release\t{release}")
        print("ensembl_species\tNOT_FOUND")
        return 2

    fasta_url = pick_ensembl_file(release, ensembl_species, "fasta")
    gtf_url = pick_ensembl_file(release, ensembl_species, "gtf")
    fasta_ok = bool(fasta_url) and (url_exists(fasta_url) if check_remote else True)
    gtf_ok = bool(gtf_url) and (url_exists(gtf_url) if check_remote else True)
    motif2tf_ok, motif2tf_strategy, motif2tf_ref = motif2tf_plan(ensembl_species)
    source = motif2tf_ref or "generated table, user table, or explicit reference"

    print(
        "\t".join(
            [
                "no",
                "query",
                "ensembl_species",
                "release",
                "assembly",
                "fasta",
                "gtf",
                "motif_collection",
                "motif2tf",
                "motif2tf_strategy",
                "motif2tf_reference",
            ]
        )
    )
    print(
        "\t".join(
            [
                "1",
                query,
                ensembl_species,
                release,
                infer_assembly_from_url(fasta_url or gtf_url, release),
                yes_no(fasta_ok),
                yes_no(gtf_ok),
                f"yes:{AERTS_MOTIF_COLLECTION}",
                motif2tf_ok,
                motif2tf_strategy,
                source,
            ]
        )
    )
    return 0


def default_list_output(release: str) -> Path:
    project_dir = os.environ.get("PROJECT_DIR", "")
    base = Path(project_dir).expanduser().resolve() if project_dir else Path.cwd()
    return base / f"ensembl_release_{release}_organism_resources.tsv"


def write_list(release: str, output: str) -> int:
    species = listed_ensembl_species(release)
    out = Path(output).expanduser() if output else default_list_output(release)
    if not out.is_absolute():
        out = (Path.cwd() / out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["no\tensembl_species\trelease\tmotif2tf\tmotif2tf_strategy\tmotif2tf_reference"]
    for idx, item in enumerate(species, start=1):
        motif2tf_ok, motif2tf_strategy, motif2tf_ref = motif2tf_plan(item)
        lines.append(
            "\t".join(
                [
                    str(idx),
                    item,
                    release,
                    motif2tf_ok,
                    motif2tf_strategy,
                    motif2tf_ref or "-",
                ]
            )
        )
    out.write_text("\n".join(lines) + "\n")
    print(f"WROTE\t{out}")
    print(f"ensembl_release\t{release}")
    print(f"species_count\t{len(species)}")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "List Ensembl release species with both FASTA and GTF, and report "
            "the motif2TF strategy available to this workflow."
        )
    )
    parser.add_argument("--organism", default=os.environ.get("ORGANISM", ""), help="Ensembl species name or common alias, for example mus_musculus or mouse.")
    parser.add_argument("--release", default=os.environ.get("ENSEMBL_RELEASE", "115"), help="Ensembl release, default: 115.")
    parser.add_argument("--list", action="store_true", help="Write a TSV listing species detected in both Ensembl release FASTA and GTF directories.")
    parser.add_argument("--output", default="", help="TSV output path for --list. Default: $PROJECT_DIR/ensembl_release_<release>_organism_resources.tsv or current directory.")
    parser.add_argument("--check-remote", action="store_true", help="HEAD-check selected FASTA/GTF URLs instead of only resolving them.")
    args = parser.parse_args()

    if args.list:
        raise SystemExit(write_list(args.release, args.output))
    if not args.organism:
        raise SystemExit("ERROR: set --organism or ORGANISM, or use --list.")
    raise SystemExit(print_one(args.organism, args.release, args.check_remote))


if __name__ == "__main__":
    main()
