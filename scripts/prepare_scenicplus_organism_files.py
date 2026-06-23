#!/usr/bin/env python
"""Install organism-specific resource files and write SCENIC+ config snippet."""
from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path


CONFIG = {
    "human": ("hsapiens", "homo_sapiens"),
    "mouse": ("mmusculus", "mus_musculus"),
    "cyno": ("custom", "custom"),
    "rat": ("custom", "custom"),
    "rabbit": ("custom", "custom"),
    "chicken": ("custom", "custom"),
    "zebrafish": ("custom", "custom"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", choices=sorted(CONFIG), default=os.environ.get("ORGANISM"))
    parser.add_argument("--resources-dir", default=None, help="Default: $PROJECT_DIR/resources")
    parser.add_argument("--scenicplus-dir", default=None, help="Default: $PROJECT_DIR/work/scenicplus")
    parser.add_argument("--out-config", default=None, help="Default: $PROJECT_DIR/work/scenicplus/organism_config.yaml")
    return parser.parse_args()


def project_dir() -> Path:
    return Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()


def resolve_project_path(path_value: str | None, default_rel: str, base: Path) -> Path:
    path = Path(path_value).expanduser() if path_value else base / default_rel
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def main() -> None:
    args = parse_args()
    if not args.organism:
        raise SystemExit("ERROR: --organism is required unless ORGANISM is set in project_env.sh.")
    pdir = project_dir()
    resources = resolve_project_path(args.resources_dir, "resources", pdir)
    scenicplus_dir = resolve_project_path(args.scenicplus_dir, "work/scenicplus", pdir)
    scenicplus_dir.mkdir(parents=True, exist_ok=True)

    species_dir = resources / args.organism
    prefix = args.organism
    annotation = species_dir / f"{prefix}.ucsc.standard.genome_annotation.tsv"
    chromsizes = species_dir / f"{prefix}.ucsc.standard.chromsizes.tsv"
    if not annotation.exists():
        raise FileNotFoundError(annotation)
    if not chromsizes.exists():
        raise FileNotFoundError(chromsizes)

    shutil.copy2(annotation, scenicplus_dir / "genome_annotation.tsv")
    shutil.copy2(chromsizes, scenicplus_dir / "chromsizes.tsv")

    data_species, motif_species = CONFIG[args.organism]
    config = f"""params_data_preparation:
  species: "{data_species}"
  biomart_host: "http://www.ensembl.org/"
params_motif_enrichment:
  species: "{motif_species}"
  annotation_version: "v10nr_clust"
  annotations_to_use: "Direct_annot Orthology_annot"
"""
    out_config = resolve_project_path(args.out_config, "work/scenicplus/organism_config.yaml", pdir)
    out_config.parent.mkdir(parents=True, exist_ok=True)
    out_config.write_text(config)
    print(f"WROTE {out_config}")
    print(f"COPIED {annotation} -> {scenicplus_dir / 'genome_annotation.tsv'}")
    print(f"COPIED {chromsizes} -> {scenicplus_dir / 'chromsizes.tsv'}")


if __name__ == "__main__":
    main()
