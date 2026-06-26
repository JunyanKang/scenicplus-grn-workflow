"""Shared Ensembl organism and motif2TF resource helpers."""
from __future__ import annotations

import html.parser
import re
import urllib.error
import urllib.request
from pathlib import Path


ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
AERTS_MOTIF_COLLECTION = "https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/"

ALIASES = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "zebrafish": "danio_rerio",
    "chicken": "gallus_gallus",
    "pig": "sus_scrofa",
    "cow": "bos_taurus",
    "dog": "canis_lupus_familiaris",
    "rabbit": "oryctolagus_cuniculus",
    "rhesus": "macaca_mulatta",
    "cyno": "macaca_fascicularis",
    "cynomolgus": "macaca_fascicularis",
    "fly": "drosophila_melanogaster",
}

DIRECT_MOTIF2TF_BY_SPECIES = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
    "gallus_gallus": "chicken",
    "drosophila_melanogaster": "fly",
}

DEFAULT_REF_BY_SPECIES = {
    "macaca_fascicularis": "human",
    "macaca_mulatta": "human",
    "oryctolagus_cuniculus": "human",
    "sus_scrofa": "human",
    "bos_taurus": "human",
    "canis_lupus_familiaris": "human",
    "rattus_norvegicus": "mouse",
}


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.links.append(value)


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def fetch_links(url: str) -> list[str]:
    with urllib.request.urlopen(url, timeout=30) as handle:
        text = handle.read().decode("utf-8", "replace")
    parser = LinkParser()
    parser.feed(text)
    return parser.links


def species_from_release_dir(release: str, resource: str) -> list[str]:
    links = fetch_links(f"{ENSEMBL_FTP}/release-{release}/{resource}/")
    species = []
    for link in links:
        clean = link.strip("/")
        if clean and clean not in {"ancestral_alleles", "compara", "multi"} and not clean.startswith("?"):
            species.append(clean)
    return sorted(set(species))


def listed_ensembl_species(release: str) -> list[str]:
    fasta_species = set(species_from_release_dir(release, "fasta"))
    gtf_species = set(species_from_release_dir(release, "gtf"))
    return sorted(fasta_species.intersection(gtf_species))


def resolve_organism(query: str, release: str | None = None, listed_species: list[str] | None = None) -> str:
    q = norm(query)
    if q in ALIASES:
        return ALIASES[q]
    if listed_species is None and release:
        try:
            listed_species = listed_ensembl_species(release)
        except urllib.error.URLError:
            listed_species = None
    if listed_species:
        by_norm = {norm(item): item for item in listed_species}
        if q in by_norm:
            return by_norm[q]
    return q


def pick_ensembl_file(release: str, species: str, kind: str) -> str | None:
    if kind == "fasta":
        base = f"{ENSEMBL_FTP}/release-{release}/fasta/{species}/dna/"
        suffixes = [".dna.primary_assembly.fa.gz", ".dna.toplevel.fa.gz"]
    elif kind == "gtf":
        base = f"{ENSEMBL_FTP}/release-{release}/gtf/{species}/"
        suffixes = [f".{release}.gtf.gz"]
    else:
        raise ValueError(f"Unsupported Ensembl resource kind: {kind}")
    try:
        links = fetch_links(base)
    except urllib.error.URLError:
        return None
    candidates = [x for x in links if any(x.endswith(suffix) for suffix in suffixes)]
    if kind == "fasta":
        candidates = sorted(candidates, key=lambda x: (".primary_assembly." not in x, x))
    else:
        candidates = sorted(candidates)
    return base + candidates[0] if candidates else None


def infer_assembly_from_url(url: str | None, release: str) -> str:
    if not url:
        return "detected_from_ensembl"
    name = Path(url).name
    name = re.sub(r"\.(dna\.primary_assembly|dna\.toplevel)\.fa\.gz$", "", name)
    name = re.sub(rf"\.{re.escape(str(release))}\.gtf\.gz$", "", name)
    parts = name.split(".")
    if len(parts) >= 2:
        return ".".join(parts[1:])
    return "detected_from_ensembl"


def infer_biomart_orthology_prefix(species: str) -> str:
    parts = species.split("_")
    if len(parts) < 2:
        return species.replace("_", "")
    return parts[0][0] + "".join(parts[1:])


def motif2tf_plan(species: str, explicit_ref: str = "auto", has_user_table: bool = False) -> tuple[str, str, str | None]:
    if has_user_table:
        return "yes", "user_supplied_species_specific_table", None
    if species in DIRECT_MOTIF2TF_BY_SPECIES:
        ref = DIRECT_MOTIF2TF_BY_SPECIES[species]
        return "yes", "direct", ref
    if explicit_ref != "auto":
        return "yes", f"map_from_{explicit_ref}", explicit_ref
    if species in DEFAULT_REF_BY_SPECIES:
        ref = DEFAULT_REF_BY_SPECIES[species]
        return "yes", f"map_from_{ref}", ref
    return "conditional", "generate_or_ref_or_user_table", None
