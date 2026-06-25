from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import gzip
import os
import shutil
import subprocess
import sys

import pandas as pd
import pysam


PROJECT = Path(os.environ.get("PROJECT_DIR", ".")).expanduser().resolve()
INPUTS = PROJECT / "inputs"
WORK = PROJECT / "work" / "metacell_fragments"
RESULTS = PROJECT / "results" / "metacells"


def write_project_setting(key: str, value: str) -> None:
    config = PROJECT / "scenicplus_project.env"
    lines = config.read_text().splitlines() if config.exists() else []
    assignment = f"{key}={value!r}"
    for i, line in enumerate(lines):
        if line.startswith(f"{key}="):
            lines[i] = assignment
            break
    else:
        lines.append(assignment)
    config.write_text("\n".join(lines) + "\n")


def resolve(path_value: str) -> Path:
    path = Path(str(path_value)).expanduser()
    if not path.is_absolute():
        path = PROJECT / path
    return path.resolve()


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def main() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(
            "Usage: reassign_fragments_to_metacells.py\n\n"
            "Rewrite single-cell ATAC fragments to metacell barcodes using "
            "inputs/metacell_membership.tsv and inputs/sample_sheet.tsv. "
            "Writes work/metacell_fragments/*.tsv.gz plus tabix indexes, then "
            "updates inputs/sample_sheet.tsv to point to metacell fragment files."
        )
        return
    WORK.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    sample_sheet_path = require(INPUTS / "sample_sheet.tsv")
    membership_path = require(INPUTS / "metacell_membership.tsv")
    sample_sheet = pd.read_csv(sample_sheet_path, sep="\t", dtype=str).fillna("")
    membership = pd.read_csv(membership_path, sep="\t", dtype=str).fillna("")
    required_membership = {"sample_id", "barcode", "metacell_barcode"}
    if not required_membership.issubset(membership.columns):
        raise ValueError(f"{membership_path} must contain columns: {sorted(required_membership)}")
    if "fragments_tsv_gz" not in sample_sheet.columns:
        raise ValueError(f"{sample_sheet_path} must contain fragments_tsv_gz")
    if not (INPUTS / "sample_sheet.single_cell_fragments.tsv").exists():
        sample_sheet.to_csv(INPUTS / "sample_sheet.single_cell_fragments.tsv", sep="\t", index=False)

    bgzip = shutil.which("bgzip")
    if bgzip is None:
        raise FileNotFoundError("bgzip was not found in PATH")

    rows = []
    updated = sample_sheet.copy()
    for idx, row in sample_sheet.iterrows():
        sample_id = str(row["sample_id"])
        fragments_in = require(resolve(row["fragments_tsv_gz"]))
        sample_map = defaultdict(list)
        sub = membership.loc[membership["sample_id"].astype(str) == sample_id]
        for item in sub.itertuples(index=False):
            sample_map[str(item.barcode)].append(str(item.metacell_barcode))
        if not sample_map:
            raise ValueError(f"No metacell memberships for sample {sample_id}")

        tmp_tsv = WORK / f"{sample_id}.metacell.fragments.tsv"
        out_gz = WORK / f"{sample_id}.metacell.fragments.tsv.gz"
        n_in = 0
        n_written = 0
        n_used_cells: set[str] = set()
        with open_text(fragments_in) as fin, open(tmp_tsv, "wt") as fout:
            for line in fin:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    raise ValueError(f"Malformed fragment line in {fragments_in}: {line[:100]}")
                n_in += 1
                barcode = fields[3]
                targets = sample_map.get(barcode)
                if targets is None:
                    continue
                n_used_cells.add(barcode)
                for metacell_barcode in targets:
                    fields[3] = metacell_barcode
                    fout.write("\t".join(fields) + "\n")
                    n_written += 1
        with open(out_gz, "wb") as out_handle:
            subprocess.run([bgzip, "-f", "-c", str(tmp_tsv)], stdout=out_handle, check=True)
        tmp_tsv.unlink()
        pysam.tabix_index(str(out_gz), preset="bed", force=True)
        updated.loc[idx, "fragments_tsv_gz"] = str(out_gz)
        rows.append(
            {
                "sample_id": sample_id,
                "input_fragments": str(fragments_in),
                "metacell_fragments": str(out_gz),
                "input_fragment_lines": n_in,
                "written_fragment_lines": n_written,
                "source_barcodes_used": len(n_used_cells),
                "metacell_barcodes": len(set(sum(sample_map.values(), []))),
            }
        )
    updated.to_csv(sample_sheet_path, sep="\t", index=False)
    pd.DataFrame(rows).to_csv(RESULTS / "metacell_fragment_reassignment.tsv", sep="\t", index=False)
    write_project_setting("ACTIVE_SAMPLE_SHEET", str(sample_sheet_path))
    write_project_setting("METACELL_FRAGMENT_REASSIGNMENT", str(RESULTS / "metacell_fragment_reassignment.tsv"))
    print("WROTE", RESULTS / "metacell_fragment_reassignment.tsv")
    print("UPDATED", sample_sheet_path)


if __name__ == "__main__":
    main()
