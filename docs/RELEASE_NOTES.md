# Release Notes

## scenicplus-grn-installer v0.1.5

This release packages a reproducible SCENIC+/GRN conda installer and a strict matched snRNA+snATAC workflow for annotated scMultiome projects.

Main deliverables:

- Dedicated conda environment installer for `scenicplus-grn`.
- Offline-capable bundled source archives for unstable GitHub access.
- Version-locked SCENIC+, pycisTopic, pycistarget, create_cisTarget_databases and Cluster-Buster layer.
- MALLET-enabled pycisTopic topic modeling path.
- hdWGCNA-based metacell workflow support.
- Step-by-step workflow scripts for:
  - annotated object inspection/export,
  - metacell aggregation,
  - ATAC fragment/peak standardization,
  - pycisTopic object/model/region set/DAR generation,
  - custom cisTarget database construction,
  - SCENIC+ Snakemake config/preflight/run,
  - sample-level eRegulon AUC statistics,
  - PDF figure generation with source tables.
- English and Chinese documentation.
- Automatic resource planning for CPU/load/memory-sensitive steps.
- Resumable partial custom cisTarget database construction for large region-by-motif matrices, using the official partial, combine and scores-to-rankings scripts from `create_cisTarget_databases`.
- Release-package cleanup: offline sources are distributed as `archives/vendor.tar.gz` and expanded at runtime to hidden `.vendor/`.
- Scientific output-tier audit: postprocessing now separates QC, chromatin topics, region sets/DARs, motif/cisTarget, eRegulon activity, condition effects and mechanism-view readiness.
- Cleaner installer root layout: helper entry points live in `bin/`, environment recipes and locks in `config/`, and long workflow/version documents in `docs/`.
- Hardened environment checks: writable numba/matplotlib cache directories are set during checks to prevent copied-environment `scanpy` cache failures.
- macOS arm64 recipe now installs conda-provided `samtools`/`htslib`; Linux x86_64 recipe explicitly installs `bzip2`.
- Custom cisTarget resource planning now defaults to more realistic `cbust` worker memory and documents sequential partial execution as the safe default.
- Custom cisTarget database construction now emits periodic heartbeat records during long Cluster-Buster scans. The heartbeat reports elapsed time, child-process activity and watched output-file status, so a quiet log is easier to distinguish from a failed run.
- The heartbeat interval is controlled by `heartbeat_seconds` in `inputs/cistarget_db_params.tsv`; the default is 600 seconds.
- Environment self-checks now fail if `samtools`, `tabix`, `bgzip` or other core workflow tools resolve outside the active conda environment.

The installer is intended to be unpacked anywhere, copied under a conda-style root when requested, and run with:

```bash
bash install.sh
```

The workflow scripts are installed to:

```text
$CONDA_PREFIX/share/scenicplus-grn
```
