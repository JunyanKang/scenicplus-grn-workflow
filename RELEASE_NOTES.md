# Release Notes

## scenicplus-grn-installer v0.1.2

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

The installer is intended to be unpacked anywhere, copied under a conda-style root when requested, and run with:

```bash
bash install.sh
```

The workflow scripts are installed to:

```text
$CONDA_PREFIX/share/scenicplus-grn
```
