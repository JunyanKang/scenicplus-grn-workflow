# Release Notes

## scenicplus-grn-installer v0.1.39

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
- Custom cisTarget FASTA generation now preserves plain UCSC-style
  `chr:start-end` region names for pycistarget compatibility.
- SCENIC+ motif enrichment now has stage-specific `dem_n_cpu` and `ctx_n_cpu`
  controls with dynamic defaults based on detected memory and custom DB size.
  This prevents laptop runs from launching too many DEM/cisTarget workers while
  still allowing larger servers to use more resources.
- Installer defaults no longer assume project-specific condition labels, color
  keys, cell labels or example marker genes; project-specific priorities belong
  in each project directory.
- The split motif-enrichment runner now has `--mode status`, a hard completion
  gate that detects missing, invalid or half-written chunk outputs before
  students continue to preflight and Snakemake.
- The step-by-step guides have been rewritten as release-quality student-facing
  manuals with runnable command numbering and stage-level biological rationale.
- The pycisTopic workflow command now manages its own timestamped log under
  `$PROJECT_DIR/logs/`, so the documented command no longer needs shell `tee`.
- Bounded long-running loops now show terminal progress indicators: custom
  cisTarget partial parts, split motif-enrichment chunks and pycisTopic stages.
- Student-editable placeholders in the README and step-by-step guides are now
  marked with shell comments, including the `corrected_metadata_column` export
  override example.
- The README now distinguishes source checkout layout from release-package
  layout: `archives/vendor.tar.gz` is distributed as a release artifact, not as
  a Git-tracked repository file.
- README installation content is shorter: installed script path details and the
  project-initialization example were removed from README and left to the
  step-by-step workflow guides.
- GitHub-default documentation is now Chinese: `README.md` and
  `docs/SCENICPLUS_STEP_BY_STEP.md` contain the Chinese guide, while English is
  preserved as `.en.md`.
- The installer code is now explicitly licensed under MIT, with bundled
  third-party source archives documented in `THIRD_PARTY_NOTICES.md`.
- Duplicate `.zh-CN.md` documentation files were removed because Chinese is now
  the default documentation language.
- The default Chinese README no longer displays an English-documentation
  pointer near the top of the page.
- R plotting wrappers now handle `--help` without requiring project result
  tables.
- Generated SCENIC+ Snakefiles are patched to honor the stage-specific motif
  enrichment thread settings.
- Large custom cisTarget databases can now use
  `run_scenicplus_motif_enrichment_split.py`, which runs DEM and cisTarget
  enrichment in independent region-set-family chunks and combines them through
  SCENIC+ `prepare_menr`. The chunk runner is resumable and dynamically chooses
  safe chunk concurrency from detected memory and DB size.
- Postprocess condition statistics now handle single-condition,
  two-condition and multi-condition projects explicitly, and optional
  `priority_eregulons` labels keep project-specific volcano priorities in the
  project parameter file rather than in installer defaults.

The installer is intended to be unpacked anywhere, copied under a conda-style root when requested, and run with:

```bash
bash install.sh
```

The workflow scripts are installed to:

```text
$CONDA_PREFIX/share/scenicplus-grn
```
