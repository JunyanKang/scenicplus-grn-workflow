# Changelog

## 0.1.5 - 2026-06-23

Adds long-running custom cisTarget heartbeat logging.

Highlights:

- `build_custom_cistarget_db.py` now emits periodic `HEARTBEAT` records while official `create_cistarget_motif_databases.py` partial, combine or rankings-conversion commands are still running.
- Heartbeats report elapsed time, main process PID, child process summary and watched output-file status without changing the official cisTarget computation path.
- Adds `heartbeat_seconds` to `inputs/cistarget_db_params.tsv`; the default interval is 600 seconds to keep server logs quiet while still proving that long scans are alive.
- Environment self-checks now fail if core workflow tools such as `samtools`, `tabix` or `bgzip` resolve outside the active conda environment.
- Updates English and Chinese step-by-step guides to explain why Cluster-Buster scans can remain silent until workers return.

## 0.1.4 - 2026-06-23

Hardens environment checks and resource-aware custom cisTarget execution.

Highlights:

- Sets writable `NUMBA_CACHE_DIR` and `MPLCONFIGDIR` during environment checks to avoid `scanpy`/numba cache failures on copied or restricted environments.
- Adds conda-provided `samtools`/`htslib` to the macOS arm64 recipe so `samtools`, `tabix` and `bgzip` resolve from the SCENIC+ environment.
- Adds explicit `bzip2` to the Linux x86_64 recipe for CentOS7/glibc 2.17 installs.
- Simplifies README files by keeping project workflow commands in the step-by-step guides.
- Updates custom cisTarget auto resource planning to favor more workers inside one active partial instead of parallel partial jobs that duplicate region-by-motif memory pressure.

## 0.1.3 - 2026-06-23

Normalizes the installer root layout.

Highlights:

- Keeps the root directory focused on `install.sh`, `README*`, `VERSION` and major functional directories.
- Moves helper entry points to `bin/`.
- Moves conda recipes, pip constraints, Snakemake template and locks to `config/`.
- Moves long workflow/version documents to `docs/`.
- Updates installer copy logic, workflow asset checks and README/step-by-step paths for the new layout.

## 0.1.2 - 2026-06-23

Normalizes the installer package layout and adds science-layered output auditing.

Highlights:

- Replaces root-level uncompressed `vendor/` distribution with `archives/vendor.tar.gz`.
- Expands bundled archives at runtime to hidden `.vendor/` cache so restricted-network installs still have local GitHub and MALLET fallbacks.
- Clarifies the boundary between user-facing `scripts/` and internal `modules/`.
- Adds `audit_scenicplus_output_tiers.py` and wires it into `run_scenicplus_postprocess.py --task audit`.
- Reorganizes SCENIC+ output guidance by scientific evidence layer: QC, chromatin topics, region sets/DARs, motif/cisTarget, eRegulon activity, condition effects and mechanism views.

## 0.1.1 - 2026-06-23

Adds resumable partial custom cisTarget database construction for large projects.

Highlights:

- Uses official `create_cisTarget_databases` partial scoring, combine, and scores-to-rankings conversion scripts.
- Automatically enables partial mode for large region-by-motif matrices.
- Skips existing non-empty partial score files on rerun, so interrupted database builds can resume.
- Lowers the Cluster-Buster worker memory estimate while preserving realtime CPU/load and available-memory checks.

## 0.1.0 - 2026-06-23

Initial version-controlled release of the reusable `scenicplus-grn-installer`.

Highlights:

- Provides macOS arm64 and Linux x86_64 conda recipes for a dedicated SCENIC+/GRN environment.
- Bundles pinned GitHub source archives for SCENIC+, pycisTopic, pycistarget, create_cisTarget_databases, Cluster-Buster, hdWGCNA and related dependencies.
- Adds MALLET support for large pycisTopic topic modeling.
- Adds AutoZyme as a no-dependency overlay that does not perturb pinned package versions.
- Installs workflow scripts under `$CONDA_PREFIX/share/scenicplus-grn`.
- Provides matched snRNA+snATAC step-by-step workflow scripts for annotated scMultiome objects plus original ATAC fragments and peaks.
- Adds Chinese documentation files:
  - `README.zh-CN.md`
  - `SCENICPLUS_STEP_BY_STEP.zh-CN.md`
- Adds automatic per-step resource detection for pycisTopic and custom cisTarget database construction.
- Adds project-local `MPLCONFIGDIR` and `NUMBA_CACHE_DIR` handling to avoid scanpy/numba cache failures on servers or restricted home directories.
