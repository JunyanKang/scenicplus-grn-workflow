# Changelog

## 0.1.34 - 2026-06-26

Adds terminal progress indicators where the workflow has real bounded work
units.

Highlights:

- Adds `tqdm` progress bars for custom cisTarget partial database parts.
- Adds `tqdm` progress bars for split DEM/cisTarget region-set chunks.
- Adds a stage-level pycisTopic workflow progress indicator while retaining the
  timestamped run log.
- Keeps heartbeat logging for long single-process jobs where a true percentage
  is not available.

## 0.1.33 - 2026-06-26

Simplifies pycisTopic workflow logging.

Highlights:

- `spgrn-run-pycistopic-workflow` now writes a timestamped log under
  `$PROJECT_DIR/logs/` while still streaming output to the terminal.
- Step-by-step commands no longer require manual shell `tee` for pycisTopic.

## 0.1.32 - 2026-06-26

Rewrites the step-by-step workflow guides as release-quality student-facing
manuals.

Highlights:

- Keeps numbered runnable commands while removing development-note style text.
- Adds concise biological and bioinformatics rationale for each workflow stage.
- Keeps the English and Chinese guides aligned on the same 0-11 structure.
- Rechecks documented `spgrn-*` commands against the installed command surface.

## 0.1.31 - 2026-06-26

Hardens the optional split motif-enrichment path.

Highlights:

- Adds `spgrn-run-scenicplus-motif-enrichment-split --mode status` to validate
  chunk HDF5 files, DEM empty-result diagnostics and `prepare_menr` outputs.
- Updates the step-by-step guides so students use the status gate before Step 9
  instead of judging completion from file size alone.

## 0.1.30 - 2026-06-25

Separates installer-level workflow logic from project-specific analysis choices.

Highlights:

- Removes hard-coded condition-order preferences from generic plotting and source-table extraction.
- Stops mapping legacy project-specific color keys onto generic condition colors.
- Replaces project-specific changelog examples with generic wording.
- Keeps priority eRegulon labels as an explicit project-level input instead of an installer default.

## 0.1.29 - 2026-06-25

Tightens condition-statistics semantics and keeps project-specific figure priorities outside installer defaults.

Highlights:

- Adds optional `priority_eregulons` support for R-rendered condition volcano labels without embedding project-specific genes or cell labels.
- Clarifies automatic condition statistics: single-condition projects produce descriptive sample means, two-condition projects report comparison-minus-reference deltas, and multi-condition projects use an omnibus sample-level test.
- Replaces stale project-specific default colors in legacy Python plotting paths with generic fallback colors.
- Synchronizes English and Chinese step-by-step wording for condition statistics and priority volcano labels.

## 0.1.28 - 2026-06-25

Tightens SCENIC+ postprocess ordering, style-parameter documentation and eRegulon sign display.

Highlights:

- Keeps Step 11 as the single postprocess stage: audit first, 01-08 eRegulon PDFs/source tables second, and 09+ condition-statistics PDFs/tables last.
- Documents `spgrn-run-scenicplus-postprocess --task` and `--layer` values directly in the step-by-step guides and in the generated R plotting style document.
- Adds a compact output-location table instead of a long static PDF/source-file list.
- Generates `plot_style_parameters.md` from the R renderers and keeps `plot_style_parameters.tsv` as a clean two-column `parameter/value` file.
- Drops stale project-specific style keys when rewriting the style TSV while still mapping legacy color keys onto the generic style parameters.
- Preserves SCENIC+ eRegulon signs in display labels, so TF-level RNA expression is not conflated with target gene-set AUCell activity.
- Adds clearer command-line help for the postprocess wrapper.

## 0.1.27 - 2026-06-25

Switches SCENIC+ postprocess figures to an R rendering layer and fixes true condition-resolved plotting.

Highlights:

- Extracts SCENIC+ figure source tables with Python, then renders all PDF figures with R/ggplot2/patchwork.
- Makes `_condition_` heatmap and dot-heatmap outputs true `cell_label_condition` joint-group displays, ordered as cell label followed by condition.
- Makes condition UMAP files true condition facets and merges the eRegulon activity embedding into one PDF with all-cells plus condition panels.
- Stops generating condition-suffixed model-level region-gene, target-region-overlap and network PDFs.
- Adds `plot_style_parameters.tsv` for rerunning figures after editing colors, font sizes, point sizes, line widths and panel layout.
- Defines condition statistics automatically for single-condition, two-condition and multi-condition projects; two-condition deltas are recorded with the resolved `contrast`, and reference/comparison overrides remain available when needed.
- Clusters heatmap eRegulon rows from the plotted numeric matrix while preserving explicit `cell_label_condition` column order for condition heatmaps.
- Labels target genes in the compact TF-target network PDF.
- Merges postprocess documentation so 01-08 standard figures and 09+ condition-statistics figures are generated through one ordered `--task all` command.
- Automatically writes `plot_style_parameters.md` next to `plot_style_parameters.tsv` to document R plotting parameters and the one-command redraw workflow.
- Updates self-check requirements for the new extractor and R renderers, and verifies `samtools`, `tabix` and `bgzip` resolve from the conda environment.

## 0.1.26 - 2026-06-25

Removes metacell-stage ATAC doublet placeholder outputs.

Highlights:

- Stops writing `atac_doublets.tsv` and `atac_doublet_diagnostics.pdf` for metacell-based pycisTopic runs.
- Removes ATAC doublet files from pycisTopic completion checks, final checks, manifests and output-tier audit requirements.
- Removes unused doublet parameters from new pycisTopic parameter files.
- Documents doublet filtering as an upstream single-cell QC prerequisite rather than a current metacell-stage analysis output.

## 0.1.25 - 2026-06-25

Flattens SCENIC+ direct/extended figure outputs into one figure directory.

Highlights:

- Writes direct and extended outputs into `results/scenicplus_figures` without layer subdirectories.
- Appends `_direct` or `_extended` to figure and source-data filenames to avoid collisions.
- Adds condition-grouped counterparts for the main eRegulon figure series using `_condition_direct` and `_condition_extended` suffixes.
- Updates postprocess defaults, audit checks, and step-by-step output paths.

## 0.1.24 - 2026-06-25

Clusters numeric axes in SCENIC+ heatmap-style figures.

Highlights:

- Applies hierarchical clustering to heatmap rows and columns before plotting.
- Applies the same numeric ordering to eRegulon dot heatmaps.
- Clusters condition sample heatmaps and cell-label condition-effect heatmaps.
- Enlarges dot-heatmap bottom margin to avoid clipped rotated labels.

## 0.1.23 - 2026-06-25

Systematically numbers SCENIC+ figure PDFs within each layer directory.

Highlights:

- Names all direct and extended figure PDFs with directory-level numeric prefixes.
- Uses automatic numbering for condition-level PDFs, including the number of cell-label-stratified comparisons actually present.
- Updates audit checks and step-by-step output paths to the numbered PDF names.

## 0.1.22 - 2026-06-25

Moves condition-level AUC statistics into layer figure directories.

Highlights:

- Writes condition-level tables and PDFs under `results/scenicplus_figures/direct` and `results/scenicplus_figures/extended`.
- Splits the former multi-page condition statistics PDF into one-panel vector PDFs with numbered filenames.
- Adds cell-label-stratified condition AUC statistics and per-label volcano PDFs when both conditions are present in a label.
- Updates postprocess defaults, audit checks, and step-by-step output paths.

## 0.1.21 - 2026-06-25

Fixes SCENIC+ postprocess UMAP plotting and condition-level display labels.

Highlights:

- Preserves annotated-object UMAP coordinates during metacell metadata export.
- Records postprocess `umap_x` and `umap_y` from the selected reduction, so
  eRegulon AUC UMAP plots no longer fall back to a placeholder page.
- Uses manuscript-readable eRegulon labels such as `Sox5 (62 targets)` in
  figures while keeping raw SCENIC+ names in source tables.
- Adds display labels and an explicit condition contrast column to
  condition-level eRegulon AUC statistics.
- Widens fixed UMAP panel layout to avoid colorbar and label clipping.

## 0.1.20 - 2026-06-25

Reduces release archive size.

Highlights:

- Release tarballs include `archives/vendor.tar.gz` only, not the extracted
  `.vendor/` directory.
- Excludes Python bytecode caches and macOS `.DS_Store` files from release
  tarballs.

## 0.1.19 - 2026-06-25

Aligns the step-by-step guide with a successful end-to-end project run.

Highlights:

- Defaults pycisTopic topic modeling to the installed MALLET backend.
- Removes manual file-inspection snippets and over-detailed auto-resource
  tuning text from the step-by-step guides.
- Keeps the student-facing workflow focused on required commands, biological
  purpose, and output interpretation.

## 0.1.18 - 2026-06-25

Shortens user-facing workflow commands.

Highlights:

- Installs environment-local `spgrn-*` wrappers instead of long
  `scenicplus-grn-*` wrappers.
- Removes old long wrappers during installation so each environment exposes one
  command style.
- Updates English and Chinese README and step-by-step guides to use the short
  wrapper names.

## 0.1.7 - 2026-06-24

Adds resource-aware split motif enrichment for large SCENIC+ custom cisTarget databases.

Highlights:

- Adds `run_scenicplus_motif_enrichment_split.py`, which runs DEM and
  cisTarget motif enrichment as independent region-set-family chunks and then
  calls SCENIC+ `prepare_menr` with all resulting HDF5 files.
- The split workflow is resumable: existing non-empty chunk HDF5 outputs are
  skipped unless `--force` is used.
- Adds automatic `--max-parallel-chunks auto` planning based on current
  available memory, total memory, DB size and configured SCENIC+ CPU limits.
- Updates English and Chinese step-by-step guides to recommend the split
  motif-enrichment route for large custom ranking/scores databases.
- Adds one-command reports for annotated-object pre-export review and
  pycisTopic completion checks, replacing manual `cat`/`ls` audit snippets in
  the step-by-step guides.
- Split motif enrichment now clamps BLAS/OpenMP thread environment variables to
  the resolved stage worker count for more predictable laptop/server resource
  use.
- DEM chunks that produce no HDF5 under formal thresholds now auto-write an
  `.empty.tsv` marker, launch a relaxed-threshold diagnostic report automatically,
  and continue if at least one DEM/CTX HDF5 exists; all-DEM-empty runs now stop
  with a hard error and a pointer to the diagnostic report directory.

## 0.1.6 - 2026-06-24

Fixes large custom cisTarget database handoff into SCENIC+ motif enrichment.

Highlights:

- `build_custom_cistarget_db.py` now writes consensus FASTA headers as plain
  UCSC-style `chr:start-end` region names, avoiding duplicated
  `chr:start-end::chr:start-end` labels that pycistarget cannot parse.
- Adds `dem_n_cpu` and `ctx_n_cpu` to SCENIC+ config generation so DEM and
  cisTarget motif enrichment can use stage-specific worker counts.
- Automatically estimates safe motif-enrichment workers from detected machine
  memory and custom DB size, scaling down on laptops and allowing more workers
  on larger servers.
- Patches generated SCENIC+ Snakefiles to use the stage-specific motif
  enrichment thread settings.
- Restores complete eRegulon inference defaults in the packaged config
  template.

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
