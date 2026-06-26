# SCENIC+ Step-by-Step Workflow For Annotated Matched snRNA + snATAC Data

This guide is for matched snRNA + snATAC / 10x Multiome-style data after single-cell QC, integration, clustering and annotation have already been completed. It follows the official SCENIC+ sequence:
annotated single-cell object → ATAC fragment/peak standardization → pycisTopic cisTopic model → custom cisTarget database → SCENIC+ Snakemake inference.
It starts from an annotated object plus the original ATAC fragment and peak files and ends with the official SCENIC+ Snakemake workflow.

Official references:

```text
SCENIC+ documentation: https://scenicplus.readthedocs.io/en/latest/
SCENIC+ running tutorial: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html
pycisTopic tutorial: https://pycistopic.readthedocs.io/en/latest/tutorials.html
SCENIC+ GitHub: https://github.com/aertslab/scenicplus
Aerts cisTarget resources: https://resources.aertslab.org/cistarget/
```

## 0. Initialize Environment And Project

This step validates the installed conda environment, creates the project folder tree, updates the base project keys in `$PROJECT_DIR/scenicplus_project.env`, and writes `$PROJECT_DIR/project_env.sh`.

Before entering the command chain, run the installation/runtime integrity checks and confirm the execution mode matches the "matched scRNA+snATAC" path in official SCENIC+ docs. A workflow should only proceed if both checks pass.

Parameter reference:

```text
CONDA_ROOT                 Absolute path to the conda, miniforge,
                           miniconda, mambaforge or anaconda root directory.
ENV_NAME                   Dedicated SCENIC+ environment name.
PROJECT_DIR                Absolute path to the SCENIC+ analysis
                           root. Workflow-created inputs, work files, logs,
                           resources and results are written under this path.
ORGANISM                   One supported organism key.
AUTOZYME                   on or off.
ENSEMBL_RELEASE            Ensembl release used for genome resources.
ANNOTATED_OBJECT           Absolute path to the active annotated
                           Seurat RDS/QS or AnnData h5ad object.
CELL_LABEL_COLUMN          Metadata column in ANNOTATED_OBJECT used
                           as the SCENIC+ cell-group label.
ATAC_INPUT_LAYOUT          One supported ATAC input layout key.
ATAC_DATA_ROOT             Absolute path to the ATAC input root for
                           the selected ATAC_INPUT_LAYOUT.
```

Supported `ORGANISM` values:

```text
human      Homo sapiens GRCh38, chr1-chr22, chrX, chrY
mouse      Mus musculus GRCm39, chr1-chr19, chrX, chrY
cyno       Macaca fascicularis 6.0, chr1-chr20, chrX
rat        Rattus norvegicus GRCr8, chr1-chr20, chrX, chrY
rabbit     Oryctolagus cuniculus OryCun2.0, chr1-chr21, chrX
chicken    Gallus gallus GRCg7b, chr1-chr39, chrZ, chrW
zebrafish  Danio rerio GRCz11, chr1-chr25
```

Supported `ATAC_INPUT_LAYOUT` values:

```text
split_ge_arc
  ATAC_DATA_ROOT/
  |-- fragments/
  |   |-- sample_1_arc/
  |   |   |-- sample_1_A_fragments.tsv.gz
  |   |   `-- sample_1_A_fragments.tsv.gz.tbi
  |   `-- sample_2_arc/
  |       |-- sample_2_A_fragments.tsv.gz
  |       `-- sample_2_A_fragments.tsv.gz.tbi
  `-- bed/
      |-- sample_1_arc/
      |   `-- sample_1_A_peaks.bed
      `-- sample_2_arc/
          `-- sample_2_A_peaks.bed

cellranger_outs
  ATAC_DATA_ROOT/
  |-- sample_1/
  |   `-- outs/
  |       |-- fragments.tsv.gz
  |       |-- fragments.tsv.gz.tbi
  |       `-- peaks.bed
  `-- sample_2/
      `-- outs/
          |-- fragments.tsv.gz
          |-- fragments.tsv.gz.tbi
          `-- peaks.bed
```

First, enter the project parameters in the terminal:

```bash
export CONDA_ROOT=/absolute/path/to/conda
export ENV_NAME=scenicplus-grn
export PROJECT_DIR=/absolute/path/to/grn_project/scenicplus_analysis
export ORGANISM=mouse
export AUTOZYME=on
export ENSEMBL_RELEASE=115
export ANNOTATED_OBJECT=/absolute/path/to/active_annotated_multiome_object.rds
export CELL_LABEL_COLUMN=cell_annotation
export ATAC_INPUT_LAYOUT=split_ge_arc
export ATAC_DATA_ROOT=/absolute/path/to/atac_input_root
```

Pre-execution review gate:

```text
Confirm this is the intended analysis scope:
- An annotated scMultiome object is available as Seurat RDS/QS or AnnData h5ad.
- Original ATAC fragments and ATAC peak files are available for each sample.
- The intended analysis is matched scRNA+snATAC SCENIC+, not unmatched pseudo-integration.
- CONDA_ROOT, ENV_NAME, PROJECT_DIR, ANNOTATED_OBJECT and ATAC_DATA_ROOT are absolute paths.
- PROJECT_DIR is a dedicated SCENIC+ analysis root.
```

Run the installed environment and workflow integrity check:

```bash
mkdir -p "$PROJECT_DIR/logs"
"$CONDA_ROOT/envs/$ENV_NAME/bin/spgrn-check" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  2>&1 | tee "$PROJECT_DIR/logs/pre_step0_check_environment.log"
```

Then run the one-step initializer. It validates the parameter values, checks the conda environment, updates the project settings file, and initializes the project runtime files:

```bash
"$CONDA_ROOT/envs/$ENV_NAME/bin/spgrn-initialize"
```

Before running downstream steps, confirm every parameter here once (absolute paths, organism, layout, object path) and keep this file as the source of truth.

After the initializer reports `PROJECT INITIALIZATION OK`, load the project runtime variables:

```bash
source "$PROJECT_DIR/project_env.sh"
```

## 1. Prepare Organism Resources

Use the installed resource-preparation script instead of hand-writing genome download, chromosome allowlist, UCSC conversion, annotation-table, motif collection and motif2TF commands. The script is resumable, writes a log, writes `resources/resource_status.tsv`, and records all file paths and checksums in `resources/resource_manifest.json`.

All prepared genome resources are converted to UCSC chromosome style. The script keeps only standard primary chromosomes and filters random, unplaced, alt, haplotype and mitochondrial records from FASTA, GTF, chromsizes and the SCENIC+ genome annotation table. Use the same standard-chromosome policy later for fragments, peaks, consensus peaks, region_sets BED files and cisTarget DB regions.

For supported organisms without an Aerts public motif2TF table, the script downloads the official human HGNC v10 motif2TF table and maps its TF gene names to the target organism through Ensembl BioMart one-to-one orthology. This is not de novo motif discovery. The cached orthology table and the mapping audit are recorded in `resources/resource_manifest.json`.

Prepare the project organism recorded in `$PROJECT_DIR/project_env.sh`:

```bash
spgrn-prepare-official-resources
```

If UCSC-standard genome resources already exist and only the cisTarget motif collection or motif2TF table is missing, use:

```bash
spgrn-prepare-official-resources --motifs-only
```

The `--motifs-only` command prepares the motif collection, species motif2TF table and project `inputs/cistarget_db/motif_annotations.tbl` without rebuilding UCSC FASTA, GTF, chromsizes or genome annotation files.

After either preparation command succeeds, reload `project_env.sh`. The resource-preparation script updates resource-derived values such as `MACS_GENOME_SIZE`.

```bash
source "$PROJECT_DIR/project_env.sh"
```

Use `--mode check` to check existing resources without rebuilding. Use `--mode status` to print and log a resource status summary.

```bash
spgrn-prepare-official-resources --mode check

spgrn-prepare-official-resources --mode status
```

Prepare all supported organisms only when building a reusable shared resource cache. In `--organism all` mode, the script prepares organism resource caches but keeps the project organism selected in Step 0 as the active analysis organism.

```bash
spgrn-prepare-official-resources \
  --organism all \
  --mode prepare
```

## 2. Inspect And Export The Annotated Object

SCENIC+ needs RNA counts and metadata for the active cells. Start from an annotated scMultiome object, such as Seurat RDS/QS or AnnData h5ad.

Run all commands after sourcing the project environment:

```bash
source "$PROJECT_DIR/project_env.sh"
```

The annotated object must contain:

```text
RNA raw-count assay or layer inside the annotated object
sample column
condition column
cell-type or cell-state annotation column
barcodes that can be matched to the fragment files
```

First inspect the annotated object. This object should already contain the cells intended for the SCENIC+ run and an embedding computed for those cells. For example, if the analysis is restricted to a selected tissue, lineage or developmental compartment, provide the annotated object produced after that subset was reprocessed, not the whole-atlas object.

The unified inspector detects the file format from the suffix (`.rds`, `.qs`, or `.h5ad`), then detects candidate sample, condition, label, barcode, RNA assay, raw-count layer and metacell embedding fields where available.

```bash
spgrn-inspect-annotated-object
```

Inspect the active analysis object before choosing `CELL_LABEL_COLUMN`.

Review the report and generated parameter table with one command:

```bash
spgrn-review-annotated-object-inspection
```

The review writes `results/annotated_object/annotated_object_pre_export_review.md` and `.tsv`, then fails if required export fields are missing.

Confirm `cell_label_column`, `assay`, `layer` and `reduction`. `cell_label_column` is the user-facing grouping choice for this step. During export, this source column is written into `inputs/cell_metadata.tsv` as the fixed downstream column `cell_label`. The other fields are detected from the object and recorded because downstream scripts need them: `assay` and `layer` define the RNA count matrix exported from the annotated object to `inputs/gex.h5ad`; sample and condition columns define biological grouping; barcode fields connect the annotated object to ATAC fragments; `reduction` becomes the coordinate system for Step 3 metacell aggregation.

For Seurat input, the generated parameter table is `inputs/annotated_object_params.tsv`; for h5ad input, it is `inputs/annotated_h5ad_params.tsv`.

If any detected field is wrong, edit the generated parameter table before export.

If `CELL_LABEL_COLUMN` is not the intended grouping, pass the correct field at export time:

```bash
spgrn-export-annotated-object --cell-label-column corrected_metadata_column
```

Export the active RNA and metadata:

```bash
spgrn-export-annotated-object
```

Export keeps all cells in the provided active object. For a selected-population analysis, provide the selected, reprocessed object as `ANNOTATED_OBJECT`; do not subset a whole-atlas object at this step.

The export command records `ANNOTATED_OBJECT`, `ANNOTATED_OBJECT_FORMAT`, `CELL_LABEL_COLUMN`, `ACTIVE_GEX_H5AD` and `ACTIVE_CELL_METADATA` in `$PROJECT_DIR/scenicplus_project.env`.

`h5Seurat`, `loom`, and other formats should be converted to Seurat RDS/QS or AnnData h5ad before this step.

Step 2 writes:

```text
inputs/cell_metadata.tsv
inputs/gex.h5ad
results/annotated_object/cells_by_sample_and_label.tsv
results/annotated_object/annotated_object_summary.tsv
inputs/grn_label_summary.tsv
```

Downstream pycisTopic must be rebuilt from the active cells and their original fragments.

Review `inputs/grn_label_summary.tsv` before running pycisTopic.

Create or update the workflow parameter tables used by the following steps:

```bash
spgrn-setup-workflow-params
```

The script writes default project parameter tables under `inputs/` and records their paths in `$PROJECT_DIR/scenicplus_project.env`. Existing scalar parameter values are kept unless `--force` is used. Sample-specific tables are refreshed when their `sample_id` values no longer match the active sample sheet. To change a value, either edit the corresponding `inputs/*_params.tsv` file or use a command-line override:

```bash
spgrn-setup-workflow-params \
  --section pycistopic \
  --set pycistopic.n_iter=300
```

## 3. Metacell Aggregation

For large matched snRNA+snATAC projects, create metacells within `sample_id × cell_label` before pycisTopic and SCENIC+. This follows the metacell-first strategy used to reduce memory load while preserving annotated cell states.

The installed metacell implementation uses hdWGCNA and requires a Seurat RDS/QS object. If Step 2 started from h5ad, provide a matched Seurat RDS/QS in `inputs/metacell_params.tsv` before running this step.

Step 3 initializes `assay`, `layer` and `reduction` from Step 2. For matched multiome analyses, use a WNN UMAP computed for the active analysis cells. Use RNA UMAP only when the ATAC modality is unavailable or deliberately excluded.

```bash
spgrn-setup-workflow-params --section metacell
spgrn-prepare-metacell-inputs-from-seurat
spgrn-make-metacell-gex-h5ad
```

Metacell aggregation is required in this workflow before pycisTopic and SCENIC+.


Default metacell parameters are stored in:

```text
inputs/metacell_params.tsv
```

To change them, edit that file or rerun the setup script with `--set`, for example:

```bash
spgrn-setup-workflow-params \
  --section metacell \
  --set metacell.k=12 \
  --set metacell.max_shared=2
```

After the first two commands, `inputs/cell_metadata.tsv` and `inputs/gex.h5ad` refer to metacells. The original single-cell metadata is retained as:

```text
inputs/cell_metadata.single_cell.tsv
inputs/metacell_membership.tsv
results/metacells/metacell_summary.tsv
```

## 4. Prepare Active ATAC Fragments And Peaks

Create `$PROJECT_DIR/inputs/sample_sheet.tsv` from the ATAC file layout selected in Step 0.

If `ATAC_INPUT_LAYOUT` and `ATAC_DATA_ROOT` were set in Step 0, register them with:

```bash
spgrn-set-atac-input-params
```

Then generate and validate the active ATAC sample sheet:

```bash
spgrn-make-sample-sheet-from-atac-inputs
spgrn-validate-and-prepare-sample-sheet
```

`make_sample_sheet_from_atac_inputs.py` dispatches layout parsing based on `inputs/atac_input_params.tsv` produced in the previous command.


The ATAC input parameters are stored in `inputs/atac_input_params.tsv` and recorded in `$PROJECT_DIR/scenicplus_project.env`. The generated sample sheet is:

```text
inputs/sample_sheet.tsv
```

Review the `condition` column before continuing.

Run one command to standardize ATAC peaks and fragments, create a tabix index for each standardized fragment file, and update the active sample sheet to point to the standardized inputs:

```bash
spgrn-standardize-atac-inputs
```

Reassign the standardized single-cell fragments to the metacell barcodes produced in Step 3:

```bash
spgrn-reassign-fragments-to-metacells
```

Validate the updated active sample sheet after fragment reassignment:

```bash
spgrn-validate-and-prepare-sample-sheet
```

## 5. Create `inputs/cistopic_obj.pkl` With pycisTopic

Build the ATAC object from UCSC-standard fragments and a project-specific consensus peak universe. This step starts from `inputs/sample_sheet.tsv` and `inputs/cell_metadata.tsv`. The intended cell-name convention is `barcode-sample_id` on both the RNA and cisTopic sides.

`cistopic_obj.pkl`, consensus peaks, topic region sets and DAR region sets are specific to the active cells or metacells and active `cell_label` grouping. If the analysis changes from whole atlas to a selected population, or from cells to metacells, rebuild this step from the selected barcodes and original fragments. Subsetting an old whole-atlas cisTopic object is only appropriate for inspecting or plotting the already inferred atlas-level model.

Create or update QC, peak-calling, DAR and topic-model parameter files. The `chromsizes` and `genome_size` values come from the public resources prepared for `ORGANISM`; sample-level QC rows are generated from the active `inputs/sample_sheet.tsv`.

The pycisTopic pseudobulk step can be slow and memory-heavy. Keep `resume_pseudobulk=1` unless there is a reason to rebuild all pseudobulk BED files from scratch. On rerun, the workflow reuses existing non-empty gzip-valid BED files and rebuilds only missing or corrupted label files. A run interrupted during gzip writing may leave a partial `.bed.gz`; the workflow tests gzip integrity before reuse.

```bash
spgrn-setup-workflow-params --section pycistopic
```

This writes or updates:

```text
inputs/atac_qc_thresholds.tsv
inputs/topic_model_grid.tsv
inputs/pycistopic_params.tsv
```

To change the topic grid or pycisTopic parameters, use explicit overrides:

```bash
spgrn-setup-workflow-params \
  --section pycistopic \
  --set topic_grid.n_topics=10,20,30 \
  --set pycistopic.n_iter=300
```

The default topic-model backend is MALLET, because it is faster and more stable for large peak-by-cell/metacell matrices. The workflow records the active MALLET runtime in `results/pycistopic/model_selection/mallet_runtime.tsv`.

If switching between CGS and MALLET, keep all topic models from the same backend. Archive or remove interrupted `Topic*.pkl` files before rerunning topic modeling.

Because Step 3 produces metacell inputs, the default `analysis_unit` is `metacell`. Set `pycistopic.analysis_unit=cell` only for a cell-level run.

Run:

```bash
spgrn-run-pycistopic-workflow 2>&1 | tee "$PROJECT_DIR/logs/run_pycistopic_workflow.log"
```

If a pseudobulk run was interrupted, validate existing pseudobulk files before continuing. Remove only files reported as invalid; valid files will be reused.

```bash
spgrn-validate-pseudobulk-files
```

Check that pycisTopic produced all files required by SCENIC+:

```bash
spgrn-check-pycistopic-completion
```

The check writes `results/pycistopic/qc/pycistopic_completion_check.md` and `.tsv`, then fails if any required output is missing or empty.

Check that the cisTopic object and active RNA matrix use matching cell IDs:

```bash
spgrn-check-cistopic-cell-names
```

## 6. Create `inputs/region_sets/`

Required structure produced by Step 5:

```text
inputs/region_sets/
  Topics_otsu/Topic1.bed
  Topics_top_3k/Topic1.bed
  DARs_cell_label/CellLabelA_VS_rest.bed
```

Validate and standardize all region-set BED files:

```bash
spgrn-standardize-region-sets
```

The script writes:

```text
results/pycistopic/qc/region_set_standardization.tsv
```

## 7. Build The Custom cisTarget Region Database

The custom cisTarget database must be built from the same UCSC-standard consensus region universe used by pycisTopic and region sets. Create or update the cisTarget database parameters, then build the database. This step is required for final SCENIC+ inference:

```bash
spgrn-setup-workflow-params --section cistarget
spgrn-build-custom-cistarget-db
```

Parameters are stored in:

```text
inputs/cistarget_db_params.tsv
```

The script manages resource use and resumable execution. The resolved plan is
written before the long motif-scanning step starts:

```text
results/cistarget_db/custom_cistarget_resource_plan.tsv
```

For large region-by-motif matrices, the wrapper uses the official partial
cisTarget workflow. It runs `create_cistarget_motif_databases.py --partial`,
skips existing non-empty part files on rerun, combines partial score databases
with the official combine script, then converts motif-vs-region scores to the
final regions-vs-motifs rankings database with the official conversion script.

Long Cluster-Buster scans may not print motif-level progress until one or more
workers return. The wrapper writes a lightweight heartbeat while the official
script is still running. The interval is controlled by `heartbeat_seconds` in
`inputs/cistarget_db_params.tsv` and defaults to 600 seconds. Set it higher for
quiet server logs or lower when actively debugging.

Partial jobs are run sequentially by default to avoid duplicating the large
region-by-motif score structure across concurrent partial processes.

Expected:

```text
inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather
inputs/cistarget_db/custom.regions_vs_motifs.scores.feather
results/cistarget_db/custom_cistarget_db_manifest.tsv
```

## 8. Initialize SCENIC+ And Generate Config

SCENIC+ is run on all active metacells together to build one shared eRegulon universe. If the project contains multiple conditions, keep them together at this step; condition-specific effects are tested later on the shared eRegulon AUC matrix.

```bash
spgrn-initialize-scenicplus-snakemake
```

This initializes the SCENIC+ Snakemake directory, copies the organism-specific annotation files, writes `work/scenicplus/organism_config.yaml`, creates or updates `inputs/scenicplus_config_params.tsv`, and generates the final Snakemake config:

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

Parameters usually worth reviewing:

```text
seed                      Reproducibility seed.
search_space_upstream     Gene search-space upstream window.
search_space_downstream   Gene search-space downstream window.
search_space_extend_tss   TSS extension window.
dem_motif_hit_thr         DEM motif hit threshold; keep 3.0 unless justified.
ctx_nes_threshold         cisTarget NES threshold; keep 3.0 unless justified.
rho_threshold             Region-gene correlation threshold.
min_target_genes          Minimum target genes retained per eRegulon.
```

8.2-Optional split motif enrichment for very large custom databases:

Very large ranking and score databases can still make a single motif-enrichment
process retain too much memory while it iterates over DARs and topic-derived
region sets. In that case, run motif enrichment as independent region-set-family
chunks before the formal Snakemake run:

```bash
spgrn-run-scenicplus-motif-enrichment-split --mode both
spgrn-run-scenicplus-motif-enrichment-split --mode status
```

This command creates `inputs/region_sets_split/`, runs DEM and cisTarget on each
region-set family as separate SCENIC+ processes, then calls `prepare_menr` with
all resulting HDF5 files. `--mode status` is the required completion gate before
Step 9. It verifies valid HDF5 signatures for completed chunk outputs, records
DEM chunks that are formally empty only when their diagnostic report exists, and
checks that `prepare_menr` produced valid `cistromes_direct.h5ad`,
`cistromes_extended.h5ad` and non-empty `tf_names.txt`.

It writes a resource plan, chunk manifest and completion reports here:

```text
work/scenicplus/motif_enrichment_split/motif_enrichment_split_resource_plan.tsv
work/scenicplus/motif_enrichment_split/motif_enrichment_split_chunks.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.md
```

The command is resumable: completed non-empty chunk HDF5 files are skipped unless
`--force` is used.

If a DEM chunk completes but produces no HDF5 under the configured thresholds,
the runner records an `.empty.tsv` marker, automatically launches a
relaxed-threshold diagnostic for that chunk, and lets `--mode status` decide
whether the split run is safe to use.

The diagnostic output is written to:

```text
work/scenicplus/motif_enrichment_split/ (chunk-relaxed DEM outputs and logs)
results/scenicplus_diagnostics/*_relaxed_threshold_diagnostic.tsv/md
```

If all DEM chunks are empty (i.e., no formal-threshold DEM HDF5 is generated),
the workflow stops with a hard error and points to the full diagnostic reports.
If `--mode status` fails for any missing, invalid or half-written output, do not
continue to Step 9. Rerun the split command; use `--force` only when replacing
old partial outputs intentionally.

If these values need to change, edit `$PROJECT_DIR/inputs/scenicplus_config_params.tsv` or rerun the setup script with `--set`, then rerun the generator. If paths or advanced SCENIC+ options need to change, edit the generated file:

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

Cell-name overlap, UCSC chromosome consistency and core SCENIC+ config invariants are checked by Step 9 before Snakemake is allowed to run.

If Step 9 fails, do not run Step 10.


## 9. Preflight Checks Before SCENIC+

Run all checks before the Snakemake dry run. These checks should fail loudly but use project-configurable thresholds.

```bash
spgrn-setup-workflow-params --section preflight
spgrn-preflight-scenicplus-inputs 2>&1 | tee "$PROJECT_DIR/logs/preflight_scenicplus_inputs.log"
```

The thresholds are recorded in `inputs/preflight_thresholds.tsv` and can be adjusted if your data are noisy, without changing output names.


Record software versions and resource checksums:

```bash
spgrn-record-scenicplus-provenance
```

## 10. Dry Run, Run, And Stability Record

```bash
spgrn-setup-workflow-params --section snakemake
spgrn-run-scenicplus-snakemake --mode dryrun
```

```bash
spgrn-run-scenicplus-snakemake --mode run
```

Key parameter file: `inputs/snakemake_params.tsv`.

Main outputs:

```text
results/scenicplus/eRegulons_direct.tsv
results/scenicplus/eRegulons_extended.tsv
results/scenicplus/AUCell_direct.h5mu
results/scenicplus/AUCell_extended.h5mu
results/scenicplus/scplus_mdata.h5mu
```

For formal analyses, keep two independent inference records to assess high-confidence edge stability. Archive the first completed run:

```bash
mkdir -p results/scenicplus_stability/run1
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run1/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run1/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run1/
```

Second independent inference:

```bash
spgrn-run-scenicplus-snakemake --mode dryrun --rerun-inference
spgrn-run-scenicplus-snakemake --mode run --rerun-inference
```

`--rerun-inference` reruns the SCENIC+ inference chain without rebuilding genome resources. Archive the second completed run:

```bash
mkdir -p results/scenicplus_stability/run2
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run2/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run2/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run2/
```

Compare edge stability:

```bash
spgrn-compare-scenicplus-stability \
  --run-a results/scenicplus_stability/run1 \
  --run-b results/scenicplus_stability/run2 \
  --out results/scenicplus_stability/stability_summary.tsv
```

Report `results/scenicplus_stability/stability_summary.tsv` with the final GRN results.

## 11. SCENIC+ Postprocess Figures And Condition Statistics

After Step 10, run postprocessing as one ordered block. This stage converts the shared SCENIC+ eRegulon universe into source tables, vector PDFs, an output-tier audit and sample-level condition statistics. Keep the inference step and the condition test conceptually separate: SCENIC+ is inferred once on all active metacells, then eRegulon AUCell activity is compared across biological samples and conditions.

11.1-Create or update postprocessing parameters:

```bash
spgrn-setup-workflow-params --section postprocess
```

11.2-Generate the output-tier audit, all standard PDFs, source tables and condition statistics:

```bash
spgrn-run-scenicplus-postprocess --task all --layer all
```

Step 11 is kept as one postprocess stage because the audit, figure source tables, PDFs and condition statistics all describe the same completed SCENIC+ inference result.

Postprocess command options:

```text
--task audit    Output-tier audit only.
--task figures  Regenerate source tables and 01-08 eRegulon PDFs.
--task stats    Regenerate 09+ condition-statistics tables and PDFs.
--task all      Generate the complete postprocess output set.

--layer direct    Use direct eRegulons only.
--layer extended  Use extended eRegulons only.
--layer all       Run both direct and extended eRegulons.
```

`--task` controls the analysis block and `--layer` controls the eRegulon layer. The installed wrapper does not target a single named PDF; rerun the relevant task/layer, or edit the source table and custom-draw a focused manuscript panel separately.

Main output locations:

```text
results/scenicplus_output_tiers/       Output-tier audit TSV/PDF for checking which result families are interpretable.
results/scenicplus_figures/01-05_*     eRegulon activity, specificity, dot heatmap and embedding PDFs.
results/scenicplus_figures/06-08_*     Region-gene, target-region overlap and TF-target network PDFs.
results/scenicplus_figures/09-15_*     Sample-level condition statistics PDFs.
results/scenicplus_figures/source_*    Figure source tables for custom manuscript redraws.
results/scenicplus_figures/*stats*.tsv Condition-statistics source tables with the resolved contrast direction.
logs/                                  Postprocess logs for each extraction, R render and condition-statistics step.
```

R plotting style is controlled by:

```text
results/scenicplus_figures/plot_style_parameters.tsv
results/scenicplus_figures/plot_style_parameters.md
```

`plot_style_parameters.tsv` is the editable file. It contains only `parameter` and `value`, so it can be reused across projects without project-specific labels. The Markdown file is generated with the parameter dictionary and the one-command redraw instruction. After editing colors, font sizes, line widths, point sizes, alpha values or panel layout values, rerun 11.2 to redraw all PDFs.

Optional priority labels for condition volcano plots can be supplied through `priority_eregulons` in `inputs/postprocess_params.tsv`. Use a TSV with `eregulon` or `display_label`; `cell_label` and `layer` columns are respected when present.

Condition statistics are automatic. One condition gives descriptive sample means only; two conditions report `delta_mean_auc` as comparison minus reference and record the resolved direction in `contrast`; more than two conditions use an omnibus sample-level test without a single delta. Set `reference_condition` or `comparison_condition` in `inputs/postprocess_params.tsv` only when the automatic two-condition reference is not the intended baseline.

Scientific output layers:

```text
0. Output-tier audit and input confidence: verifies which SCENIC+ outputs are present and which output families can be interpreted.
1. eRegulon activity and specificity: heatmap/dot/embedding views for cell-state regulatory programs.
2. Condition-resolved eRegulon activity: the same eRegulon universe split by cell label and condition.
3. Region-gene and target-region structure: model-level enhancer-target and region-overlap summaries.
4. TF-target network overview: compact network view for selected TF-target links.
5. Sample-level condition statistics: condition effects from biological-sample mean AUCell scores.
```

Interpretation rules:

```text
AUC heatmap: use for group-level eRegulon activity patterns; rows and columns are clustered from the plotted numeric matrix.
Condition AUC heatmap: compare condition shifts within each cell label; columns are ordered as cell_label_condition and rows are clustered.
Specificity heatmap: use for prioritizing cell-state-specific eRegulons; rows and columns are clustered from the plotted numeric matrix.
Dot heatmap: color encodes relative AUC; dot size encodes the fraction of active cells.
UMAP/activity embedding: use for localization and state separation, not as standalone statistical evidence.
Region-gene and overlap PDFs: model-level structural views; interpret condition effects from AUCell statistics, not from these structural plots alone.
Network PDF: use as a compact overview; use source edge tables for final focused network panels.
Condition volcano and condition heatmap: interpret sample-level effect size first; report FDR only when there are enough independent biological samples.
eRegulon signs: labels retain SCENIC+ signs such as +/+, -/+ and -/-. Do not collapse these into TF RNA expression; eRegulon AUC is target gene-set activity.
```

For a main figure, redraw selected eRegulons from the source tables rather than showing every detected regulon. Keep direct and extended eRegulons separated unless the figure explicitly states how they were merged.
