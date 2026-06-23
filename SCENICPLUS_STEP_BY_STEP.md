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
bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/check_environment.sh" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  2>&1 | tee "$PROJECT_DIR/logs/pre_step0_check_environment.log"
```

Then run the one-step initializer. It validates the parameter values, checks the conda environment, updates the project settings file, and initializes the project runtime files:

```bash
bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/initialize_scenicplus_project.sh"
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
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py
```

If UCSC-standard genome resources already exist and only the cisTarget motif collection or motif2TF table is missing, use:

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py --motifs-only
```

The `--motifs-only` command prepares the motif collection, species motif2TF table and project `inputs/cistarget_db/motif_annotations.tbl` without rebuilding UCSC FASTA, GTF, chromsizes or genome annotation files.

After either preparation command succeeds, reload `project_env.sh`. The resource-preparation script updates resource-derived values such as `MACS_GENOME_SIZE`.

```bash
source "$PROJECT_DIR/project_env.sh"
```

Use `--mode check` to check existing resources without rebuilding. Use `--mode status` to print and log a resource status summary.

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py --mode check

python $SCENICPLUS_HOME/scripts/prepare_official_resources.py --mode status
```

Prepare all supported organisms only when building a reusable shared resource cache. In `--organism all` mode, the script prepares organism resource caches but keeps the project organism selected in Step 0 as the active analysis organism.

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py \
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
python $SCENICPLUS_HOME/scripts/inspect_annotated_object.py
```

Inspect the active analysis object before choosing `CELL_LABEL_COLUMN`.

Review the report and the generated parameter table:

```bash
cat "$PROJECT_DIR/results/annotated_object/annotated_object_inspection_report.md"
if [[ -f "$PROJECT_DIR/inputs/annotated_object_params.tsv" ]]; then
  cat "$PROJECT_DIR/inputs/annotated_object_params.tsv"
elif [[ -f "$PROJECT_DIR/inputs/annotated_h5ad_params.tsv" ]]; then
  cat "$PROJECT_DIR/inputs/annotated_h5ad_params.tsv"
else
  echo "Missing annotated_object params file" >&2
  exit 1
fi
```

Additional previews are written under `results/annotated_object/`.

Confirm `cell_label_column`, `assay`, `layer` and `reduction`. `cell_label_column` is the user-facing grouping choice for this step. During export, this source column is written into `inputs/cell_metadata.tsv` as the fixed downstream column `cell_label`. The other fields are detected from the object and recorded because downstream scripts need them: `assay` and `layer` define the RNA count matrix exported from the annotated object to `inputs/gex.h5ad`; sample and condition columns define biological grouping; barcode fields connect the annotated object to ATAC fragments; `reduction` becomes the coordinate system for Step 3 metacell aggregation.

For Seurat input, the generated parameter table is `inputs/annotated_object_params.tsv`; for h5ad input, it is `inputs/annotated_h5ad_params.tsv`.

If any detected field is wrong, edit the generated parameter table before export.

If `CELL_LABEL_COLUMN` is not the intended grouping, pass the correct field at export time:

```bash
python $SCENICPLUS_HOME/scripts/export_annotated_object.py --cell-label-column corrected_metadata_column
```

Export the active RNA and metadata:

```bash
python $SCENICPLUS_HOME/scripts/export_annotated_object.py
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
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py
```

The script writes default project parameter tables under `inputs/` and records their paths in `$PROJECT_DIR/scenicplus_project.env`. Existing scalar parameter values are kept unless `--force` is used. Sample-specific tables are refreshed when their `sample_id` values no longer match the active sample sheet. To change a value, either edit the corresponding `inputs/*_params.tsv` file or use a command-line override:

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section pycistopic \
  --set pycistopic.n_iter=300
```

## 3. Metacell Aggregation

For large matched snRNA+snATAC projects, create metacells within `sample_id × cell_label` before pycisTopic and SCENIC+. This follows the metacell-first strategy used to reduce memory load while preserving annotated cell states.

The installed metacell implementation uses hdWGCNA and requires a Seurat RDS/QS object. If Step 2 started from h5ad, provide a matched Seurat RDS/QS in `inputs/metacell_params.tsv` before running this step.

Step 3 initializes `assay`, `layer` and `reduction` from Step 2. For matched multiome analyses, use a WNN UMAP computed for the active analysis cells. Use RNA UMAP only when the ATAC modality is unavailable or deliberately excluded.

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section metacell
Rscript $SCENICPLUS_HOME/scripts/prepare_metacell_inputs_from_seurat.R
python $SCENICPLUS_HOME/scripts/make_metacell_gex_h5ad.py
```

Metacell aggregation is required in this workflow before pycisTopic and SCENIC+.


Default metacell parameters are stored in:

```text
inputs/metacell_params.tsv
```

To change them, edit that file or rerun the setup script with `--set`, for example:

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
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
python $SCENICPLUS_HOME/scripts/set_atac_input_params.py
```

Then generate and validate the active ATAC sample sheet:

```bash
python $SCENICPLUS_HOME/scripts/make_sample_sheet_from_atac_inputs.py
python $SCENICPLUS_HOME/scripts/validate_and_prepare_sample_sheet.py
```

`make_sample_sheet_from_atac_inputs.py` dispatches layout parsing based on `inputs/atac_input_params.tsv` produced in the previous command.


The ATAC input parameters are stored in `inputs/atac_input_params.tsv` and recorded in `$PROJECT_DIR/scenicplus_project.env`. The generated sample sheet is:

```text
inputs/sample_sheet.tsv
```

Review the `condition` column before continuing.

Run one command to standardize ATAC peaks and fragments, create a tabix index for each standardized fragment file, and update the active sample sheet to point to the standardized inputs:

```bash
python $SCENICPLUS_HOME/scripts/standardize_atac_inputs.py
```

Reassign the standardized single-cell fragments to the metacell barcodes produced in Step 3:

```bash
python $SCENICPLUS_HOME/scripts/reassign_fragments_to_metacells.py
```

Validate the updated active sample sheet after fragment reassignment:

```bash
python $SCENICPLUS_HOME/scripts/validate_and_prepare_sample_sheet.py
```

## 5. Create `inputs/cistopic_obj.pkl` With pycisTopic

Build the ATAC object from UCSC-standard fragments and a project-specific consensus peak universe. This step starts from `inputs/sample_sheet.tsv` and `inputs/cell_metadata.tsv`. The intended cell-name convention is `barcode-sample_id` on both the RNA and cisTopic sides.

`cistopic_obj.pkl`, consensus peaks, topic region sets and DAR region sets are specific to the active cells or metacells and active `cell_label` grouping. If the analysis changes from whole atlas to a selected population, or from cells to metacells, rebuild this step from the selected barcodes and original fragments. Subsetting an old whole-atlas cisTopic object is only appropriate for inspecting or plotting the already inferred atlas-level model.

Create or update QC, peak-calling, DAR, doublet and topic-model parameter files. The `chromsizes` and `genome_size` values come from the public resources prepared for `ORGANISM`; sample-level QC rows are generated from the active `inputs/sample_sheet.tsv`.

The pycisTopic pseudobulk step can be slow and memory-heavy. Keep `resume_pseudobulk=1` unless there is a reason to rebuild all pseudobulk BED files from scratch. On rerun, the workflow reuses existing non-empty gzip-valid BED files and rebuilds only missing or corrupted label files. A run interrupted during gzip writing may leave a partial `.bed.gz`; the workflow tests gzip integrity before reuse.

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section pycistopic
```

This writes or updates:

```text
inputs/atac_qc_thresholds.tsv
inputs/topic_model_grid.tsv
inputs/pycistopic_params.tsv
```

To change the topic grid or pycisTopic parameters, use explicit overrides:

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section pycistopic \
  --set topic_grid.n_topics=10,20,30 \
  --set pycistopic.n_iter=300
```

The default topic-model backend is `cgs`, which is dependency-light and fully contained in pycisTopic. For large peak-by-cell/metacell matrices, use the pycisTopic MALLET backend as an explicit backend choice.

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section pycistopic \
  --set pycistopic.lda_backend=mallet \
  --set pycistopic.mallet_path="$CONDA_ENV_PREFIX/bin/mallet" \
  --set pycistopic.mallet_memory_gb=auto
```

If MALLET is used, keep all topic models from the same backend. Archive or remove interrupted CGS `Topic*.pkl` files before rerunning with MALLET. `mallet_memory_gb=auto` detects available memory and sets Java heap through `MALLET_MEMORY`. The workflow records the active MALLET runtime in `results/pycistopic/model_selection/mallet_runtime.tsv`.

`topic_n_cpu=auto` is resolved separately from pseudobulk and MACS2 workers. Set an explicit value, for example `3`, when MALLET topic modeling is CPU-bound and memory is stable.

Because Step 3 produces metacell inputs, the default `analysis_unit` is `metacell`. Set `pycistopic.analysis_unit=cell` only for a cell-level run.

Run:

```bash
python $SCENICPLUS_HOME/scripts/run_pycistopic_workflow.py 2>&1 | tee "$PROJECT_DIR/logs/run_pycistopic_workflow.log"
```

Inspect the automatic worker plan and pseudobulk resume plan:

```bash
cat results/pycistopic/qc/parallelism_plan.tsv
cat results/pycistopic/qc/pseudobulk_resume_plan.tsv
```

If a pseudobulk run was interrupted, validate existing pseudobulk files before continuing. Remove only files reported as invalid; valid files will be reused.

```bash
python $SCENICPLUS_HOME/scripts/validate_pseudobulk_files.py
```

The workflow is complete only when these outputs exist (all are required for SCENIC+):

```text
inputs/cistopic_obj.pkl
work/pycistopic/consensus_peaks.bed
inputs/region_sets/Topics_otsu/*.bed
inputs/region_sets/Topics_top_3k/*.bed
inputs/region_sets/DARs_cell_label/*.bed
results/pycistopic/consensus_peaks.tsv
results/pycistopic/qc/cell_qc_metrics.tsv
results/pycistopic/qc/cistopic_cell_qc.pdf
results/pycistopic/qc/consensus_peak_qc.pdf
results/pycistopic/qc/atac_doublets.tsv
results/pycistopic/qc/atac_doublet_diagnostics.pdf
results/pycistopic/qc/cells_by_sample_and_label.tsv
results/pycistopic/qc/parallelism_plan.tsv
results/pycistopic/qc/pseudobulk_resume_plan.tsv
results/pycistopic/dar/*.tsv
results/pycistopic/dar/dar_summary.pdf
results/pycistopic/model_selection/topic_model_metrics.tsv
results/pycistopic/model_selection/topic_model_metrics.pdf
results/pycistopic/model_selection/topic_qc_metrics.tsv
results/pycistopic/model_selection/topic_qc_metrics.pdf
results/pycistopic/model_selection/topic_otsu_thresholds.pdf
results/pycistopic/model_selection/topic_region_set_summary.tsv
results/pycistopic/model_selection/topic_region_set_summary.pdf
results/pycistopic/model_selection/selected_model.txt
results/pycistopic/pycistopic_manifest.json
```

Check that the cisTopic object and active RNA matrix use matching cell IDs:

```bash
python $SCENICPLUS_HOME/scripts/check_cistopic_cell_names.py
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
python $SCENICPLUS_HOME/scripts/standardize_region_sets.py
```

The script writes:

```text
results/pycistopic/qc/region_set_standardization.tsv
```

## 7. Build The Custom cisTarget Region Database

The custom cisTarget database must be built from the same UCSC-standard consensus region universe used by pycisTopic and region sets. Create or update the cisTarget database parameters, then build the database. This step is required for final SCENIC+ inference:

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section cistarget
python $SCENICPLUS_HOME/scripts/build_custom_cistarget_db.py
```

Parameters are stored in:

```text
inputs/cistarget_db_params.tsv
```

For large region universes, keep `n_cpu=auto` unless a scheduler requires a fixed
thread count. At launch time, `build_custom_cistarget_db.py` inspects CPU count,
current load and available memory, then resolves a conservative worker count for
`create_cistarget_motif_databases.py -t`. The resolved plan is written before
the long motif-scanning step starts:

```text
results/cistarget_db/custom_cistarget_resource_plan.tsv
```

For large region-by-motif matrices, `use_partial=auto` enables the official
partial cisTarget workflow automatically. The wrapper runs
`create_cistarget_motif_databases.py --partial`, skips existing non-empty part
files on rerun, combines partial score databases with the official combine
script, then converts motif-vs-region scores to the final
regions-vs-motifs rankings database with the official conversion script. The
resource plan records `use_partial`, `partial_n_parts` and the estimated full
float32 score matrix size.

Expected:

```text
inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather
inputs/cistarget_db/custom.regions_vs_motifs.scores.feather
results/cistarget_db/custom_cistarget_db_manifest.tsv
```

## 8. Initialize SCENIC+ And Generate Config

SCENIC+ is run on all active metacells together to build one shared eRegulon universe. If the project contains multiple conditions, keep them together at this step; condition-specific effects are tested later on the shared eRegulon AUC matrix.

```bash
python $SCENICPLUS_HOME/scripts/initialize_scenicplus_snakemake.py
```

This initializes the SCENIC+ Snakemake directory, copies the organism-specific annotation files, writes `work/scenicplus/organism_config.yaml`, creates or updates `inputs/scenicplus_config_params.tsv`, and generates the final Snakemake config:

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

Important parameters to review:

```text
n_cpu                     CPU threads used by SCENIC+.
seed                      Reproducibility seed.
nr_cells_per_metacells    Used by SCENIC+ for non-multiome linkage. In this
                          workflow `is_multiome` stays true because RNA and
                          ATAC cell IDs are matched; precomputed metacell IDs
                          from Step 3 are treated as the matched units.
search_space_upstream     Gene search-space upstream window.
search_space_downstream   Gene search-space downstream window.
search_space_extend_tss   TSS extension window.
dem_motif_hit_thr         DEM motif hit threshold; keep 3.0 unless justified.
ctx_nes_threshold         cisTarget NES threshold; keep 3.0 unless justified.
rho_threshold             Region-gene correlation threshold.
min_target_genes          Minimum target genes retained per eRegulon.
```

If these values need to change, edit `$PROJECT_DIR/inputs/scenicplus_config_params.tsv` or rerun the setup script with `--set`, then rerun the generator. If paths or advanced SCENIC+ options need to change, edit the generated file:

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

Cell-name overlap, UCSC chromosome consistency and core SCENIC+ config invariants are checked by Step 9 before Snakemake is allowed to run.

If Step 9 fails, do not run Step 10.


## 9. Preflight Checks Before SCENIC+

Run all checks before the Snakemake dry run. These checks should fail loudly but use project-configurable thresholds.

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section preflight
python $SCENICPLUS_HOME/scripts/preflight_scenicplus_inputs.py 2>&1 | tee "$PROJECT_DIR/logs/preflight_scenicplus_inputs.log"
```

The thresholds are recorded in `inputs/preflight_thresholds.tsv` and can be adjusted if your data are noisy, without changing output names.


Record software versions and resource checksums:

```bash
python $SCENICPLUS_HOME/scripts/record_scenicplus_provenance.py
```

## 10. Dry Run, Run, And Stability Record

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section snakemake
python $SCENICPLUS_HOME/scripts/run_scenicplus_snakemake.py --mode dryrun
```

```bash
python $SCENICPLUS_HOME/scripts/run_scenicplus_snakemake.py --mode run
```

Run parameters are stored in:

```text
inputs/snakemake_params.tsv
```

Expected outputs:

```text
results/scenicplus/eRegulons_direct.tsv
results/scenicplus/eRegulons_extended.tsv
results/scenicplus/AUCell_direct.h5mu
results/scenicplus/AUCell_extended.h5mu
results/scenicplus/scplus_mdata.h5mu
```

For formal analyses, run at least two independent SCENIC+ executions or one standard execution plus a deterministic patched execution, then compare the high-confidence edges. Copy each completed result directory before starting the next run:

```bash
mkdir -p results/scenicplus_stability/run1 results/scenicplus_stability/run2
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run1/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run1/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run1/

# After changing output paths, seed, or using a duplicated project directory for run2:
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run2/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run2/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run2/
```

Run:

```bash
python $SCENICPLUS_HOME/scripts/compare_scenicplus_stability.py \
  --run-a results/scenicplus_stability/run1 \
  --run-b results/scenicplus_stability/run2 \
  --out results/scenicplus_stability/stability_summary.tsv
```

The global SCENIC+ seed controls part of the workflow, but region-to-gene GBM can still be sensitive to implementation details and dependency versions. Report the stability summary with the final GRN results.

## 11. Condition-Level eRegulon AUC Statistics

Run differential eRegulon activity testing at the sample level on the shared eRegulon universe from Step 10. Do not compare separately inferred condition-specific eRegulon lists as the primary analysis. This script reads the SCENIC+ AUCell `.h5mu`, joins `inputs/cell_metadata.tsv`, averages AUC per sample, and tests condition differences across sample-level means. It always writes the statistical table and a PDF diagnostic page.

Create or update the postprocessing parameters, then run for direct and extended eRegulons:

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section postprocess
python $SCENICPLUS_HOME/scripts/run_scenicplus_postprocess.py --task stats --layer all
```

Required outputs:

```text
results/scenicplus_stats/auc_by_condition_direct/sample_mean_auc.tsv
results/scenicplus_stats/auc_by_condition_direct/condition_eregulon_auc_statistics.tsv
results/scenicplus_stats/auc_by_condition_direct/condition_eregulon_auc_statistics.pdf
results/scenicplus_stats/auc_by_condition_extended/sample_mean_auc.tsv
results/scenicplus_stats/auc_by_condition_extended/condition_eregulon_auc_statistics.tsv
results/scenicplus_stats/auc_by_condition_extended/condition_eregulon_auc_statistics.pdf
```

Report `delta_mean_auc` and FDR only when there are enough independent samples. If a condition has fewer than two samples, keep the sample-mean table and PDF, but do not interpret per-cell P values as biological replication.

## 12. SCENIC+ Result Figures

A complete SCENIC+ analysis should not stop at tables, but the figures should be organized by scientific evidence layer rather than by plotting convenience. First audit which evidence layers are ready:

```bash
python $SCENICPLUS_HOME/scripts/run_scenicplus_postprocess.py --task audit
```

The audit writes:

```text
results/scenicplus_output_tiers/scenicplus_output_tier_audit.tsv
results/scenicplus_output_tiers/scenicplus_output_tier_audit.pdf
```

Scientific output tiers:

```text
0. input/QC confidence: active RNA/ATAC cells, metadata, fragments, doublet and ATAC QC
1. chromatin topics: major accessibility states and topic model quality
2. region sets and DARs: topic-region programs and differential accessible regions
3. motif/cisTarget evidence: project-specific motif enrichment database and motif2TF mapping
4. eRegulon activity: TF-region-gene eRegulons and cell-state activity patterns
5. condition effects: sample-level differential eRegulon activity across conditions
6. mechanism views: focused TF-target networks and locus/coverage views when the required inputs exist
```

After the audit, generate the standard eRegulon figure set for each direct and extended eRegulon layer. The installer provides `plot_scenicplus_publication_outputs.py`, which writes vector PDFs and source data tables for:

```text
1. eRegulon AUC heatmap across cell states or labels
2. RSS-style eRegulon specificity heatmap
3. heatmap-dotplot: color = group AUC z-score, size = active-cell fraction
4. UMAP overlays for high-priority eRegulons
5. compact TF-target network view
```

The visual style is intentionally clean and journal-ready: white background, Arial/Helvetica fonts, low-saturation colors, small but readable labels, no heavy grid, and all panels exported as PDF. Each figure has a matching source table.

Run direct and extended eRegulon figure generation:

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section postprocess
python $SCENICPLUS_HOME/scripts/run_scenicplus_postprocess.py --task figures --layer all
```

Postprocessing parameters are stored in:

```text
inputs/postprocess_params.tsv
```

Required outputs for each layer:

```text
results/scenicplus_figures/<direct_or_extended>/eregulon_auc_heatmap.pdf
results/scenicplus_figures/<direct_or_extended>/eregulon_specificity_heatmap.pdf
results/scenicplus_figures/<direct_or_extended>/eregulon_dot_heatmap.pdf
results/scenicplus_figures/<direct_or_extended>/eregulon_auc_umap.pdf
results/scenicplus_figures/<direct_or_extended>/tf_target_network.pdf
results/scenicplus_figures/<direct_or_extended>/source_auc_mean_by_group.tsv
results/scenicplus_figures/<direct_or_extended>/source_auc_heatmap_zscore.tsv
results/scenicplus_figures/<direct_or_extended>/source_rss_specificity.tsv
results/scenicplus_figures/<direct_or_extended>/source_dot_heatmap.tsv
results/scenicplus_figures/<direct_or_extended>/source_tf_target_network_edges.tsv
results/scenicplus_figures/<direct_or_extended>/source_selected_top_eregulons.tsv
```

Interpretation rules:

```text
AUC heatmap: use for group-level eRegulon activity patterns.
Specificity heatmap: use for prioritizing cell-state-specific eRegulons.
Dot heatmap: color encodes relative AUC; dot size encodes the fraction of active cells.
UMAP overlay: use for localization of selected eRegulons, not as the only statistical evidence.
Network PDF: use as a compact overview; use the source edge TSV for final custom network layouts.
```

For a main figure, redraw selected eRegulons from the source tables rather than showing every detected regulon. Keep direct and extended eRegulons separated unless the figure explicitly states how they were merged.
