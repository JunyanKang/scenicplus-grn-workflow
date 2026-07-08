# SCENIC+ Matched snRNA+snATAC Step-by-Step Workflow

This guide is for annotated matched snRNA+snATAC / scMultiome data. The starting
point is an active annotated object and matched ATAC fragments/peaks from the
same samples. The endpoint is a SCENIC+ eRegulon analysis with AUCell activity,
condition statistics, source tables and vector PDF outputs.

Core rules:

- Infer one shared eRegulon universe from all active cells or metacells.
- Compare conditions after inference on the shared eRegulon AUC matrix.
- For selected tissue, lineage or cell-state analyses, provide a selected and
  reprocessed active annotated object.
- RNA counts, metadata, ATAC fragments, region sets and custom cisTarget
  databases must refer to the same active cells/metacells and the same
  UCSC-style genome coordinate system.

Official resources:

```text
SCENIC+ documentation: https://scenicplus.readthedocs.io/en/latest/
SCENIC+ running tutorial: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html
pycisTopic documentation: https://pycistopic.readthedocs.io/en/latest/
SCENIC+ GitHub: https://github.com/aertslab/scenicplus
Aerts cisTarget resources: https://resources.aertslab.org/cistarget/
```

## 0. Initialize Environment And Project

Purpose: define the conda environment and project directory. These values are
recorded in `$PROJECT_DIR/scenicplus_project.env` and
`$PROJECT_DIR/project_env.sh`; later steps append organism, annotated object,
cell-label and ATAC input settings when they become relevant.

Bioinformatics rationale: first fix the runtime and output directory. Each
later step writes its own parameters only when the relevant biological input is
known, reducing the chance of carrying a wrong object or ATAC directory through
the full workflow.

0.1-Enter project parameters:

```bash
# Replace every value in this block before running.
# CONDA_ROOT is the real conda/miniforge/miniconda/mambaforge/anaconda root.
export CONDA_ROOT=/absolute/path/to/conda
# ENV_NAME is the SCENIC+ environment created by the workflow installation.
export ENV_NAME=scenicplus-grn
# Add installed workflow commands to this terminal PATH.
export PATH="$CONDA_ROOT/envs/$ENV_NAME/bin:$PATH"
# PROJECT_DIR is the dedicated SCENIC+ analysis root; inputs/work/logs/results are written there.
export PROJECT_DIR=/absolute/path/to/grn_project/scenicplus_analysis
# AUTOZYME is on or off.
export AUTOZYME=on
```

0.2-Check the installed environment:

```bash
spgrn-check
```

0.3-Initialize the project:

```bash
spgrn-initialize
```

0.4-Activate the environment and load project variables:

```bash
source "$CONDA_ROOT/bin/activate" "$ENV_NAME"
source "$PROJECT_DIR/project_env.sh"
```

## 1. Prepare Organism Resources

Purpose: confirm that the study organism has both FASTA and GTF in the
selected Ensembl release, then prepare the motif collection and motif2TF table
needed by SCENIC+. FASTA/GTF come from Ensembl; the remaining SCENIC+ standard
resources are derived and audited by the workflow.

1.0-Select organism and Ensembl release, then write them to the project config:

```bash
# ENSEMBL_RELEASE defaults to 115; use --release first to query another release and write the available-species TSV.
export ENSEMBL_RELEASE=115
spgrn-query-organism-resources --list --release "$ENSEMBL_RELEASE"

# ORGANISM must be edited after checking the query output. Ensembl species names are preferred.
export ORGANISM=mus_musculus

# Query the selected organism.
spgrn-query-organism-resources --organism "$ORGANISM" --release "$ENSEMBL_RELEASE"

# Write ORGANISM and ENSEMBL_RELEASE to scenicplus_project.env and project_env.sh.
spgrn-initialize
source "$PROJECT_DIR/project_env.sh"
```

`--list` writes `ensembl_release_<release>_organism_resources.tsv` with only Ensembl species that have both FASTA and GTF and reports the
motif2TF preparation strategy. Use `--organism` to inspect the assembly,
FASTA/GTF availability, motif collection and motif2TF strategy for one species.

1.1-Prepare resources. Choose one motif2TF route:

- direct: human, mouse, fly and chicken have Aerts v10 direct motif2TF tables.

```bash
spgrn-prepare-official-resources
```

- mapping: other species can be explicitly mapped with `--ref human|mouse|fly|chicken` through Ensembl BioMart orthology and audit files.

```bash
spgrn-prepare-official-resources --ref human
```

For organisms with extensive paralog expansion, use the paralog-aware policy
and review the mapping audit:

```bash
spgrn-prepare-official-resources --ref human --orthology-policy paralog-aware
```

- generated: for organisms without a direct table or safe default mapping, build an audited symbol-evidence table from Aerts direct motif evidence and target gene symbols.

```bash
spgrn-prepare-official-resources --generate-motif2tf
```

Generated tables strictly match target annotation gene symbols against Aerts
direct human/mouse/fly/chicken motif2TF evidence. Review
`resources/<organism>/*_motif2tf_generated_symbol_audit.tsv` before treating
the table as biological evidence.

- user table: if a manually reviewed species-specific motif2TF table already exists, provide that table.

```bash
# MOTIF2TF_TABLE must be edited to the real motif_annotations.tbl path.
export MOTIF2TF_TABLE=/absolute/path/to/motif_annotations.tbl
spgrn-prepare-official-resources
```

1.2-Reload resource-derived variables:

```bash
source "$PROJECT_DIR/project_env.sh"
```

1.3-Check resource status:

```bash
spgrn-prepare-official-resources --mode status
```

If only motif resources are missing:

```bash
spgrn-prepare-official-resources --motifs-only
```

Main outputs:

```text
resources/resource_manifest.json
resources/resource_status.tsv
resources/motif2tf/motif_annotations.<organism>.<strategy>[.<reference>].tbl
resources/motif2tf/motif_annotations.<organism>.active.tsv
inputs/cistarget_db/motif_annotations.tbl
resources/<organism>/*_chromosome_audit.tsv
resources/<organism>/*_motif2tf_*_audit.tsv
```

The strategy-named table under `resources/motif2tf/` is the canonical copy.
`inputs/cistarget_db/motif_annotations.tbl` remains the standard SCENIC+ input
file for the current active species.

## 2. Inspect And Export The Active Annotated Object

Purpose: export RNA counts and metadata for the active cells/metacells. This
step does not annotate cells; it confirms the assay/layer, sample, condition,
cell label, barcode and embedding fields used downstream.

Bioinformatics rationale: `CELL_LABEL_COLUMN` controls metacell aggregation,
DAR calling, region-set naming and eRegulon AUC summaries. It should represent
the main cell type, state or developmental stage used for the GRN analysis.

2.1-Set and inspect the annotated object:

```bash
# ANNOTATED_OBJECT is the active annotated Seurat RDS/QS or AnnData h5ad object.
export ANNOTATED_OBJECT=/absolute/path/to/active_annotated_multiome_object.rds
spgrn-initialize
source "$PROJECT_DIR/project_env.sh"
spgrn-inspect-annotated-object
```

2.2-Review the pre-export report:

```bash
spgrn-review-annotated-object-inspection
```

If automatic field detection is wrong, edit:

```text
inputs/annotated_object_params.tsv
inputs/annotated_h5ad_params.tsv
```

2.3-Set cell label and export active RNA and metadata:

```bash
# CELL_LABEL_COLUMN is the metadata column used as the SCENIC+ cell label.
export CELL_LABEL_COLUMN=cell_annotation
spgrn-initialize
source "$PROJECT_DIR/project_env.sh"
spgrn-export-annotated-object
```

To override the selected cell-label column:

```bash
# corrected_metadata_column must be replaced with a real metadata column name.
spgrn-export-annotated-object --cell-label-column corrected_metadata_column
```

2.4-Create workflow parameter tables:

```bash
spgrn-setup-workflow-params
```

Main outputs:

```text
inputs/gex.h5ad
inputs/cell_metadata.tsv
inputs/grn_label_summary.tsv
results/annotated_object/annotated_object_summary.tsv
results/annotated_object/annotated_object_pre_export_review.md
```

Review `inputs/grn_label_summary.tsv` before continuing.

## 3. Build Metacells

Purpose: aggregate similar cells within the `sample_id × cell_label` framework
to stabilize RNA/ATAC signal and reduce SCENIC+ compute cost.

Bioinformatics rationale: a metacell is a statistical unit, not a new cell
type. Build metacells within comparable sample and label groups to avoid
averaging distinct lineages or states.

This workflow treats metacells as required. The installed implementation uses
hdWGCNA and needs a Seurat RDS/QS object. If Step 2 started from h5ad, provide
a matched Seurat object in `inputs/metacell_params.tsv`.

3.1-Create or update metacell parameters:

```bash
spgrn-setup-workflow-params --section metacell
```

3.2-Generate metacell membership and metadata:

```bash
spgrn-prepare-metacell-inputs-from-seurat
```

3.3-Generate metacell RNA h5ad:

```bash
spgrn-make-metacell-gex-h5ad
```

Main outputs:

```text
inputs/metacell_params.tsv
inputs/metacell_membership.tsv
inputs/cell_metadata.single_cell.tsv
inputs/cell_metadata.tsv
inputs/gex.h5ad
results/metacells/metacell_summary.tsv
```

For matched multiome data, use a WNN UMAP computed on the active cells when
available. Use RNA UMAP only if ATAC is unavailable or intentionally excluded.

## 4. Prepare Active ATAC Fragments And Peaks

Purpose: create the active sample sheet, standardize peaks/fragments and
rewrite single-cell fragments to metacell barcodes.

Bioinformatics rationale: cells exported from the RNA object must have matched
ATAC reads in the fragment files. This step aligns sample IDs, barcodes,
fragment indexes and peak coordinates.

Supported ATAC input layouts:

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

4.1-Register ATAC input parameters:

```bash
# ATAC_INPUT_LAYOUT is one supported ATAC layout and must match ATAC_DATA_ROOT.
export ATAC_INPUT_LAYOUT=split_ge_arc
# ATAC_DATA_ROOT is the root directory containing ATAC fragments and peaks.
export ATAC_DATA_ROOT=/absolute/path/to/atac_input_root
spgrn-set-atac-input-params
```

4.2-Generate the sample sheet:

```bash
spgrn-make-sample-sheet-from-atac-inputs
```

4.3-Validate the sample sheet:

```bash
spgrn-validate-and-prepare-sample-sheet
```

4.4-Standardize ATAC inputs:

```bash
spgrn-standardize-atac-inputs
```

4.5-Reassign fragments to metacell barcodes:

```bash
spgrn-reassign-fragments-to-metacells
```

4.6-Validate the sample sheet again:

```bash
spgrn-validate-and-prepare-sample-sheet
```

Main outputs:

```text
inputs/atac_input_params.tsv
inputs/sample_sheet.tsv
work/standard_peaks/
work/metacell_fragments/
results/metacells/metacell_fragment_reassignment.tsv
```

Check sample, condition, fragment path and peak path columns in
`inputs/sample_sheet.tsv` before continuing. Standardized peak paths and
metacell fragment paths are written back into `inputs/sample_sheet.tsv`.

## 5. Run The pycisTopic Workflow

Purpose: build the cisTopic object from ATAC fragments and consensus peaks,
learn co-accessible region topics, and export topics, DARs and region sets.

Bioinformatics rationale: pycisTopic describes chromatin accessibility
programs. It is not the final GRN; it supplies topic-derived and DAR-derived
region sets for motif enrichment and TF-region-gene linking.

`cistopic_obj.pkl`, consensus peaks, topic region sets and DAR region sets are
specific to the active cells/metacells and active `cell_label`. Rebuild them
from original fragments if the analysis scope changes.

5.1-Create or update pycisTopic parameters:

```bash
spgrn-setup-workflow-params --section pycistopic
```

5.2-Run pycisTopic:

```bash
spgrn-run-pycistopic-workflow
```

The command writes its own timestamped log under `$PROJECT_DIR/logs/`.

5.3-If a previous run was interrupted, validate pseudobulk files:

```bash
spgrn-validate-pseudobulk-files
```

5.4-Check pycisTopic completion:

```bash
spgrn-check-pycistopic-completion
```

5.5-Check RNA and cisTopic cell IDs:

```bash
spgrn-check-cistopic-cell-names
```

Main outputs:

```text
inputs/cistopic_obj.pkl
inputs/region_sets/
results/pycistopic/qc/pycistopic_completion_check.md
results/pycistopic/model_selection/
```

The default topic-model backend is MALLET. If switching between CGS and MALLET,
keep all topic models in one run from the same backend.

## 6. Standardize Region Sets

Purpose: ensure topic regions and DAR regions use one UCSC-style BED format
before motif enrichment.

Bioinformatics rationale: region sets define motif-enrichment foregrounds.
Wrong chromosome style or malformed BED records break the background model.

6.1-Standardize region sets:

```bash
spgrn-standardize-region-sets
```

Expected structure:

```text
inputs/region_sets/
  Topics_otsu/
  Topics_top_3k/
  DARs_cell_label/
```

Main output:

```text
results/pycistopic/qc/region_set_standardization.tsv
```

## 7. Build The Custom cisTarget Database

Purpose: build motif ranking and score databases from the project consensus
region universe.

Bioinformatics rationale: public cisTarget databases may not match the regions
detectable in this project. A custom database keeps motif enrichment on the
same ATAC peak universe used for pycisTopic. Rankings support cisTarget
enrichment; scores support DEM enrichment.

7.1-Create or update cisTarget parameters:

```bash
spgrn-setup-workflow-params --section cistarget
```

7.2-Build the custom database:

```bash
spgrn-build-custom-cistarget-db
```

Main outputs:

```text
inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather
inputs/cistarget_db/custom.regions_vs_motifs.scores.feather
results/cistarget_db/custom_cistarget_resource_plan.tsv
results/cistarget_db/custom_cistarget_db_manifest.tsv
```

Large databases use the official partial/combine/convert path and can resume
from completed parts. Heartbeat logs report long motif-scanning progress.

## 8. Initialize SCENIC+ Snakemake

Purpose: write RNA, ATAC, region sets, custom cisTarget databases and genome
annotation into the SCENIC+ Snakemake config.

Bioinformatics rationale: SCENIC+ should infer one shared eRegulon universe on
all active metacells. Multiple conditions stay together here; condition effects
are tested in Step 11.

8.1-Initialize the config:

```bash
spgrn-initialize-scenicplus-snakemake
```

Main outputs:

```text
work/scenicplus/Snakemake/config/config.yaml
work/scenicplus/organism_config.yaml
inputs/scenicplus_config_params.tsv
```

Parameters usually worth reviewing:

```text
seed
search_space_upstream
search_space_downstream
search_space_extend_tss
dem_motif_hit_thr
ctx_nes_threshold
rho_threshold
min_target_genes
```

8.2-Optional split motif enrichment for very large custom databases:

```bash
spgrn-run-scenicplus-motif-enrichment-split --mode both
spgrn-run-scenicplus-motif-enrichment-split --mode status
```

`--mode both` runs DEM and cisTarget motif enrichment by region-set family and
then runs `prepare_menr`. `--mode status` is the completion gate before Step 9:
it checks chunk HDF5 signatures, DEM empty diagnostics, `cistromes_direct.h5ad`,
`cistromes_extended.h5ad` and `tf_names.txt`.

Main outputs:

```text
work/scenicplus/motif_enrichment_split/motif_enrichment_split_resource_plan.tsv
work/scenicplus/motif_enrichment_split/motif_enrichment_split_chunks.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.md
```

If `--mode status` fails, do not continue to Step 9. Rerun the split command;
use `--force` only when intentionally replacing old partial outputs.

## 9. SCENIC+ Preflight Checks

Purpose: check cell-name overlap, UCSC chromosome consistency, motif
annotations, region names and core config invariants before the Snakemake dry
run.

Bioinformatics rationale: preflight failures usually indicate inconsistent
input definitions, not poor biology. Running Snakemake after a failed preflight
wastes time and produces hard-to-interpret output.

9.1-Create or update preflight parameters:

```bash
spgrn-setup-workflow-params --section preflight
```

9.2-Run preflight:

```bash
spgrn-preflight-scenicplus-inputs
```

9.3-Record software versions and resource checksums:

```bash
spgrn-record-scenicplus-provenance
```

Main outputs:

```text
inputs/preflight_thresholds.tsv
logs/preflight_scenicplus_inputs_*.log
results/provenance/
```

## 10. Run SCENIC+ Inference And Record Stability

Purpose: run SCENIC+ inference to build TF-to-gene, region-to-gene and
TF-region-gene eRegulons, then compute AUCell activity.

Bioinformatics rationale: dry run checks the DAG and file dependencies. The
formal run creates the GRN. For formal analyses, keep two independent
inference records to evaluate high-confidence edge stability.

10.1-Create or update Snakemake parameters:

```bash
spgrn-setup-workflow-params --section snakemake
```

10.2-Run dry run:

```bash
spgrn-run-scenicplus-snakemake --mode dryrun
```

10.3-Run inference:

```bash
spgrn-run-scenicplus-snakemake --mode run
```

Main SCENIC+ internal rules:

```text
prepare_GEX_ACC_multiome
get_search_space
motif_enrichment_dem
motif_enrichment_cistarget
prepare_menr
TF_to_gene
region_to_gene
eGRN_direct / eGRN_extended
AUCell_direct / AUCell_extended
scplus_mudata
```

10.4-Archive run1:

```bash
mkdir -p results/scenicplus_stability/run1
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run1/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run1/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run1/
```

10.5-Run the second independent inference dry run:

```bash
spgrn-run-scenicplus-snakemake --mode dryrun --rerun-inference
```

10.6-Run the second independent inference:

```bash
spgrn-run-scenicplus-snakemake --mode run --rerun-inference
```

10.7-Archive run2:

```bash
mkdir -p results/scenicplus_stability/run2
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run2/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run2/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run2/
```

10.8-Compare stability:

```bash
spgrn-compare-scenicplus-stability \
  --run-a results/scenicplus_stability/run1 \
  --run-b results/scenicplus_stability/run2 \
  --out results/scenicplus_stability/stability_summary.tsv
```

Main outputs:

```text
results/scenicplus/eRegulons_direct.tsv
results/scenicplus/eRegulons_extended.tsv
results/scenicplus/AUCell_direct.h5mu
results/scenicplus/AUCell_extended.h5mu
results/scenicplus/scplus_mdata.h5mu
results/scenicplus_stability/stability_summary.tsv
```

## 11. Postprocess, Figures And Condition Statistics

Purpose: convert completed SCENIC+ inference into an output-tier audit, source
tables, vector PDFs and sample-level condition statistics.

Bioinformatics rationale: network inference and condition testing are separate.
SCENIC+ infers eRegulons once on all active metacells; condition statistics are
computed later from the shared AUCell activity matrix across biological samples.

11.1-Create or update postprocess parameters:

```bash
spgrn-setup-workflow-params --section postprocess
```

11.2-Generate complete postprocess outputs:

```bash
spgrn-run-scenicplus-postprocess --task all --layer all
```

Options:

```text
--task audit      Output-tier audit only.
--task figures    Source tables and 01-08 eRegulon PDFs.
--task stats      09+ condition-statistics tables and PDFs.
--task all        Complete postprocess output set.

--layer direct    Direct eRegulons.
--layer extended  Extended eRegulons.
--layer all       Both direct and extended layers.
```

`regulon_sign_filter=tf_positive` in `inputs/postprocess_params.tsv` is the
default manuscript-facing setting and keeps only `+/+` eRegulons. Set it to
`all` only when exploratory signed regulon output is needed.

Main outputs:

```text
results/scenicplus_output_tiers/
results/scenicplus_figures/01-08_*.pdf
results/scenicplus_figures/09-15_*.pdf
results/scenicplus_figures/source_*.tsv
results/scenicplus_figures/*stats*.tsv
results/scenicplus_figures/plot_style_parameters.tsv
results/scenicplus_figures/plot_style_parameters.md
```

Condition statistics are automatic:

```text
Single condition: descriptive sample means only.
Two conditions: delta_mean_auc = comparison - reference, with direction recorded in contrast.
More than two conditions: sample-level omnibus test, no forced single delta.
```

Set `reference_condition` or `comparison_condition` in
`inputs/postprocess_params.tsv` only when the automatic two-condition reference
is not the intended baseline.

Interpretation rules:

```text
AUC heatmap: group-level eRegulon activity; rows and columns are clustered.
Condition heatmap: condition shifts within cell labels; columns are cell_label_condition.
Specificity heatmap: prioritize cell-state-specific eRegulons.
Dot heatmap: color is relative AUC; dot size is active fraction.
Embedding: localizes eRegulon activity but is not standalone statistical evidence.
Region-gene / overlap: model structure, not direct condition effect.
Network: TF-target overview; redraw focused networks from source tables.
Volcano: prioritize sample-level effect size; interpret FDR carefully with few samples.
eRegulon signs: `+/+`, `-/+` and `-/-` are SCENIC+ regulon signs, not TF RNA
expression. Main figures and condition statistics keep TF-positive `+/+`
eRegulons by default. To explore all signed regulons, set
`regulon_sign_filter` to `all` in `inputs/postprocess_params.tsv`.
```

For manuscripts, redraw a focused set of biologically relevant eRegulons from
the source tables. Keep direct and extended eRegulons separate unless the figure
legend explicitly states the merge rule.
