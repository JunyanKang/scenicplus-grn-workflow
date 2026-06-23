# SCENIC+ 注释后 matched snRNA + snATAC 数据逐步分析流程

本指南适用于已经完成单细胞 QC、整合、聚类和注释的 matched snRNA + snATAC / 10x Multiome 风格数据。流程遵循 SCENIC+ 官方分析顺序：
annotated single-cell object → ATAC fragment/peak standardization → pycisTopic cisTopic model → custom cisTarget database → SCENIC+ Snakemake inference。

本流程从一个已经注释好的 scMultiome 对象，以及原始 ATAC fragment 和 peak 文件开始，最后进入官方 SCENIC+ Snakemake 工作流。

官方资源：

```text
SCENIC+ documentation: https://scenicplus.readthedocs.io/en/latest/
SCENIC+ running tutorial: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html
pycisTopic tutorial: https://pycistopic.readthedocs.io/en/latest/tutorials.html
SCENIC+ GitHub: https://github.com/aertslab/scenicplus
Aerts cisTarget resources: https://resources.aertslab.org/cistarget/
```

## 0. Initialize Environment And Project

这一步会验证已安装的 conda 环境，创建项目目录结构，更新 `$PROJECT_DIR/scenicplus_project.env` 中的基础项目参数，并写出 `$PROJECT_DIR/project_env.sh`。

在正式执行命令链前，应确认安装和运行时完整性检查通过，并确认当前分析属于官方 SCENIC+ 文档中的 matched scRNA+snATAC 路径。

参数说明：

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

支持的 `ORGANISM`：

```text
human      Homo sapiens GRCh38, chr1-chr22, chrX, chrY
mouse      Mus musculus GRCm39, chr1-chr19, chrX, chrY
cyno       Macaca fascicularis 6.0, chr1-chr20, chrX
rat        Rattus norvegicus GRCr8, chr1-chr20, chrX, chrY
rabbit     Oryctolagus cuniculus OryCun2.0, chr1-chr21, chrX
chicken    Gallus gallus GRCg7b, chr1-chr39, chrZ, chrW
zebrafish  Danio rerio GRCz11, chr1-chr25
```

支持的 `ATAC_INPUT_LAYOUT`：

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

首先在终端输入项目参数：

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

执行前确认：

```text
Confirm this is the intended analysis scope:
- An annotated scMultiome object is available as Seurat RDS/QS or AnnData h5ad.
- Original ATAC fragments and ATAC peak files are available for each sample.
- The intended analysis is matched scRNA+snATAC SCENIC+, not unmatched pseudo-integration.
- CONDA_ROOT, ENV_NAME, PROJECT_DIR, ANNOTATED_OBJECT and ATAC_DATA_ROOT are absolute paths.
- PROJECT_DIR is a dedicated SCENIC+ analysis root.
```

运行已安装环境和 workflow 完整性检查：

```bash
mkdir -p "$PROJECT_DIR/logs"
bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/check_environment.sh" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  2>&1 | tee "$PROJECT_DIR/logs/pre_step0_check_environment.log"
```

然后运行一步式初始化。它会验证参数值、检查 conda 环境、更新项目设置文件，并初始化项目运行文件：

```bash
bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/initialize_scenicplus_project.sh"
```

下游步骤开始前，请再次确认这里的每个参数，尤其是绝对路径、物种、ATAC layout 和对象路径。`scenicplus_project.env` 是后续项目参数的记录来源。

初始化成功后加载运行变量：

```bash
source "$PROJECT_DIR/project_env.sh"
```

## 1. Prepare Organism Resources

使用已安装的资源准备脚本，而不是手写 genome 下载、chromosome allowlist、UCSC 转换、annotation table、motif collection 和 motif2TF 命令。脚本支持断点续跑，会写日志、写出 `resources/resource_status.tsv`，并把所有文件路径和 checksum 记录到 `resources/resource_manifest.json`。

所有 genome resources 都会转换到 UCSC chromosome style。脚本只保留标准 primary chromosomes，并从 FASTA、GTF、chromsizes 和 SCENIC+ genome annotation table 中过滤 random、unplaced、alt、haplotype 和 mitochondrial records。后续 fragments、peaks、consensus peaks、region_sets BED 和 cisTarget DB regions 使用同一套标准染色体策略。

对于没有 Aerts public motif2TF table 的支持物种，脚本会下载 official human HGNC v10 motif2TF table，并通过 Ensembl BioMart one-to-one orthology 把 TF gene names 映射到目标物种。这不是从头 motif discovery。cached orthology table 和 mapping audit 会记录到 `resources/resource_manifest.json`。

准备 Step 0 中选定的项目物种：

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py
```

如果 UCSC-standard genome resources 已经存在，只缺 cisTarget motif collection 或 motif2TF table：

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py --motifs-only
```

`--motifs-only` 会准备 motif collection、species motif2TF table 和项目级 `inputs/cistarget_db/motif_annotations.tbl`，不会重建 UCSC FASTA、GTF、chromsizes 或 genome annotation files。

资源准备成功后重新加载 `project_env.sh`。资源准备脚本会更新 `MACS_GENOME_SIZE` 等资源派生变量。

```bash
source "$PROJECT_DIR/project_env.sh"
```

检查已有资源而不重建：

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py --mode check

python $SCENICPLUS_HOME/scripts/prepare_official_resources.py --mode status
```

只有在构建可复用共享资源缓存时，才准备所有支持物种：

```bash
python $SCENICPLUS_HOME/scripts/prepare_official_resources.py \
  --organism all \
  --mode prepare
```

## 2. Inspect And Export The Annotated Object

SCENIC+ 需要 active cells 的 RNA counts 和 metadata。从已经注释好的 scMultiome object 开始，例如 Seurat RDS/QS 或 AnnData h5ad。

所有命令先加载项目环境：

```bash
source "$PROJECT_DIR/project_env.sh"
```

annotated object 需要包含：

```text
RNA raw-count assay or layer inside the annotated object
sample column
condition column
cell-type or cell-state annotation column
barcodes that can be matched to the fragment files
```

先检查 annotated object。这个对象应该已经只包含本次 SCENIC+ 需要分析的细胞，并且 embedding 应该是在这些细胞上计算得到的。例如，如果只分析某个组织、谱系或发育 compartment，应提供 subset 后重新处理得到的对象，而不是 whole-atlas object。

统一 inspector 会根据后缀自动识别 `.rds`、`.qs` 或 `.h5ad`，并检测可用的 sample、condition、label、barcode、RNA assay、raw-count layer 和 metacell embedding 字段。

```bash
python $SCENICPLUS_HOME/scripts/inspect_annotated_object.py
```

导出前检查报告和参数表：

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

额外预览文件会写到 `results/annotated_object/`。

确认 `cell_label_column`、`assay`、`layer` 和 `reduction`。`cell_label_column` 是用户在这一步选择的分组字段；导出时这个源字段会统一写成下游固定列 `cell_label`。`assay` 和 `layer` 决定导出到 `inputs/gex.h5ad` 的 RNA count matrix；sample 和 condition 列决定生物学分组；barcode 字段用于连接 annotated object 和 ATAC fragments；`reduction` 用于 Step 3 metacell aggregation。

Seurat 输入会生成 `inputs/annotated_object_params.tsv`；h5ad 输入会生成 `inputs/annotated_h5ad_params.tsv`。

如果自动检测字段不对，在导出前编辑对应参数表。

如果 `CELL_LABEL_COLUMN` 不是目标分组字段，可以在导出时指定：

```bash
python $SCENICPLUS_HOME/scripts/export_annotated_object.py --cell-label-column corrected_metadata_column
```

导出 active RNA 和 metadata：

```bash
python $SCENICPLUS_HOME/scripts/export_annotated_object.py
```

导出会保留输入 active object 中的所有细胞。做 selected-population analysis 时，应直接把已筛选并重新处理的对象作为 `ANNOTATED_OBJECT`，不要在这一步从 whole-atlas object 临时 subset。

导出命令会把 `ANNOTATED_OBJECT`、`ANNOTATED_OBJECT_FORMAT`、`CELL_LABEL_COLUMN`、`ACTIVE_GEX_H5AD` 和 `ACTIVE_CELL_METADATA` 记录到 `$PROJECT_DIR/scenicplus_project.env`。

`h5Seurat`、`loom` 和其他格式应先转换成 Seurat RDS/QS 或 AnnData h5ad。

Step 2 写出：

```text
inputs/cell_metadata.tsv
inputs/gex.h5ad
results/annotated_object/cells_by_sample_and_label.tsv
results/annotated_object/annotated_object_summary.tsv
inputs/grn_label_summary.tsv
```

下游 pycisTopic 必须基于 active cells 和原始 fragments 重建。

运行 pycisTopic 前检查 `inputs/grn_label_summary.tsv`。

创建或更新后续步骤使用的参数表：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py
```

脚本会在 `inputs/` 下写默认项目参数表，并把路径记录到 `$PROJECT_DIR/scenicplus_project.env`。除非使用 `--force`，已有 scalar 参数值会保留。sample-specific tables 会在 `sample_id` 不再匹配 active sample sheet 时刷新。修改参数可以编辑对应 `inputs/*_params.tsv`，或使用命令行覆盖：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section pycistopic \
  --set pycistopic.n_iter=300
```

## 3. Metacell Aggregation

对于较大的 matched snRNA+snATAC 项目，在 pycisTopic 和 SCENIC+ 之前按 `sample_id × cell_label` 创建 metacells。这符合 metacell-first 策略，可降低内存压力，同时保留已注释细胞状态。

安装包中的 metacell 实现使用 hdWGCNA，需要 Seurat RDS/QS object。如果 Step 2 从 h5ad 开始，应在运行本步骤前，在 `inputs/metacell_params.tsv` 中提供匹配的 Seurat RDS/QS。

Step 3 会从 Step 2 初始化 `assay`、`layer` 和 `reduction`。matched multiome 分析中，应使用在 active analysis cells 上计算的 WNN UMAP。只有 ATAC modality 不可用或明确排除时才使用 RNA UMAP。

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section metacell
Rscript $SCENICPLUS_HOME/scripts/prepare_metacell_inputs_from_seurat.R
python $SCENICPLUS_HOME/scripts/make_metacell_gex_h5ad.py
```

本 workflow 中，metacell aggregation 是 pycisTopic 和 SCENIC+ 前的必经步骤。

默认 metacell 参数：

```text
inputs/metacell_params.tsv
```

修改参数可编辑该文件，或用 `--set`：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section metacell \
  --set metacell.k=12 \
  --set metacell.max_shared=2
```

前两个命令完成后，`inputs/cell_metadata.tsv` 和 `inputs/gex.h5ad` 会指向 metacells。原始 single-cell metadata 会保留为：

```text
inputs/cell_metadata.single_cell.tsv
inputs/metacell_membership.tsv
results/metacells/metacell_summary.tsv
```

## 4. Prepare Active ATAC Fragments And Peaks

根据 Step 0 选择的 ATAC 文件布局创建 `$PROJECT_DIR/inputs/sample_sheet.tsv`。

如果 Step 0 已设置 `ATAC_INPUT_LAYOUT` 和 `ATAC_DATA_ROOT`，先注册：

```bash
python $SCENICPLUS_HOME/scripts/set_atac_input_params.py
```

生成并验证 active ATAC sample sheet：

```bash
python $SCENICPLUS_HOME/scripts/make_sample_sheet_from_atac_inputs.py
python $SCENICPLUS_HOME/scripts/validate_and_prepare_sample_sheet.py
```

`make_sample_sheet_from_atac_inputs.py` 会根据上一条命令产生的 `inputs/atac_input_params.tsv` 自动分发到对应 layout parser。

ATAC 输入参数保存在 `inputs/atac_input_params.tsv`，并记录到 `$PROJECT_DIR/scenicplus_project.env`。生成的 sample sheet 是：

```text
inputs/sample_sheet.tsv
```

继续前检查 `condition` 列。

用一个命令标准化 ATAC peaks 和 fragments，为每个标准化 fragment 文件建立 tabix index，并更新 active sample sheet 指向标准化输入：

```bash
python $SCENICPLUS_HOME/scripts/standardize_atac_inputs.py
```

把标准化后的 single-cell fragments 重写为 Step 3 产生的 metacell barcodes：

```bash
python $SCENICPLUS_HOME/scripts/reassign_fragments_to_metacells.py
```

fragment reassignment 后再次验证 active sample sheet：

```bash
python $SCENICPLUS_HOME/scripts/validate_and_prepare_sample_sheet.py
```

## 5. Create `inputs/cistopic_obj.pkl` With pycisTopic

基于 UCSC-standard fragments 和项目特异 consensus peak universe 构建 ATAC object。本步骤从 `inputs/sample_sheet.tsv` 和 `inputs/cell_metadata.tsv` 开始。RNA 和 cisTopic 两侧的 cell-name convention 应为 `barcode-sample_id`。

`cistopic_obj.pkl`、consensus peaks、topic region sets 和 DAR region sets 都特异于 active cells/metacells 和 active `cell_label` grouping。如果分析范围从 whole atlas 变为 selected population，或从 cells 变为 metacells，必须从 selected barcodes 和原始 fragments 重新构建。直接 subset 旧 whole-atlas cisTopic object 只适合检查或绘制既有 atlas-level model，不适合作为新的 GRN 推断输入。

创建或更新 QC、peak-calling、DAR、doublet 和 topic-model 参数文件。`chromsizes` 和 `genome_size` 来自 `ORGANISM` 对应公共资源；sample-level QC rows 根据 active `inputs/sample_sheet.tsv` 生成。

pycisTopic pseudobulk 可能很慢且占内存。除非需要从头重建所有 pseudobulk BED，保持 `resume_pseudobulk=1`。重跑时，workflow 会复用非空且 gzip-valid 的 BED 文件，只重建缺失或损坏的 label 文件。

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section pycistopic
```

写出或更新：

```text
inputs/atac_qc_thresholds.tsv
inputs/topic_model_grid.tsv
inputs/pycistopic_params.tsv
```

修改 topic grid 或 pycisTopic 参数：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section pycistopic \
  --set topic_grid.n_topics=10,20,30 \
  --set pycistopic.n_iter=300
```

默认 topic-model backend 是 `cgs`。对于较大的 peak-by-cell/metacell matrix，建议显式使用 pycisTopic MALLET backend：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py \
  --section pycistopic \
  --set pycistopic.lda_backend=mallet \
  --set pycistopic.mallet_path="$CONDA_ENV_PREFIX/bin/mallet" \
  --set pycistopic.mallet_memory_gb=auto
```

使用 MALLET 时，同一轮 topic models 应全部来自同一 backend。若从 CGS 切到 MALLET，应 archive 或移除中断的 CGS `Topic*.pkl` 后再重跑。`mallet_memory_gb=auto` 会检测可用内存并通过 `MALLET_MEMORY` 设置 Java heap。运行时信息记录在 `results/pycistopic/model_selection/mallet_runtime.tsv`。

`topic_n_cpu=auto` 会独立于 pseudobulk 和 MACS2 workers 解析。MALLET topic modeling CPU-bound 且内存稳定时，可显式设置如 `3`。

由于 Step 3 产生 metacell inputs，默认 `analysis_unit` 为 `metacell`。只有做 cell-level run 时才设置 `pycistopic.analysis_unit=cell`。

运行：

```bash
python $SCENICPLUS_HOME/scripts/run_pycistopic_workflow.py 2>&1 | tee "$PROJECT_DIR/logs/run_pycistopic_workflow.log"
```

查看自动 worker plan 和 pseudobulk resume plan：

```bash
cat results/pycistopic/qc/parallelism_plan.tsv
cat results/pycistopic/qc/pseudobulk_resume_plan.tsv
```

如果 pseudobulk 曾中断，继续前验证已有 pseudobulk 文件。只删除报告为 invalid 的文件；valid 文件会复用。

```bash
python $SCENICPLUS_HOME/scripts/validate_pseudobulk_files.py
```

SCENIC+ 所需的 pycisTopic 完成标志：

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

检查 cisTopic object 和 active RNA matrix 的 cell IDs 是否匹配：

```bash
python $SCENICPLUS_HOME/scripts/check_cistopic_cell_names.py
```

## 6. Create `inputs/region_sets/`

Step 5 产生的必要结构：

```text
inputs/region_sets/
  Topics_otsu/Topic1.bed
  Topics_top_3k/Topic1.bed
  DARs_cell_label/CellLabelA_VS_rest.bed
```

验证并标准化所有 region-set BED 文件：

```bash
python $SCENICPLUS_HOME/scripts/standardize_region_sets.py
```

脚本写出：

```text
results/pycistopic/qc/region_set_standardization.tsv
```

## 7. Build The Custom cisTarget Region Database

custom cisTarget database 必须基于 pycisTopic 和 region sets 使用的同一个 UCSC-standard consensus region universe。先创建或更新 cisTarget database 参数，再构建数据库。这一步是最终 SCENIC+ inference 的必需输入：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section cistarget
python $SCENICPLUS_HOME/scripts/build_custom_cistarget_db.py
```

参数保存于：

```text
inputs/cistarget_db_params.tsv
```

对于大型 region universe，除非调度系统要求固定线程数，否则保持 `n_cpu=auto`。启动时，`build_custom_cistarget_db.py` 会检查 CPU count、当前负载和可用内存，为 `create_cistarget_motif_databases.py -t` 解析一个保守 worker count。长时间 motif-scanning 开始前会写出：

```text
results/cistarget_db/custom_cistarget_resource_plan.tsv
```

对于大型 region-by-motif matrix，`use_partial=auto` 会自动启用官方 partial
cisTarget workflow。wrapper 会运行
`create_cistarget_motif_databases.py --partial`，重跑时跳过已经存在且非空的
part 文件；所有 part 完成后，用官方 combine 脚本合并 partial score
databases，再用官方 conversion 脚本把 motif-vs-region scores 转换为最终的
regions-vs-motifs rankings database。resource plan 会记录 `use_partial`、
`partial_n_parts` 和估算的 full float32 score matrix 大小。

预期输出：

```text
inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather
inputs/cistarget_db/custom.regions_vs_motifs.scores.feather
results/cistarget_db/custom_cistarget_db_manifest.tsv
```

## 8. Initialize SCENIC+ And Generate Config

SCENIC+ 应在所有 active metacells 上一起运行，以建立一个共享 eRegulon universe。如果项目包含多个 condition，应在这一步合并运行；condition-specific effects 后续在共享 eRegulon AUC matrix 上测试。

```bash
python $SCENICPLUS_HOME/scripts/initialize_scenicplus_snakemake.py
```

这会初始化 SCENIC+ Snakemake 目录，复制物种特异 annotation 文件，写出 `work/scenicplus/organism_config.yaml`，创建或更新 `inputs/scenicplus_config_params.tsv`，并生成最终 Snakemake config：

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

需要检查的重要参数：

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

如需修改这些值，编辑 `$PROJECT_DIR/inputs/scenicplus_config_params.tsv`，或用 `--set` 重跑 setup script，然后重新生成 config。若需修改路径或高级 SCENIC+ 选项，编辑生成文件：

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

cell-name overlap、UCSC chromosome consistency 和核心 SCENIC+ config invariants 会在 Step 9 检查。Step 9 失败时不要运行 Step 10。

## 9. Preflight Checks Before SCENIC+

Snakemake dry run 前运行所有检查。这些检查应明确失败，但阈值可通过项目参数调整。

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section preflight
python $SCENICPLUS_HOME/scripts/preflight_scenicplus_inputs.py 2>&1 | tee "$PROJECT_DIR/logs/preflight_scenicplus_inputs.log"
```

阈值记录在 `inputs/preflight_thresholds.tsv`，数据较 noisy 时可调整阈值，但不要改变输出文件名。

记录软件版本和资源 checksum：

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

运行参数保存于：

```text
inputs/snakemake_params.tsv
```

预期输出：

```text
results/scenicplus/eRegulons_direct.tsv
results/scenicplus/eRegulons_extended.tsv
results/scenicplus/AUCell_direct.h5mu
results/scenicplus/AUCell_extended.h5mu
results/scenicplus/scplus_mdata.h5mu
```

正式分析中，建议至少运行两次独立 SCENIC+，或一次标准执行加一次 deterministic patched execution，然后比较 high-confidence edges。开始第二次前复制已完成结果目录：

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

运行：

```bash
python $SCENICPLUS_HOME/scripts/compare_scenicplus_stability.py \
  --run-a results/scenicplus_stability/run1 \
  --run-b results/scenicplus_stability/run2 \
  --out results/scenicplus_stability/stability_summary.tsv
```

全局 SCENIC+ seed 控制部分 workflow，但 region-to-gene GBM 仍可能受实现细节和依赖版本影响。最终 GRN 结果应报告 stability summary。

## 11. Condition-Level eRegulon AUC Statistics

在 Step 10 共享 eRegulon universe 上，用 sample-level AUC 做 differential eRegulon activity testing。不要把分别推断的 condition-specific eRegulon lists 作为主要比较。这一脚本读取 SCENIC+ AUCell `.h5mu`，连接 `inputs/cell_metadata.tsv`，按 sample 平均 AUC，并在 sample-level means 上测试 condition differences。它总会写统计表和 PDF diagnostic page。

创建或更新 postprocessing 参数，并对 direct 和 extended eRegulons 运行：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section postprocess
python $SCENICPLUS_HOME/scripts/run_scenicplus_postprocess.py --task stats --layer all
```

必要输出：

```text
results/scenicplus_stats/auc_by_condition_direct/sample_mean_auc.tsv
results/scenicplus_stats/auc_by_condition_direct/condition_eregulon_auc_statistics.tsv
results/scenicplus_stats/auc_by_condition_direct/condition_eregulon_auc_statistics.pdf
results/scenicplus_stats/auc_by_condition_extended/sample_mean_auc.tsv
results/scenicplus_stats/auc_by_condition_extended/condition_eregulon_auc_statistics.tsv
results/scenicplus_stats/auc_by_condition_extended/condition_eregulon_auc_statistics.pdf
```

只有在有足够独立样本时才报告 `delta_mean_auc` 和 FDR。如果某个 condition 少于两个样本，保留 sample-mean table 和 PDF，但不要把 per-cell P values 当作生物学重复解释。

## 12. SCENIC+ Result Figures

完整 SCENIC+ 分析不应只停在表格，但图形输出应按科学证据层级组织，而不是按绘图函数是否容易调用来组织。先检查每一层证据是否已经具备：

```bash
python $SCENICPLUS_HOME/scripts/run_scenicplus_postprocess.py --task audit
```

审计输出：

```text
results/scenicplus_output_tiers/scenicplus_output_tier_audit.tsv
results/scenicplus_output_tiers/scenicplus_output_tier_audit.pdf
```

科学输出层级：

```text
0. input/QC confidence: active RNA/ATAC cells, metadata, fragments, doublet and ATAC QC
1. chromatin topics: major accessibility states and topic model quality
2. region sets and DARs: topic-region programs and differential accessible regions
3. motif/cisTarget evidence: project-specific motif enrichment database and motif2TF mapping
4. eRegulon activity: TF-region-gene eRegulons and cell-state activity patterns
5. condition effects: sample-level differential eRegulon activity across conditions
6. mechanism views: focused TF-target networks and locus/coverage views when the required inputs exist
```

完成审计后，为 direct 和 extended eRegulon layer 各生成标准 eRegulon 图和 source data。安装器提供 `plot_scenicplus_publication_outputs.py`，会写出 vector PDFs 和对应 source tables：

```text
1. eRegulon AUC heatmap across cell states or labels
2. RSS-style eRegulon specificity heatmap
3. heatmap-dotplot: color = group AUC z-score, size = active-cell fraction
4. UMAP overlays for high-priority eRegulons
5. compact TF-target network view
```

视觉风格为白底、Arial/Helvetica 字体、低饱和色、小而清晰的标签、无重网格，并全部导出为 PDF。每张图都有对应 source table。

生成 direct 和 extended eRegulon figures：

```bash
python $SCENICPLUS_HOME/scripts/setup_workflow_params.py --section postprocess
python $SCENICPLUS_HOME/scripts/run_scenicplus_postprocess.py --task figures --layer all
```

Postprocessing 参数：

```text
inputs/postprocess_params.tsv
```

每个 layer 的必要输出：

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

解释规则：

```text
AUC heatmap: use for group-level eRegulon activity patterns.
Specificity heatmap: use for prioritizing cell-state-specific eRegulons.
Dot heatmap: color encodes relative AUC; dot size encodes the fraction of active cells.
UMAP overlay: use for localization of selected eRegulons, not as the only statistical evidence.
Network PDF: use as a compact overview; use the source edge TSV for final custom network layouts.
```

正式主图建议从 source tables 中重绘少数关键 eRegulons，而不是展示全部 detected regulons。direct 和 extended eRegulons 应分开呈现，除非图中明确说明如何合并。
