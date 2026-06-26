# SCENIC+ 注释后 matched snRNA + snATAC 数据逐步分析流程

本指南适用于已经完成单细胞 QC、整合、聚类和注释的 matched snRNA + snATAC / 10x Multiome 风格数据。流程遵循 SCENIC+ 官方分析顺序：
annotated single-cell object → ATAC fragment/peak standardization → pycisTopic cisTopic model → custom cisTarget database → SCENIC+ Snakemake inference。

本流程从一个已经注释好的 scMultiome 对象，以及原始 ATAC fragment 和 peak 文件开始，最后进入官方 SCENIC+ Snakemake 工作流。

分析逻辑：

```text
1. 先固定分析对象和细胞标签，避免后续 GRN 混入错误细胞范围。
2. 再用同一批 active cells 同步准备 RNA counts、metadata、ATAC fragments 和 peaks。
3. pycisTopic 用 ATAC co-accessibility 描述染色质状态，输出 topics、DARs 和 region sets。
4. custom cisTarget database 把项目自己的 peak universe 转成 motif ranking/score 背景。
5. SCENIC+ 在同一套细胞和 region universe 上连接 TF、enhancer-region 和 target gene。
6. 最后用 shared eRegulon universe 比较细胞状态和 condition effects。
```

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

生物信息学目的：先固定物种、分析对象、细胞标签和 ATAC 输入布局。SCENIC+ 的 TF-region-gene 连接高度依赖这些基础定义；如果这里混入非目标细胞或错误 genome build，后续 eRegulon 不能解释为目标组织/状态的调控网络。

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

0.1-在终端输入项目参数：

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

0.2-运行已安装环境和 workflow 完整性检查：

```bash
mkdir -p "$PROJECT_DIR/logs"
"$CONDA_ROOT/envs/$ENV_NAME/bin/spgrn-check" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  2>&1 | tee "$PROJECT_DIR/logs/pre_step0_check_environment.log"
```

0.3-运行一步式初始化。它会验证参数值、检查 conda 环境、更新项目设置文件，并初始化项目运行文件：

```bash
"$CONDA_ROOT/envs/$ENV_NAME/bin/spgrn-initialize"
```

下游步骤开始前，请再次确认这里的每个参数，尤其是绝对路径、物种、ATAC layout 和对象路径。`scenicplus_project.env` 是后续项目参数的记录来源。

0.4-加载项目运行变量：

```bash
source "$PROJECT_DIR/project_env.sh"
```

## 1. Prepare Organism Resources

分析目的：准备统一的 UCSC 染色体格式、gene annotation、chromsizes、motif collection 和 motif2TF 表。这里定义的是全流程的坐标系统和 motif-to-TF 解释框架；后续 peak、region set、cisTarget database 和 SCENIC+ search space 都必须使用同一套坐标和注释。

使用已安装的资源准备脚本，而不是手写 genome 下载、chromosome allowlist、UCSC 转换、annotation table、motif collection 和 motif2TF 命令。脚本支持断点续跑，会写日志、写出 `resources/resource_status.tsv`，并把所有文件路径和 checksum 记录到 `resources/resource_manifest.json`。

所有 genome resources 都会转换到 UCSC chromosome style。脚本只保留标准 primary chromosomes，并从 FASTA、GTF、chromsizes 和 SCENIC+ genome annotation table 中过滤 random、unplaced、alt、haplotype 和 mitochondrial records。后续 fragments、peaks、consensus peaks、region_sets BED 和 cisTarget DB regions 使用同一套标准染色体策略。

对于没有 Aerts public motif2TF table 的支持物种，脚本会下载 official human HGNC v10 motif2TF table，并通过 Ensembl BioMart one-to-one orthology 把 TF gene names 映射到目标物种。这不是从头 motif discovery。cached orthology table 和 mapping audit 会记录到 `resources/resource_manifest.json`。

1.1-准备 Step 0 中选定的项目物种：

```bash
spgrn-prepare-official-resources
```

如果 UCSC-standard genome resources 已经存在，只缺 cisTarget motif collection 或 motif2TF table：

```bash
spgrn-prepare-official-resources --motifs-only
```

`--motifs-only` 会准备 motif collection、species motif2TF table 和项目级 `inputs/cistarget_db/motif_annotations.tbl`，不会重建 UCSC FASTA、GTF、chromsizes 或 genome annotation files。

资源准备成功后，重新加载 `project_env.sh`。资源准备脚本会更新 `MACS_GENOME_SIZE` 等资源派生变量。

1.2-重新加载项目运行变量：

```bash
source "$PROJECT_DIR/project_env.sh"
```

检查已有资源而不重建：

```bash
spgrn-prepare-official-resources --mode check

spgrn-prepare-official-resources --mode status
```

只有在构建可复用共享资源缓存时，才准备所有支持物种：

```bash
spgrn-prepare-official-resources \
  --organism all \
  --mode prepare
```

## 2. Inspect And Export The Annotated Object

分析目的：从已经注释好的 scMultiome 对象导出 active cells 的 RNA count matrix 和 metadata。这里不重新做细胞注释；它的作用是确认对象中哪些 metadata column、assay/layer 和 embedding 会被后续 SCENIC+ 使用。

解释重点：`CELL_LABEL_COLUMN` 应代表后续 GRN 要比较的主要细胞状态或细胞类型。这个标签会进入 metacell aggregation、DAR calling、region-set 命名和 eRegulon AUC 汇总。

SCENIC+ 需要 active cells 的 RNA counts 和 metadata。从已经注释好的 scMultiome object 开始，例如 Seurat RDS/QS 或 AnnData h5ad。

2.1-加载项目运行变量：

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

2.2-检查 annotated object：

```bash
spgrn-inspect-annotated-object
```

2.3-检查导出前报告和参数表：

```bash
spgrn-review-annotated-object-inspection
```

该脚本会写出 `results/annotated_object/annotated_object_pre_export_review.md` 和 `.tsv`，如果必要导出字段缺失会直接失败。

确认 `cell_label_column`、`assay`、`layer` 和 `reduction`。`cell_label_column` 是用户在这一步选择的分组字段；导出时这个源字段会统一写成下游固定列 `cell_label`。`assay` 和 `layer` 决定导出到 `inputs/gex.h5ad` 的 RNA count matrix；sample 和 condition 列决定生物学分组；barcode 字段用于连接 annotated object 和 ATAC fragments；`reduction` 用于 Step 3 metacell aggregation。

Seurat 输入会生成 `inputs/annotated_object_params.tsv`；h5ad 输入会生成 `inputs/annotated_h5ad_params.tsv`。

如果自动检测字段不对，在导出前编辑对应参数表。

如果 `CELL_LABEL_COLUMN` 不是目标分组字段，可以在导出时指定：

```bash
spgrn-export-annotated-object --cell-label-column corrected_metadata_column
```

2.4-导出 active RNA 和 metadata：

```bash
spgrn-export-annotated-object
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

2.5-创建或更新后续步骤使用的参数表：

```bash
spgrn-setup-workflow-params
```

脚本会在 `inputs/` 下写默认项目参数表，并把路径记录到 `$PROJECT_DIR/scenicplus_project.env`。除非使用 `--force`，已有 scalar 参数值会保留。sample-specific tables 会在 `sample_id` 不再匹配 active sample sheet 时刷新。修改参数可以编辑对应 `inputs/*_params.tsv`，或使用命令行覆盖：

```bash
spgrn-setup-workflow-params \
  --section pycistopic \
  --set pycistopic.n_iter=300
```

## 3. Metacell Aggregation

分析目的：把相近细胞聚合成 metacells，提高 RNA/ATAC 信号稳定性，并降低 SCENIC+ 推断的计算量。metacell 不是新的生物学细胞类型；它是保留原始细胞结构的统计单位。

解释重点：metacells 应在同一 cell-label 框架内生成，避免把不同发育状态或不同谱系强行平均。每个 metacell 的 sample、condition 和 cell label 必须可追溯。

对于较大的 matched snRNA+snATAC 项目，在 pycisTopic 和 SCENIC+ 之前按 `sample_id × cell_label` 创建 metacells。这符合 metacell-first 策略，可降低内存压力，同时保留已注释细胞状态。

安装包中的 metacell 实现使用 hdWGCNA，需要 Seurat RDS/QS object。如果 Step 2 从 h5ad 开始，应在运行本步骤前，在 `inputs/metacell_params.tsv` 中提供匹配的 Seurat RDS/QS。

Step 3 会从 Step 2 初始化 `assay`、`layer` 和 `reduction`。matched multiome 分析中，应使用在 active analysis cells 上计算的 WNN UMAP。只有 ATAC modality 不可用或明确排除时才使用 RNA UMAP。

3.1-创建或更新 metacell 参数：

```bash
spgrn-setup-workflow-params --section metacell
```

3.2-从 Seurat 对象生成 metacell 输入：

```bash
spgrn-prepare-metacell-inputs-from-seurat
```

3.3-生成 metacell `gex.h5ad`：

```bash
spgrn-make-metacell-gex-h5ad
```

本 workflow 中，metacell aggregation 是 pycisTopic 和 SCENIC+ 前的必经步骤。

默认 metacell 参数：

```text
inputs/metacell_params.tsv
```

修改参数可编辑该文件，或用 `--set`：

```bash
spgrn-setup-workflow-params \
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

分析目的：把原始 ATAC fragments 限定到 active cells/metacells，并统一 peak/fragment 坐标格式。SCENIC+ 的 enhancer evidence 来自这些 ATAC regions；RNA 对象中有细胞不等于 ATAC fragment 中也有可用 reads，因此这里必须做 barcode 和 sample 对齐。

解释重点：fragments 决定实际可用的 accessibility signal，peaks 决定初始候选 region universe。后续 consensus peaks、DARs 和 cisTarget database 都从这里延伸。

根据 Step 0 选择的 ATAC 文件布局创建 `$PROJECT_DIR/inputs/sample_sheet.tsv`。

4.1-注册 ATAC 输入参数：

```bash
spgrn-set-atac-input-params
```

4.2-生成 active ATAC sample sheet：

```bash
spgrn-make-sample-sheet-from-atac-inputs
```

4.3-验证 active ATAC sample sheet：

```bash
spgrn-validate-and-prepare-sample-sheet
```

`make_sample_sheet_from_atac_inputs.py` 会根据上一条命令产生的 `inputs/atac_input_params.tsv` 自动分发到对应 layout parser。

ATAC 输入参数保存在 `inputs/atac_input_params.tsv`，并记录到 `$PROJECT_DIR/scenicplus_project.env`。生成的 sample sheet 是：

```text
inputs/sample_sheet.tsv
```

继续前检查 `condition` 列。

4.4-标准化 ATAC peaks 和 fragments：

```bash
spgrn-standardize-atac-inputs
```

4.5-把 single-cell fragments 重写为 metacell barcodes：

```bash
spgrn-reassign-fragments-to-metacells
```

4.6-再次验证 active sample sheet：

```bash
spgrn-validate-and-prepare-sample-sheet
```

## 5. Create `inputs/cistopic_obj.pkl` With pycisTopic

分析目的：用 pycisTopic 从 ATAC peak-by-cell/metacell matrix 中学习 co-accessible region topics。topics 可理解为染色质可及性程序，通常对应细胞状态、发育轴或调控模块。

解释重点：topic model 不是最终 GRN。它提供 region sets 和 accessibility programs，供后续 motif enrichment 和 SCENIC+ 连接 TF-region-gene。

基于 UCSC-standard fragments 和项目特异 consensus peak universe 构建 ATAC object。本步骤从 `inputs/sample_sheet.tsv` 和 `inputs/cell_metadata.tsv` 开始。RNA 和 cisTopic 两侧的 cell-name convention 应为 `barcode-sample_id`。

`cistopic_obj.pkl`、consensus peaks、topic region sets 和 DAR region sets 都特异于 active cells/metacells 和 active `cell_label` grouping。如果分析范围从 whole atlas 变为 selected population，或从 cells 变为 metacells，必须从 selected barcodes 和原始 fragments 重新构建。直接 subset 旧 whole-atlas cisTopic object 只适合检查或绘制既有 atlas-level model，不适合作为新的 GRN 推断输入。

创建或更新 QC、peak-calling、DAR 和 topic-model 参数文件。`chromsizes` 和 `genome_size` 来自 `ORGANISM` 对应公共资源；sample-level QC rows 根据 active `inputs/sample_sheet.tsv` 生成。

pycisTopic pseudobulk 可能很慢且占内存。除非需要从头重建所有 pseudobulk BED，保持 `resume_pseudobulk=1`。重跑时，workflow 会复用非空且 gzip-valid 的 BED 文件，只重建缺失或损坏的 label 文件。

5.1-创建或更新 pycisTopic 参数：

```bash
spgrn-setup-workflow-params --section pycistopic
```

写出或更新：

```text
inputs/atac_qc_thresholds.tsv
inputs/topic_model_grid.tsv
inputs/pycistopic_params.tsv
```

修改 topic grid 或 pycisTopic 参数：

```bash
spgrn-setup-workflow-params \
  --section pycistopic \
  --set topic_grid.n_topics=10,20,30 \
  --set pycistopic.n_iter=300
```

默认 topic-model backend 是 MALLET，更适合较大的 peak-by-cell/metacell matrix。workflow 会把实际 MALLET runtime 写入 `results/pycistopic/model_selection/mallet_runtime.tsv`。

如果在 CGS 和 MALLET 之间切换，同一轮 topic models 应全部来自同一 backend。切换前应 archive 或移除中断的 `Topic*.pkl` 后再重跑。

由于 Step 3 产生 metacell inputs，默认 `analysis_unit` 为 `metacell`。只有做 cell-level run 时才设置 `pycistopic.analysis_unit=cell`。

5.2-运行 pycisTopic workflow：

```bash
spgrn-run-pycistopic-workflow 2>&1 | tee "$PROJECT_DIR/logs/run_pycistopic_workflow.log"
```

如果 pseudobulk 曾中断，继续前验证已有 pseudobulk 文件。只删除报告为 invalid 的文件；valid 文件会复用。

```bash
spgrn-validate-pseudobulk-files
```

5.3-检查 pycisTopic 是否已经生成 SCENIC+ 所需的全部输入：

```bash
spgrn-check-pycistopic-completion
```

该脚本会写出 `results/pycistopic/qc/pycistopic_completion_check.md` 和 `.tsv`，如果任何必需输出缺失或为空会直接失败。

5.4-检查 cisTopic object 和 active RNA matrix 的 cell IDs 是否匹配：

```bash
spgrn-check-cistopic-cell-names
```

## 6. Create `inputs/region_sets/`

分析目的：生成进入 motif enrichment 的 region sets，包括 topic-derived regions 和 cell-label/DAR-derived regions。topic sets 捕捉共可及性程序；DAR sets 捕捉细胞状态间 differential accessibility。

解释重点：SCENIC+ 的 motif enrichment 不应只依赖一种 region set。topic、DAR 和 top-region sets 互相补充，可以区分广泛 accessibility state 与更具体的 cell-state regulatory elements。

Step 5 产生的必要结构：

```text
inputs/region_sets/
  Topics_otsu/Topic1.bed
  Topics_top_3k/Topic1.bed
  DARs_cell_label/CellLabelA_VS_rest.bed
```

6.1-验证并标准化所有 region-set BED 文件：

```bash
spgrn-standardize-region-sets
```

脚本写出：

```text
results/pycistopic/qc/region_set_standardization.tsv
```

## 7. Build The Custom cisTarget Region Database

分析目的：用当前项目的 consensus peak/region universe 构建 motif ranking 和 score database。这样 motif enrichment 的背景与本项目真实可检测的 ATAC regions 一致，而不是套用不匹配的公共 region universe。

解释重点：rankings database 支持 cisTarget enrichment，scores database 支持 DEM enrichment。两者服务于不同 motif enrichment 统计，但都用于后续 eRegulon 构建。

custom cisTarget database 必须基于 pycisTopic 和 region sets 使用的同一个 UCSC-standard consensus region universe。这一步是最终 SCENIC+ inference 的必需输入。

7.1-创建或更新 cisTarget database 参数：

```bash
spgrn-setup-workflow-params --section cistarget
```

7.2-构建 custom cisTarget database：

```bash
spgrn-build-custom-cistarget-db
```

参数保存于：

```text
inputs/cistarget_db_params.tsv
```

脚本会自动管理资源使用和断点续跑。长时间 motif-scanning 的运行记录写在：

```text
results/cistarget_db/custom_cistarget_resource_plan.tsv
```

大型 region-by-motif matrix 会使用官方 partial/combine/convert 路径生成最终 rankings 和 scores database。

主要输出：

```text
inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather
inputs/cistarget_db/custom.regions_vs_motifs.scores.feather
results/cistarget_db/custom_cistarget_db_manifest.tsv
```

## 8. Initialize SCENIC+ And Generate Config

SCENIC+ 应在所有 active metacells 上一起运行，以建立一个共享 eRegulon universe。如果项目包含多个 condition，应在这一步合并运行；condition-specific effects 后续在共享 eRegulon AUC matrix 上测试。

分析目的：把前面准备好的 RNA、ATAC、region sets、motif database 和 genome resources 写入 SCENIC+ Snakemake 配置。这里决定后续 inference 使用同一个 search space 和 eRegulon universe。

8.1-初始化 SCENIC+ Snakemake 配置：

```bash
spgrn-initialize-scenicplus-snakemake
```

这会初始化 SCENIC+ Snakemake 目录，复制物种特异 annotation 文件，写出 `work/scenicplus/organism_config.yaml`，创建或更新 `inputs/scenicplus_config_params.tsv`，并生成最终 Snakemake config：

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

通常需要按研究设计检查的参数：

```text
seed                      可重复性随机种子。
search_space_upstream     gene upstream enhancer search-space window。
search_space_downstream   gene downstream enhancer search-space window。
search_space_extend_tss   TSS extension window。
dem_motif_hit_thr         DEM motif hit threshold。
ctx_nes_threshold         cisTarget NES threshold。
rho_threshold             region-gene correlation threshold。
min_target_genes          每个 eRegulon 至少保留的 target genes 数。
```

大型 custom cisTarget database 可在正式 Snakemake run 前按 region-set family 分块运行 motif enrichment：

```bash
spgrn-run-scenicplus-motif-enrichment-split --mode both
spgrn-run-scenicplus-motif-enrichment-split --mode status
```

这个命令会创建 `inputs/region_sets_split/`，把 DEM 和 cisTarget 拆成多个独立 SCENIC+ 进程运行，然后用所有 HDF5 结果调用 `prepare_menr`。分块清单写在：
`--mode status` 是进入 Step 9 前必须通过的完成检查：它会验证 chunk HDF5 文件签名，确认正式阈值下为空的 DEM chunk 已有诊断报告，并检查 `prepare_menr` 已生成有效的 `cistromes_direct.h5ad`、`cistromes_extended.h5ad` 和非空 `tf_names.txt`。

```text
work/scenicplus/motif_enrichment_split/motif_enrichment_split_resource_plan.tsv
work/scenicplus/motif_enrichment_split/motif_enrichment_split_chunks.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.md
```

命令支持断点续跑：已完成且非空的 chunk HDF5 会自动跳过，除非使用 `--force`。

如果某个 DEM chunk 在正式阈值下完成但没生成 HDF5，runner 会写入 `.empty.tsv` marker，并自动触发该 chunk 的 relaxed-threshold diagnostic；是否可以继续由 `--mode status` 判定。

诊断报告输出到：

```text
work/scenicplus/motif_enrichment_split/（chunk 级 relaxed DEM 输出与日志）
results/scenicplus_diagnostics/*_relaxed_threshold_diagnostic.tsv/md
```

如果全部 DEM chunk 都为空（即正式阈值下没有 DEM HDF5），流程会中止，并在错误信息中提示查看诊断报告。若 `--mode status` 因缺失、无效或半写入文件失败，不要进入 Step 9；重新运行 split 命令，只有明确要替换旧半成品时才使用 `--force`。

如需修改这些值，编辑 `$PROJECT_DIR/inputs/scenicplus_config_params.tsv`，或用 `--set` 重跑 setup script，然后重新生成 config。若需修改路径或高级 SCENIC+ 选项，编辑生成文件：

```text
$PROJECT_DIR/work/scenicplus/Snakemake/config/config.yaml
```

cell-name overlap、UCSC chromosome consistency 和核心 SCENIC+ config invariants 会在 Step 9 检查。Step 9 失败时不要运行 Step 10。

## 9. Preflight Checks Before SCENIC+

Snakemake dry run 前运行所有检查。这些检查应明确失败，但阈值可通过项目参数调整。

分析目的：在正式长任务前确认 RNA cells、ATAC cells、fragment barcodes、region names、motif annotations 和 genome coordinates 彼此一致。这里失败通常说明输入定义不一致，而不是 SCENIC+ 生物学结果不好。

9.1-创建或更新 preflight 参数：

```bash
spgrn-setup-workflow-params --section preflight
```

9.2-运行 preflight checks：

```bash
spgrn-preflight-scenicplus-inputs 2>&1 | tee "$PROJECT_DIR/logs/preflight_scenicplus_inputs.log"
```

阈值记录在 `inputs/preflight_thresholds.tsv`，数据较 noisy 时可调整阈值，但不要改变输出文件名。

9.3-记录软件版本和资源 checksum：

```bash
spgrn-record-scenicplus-provenance
```

## 10. Dry Run, Run, And Stability Record

分析目的：运行 SCENIC+ inference，建立 TF-to-gene、region-to-gene 和 TF-region-gene eRegulons。dry run 只检查 DAG 和文件依赖；正式 run 才生成网络结果。

10.1-创建或更新 Snakemake 参数：

```bash
spgrn-setup-workflow-params --section snakemake
```

10.2-运行 SCENIC+ dry run：

```bash
spgrn-run-scenicplus-snakemake --mode dryrun
```

10.3-运行正式 SCENIC+ inference：

```bash
spgrn-run-scenicplus-snakemake --mode run
```

这个外层命令会由 Snakemake 自动展开 SCENIC+ inference chain。日志中常见的内部规则包括：

```text
prepare_GEX_ACC_multiome   合并 matched RNA 和 ATAC matrix
get_search_space           建立 gene 周围候选 enhancer search space
motif_enrichment_dem       用 score database 做 DEM motif enrichment
motif_enrichment_cistarget 用 ranking database 做 cisTarget enrichment
prepare_menr               整合 DEM/cisTarget motif enrichment 结果
TF_to_gene                 用 TF expression 预测 target gene expression
region_to_gene             用 region accessibility 连接 candidate target genes
eGRN_direct/extended       合并 TF-gene、region-gene 和 motif/cistrome evidence
AUCell_direct/extended     计算每个细胞/metacell 的 eRegulon activity
scplus_mudata              汇总最终 SCENIC+ MuData 输出
```

关键参数文件：`inputs/snakemake_params.tsv`。

主要输出：

```text
results/scenicplus/eRegulons_direct.tsv
results/scenicplus/eRegulons_extended.tsv
results/scenicplus/AUCell_direct.h5mu
results/scenicplus/AUCell_extended.h5mu
results/scenicplus/scplus_mdata.h5mu
```

正式分析至少保留两次独立 inference 记录，用于评估 high-confidence edges 的稳定性。

10.4-归档第一次完成结果：

```bash
mkdir -p results/scenicplus_stability/run1
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run1/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run1/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run1/
```

10.5-运行第二次独立 inference dry run：

```bash
spgrn-run-scenicplus-snakemake --mode dryrun --rerun-inference
```

10.6-运行第二次独立 inference：

```bash
spgrn-run-scenicplus-snakemake --mode run --rerun-inference
```

`--rerun-inference` 只重跑 SCENIC+ inference chain，不重建 genome resources。

10.7-第二次完成后归档 run2：

```bash
mkdir -p results/scenicplus_stability/run2
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run2/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run2/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run2/
```

10.8-比较两次 inference 的 edge 稳定性：

```bash
spgrn-compare-scenicplus-stability \
  --run-a results/scenicplus_stability/run1 \
  --run-b results/scenicplus_stability/run2 \
  --out results/scenicplus_stability/stability_summary.tsv
```

最终 GRN 结果同时报告 `results/scenicplus_stability/stability_summary.tsv`。

## 11. SCENIC+ Postprocess Figures And Condition Statistics

Step 10 结束后，把 postprocess 作为一个有顺序的整体运行。这一步把共享的 SCENIC+ eRegulon universe 转化为 source tables、矢量 PDF、输出层级审计和 sample-level condition statistics。网络推断和 condition testing 要分开理解：SCENIC+ 在所有 active metacells 上共同推断一次；condition 比较是在同一套 eRegulon AUCell activity 上按生物学样本进行。

11.1-创建或更新 postprocessing 参数：

```bash
spgrn-setup-workflow-params --section postprocess
```

11.2-生成输出层级审计、全部标准 PDF、source tables 和 condition statistics：

```bash
spgrn-run-scenicplus-postprocess --task all --layer all
```

Step 11 保持为一个 postprocess 阶段，因为 audit、figure source tables、PDF 和 condition statistics 都是在解释同一次完成的 SCENIC+ inference 结果。

Postprocess 命令参数：

```text
--task audit    只运行输出层级审计。
--task figures  重新生成 source tables 和 01-08 eRegulon PDFs。
--task stats    重新生成 09+ condition-statistics tables 和 PDFs。
--task all      生成完整 postprocess output set。

--layer direct    只处理 direct eRegulons。
--layer extended  只处理 extended eRegulons。
--layer all       同时处理 direct 和 extended eRegulons。
```

`--task` 控制分析模块，`--layer` 控制 eRegulon 层。当前 wrapper 不按单个 PDF 文件名运行；如果只想重画某类图，使用对应的 task/layer。正式论文主图建议从 source table 中单独重绘聚焦 panel。

主要输出位置：

```text
results/scenicplus_output_tiers/       输出层级审计 TSV/PDF，用来判断哪些结果家族可以解释。
results/scenicplus_figures/01-05_*     eRegulon activity、specificity、dot heatmap 和 embedding PDFs。
results/scenicplus_figures/06-08_*     Region-gene、target-region overlap 和 TF-target network PDFs。
results/scenicplus_figures/09-15_*     Sample-level condition statistics PDFs。
results/scenicplus_figures/source_*    供正式主图重绘使用的 figure source tables。
results/scenicplus_figures/*stats*.tsv condition-statistics source tables，包含实际 contrast 方向。
logs/                                  每个 extraction、R render 和 condition-statistics 步骤的日志。
```

R 绘图参数文件：

```text
results/scenicplus_figures/plot_style_parameters.tsv
results/scenicplus_figures/plot_style_parameters.md
```

`plot_style_parameters.tsv` 是可编辑文件，只包含 `parameter` 和 `value` 两列，不写入项目特异性标签。Markdown 文件会自动生成，记录参数含义和一键重跑命令。修改颜色、字号、线宽、点大小、透明度或每行 panel 数后，重新运行 11.2 即可重画全部 PDF。

Condition volcano 的优先标注可以通过 `inputs/postprocess_params.tsv` 的 `priority_eregulons` 指定。该 TSV 至少包含 `eregulon` 或 `display_label`，如包含 `cell_label` 和 `layer`，绘图时会自动按对应层级匹配。

Condition statistics 会自动判断设计：单 condition 只输出描述性 sample means；两组 condition 输出 `comparison - reference` 的 `delta_mean_auc`，并在 `contrast` 列记录方向；多 condition 输出 sample-level omnibus test，不强行给单一 delta。只有自动识别的两组 reference 不符合实验设计时，才需要在 `inputs/postprocess_params.tsv` 修改 `reference_condition` 或 `comparison_condition`。

科学输出层级：

```text
0. Output-tier audit and input confidence：确认哪些 SCENIC+ 输出存在，以及哪些输出层级可以解释。
1. eRegulon activity and specificity：展示细胞状态层面的调控程序活性和特异性。
2. Condition-resolved eRegulon activity：在同一套 eRegulon universe 上按 cell label 和 condition 拆分。
3. Region-gene and target-region structure：总结 enhancer-target link 和 target-region overlap 的模型结构。
4. TF-target network overview：展示筛选后的 TF-target link 概览。
5. Sample-level condition statistics：基于生物学样本 mean AUCell score 计算 condition effects。
```

解释规则：

```text
AUC heatmap：用于观察 group-level eRegulon activity；行和列按绘图矩阵数值聚类。
Condition AUC heatmap：用于比较每个 cell label 内的 condition shift；列按 cell_label_condition 排列，行按数值聚类。
Specificity heatmap：用于优先筛选 cell-state-specific eRegulons；行和列按绘图矩阵数值聚类。
Dot heatmap：颜色表示相对 AUC，点大小表示 active cell fraction。
UMAP/activity embedding：用于定位和状态分离展示，不能单独作为统计证据。
Region-gene and overlap PDFs：属于模型结构图；condition effects 应从 AUCell statistics 解释。
Network PDF：作为 TF-target 概览；正式聚焦网络建议用 source edge table 重绘。
Condition volcano and condition heatmap：优先解释 sample-level effect size；只有独立生物学样本足够时才报告 FDR。
eRegulon signs：标签保留 SCENIC+ 符号，例如 +/+、-/+、-/-。不要把这些符号折叠成 TF RNA expression；eRegulon AUC 反映的是 target gene-set activity。
```

正式主图建议从 source tables 中重绘少数关键 eRegulons，而不是展示全部 detected regulons。direct 和 extended eRegulons 应分开呈现，除非图中明确说明如何合并。
