# SCENIC+ matched snRNA+snATAC 逐步分析手册

本手册面向已经完成基础单细胞分析的 matched snRNA+snATAC / scMultiome
数据。起点是一个已经注释好的 active scMultiome 对象，以及同一批样本的
ATAC fragments 和 peaks。终点是 SCENIC+ eRegulon、AUCell activity、
condition statistics、source tables 和可重画的 PDF 图。

核心原则：

- SCENIC+ 应在同一批 active cells 或 metacells 上共同推断一个 shared eRegulon universe。
- condition 或处理组差异不在建库阶段分开跑，而是在同一套 eRegulon AUC 上按生物学样本比较。
- 如果只分析某个组织、谱系或细胞群，应先提供已经 subset 并重新处理的 active annotated object。
- RNA counts、metadata、ATAC fragments、region sets 和 custom cisTarget database 必须指向同一批 active cells/metacells 和同一套 UCSC-style genome coordinates。

官方资源：

```text
SCENIC+ documentation: https://scenicplus.readthedocs.io/en/latest/
SCENIC+ running tutorial: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html
pycisTopic documentation: https://pycistopic.readthedocs.io/en/latest/
SCENIC+ GitHub: https://github.com/aertslab/scenicplus
Aerts cisTarget resources: https://resources.aertslab.org/cistarget/
```

## 0. 初始化环境和项目

目的：固定本次 GRN 分析的环境、物种、active object、细胞标签列和 ATAC 输入目录。这里的定义会写入 `$PROJECT_DIR/scenicplus_project.env` 和 `$PROJECT_DIR/project_env.sh`，后续命令都从这里读取。

生物信息学逻辑：SCENIC+ 的 TF-region-gene 连接依赖 cell identity、genome build、fragment barcode 和 peak universe。Step 0 如果选错对象或物种，后面即使命令成功，结果也不能解释为目标细胞群的调控网络。

支持物种：

```text
human      Homo sapiens GRCh38, chr1-chr22, chrX, chrY
mouse      Mus musculus GRCm39, chr1-chr19, chrX, chrY
cyno       Macaca fascicularis 6.0, chr1-chr20, chrX
rat        Rattus norvegicus GRCr8, chr1-chr20, chrX, chrY
rabbit     Oryctolagus cuniculus OryCun2.0, chr1-chr21, chrX
chicken    Gallus gallus GRCg7b, chr1-chr39, chrZ, chrW
zebrafish  Danio rerio GRCz11, chr1-chr25
```

支持的 ATAC 输入布局：

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

0.1-输入项目参数：

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

参数含义：

```text
CONDA_ROOT          conda/miniforge/miniconda/mambaforge/anaconda 根目录。
ENV_NAME            安装器创建的 SCENIC+ 环境名。
PROJECT_DIR         本次 SCENIC+ 分析根目录，所有 inputs/work/logs/results 都写在这里。
ORGANISM            上方支持物种之一。
AUTOZYME            on 或 off；控制 AutoZyme runtime 是否启用。
ENSEMBL_RELEASE     准备 genome resources 使用的 Ensembl release。
ANNOTATED_OBJECT    已注释且已经限定到本次 active cells 的 Seurat RDS/QS 或 AnnData h5ad。
CELL_LABEL_COLUMN   annotated object 中作为 GRN cell label 的 metadata column。
ATAC_INPUT_LAYOUT   上方支持的 ATAC 文件布局之一。
ATAC_DATA_ROOT      ATAC fragments 和 peaks 所在根目录。
```

0.2-运行已安装环境检查：

```bash
mkdir -p "$PROJECT_DIR/logs"
"$CONDA_ROOT/envs/$ENV_NAME/bin/spgrn-check" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  2>&1 | tee "$PROJECT_DIR/logs/pre_step0_check_environment.log"
```

0.3-初始化项目：

```bash
"$CONDA_ROOT/envs/$ENV_NAME/bin/spgrn-initialize"
```

0.4-激活环境并加载项目变量：

```bash
source "$CONDA_ROOT/bin/activate" "$ENV_NAME"
source "$PROJECT_DIR/project_env.sh"
```

继续前确认：

```text
scenicplus_project.env 和 project_env.sh 已生成。
ANNOTATED_OBJECT 是本次真正要分析的 active object。
ATAC_DATA_ROOT 指向同一批样本的 fragments 和 peaks。
PROJECT_DIR 是专用分析目录，不是更大的父目录。
```

## 1. 准备物种公共资源

目的：准备 FASTA、GTF、chromsizes、SCENIC+ genome annotation、motif collection 和 motif2TF table。

生物信息学逻辑：SCENIC+ 会把 peaks、region sets、motif hits、enhancer-gene search space 和 gene annotations 放在同一个坐标系统里解释。这里统一使用 UCSC chromosome style，并只保留 standard primary chromosomes；random、unplaced、alt、haplotype 和 mitochondrial records 会被剔除。

对于没有 Aerts public motif2TF table 的支持物种，脚本会使用 human HGNC v10 motif2TF table，通过 Ensembl BioMart one-to-one orthology 映射到目标物种，并输出 mapping audit。这个过程不是从头 motif discovery。

1.1-准备 Step 0 选定物种：

```bash
spgrn-prepare-official-resources
```

1.2-重新加载资源派生变量：

```bash
source "$PROJECT_DIR/project_env.sh"
```

1.3-检查资源状态：

```bash
spgrn-prepare-official-resources --mode status
```

只补 motif collection 或 motif2TF table 时使用：

```bash
spgrn-prepare-official-resources --motifs-only
```

主要输出：

```text
resources/resource_manifest.json
resources/resource_status.tsv
inputs/cistarget_db/motif_annotations.tbl
```

## 2. 检查并导出 active annotated object

目的：从已注释 scMultiome 对象中导出 active cells/metacells 的 RNA counts 和 metadata。这里不重新注释细胞；它只确认下游 SCENIC+ 使用哪些 assay/layer、sample、condition、cell label、barcode 和 embedding。

生物信息学逻辑：`CELL_LABEL_COLUMN` 决定 metacell aggregation、DAR calling、region-set 命名和 eRegulon AUC 汇总的生物学分组。它应该代表本次 GRN 想比较的主要细胞类型、状态或发育阶段。

2.1-检查 annotated object：

```bash
spgrn-inspect-annotated-object
```

2.2-审查导出前报告：

```bash
spgrn-review-annotated-object-inspection
```

如果自动识别字段不对，编辑生成的参数表：

```text
inputs/annotated_object_params.tsv
inputs/annotated_h5ad_params.tsv
```

2.3-导出 active RNA 和 metadata：

```bash
spgrn-export-annotated-object
```

如果需要覆盖 Step 0 中的 cell label column：

```bash
spgrn-export-annotated-object --cell-label-column corrected_metadata_column
```

2.4-创建后续步骤参数表：

```bash
spgrn-setup-workflow-params
```

主要输出：

```text
inputs/gex.h5ad
inputs/cell_metadata.tsv
inputs/grn_label_summary.tsv
results/annotated_object/annotated_object_summary.tsv
results/annotated_object/annotated_object_pre_export_review.md
```

继续前检查 `inputs/grn_label_summary.tsv`。如果 cell label、sample 或 condition 明显不对，应回到 annotated object 或参数表修正。

## 3. 构建 metacells

目的：在 `sample_id × cell_label` 框架内把相近细胞聚合为 metacells，提高 RNA/ATAC 信号稳定性并降低 SCENIC+ 计算量。

生物信息学逻辑：metacell 是统计单位，不是新的细胞类型。它应该在相同 sample 和 cell label 内构建，避免把不同谱系或状态平均在一起。

本 workflow 默认 metacell 是必经步骤。安装包中的实现依赖 hdWGCNA 和 Seurat RDS/QS。如果 Step 2 从 h5ad 开始，需要在 `inputs/metacell_params.tsv` 中提供匹配的 Seurat 对象。

3.1-创建或更新 metacell 参数：

```bash
spgrn-setup-workflow-params --section metacell
```

3.2-生成 metacell membership 和 metacell metadata：

```bash
spgrn-prepare-metacell-inputs-from-seurat
```

3.3-生成 metacell RNA h5ad：

```bash
spgrn-make-metacell-gex-h5ad
```

主要输出：

```text
inputs/metacell_params.tsv
inputs/metacell_membership.tsv
inputs/cell_metadata.single_cell.tsv
inputs/cell_metadata.tsv
inputs/gex.h5ad
results/metacells/metacell_summary.tsv
```

对于 matched multiome，优先使用在 active cells 上计算的 WNN UMAP 作为 metacell 坐标。只有 ATAC modality 不可用或明确排除时才使用 RNA UMAP。

## 4. 准备 active ATAC fragments 和 peaks

目的：根据 active samples 和 active cells/metacells 生成 sample sheet，标准化 ATAC peaks/fragments，并把 single-cell fragments 重新分配到 metacell barcodes。

生物信息学逻辑：RNA 对象中的细胞必须能在 fragments 中找到可用 ATAC reads。这里确认 sample、barcode、fragment index 和 peak coordinates 一致。

4.1-注册 ATAC 输入参数：

```bash
spgrn-set-atac-input-params
```

4.2-生成 sample sheet：

```bash
spgrn-make-sample-sheet-from-atac-inputs
```

4.3-验证 sample sheet：

```bash
spgrn-validate-and-prepare-sample-sheet
```

4.4-标准化 ATAC 输入：

```bash
spgrn-standardize-atac-inputs
```

4.5-将 fragments 重写为 metacell barcodes：

```bash
spgrn-reassign-fragments-to-metacells
```

4.6-再次验证 sample sheet：

```bash
spgrn-validate-and-prepare-sample-sheet
```

主要输出：

```text
inputs/atac_input_params.tsv
inputs/sample_sheet.tsv
inputs/fragments_standardized/
inputs/peaks_standardized/
```

继续前检查 `inputs/sample_sheet.tsv` 中的 sample、condition、fragment 路径和 peak 路径。

## 5. 运行 pycisTopic workflow

目的：从 ATAC fragments 和 consensus peaks 构建 cisTopic object，学习 co-accessible region topics，并输出 topics、DARs 和 region sets。

生物信息学逻辑：pycisTopic 描述染色质可及性程序，不是最终 GRN。它为后续 motif enrichment 提供 topic-derived region sets，并为不同 cell label 的 DAR-derived region sets 提供基础。

`cistopic_obj.pkl`、consensus peaks、topic region sets 和 DAR region sets 都特异于 active cells/metacells 和 active `cell_label`。如果分析范围变化，应从原始 fragments 重新构建，不要直接 subset 旧的 whole-atlas cisTopic object 作为新 GRN 输入。

5.1-创建或更新 pycisTopic 参数：

```bash
spgrn-setup-workflow-params --section pycistopic
```

5.2-运行 pycisTopic：

```bash
spgrn-run-pycistopic-workflow 2>&1 | tee "$PROJECT_DIR/logs/run_pycistopic_workflow.log"
```

5.3-如曾中断，验证 pseudobulk 文件：

```bash
spgrn-validate-pseudobulk-files
```

5.4-检查 pycisTopic 完成度：

```bash
spgrn-check-pycistopic-completion
```

5.5-检查 RNA 和 cisTopic cell IDs：

```bash
spgrn-check-cistopic-cell-names
```

主要输出：

```text
inputs/cistopic_obj.pkl
inputs/region_sets/
results/pycistopic/qc/pycistopic_completion_check.md
results/pycistopic/model_selection/
```

默认 topic-model backend 是 MALLET。切换 CGS/MALLET 时，同一轮 topic models 应全部来自同一 backend；切换前 archive 或删除中断的 `Topic*.pkl`。

## 6. 标准化 region sets

目的：确保进入 motif enrichment 的 topic regions 和 DAR regions 都使用统一的 UCSC-style BED 格式。

生物信息学逻辑：region set 是 motif enrichment 的前景集合。格式错误、染色体命名不一致或非标准 chromosome 会导致 motif enrichment 背景不一致。

6.1-标准化 region sets：

```bash
spgrn-standardize-region-sets
```

预期结构：

```text
inputs/region_sets/
  Topics_otsu/
  Topics_top_3k/
  DARs_cell_label/
```

主要输出：

```text
results/pycistopic/qc/region_set_standardization.tsv
```

## 7. 构建 custom cisTarget database

目的：用本项目 consensus region universe 构建 motif ranking 和 score database。

生物信息学逻辑：公共 cisTarget database 的 region universe 未必匹配本项目可检测 peaks。custom database 让 motif enrichment 的背景和本项目 ATAC peak universe 一致。rankings 支持 cisTarget enrichment，scores 支持 DEM enrichment。

7.1-创建或更新 cisTarget 参数：

```bash
spgrn-setup-workflow-params --section cistarget
```

7.2-构建 custom database：

```bash
spgrn-build-custom-cistarget-db
```

主要输出：

```text
inputs/cistarget_db/custom.regions_vs_motifs.rankings.feather
inputs/cistarget_db/custom.regions_vs_motifs.scores.feather
results/cistarget_db/custom_cistarget_resource_plan.tsv
results/cistarget_db/custom_cistarget_db_manifest.tsv
```

大型数据库会使用官方 partial/combine/convert 路径，支持断点续跑。长时间 motif scanning 的 heartbeat 写入日志，用于区分“仍在运行”和“静默失败”。

## 8. 初始化 SCENIC+ Snakemake

目的：把 RNA、ATAC、region sets、custom cisTarget database 和 genome annotation 写入 SCENIC+ Snakemake 配置。

生物信息学逻辑：SCENIC+ 在所有 active metacells 上共同建立 shared eRegulon universe。多个 condition 应合并在这一步一起运行；condition-specific effects 在 Step 11 使用 shared AUC matrix 检验。

8.1-初始化配置：

```bash
spgrn-initialize-scenicplus-snakemake
```

主要输出：

```text
work/scenicplus/Snakemake/config/config.yaml
work/scenicplus/organism_config.yaml
inputs/scenicplus_config_params.tsv
```

通常需要按研究设计检查：

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

8.2-可选：大型 custom cisTarget database 分块运行 motif enrichment：

```bash
spgrn-run-scenicplus-motif-enrichment-split --mode both
spgrn-run-scenicplus-motif-enrichment-split --mode status
```

`--mode both` 会按 region-set family 拆分 DEM 和 cisTarget motif enrichment，并调用 `prepare_menr`。`--mode status` 是进入 Step 9 前的完成 gate：它检查 chunk HDF5 签名、DEM empty diagnostic、`cistromes_direct.h5ad`、`cistromes_extended.h5ad` 和 `tf_names.txt`。

主要输出：

```text
work/scenicplus/motif_enrichment_split/motif_enrichment_split_resource_plan.tsv
work/scenicplus/motif_enrichment_split/motif_enrichment_split_chunks.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.tsv
results/scenicplus_diagnostics/motif_enrichment_split_status.md
```

如果 `--mode status` 失败，不要进入 Step 9。重新运行 split；只有明确要覆盖旧半成品时才使用 `--force`。

## 9. SCENIC+ preflight checks

目的：在正式 Snakemake dry run 前检查 cell-name overlap、UCSC chromosome consistency、motif annotation、region names 和核心 config invariants。

生物信息学逻辑：preflight 失败通常说明输入定义不一致，而不是生物学结果不好。这里失败时继续跑 Snakemake 只会浪费时间并产生不可解释结果。

9.1-创建或更新 preflight 参数：

```bash
spgrn-setup-workflow-params --section preflight
```

9.2-运行 preflight：

```bash
spgrn-preflight-scenicplus-inputs 2>&1 | tee "$PROJECT_DIR/logs/preflight_scenicplus_inputs.log"
```

9.3-记录软件版本和资源 checksum：

```bash
spgrn-record-scenicplus-provenance
```

主要输出：

```text
inputs/preflight_thresholds.tsv
logs/preflight_scenicplus_inputs.log
results/provenance/
```

## 10. 运行 SCENIC+ inference 并记录稳定性

目的：运行 SCENIC+ inference，建立 TF-to-gene、region-to-gene 和 TF-region-gene eRegulons，并计算 AUCell activity。

生物信息学逻辑：dry run 只检查 DAG 和文件依赖；正式 run 才生成 eRegulon 结果。正式分析建议保留两次独立 inference 记录，用于评估 high-confidence edges 的稳定性。

10.1-创建或更新 Snakemake 参数：

```bash
spgrn-setup-workflow-params --section snakemake
```

10.2-运行 dry run：

```bash
spgrn-run-scenicplus-snakemake --mode dryrun
```

10.3-运行正式 inference：

```bash
spgrn-run-scenicplus-snakemake --mode run
```

SCENIC+ 内部主要规则：

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

10.4-归档 run1：

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

10.7-归档 run2：

```bash
mkdir -p results/scenicplus_stability/run2
cp results/scenicplus/eRegulons_direct.tsv results/scenicplus_stability/run2/
cp results/scenicplus/eRegulons_extended.tsv results/scenicplus_stability/run2/
cp work/scenicplus/region_to_gene_adj.tsv results/scenicplus_stability/run2/
```

10.8-比较稳定性：

```bash
spgrn-compare-scenicplus-stability \
  --run-a results/scenicplus_stability/run1 \
  --run-b results/scenicplus_stability/run2 \
  --out results/scenicplus_stability/stability_summary.tsv
```

主要输出：

```text
results/scenicplus/eRegulons_direct.tsv
results/scenicplus/eRegulons_extended.tsv
results/scenicplus/AUCell_direct.h5mu
results/scenicplus/AUCell_extended.h5mu
results/scenicplus/scplus_mdata.h5mu
results/scenicplus_stability/stability_summary.tsv
```

## 11. Postprocess、绘图和 condition statistics

目的：把完成的 SCENIC+ inference 转成输出层级审计、source tables、矢量 PDF 和 sample-level condition statistics。

生物信息学逻辑：网络推断和 condition testing 是两层问题。SCENIC+ 先在所有 active metacells 上共同推断 eRegulons；condition statistics 再在同一套 eRegulon AUCell activity 上按生物学样本比较。

11.1-创建或更新 postprocess 参数：

```bash
spgrn-setup-workflow-params --section postprocess
```

11.2-生成完整 postprocess 输出：

```bash
spgrn-run-scenicplus-postprocess --task all --layer all
```

参数：

```text
--task audit      只运行输出层级审计。
--task figures    重新生成 source tables 和 01-08 eRegulon PDFs。
--task stats      重新生成 09+ condition-statistics tables 和 PDFs。
--task all        生成完整 postprocess output set。

--layer direct    direct eRegulons。
--layer extended  extended eRegulons。
--layer all       同时处理 direct 和 extended。
```

主要输出：

```text
results/scenicplus_output_tiers/
results/scenicplus_figures/01-08_*.pdf
results/scenicplus_figures/09-15_*.pdf
results/scenicplus_figures/source_*.tsv
results/scenicplus_figures/*stats*.tsv
results/scenicplus_figures/plot_style_parameters.tsv
results/scenicplus_figures/plot_style_parameters.md
```

condition statistics 自动判断设计：

```text
单 condition：只输出 descriptive sample means，不做差异检验。
两组 condition：输出 delta_mean_auc = comparison - reference，并在 contrast 列记录方向。
多 condition：输出 sample-level omnibus test，不强行定义单一 delta。
```

只有自动识别的两组 reference 不符合实验设计时，才在 `inputs/postprocess_params.tsv` 修改 `reference_condition` 或 `comparison_condition`。

解释规则：

```text
AUC heatmap：查看 group-level eRegulon activity；行列按绘图矩阵聚类。
Condition heatmap：比较 cell label 内的 condition shift；列按 cell_label_condition 组织。
Specificity heatmap：优先筛选 cell-state-specific eRegulons。
Dot heatmap：颜色表示相对 AUC，点大小表示 active fraction。
Embedding：用于定位 eRegulon 活性，不单独作为统计证据。
Region-gene / overlap：解释 enhancer-target model structure，不直接代表 condition effect。
Network：展示 TF-target overview；正式主图建议用 source table 重画重点网络。
Volcano：优先解释 sample-level effect size；独立样本不足时谨慎解释 FDR。
eRegulon signs：+/+、-/+、-/- 是 SCENIC+ 的 regulon sign，不等同于 TF RNA expression。
```

正式论文图建议从 source tables 中选择少数关键 eRegulons 重绘，不要直接把全部 detected regulons 放入主图。direct 和 extended eRegulons 应分开呈现，除非图注明确说明合并规则。
