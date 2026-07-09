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

目的：固定本次 GRN 分析的 conda 环境和项目目录。这里的定义会写入 `$PROJECT_DIR/scenicplus_project.env` 和 `$PROJECT_DIR/project_env.sh`，后续步骤会逐步追加物种、annotated object、cell label 和 ATAC 输入参数。

生物信息学逻辑：先把运行环境和输出目录固定下来，后续每一步只在对应分析对象明确时写入参数，减少把错误对象或错误 ATAC 目录带入全流程的风险。

0.1-输入项目参数：

```bash
# 运行前必须替换本代码块中的所有示例值。
# CONDA_ROOT 是真实 conda/miniforge/miniconda/mambaforge/anaconda 根目录。
export CONDA_ROOT=/absolute/path/to/conda
# ENV_NAME 是 workflow 安装创建的 SCENIC+ 环境名。
export ENV_NAME=scenicplus-grn
# 将已安装环境命令加入当前终端 PATH，后续可以直接使用 spgrn-* 命令。
export PATH="$CONDA_ROOT/envs/$ENV_NAME/bin:$PATH"
# PROJECT_DIR 是本次 SCENIC+ 分析根目录，所有 inputs/work/logs/results 都写在这里。
export PROJECT_DIR=/absolute/path/to/grn_project/scenicplus_analysis
# AUTOZYME 设置为 on 或 off，控制 AutoZyme runtime 是否启用。
export AUTOZYME=on
```

0.2-运行已安装环境检查：

```bash
spgrn-check
```

0.3-初始化项目：

```bash
spgrn-initialize
```

0.4-激活环境并加载项目变量：

```bash
source "$CONDA_ROOT/bin/activate" "$ENV_NAME"
source "$PROJECT_DIR/project_env.sh"
```

## 1. 准备物种公共资源

目的：确认目标物种在指定 Ensembl release 中同时有 FASTA 和 GTF，并准备 SCENIC+ 需要的 motif collection 与 motif2TF table。FASTA/GTF 来自 Ensembl；其余 SCENIC+ 所需标准资源由 workflow 派生并审计。

1.0-选择物种和 Ensembl release，并写入项目配置：

```bash
# ENSEMBL_RELEASE 默认使用 115；如需其他版本，先用 --release 查询并输出可用物种 TSV。
export ENSEMBL_RELEASE=115
spgrn-query-organism-resources --list --release "$ENSEMBL_RELEASE"

# ORGANISM 需要按查询结果修改。建议使用 Ensembl species name；常用别名如 mouse 也可识别。
export ORGANISM=mus_musculus

# 查询本次选择的物种资源状态。
spgrn-query-organism-resources --organism "$ORGANISM" --release "$ENSEMBL_RELEASE"

# 将 ORGANISM 和 ENSEMBL_RELEASE 写入 scenicplus_project.env 与 project_env.sh。
spgrn-initialize
source "$PROJECT_DIR/project_env.sh"
```

`--list` 会写出 `ensembl_release_<release>_organism_resources.tsv`，只包含同时有 FASTA 和 GTF 的 Ensembl species，并显示 motif2TF 准备策略。`--organism` 用于查看单个物种的 assembly、FASTA/GTF 可及性、motif collection 和 motif2TF 策略。

1.1-准备资源。物种资源准备必须选择一种 motif2TF 路径：

- direct：human、mouse、fly、chicken 有 Aerts v10 direct motif2TF table。

```bash
spgrn-prepare-official-resources
```

- mapping：其他物种可显式指定 `--ref human|mouse|fly|chicken`，通过 Ensembl BioMart orthology 映射并输出 audit。

```bash
spgrn-prepare-official-resources --ref human
```

如果研究对象存在明显 paralog 扩张，使用 paralog-aware orthology 策略，并检查 mapping audit：

```bash
spgrn-prepare-official-resources --ref human --orthology-policy paralog-aware
```

- generated：没有可靠 direct table 或默认映射时，可用 Aerts direct motif evidence 与目标物种 gene symbols 生成 audited symbol-evidence table。

```bash
spgrn-prepare-official-resources --generate-motif2tf
```

generated table 使用目标物种注释中实际存在的 gene symbols，严格匹配 Aerts direct human/mouse/fly/chicken motif2TF evidence，并输出覆盖率、来源和未匹配 TF 的 audit。它适合探索性非模式物种或跨物种比较初筛；正式结论应检查 `resources/<organism>/*_motif2tf_generated_symbol_audit.tsv`。

- user table：已有人工审计过的 species-specific motif2TF table 时，直接提供该表。

```bash
# MOTIF2TF_TABLE 需要改成真实 motif_annotations.tbl 路径。
export MOTIF2TF_TABLE=/absolute/path/to/motif_annotations.tbl
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
resources/motif2tf/motif_annotations.<organism>.<strategy>[.<reference>].tbl
resources/motif2tf/motif_annotations.<organism>.active.tsv
inputs/cistarget_db/motif_annotations.tbl
resources/<organism>/*_chromosome_audit.tsv
resources/<organism>/*_motif2tf_*_audit.tsv
```

`resources/motif2tf/` 中保留带物种和策略名的 canonical motif2TF table；`inputs/cistarget_db/motif_annotations.tbl` 是当前 active species 安装给 SCENIC+ 读取的标准文件名。

## 2. 检查并导出 active annotated object

目的：从已注释 scMultiome 对象中导出 active cells/metacells 的 RNA counts 和 metadata。这里不重新注释细胞；它只确认下游 SCENIC+ 使用哪些 assay/layer、sample、condition、cell label、barcode 和 embedding。

生物信息学逻辑：`CELL_LABEL_COLUMN` 决定 metacell aggregation、DAR calling、region-set 命名和 eRegulon AUC 汇总的生物学分组。它应该代表本次 GRN 想比较的主要细胞类型、状态或发育阶段。

2.1-设置并检查 annotated object：

```bash
# ANNOTATED_OBJECT 是已注释且已经限定到本次 active cells 的 Seurat RDS/QS 或 AnnData h5ad。
export ANNOTATED_OBJECT=/absolute/path/to/active_annotated_multiome_object.rds
spgrn-initialize
source "$PROJECT_DIR/project_env.sh"
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

2.3-设置 cell label 并导出 active RNA 和 metadata：

```bash
# CELL_LABEL_COLUMN 是 annotated object 中作为 GRN cell label 的 metadata column。
export CELL_LABEL_COLUMN=cell_annotation
spgrn-initialize
source "$PROJECT_DIR/project_env.sh"
spgrn-export-annotated-object
```

如果需要覆盖已选择的 cell label column：

```bash
# corrected_metadata_column 需要手动修改为真实 metadata column 名称。
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

本 workflow 默认 metacell 是必经步骤。该步骤依赖 hdWGCNA 和 Seurat RDS/QS。如果 Step 2 从 h5ad 开始，需要在 `inputs/metacell_params.tsv` 中提供匹配的 Seurat 对象。

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

4.1-注册 ATAC 输入参数：

```bash
# ATAC_INPUT_LAYOUT 是支持的 ATAC 文件布局之一，需要和 ATAC_DATA_ROOT 目录结构匹配。
export ATAC_INPUT_LAYOUT=split_ge_arc
# ATAC_DATA_ROOT 是 ATAC fragments 和 peaks 所在根目录。
export ATAC_DATA_ROOT=/absolute/path/to/atac_input_root
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
work/standard_peaks/
work/metacell_fragments/
results/metacells/metacell_fragment_reassignment.tsv
```

继续前检查 `inputs/sample_sheet.tsv` 中的 sample、condition、fragment 路径和 peak 路径。标准化后的 peaks 以及 metacell fragments 路径会写回 `inputs/sample_sheet.tsv`。

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
spgrn-run-pycistopic-workflow
```

该命令会自动在 `$PROJECT_DIR/logs/` 下写入带时间戳的日志。

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

目的：把 RNA、ATAC、region sets、custom cisTarget database 和物种资源写入 SCENIC+ Snakemake 配置。

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

如果获得的 TF-GRN/eRegulon 很少，先按层级判断是哪一个 gate 在减少结果。SCENIC+ 不是只做 TF-gene 共表达；它要求 motif enriched regions、motif2TF annotation、TF expression、TF-to-gene links、region-to-gene links 和 TF-region-gene 三者一致性共同成立。官方参数说明见 [SCENIC+ Running tutorial: params_motif_enrichment / params_inference](https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html) 和 [SCENIC+ API: eRegulon enrichment](https://scenicplus.readthedocs.io/en/latest/api.html#eregulon-enrichment-in-cells)。

最直接影响 TF-GRN 捕获数目的参数：

| 分析 gate | 当前流程位置 | 参数 | 越严格时的影响 |
|---|---|---|---|
| topic/DAR region sets | Step 5 | `pycistopic.ntop_regions`, `pycistopic.dar_adjpval_thr`, `pycistopic.dar_log2fc_thr` | 进入 motif enrichment 的 regions 变少或偏弱，后续 TF 候选减少。 |
| cisTarget motif enrichment | Step 8.1/8.2 | `scenicplus_config.ctx_nes_threshold`, `scenicplus_config.ctx_auc_threshold`, `scenicplus_config.ctx_rank_threshold` | enriched motifs 和 cistromes 减少。 |
| DEM motif enrichment | Step 8.1/8.2 | `scenicplus_config.dem_adj_pval_thr`, `scenicplus_config.dem_log2fc_thr`, `scenicplus_config.dem_mean_fg_thr`, `scenicplus_config.dem_motif_hit_thr` | DEM 支持的 motif/TF 减少；DAR-derived evidence 尤其敏感。 |
| motif2TF annotation | Step 1/8.1 | `scenicplus_config.motif_similarity_fdr`, `scenicplus_config.orthologous_identity_threshold`, motif2TF direct/mapped/generated/user 策略 | motif enriched 但映射不到 TF，direct/extended cistromes 会少。 |
| enhancer-gene search space | Step 8.1 | `scenicplus_config.search_space_upstream`, `scenicplus_config.search_space_downstream`, `scenicplus_config.search_space_extend_tss` | 候选 region-gene pairs 变少或变宽；过窄会漏 distal enhancers，过宽会增加噪声。 |
| region-to-gene pruning | Step 8.1/10 | `scenicplus_config.quantile_thresholds_region_to_gene`, `scenicplus_config.top_n_regionTogenes_per_gene`, `scenicplus_config.top_n_regionTogenes_per_region`, `scenicplus_config.min_regions_per_gene` | 进入 eGRN 合并的 region-gene links 减少。 |
| TF-gene / region-gene concordance | Step 8.1/10 | `scenicplus_config.rho_threshold` | 相关性 gate 更严，方向一致的 TF-region-gene triples 减少。 |
| eRegulon size filter | Step 8.1/10 | `scenicplus_config.min_target_genes` | 小 target-set TF 被过滤；这是最常见的最终 eRegulon 数量 gate。 |

8.1a-可选：调整 SCENIC+ 捕获阈值。下面是探索性诊断示例，不建议不加审计地作为正式结论：

```bash
spgrn-setup-workflow-params --section scenicplus_config \
  --set scenicplus_config.min_target_genes=5 \
  --set scenicplus_config.ctx_nes_threshold=2.5 \
  --set scenicplus_config.dem_adj_pval_thr=0.10 \
  --set scenicplus_config.rho_threshold=0.03

spgrn-initialize-scenicplus-snakemake
```

如果怀疑 Step 5 产生的 topic/DAR region sets 太少，调整 pycisTopic 参数并从 Step 5 重新运行：

```bash
spgrn-setup-workflow-params --section pycistopic \
  --set pycistopic.ntop_regions=5000 \
  --set pycistopic.dar_adjpval_thr=0.10

spgrn-run-pycistopic-workflow
spgrn-standardize-region-sets
```

参数修改后的生效范围：

- Step 5 参数改变后，必须重新运行 Step 5-7，因为 region sets 和 custom cisTarget database 会改变。
- Step 8.1 中 `scenicplus_config.*` 改变后，必须重新运行 `spgrn-initialize-scenicplus-snakemake`。
- 如果已经运行过 Step 8.2 split motif enrichment，修改 DEM/cisTarget 阈值后需要用 `--force` 重新运行 Step 8.2。
- 如果只修改 `rho_threshold`、`min_target_genes` 或 region-to-gene pruning 参数，至少需要重新运行 Step 10 inference。

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
spgrn-preflight-scenicplus-inputs
```

9.3-记录软件版本和资源 checksum：

```bash
spgrn-record-scenicplus-provenance
```

主要输出：

```text
inputs/preflight_thresholds.tsv
logs/preflight_scenicplus_inputs_*.log
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

`inputs/postprocess_params.tsv` 中 `regulon_sign_filter=tf_positive` 是默认正式展示口径，只保留 `+/+` eRegulons；改为 `all` 时才输出全部 signed eRegulons。

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
eRegulon signs：`+/+`、`-/+`、`-/-` 是 SCENIC+ 的 regulon sign，不等同于 TF RNA expression。主展示和 condition statistics 默认只保留 TF-positive `+/+` eRegulons；需要探索全部 signed regulons 时，可在 `inputs/postprocess_params.tsv` 中将 `regulon_sign_filter` 改为 `all`。
```

正式论文图建议从 source tables 中选择少数关键 eRegulons 重绘，不要直接把全部 detected regulons 放入主图。direct 和 extended eRegulons 应分开呈现，除非图注明确说明合并规则。
