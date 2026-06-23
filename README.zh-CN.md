# SCENIC+ GRN Conda 安装器

这个安装包用于安装一个独立的 conda 环境，专门用于 SCENIC+/GRN 分析。设计目标是：可以把这个小型压缩包复制到工作站或 Linux 服务器，任意目录解压，然后运行 `bash install.sh` 完成安装。

## 打包内容和作用

```text
Installer entry points:
  install.sh                         Required. Main installer and bootstrap script.
  check_environment.sh               Required. Environment and workflow self-check entry point.
  install_r.R                        Required R/hdWGCNA/Seurat/Signac installation layer for the metacell workflow.

Environment recipes and pinned Python layer:
  environment-linux-64.yml           Required on Linux x86_64, compatible with glibc >= 2.17.
  environment-macos-arm64.yml        Required on Apple Silicon macOS.
  pip-constraints.txt                Required by the pinned pip supplement.

Offline/restricted-network source archives:
  archives/vendor.tar.gz             Required in release archives for robust installation when GitHub is slow or blocked.
                                      install.sh 会在运行时解压到隐藏缓存 .vendor/。
  .vendor/github/                    运行时解出的 pinned source archives；GitHub 重试失败后使用。
  .vendor/mallet/                    运行时解出的 MALLET 2.0.8 archive；默认 MALLET 安装使用。

Installed workflow assets:
  initialize_scenicplus_project.sh   User-facing one-step initializer after installation.
                                      Checks the environment, updates scenicplus_project.env, then initializes the project runtime files.
  scripts/                           Required executable workflow entry points after installation.
                                      Installed to $CONDA_PREFIX/share/scenicplus-grn/scripts/.
                                      Includes project parameter setup, raw-data sample-sheet generation,
                                      pycisTopic, cisTarget DB, SCENIC+ config, Snakemake and postprocessing wrappers.
  modules/                           Required internal helper modules imported by scripts; not user-facing commands.
                                      Installed to $CONDA_PREFIX/share/scenicplus-grn/modules/.
  scenicplus_config_template.yaml     Required Snakemake config template used by workflow generators.

Documentation and audit records:
  README.md                          User-facing installer guide.
  README.zh-CN.md                    Chinese installer guide.
  SCENICPLUS_STEP_BY_STEP.md          Strict matched snRNA+snATAC workflow guide.
  SCENICPLUS_STEP_BY_STEP.zh-CN.md    Chinese matched snRNA+snATAC workflow guide.
  VERSION                            Installer package version.
  CHANGELOG.md                       Release history.
  RELEASE_NOTES.md                   Current release notes.
  VERSION_LOCK.md                    Human-readable record of pinned analysis software versions.
  locks/environment-linux-64.solved-lock.yml
                                      Linux dry-run conda lock record with build strings for auditing/debugging.
```

## 推荐使用方式

把压缩包复制到目标机器，解压并运行安装器：

```bash
tar -xzf scenicplus-grn-installer.tar.gz
cd scenicplus-grn-installer
bash install.sh
```

脚本会先查找 conda/miniforge/miniconda/mamba/anaconda 这类 conda 根目录。也可以显式指定：

```bash
CONDA_ROOT=/path/to/conda bash install.sh
```

如果解压目录不在检测到的 conda 根目录下面，安装器会询问是否把自身复制到：

```text
$CONDA_ROOT/share/scenicplus-grn-installer
```

然后从那里继续安装。这样可复用安装包会放在 conda 安装旁边，而不是放进某个具体项目目录。

## 运行时会发生什么

`install.sh` 会依次执行：

```text
1. 启动带时间戳的安装日志。
2. 检测 conda/miniforge/miniconda/mamba 风格的根目录。
3. 请求确认检测到的根目录。
4. 必要时询问是否复制安装器到 $CONDA_ROOT/share/scenicplus-grn-installer。
5. 检查 $CONDA_ROOT、envs、pkgs 和 share 的写权限。
6. 如果 base 中没有 mamba，则先用 conda 安装 mamba。
7. 创建或更新独立环境 scenicplus-grn。
8. 安装由 conda/mamba 解析的 Python、命令行工具、基因组学工具和 R 基础层。
9. 必要时先把 `archives/vendor.tar.gz` 解压成 `.vendor/`，再安装小型 pinned pip supplement，以及固定 commit 的 SCENIC+/pycisTopic/pycistarget 源码层。每个 GitHub 源会尝试 3 次；如果 GitHub 不稳定，则使用 `.vendor/github` 中的本地归档。
10. 默认安装 MALLET 2.0.8，用于 pycisTopic MALLET LDA backend，并通过一次小型 import-file 测试验证 wrapper。
11. 安装 metacell workflow 所需的 R/hdWGCNA 层；GitHub 不可用时使用 bundled pinned source archive。
12. 运行 check_environment.sh。
13. 把安装配方复制到 $CONDA_PREFIX/share/scenicplus-grn。
14. 从已安装环境副本中运行 workflow asset checker。
```

默认会创建或更新独立环境：

```text
scenicplus-grn
```

除非显式使用 `MODE=active`，否则不会安装到当前已经激活的环境里。

## 日志

每次运行都会写日志。

默认位置：

```text
scenicplus-grn-installer/logs/install_YYYYMMDD_HHMMSS.log
```

如果解压目录不可写，日志会写到：

```text
/tmp/scenicplus-grn-installer-logs
```

也可以指定日志目录：

```bash
LOG_DIR=/path/to/logs bash install.sh
```

## 后台安装

完整安装在共享服务器上可能会跑很久。只要日志在写，不需要一直盯着屏幕。

如果有 `tmux`，推荐：

```bash
tmux new -s scenicplus-install
cd scenicplus-grn-installer
ASSUME_YES=1 CONDA_ROOT=/path/to/miniconda3 bash install.sh
```

按 `Ctrl-b` 然后 `d` 退出会话。之后重新进入：

```bash
tmux attach -t scenicplus-install
```

没有 `tmux` 时可用简单后台模式：

```bash
cd scenicplus-grn-installer
mkdir -p logs
nohup env ASSUME_YES=1 CONDA_ROOT=/path/to/miniconda3 bash install.sh \
  > logs/nohup_install_$(date +%Y%m%d_%H%M%S).out 2>&1 &
echo $! > logs/install.pid
```

查看是否仍在运行：

```bash
ps -p "$(cat logs/install.pid)"
```

追踪最新安装日志：

```bash
tail -f "$(ls -t logs/install_*.log | head -n 1)"
```

成功安装的结尾会出现：

```text
DONE: SCENIC+ environment is installed and checked.
```

如果失败，把最新的 `logs/install_*.log` 发给维护者诊断。

## 主要命令

在 Linux 上，如果 `$CONDA_ROOT/bin/mamba` 不存在，`install.sh` 会自动把 `mamba` 安装到 conda `base`。这个 bootstrap 步骤使用已有 conda solver，以兼容较旧 Miniconda；mamba 可用后，分析环境创建阶段会使用 mamba，避免经典 conda solver 长时间陷入冲突解析。

关闭自动 mamba bootstrap：

```bash
AUTO_INSTALL_MAMBA=0 bash install.sh
```

创建或更新独立环境：

```bash
bash install.sh
```

指定 conda 根目录：

```bash
CONDA_ROOT=/path/to/miniforge3 bash install.sh
```

创建不同名称的环境：

```bash
ENV_NAME=scenicplus-grn-test bash install.sh
```

修改 GitHub 重试次数：

```bash
GITHUB_TRIES=1 bash install.sh
```

完全跳过 GitHub，直接使用 bundled archives：

```bash
GITHUB_TRIES=0 bash install.sh
```

从头重建环境：

```bash
FORCE=1 bash install.sh
```

跳过 MALLET 安装：

```bash
INSTALL_MALLET=0 bash install.sh
```

默认安装 MALLET，推荐用于较大的 pycisTopic topic model。安装器会用小型 `import-file` run 验证 wrapper。若关闭 MALLET，workflow 仍可用 `pycistopic.lda_backend=cgs`。

非交互式安装：

```bash
ASSUME_YES=1 bash install.sh
```

只检查检测、迁移、权限和平台 recipe，不真正安装：

```bash
PRECHECK_ONLY=1 bash install.sh
```

不把安装器复制到 conda root 下：

```bash
RELOCATE_INSTALLER=0 bash install.sh
```

安装到当前激活的 conda 环境：

```bash
source /path/to/miniforge3/bin/activate my-existing-env
MODE=active bash install.sh
```

安装器默认拒绝安装到 `base`，除非显式允许：

```bash
ALLOW_BASE=1 MODE=active bash install.sh
```

## 安装后

在终端中输入环境和项目参数：

```bash
CONDA_ROOT=/path/to/conda
ENV_NAME=scenicplus-grn
PROJECT_DIR=/path/to/grn_project/scenicplus_analysis
ORGANISM=mouse
AUTOZYME=on
ENSEMBL_RELEASE=115
ANNOTATED_OBJECT=/path/to/active_annotated_multiome_object.rds
CELL_LABEL_COLUMN=cell_annotation
ATAC_INPUT_LAYOUT=split_ge_arc
ATAC_DATA_ROOT=/path/to/atac_input_root
```

`PROJECT_DIR` 应该是专门用于 SCENIC+ 分析的根目录，而不是更大的研究项目目录，除非你确实希望 workflow 产生的 `inputs/`、`work/`、`resources/`、`logs/` 和 `results/` 都放在那里。

运行一步式初始化。它会先检查环境，更新 `$PROJECT_DIR/scenicplus_project.env`，并初始化项目运行文件：

```bash
env \
  CONDA_ROOT="$CONDA_ROOT" \
  ENV_NAME="$ENV_NAME" \
  PROJECT_DIR="$PROJECT_DIR" \
  ORGANISM="$ORGANISM" \
  AUTOZYME="$AUTOZYME" \
  ENSEMBL_RELEASE="$ENSEMBL_RELEASE" \
  ANNOTATED_OBJECT="$ANNOTATED_OBJECT" \
  CELL_LABEL_COLUMN="$CELL_LABEL_COLUMN" \
  ATAC_INPUT_LAYOUT="$ATAC_INPUT_LAYOUT" \
  ATAC_DATA_ROOT="$ATAC_DATA_ROOT" \
  bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/initialize_scenicplus_project.sh"
```

出现 `PROJECT INITIALIZATION OK` 后，加载项目运行变量：

```bash
source "$PROJECT_DIR/project_env.sh"
```

项目设置文件写在 `$PROJECT_DIR/scenicplus_project.env`。workflow 脚本会在启动时自动检测可用内存、当前负载和 CPU 数量，并选择保守的并行 worker 数。

当前 recipe 的核心版本预期：

```text
pycisTopic 1.0.2
pycistarget 1.1
scenicplus 1.0a2
MACS2 2.2.9.1
Snakemake 7.32.4
MALLET 2.0.8
```

## 平台检测

`install.sh` 会自动选择环境文件：

```text
macOS arm64      -> environment-macos-arm64.yml
Linux x86_64     -> environment-linux-64.yml
```

Linux 下，如果没有设置 `CONDA_SUBDIR`，脚本会导出 `CONDA_SUBDIR=linux-64`。这个 recipe 面向常见 Linux x86_64 服务器，以 glibc 2.17 作为最低兼容线；通常也适用于更新 glibc 的 Rocky/Alma/CentOS Stream、Ubuntu 和 Debian。低于 glibc 2.17 的服务器不受支持。Linux ARM/aarch64 需要单独 recipe。

## 安装内容

安装器版本锁定到本地 macOS 成功 dry run。Linux recipe 是 conda-first：Python、CLI、genomics、SCENIC 相关包和 R base 由 mamba 从 conda-forge/bioconda 解析，并兼容 CentOS7/glibc 2.17。`environment-linux-64.yml` 固定顶层安装版本；`locks/environment-linux-64.solved-lock.yml` 记录完整 dry-run solve，包括传递依赖 build string，用于审计和排错。pip 只用于小型 pinned supplement 和精确 GitHub source commits。SCENIC+/pycisTopic/pycistarget source commits 在不同平台保持一致。release archive 还包含 `archives/vendor.tar.gz`；`install.sh` 会在安装前把它解压为 `.vendor/`，避免受限网络下 GitHub 失败阻断安装。本地 bundled-source 安装时，安装器会显式把 pinned package versions 传给 setuptools-scm，因此 GitHub archive tarball 不需要 `.git` 元数据也能可重复构建。完整版本摘要见 `VERSION_LOCK.md`。

Workflow 脚本和 helper modules 安装在独立 conda 环境中，不安装到全局系统目录。激活后，`$CONDA_PREFIX` 指向该环境，脚本位于：

```text
$CONDA_PREFIX/share/scenicplus-grn/scripts/
```

step-by-step 指南中为了方便会定义：

```bash
export SCENICPLUS_HOME="$CONDA_PREFIX/share/scenicplus-grn"
```

随后用 `$SCENICPLUS_HOME/scripts/<script>.py` 调用脚本。

核心 Python、命令行和基因组学层：

```text
Python                       3.10.20 macOS; 3.10.19 Linux glibc2.17 recipe
setuptools                   80.10.2
wheel                        0.47.0
MACS2                        2.2.9.1
Snakemake                    7.32.4
scanpy                       1.10.4
anndata                      0.10.8
mudata                       0.2.3
pandas                       1.5.3
numpy                        1.23.5
scipy                        1.14.1
scikit-learn                 1.7.2
pyscenic                     0.12.1
ctxcore                      0.2.0
arboreto                     0.1.6
dask/distributed             2023.4.1
ray-default                  2.53.0 macOS; 2.9.3 Linux glibc2.17 recipe
pyarrow                      23.0.0 macOS; 14.0.2 Linux glibc2.17 recipe
polars                       1.41.2 macOS; 1.35.2 Linux glibc2.17 recipe
flatbuffers                  25.12.19 macOS; 25.9.23 Linux glibc2.17 recipe
python-flatbuffers           25.9.23
leidenalg                    0.12.0 macOS; 0.10.2 Linux glibc2.17 recipe
python-igraph                1.0.0 macOS; 0.11.8 Linux glibc2.17 recipe
lxml                         6.1.1 macOS; 5.3.1 Linux glibc2.17 recipe
bokeh                        3.9.1
tornado                      6.5.7 macOS; 6.5.2 Linux glibc2.17 recipe
MALLET                       2.0.8, installed under $CONDA_PREFIX/opt/mallet-2.0.8 with $CONDA_PREFIX/bin/mallet Python wrapper
```

Aerts / SCENIC+ source layer:

```text
pycisTopic                   1.0.2, commit 219225df56b32738d82cd14532b187a1483de04f
pycistarget                  1.1, commit 5aa517604e4842539a7531c16905825dc7cb80fb
scenicplus                   1.0a2, commit e82b82f14b76618b850dfe442efc2421bb34f3b4
create_cisTarget_databases   commit 304d5dc1b15e5c923908a50a1ec291c3faaccf9c
Cluster-Buster / cbust       commit 5911cd6201b767a43316ce613afc6c9255dc3511
LoomXpy                      0.4.2, commit 61995ff10940968eac2cee8fe48300ab477a15d0
```

AutoZyme acceleration layer:

```text
AutoZyme Python              0.3.1, commit 35f91f2229eb44d82710470803865d3c15102716, installed with pip --no-deps --no-build-isolation
AutoZyme R                   0.3.1, same commit, installed from bundled source with dependencies=FALSE
r-rcppdist                   0.1.1, conda-provided for R AutoZyme linking
RcppParallel                 5.1.11-1, installed from CRAN source with dependencies=FALSE if absent
```

AutoZyme 作为 no-dependency overlay 处理，不能升级、降级或替换已固定的 Scanpy、Seurat、SCENIC+、pycisTopic 或 R 包栈。安装器默认 `INSTALL_AUTOZYME=1` 和 `INSTALL_AUTOZYME_R=1`。分析脚本运行时可用 `AUTOZYME_DISABLED=1` 关闭 patch 激活。AutoZyme 可能输出与邻近上游版本有关的 warning；这些 warning 会进入日志，不改变包版本。

GitHub 重试失败后使用的 bundled source archives 存在于 `archives/vendor.tar.gz` 内，运行时会解压为：

```text
.vendor/github/pycisTopic-219225df56b32738d82cd14532b187a1483de04f.tar.gz
.vendor/github/pycistarget-5aa517604e4842539a7531c16905825dc7cb80fb.tar.gz
.vendor/github/scenicplus-e82b82f14b76618b850dfe442efc2421bb34f3b4.tar.gz
.vendor/github/create_cisTarget_databases-304d5dc1b15e5c923908a50a1ec291c3faaccf9c.tar.gz
.vendor/github/cluster-buster-5911cd6201b767a43316ce613afc6c9255dc3511.tar.gz
.vendor/github/LoomXpy-61995ff10940968eac2cee8fe48300ab477a15d0.tar.gz
.vendor/github/hdWGCNA-afa09abb890f5be087b63e510a7346e8e1952ecc.tar.gz
.vendor/github/SHA256SUMS
```

Bundled archives 使用前会做 checksum 验证。若上传后的压缩包损坏，安装器会以 checksum error 停止，不会继续使用不完整源码树。

R、Bioconductor 和 hdWGCNA 层：

```text
R                            4.5.3 macOS; 4.4.3 Linux glibc2.17 recipe
cmake                        4.3.3 Linux glibc2.17 recipe
pkg-config                   0.29.2 Linux glibc2.17 recipe
libuv                        1.52.1 Linux glibc2.17 recipe
xz / liblzma / liblzma-devel 5.8.3 Linux glibc2.17 recipe
BiocManager                  1.30.27
remotes                      2.5.0
Seurat                       5.5.0
Signac                       1.17.1 macOS/R 4.5; 1.16.0 Linux/R 4.4
WGCNA                        1.74
R igraph                     2.1.4 Linux/R 4.4
hdWGCNA                      0.4.11, commit afa09abb890f5be087b63e510a7346e8e1952ecc
GenomicRanges                1.62.1 macOS/R 4.5; 1.58.0 Linux/R 4.4
GeneOverlap                  1.46.0 macOS/R 4.5; 1.42.0 Linux/R 4.4
UCell                        2.14.0 macOS/R 4.5; 2.10.1 Linux/R 4.4
impute                       1.84.0 macOS/R 4.5; 1.80.0 Linux/R 4.4
preprocessCore               1.72.0 macOS/R 4.5; 1.68.0 Linux/R 4.4
fs                           2.1.0 Linux/R 4.4
Hmisc                        5.2-5 Linux/R 4.4
htmlTable                    2.5.0 Linux/R 4.4
htmlwidgets                  1.6.4 Linux/R 4.4
rmarkdown                    2.31 Linux/R 4.4
```

Linux 下，`cmake`、`pkg-config`、`libuv`、`xz`、`liblzma-devel`、Seurat、Signac、WGCNA、R igraph 和重型 Seurat 依赖栈会先通过 conda 安装，再运行 R source layer。这样能避免服务器系统 CMake 太旧导致 CRAN 包如 `fs` 编译失败，也避免在共享服务器上源码编译 R igraph/Seurat。R igraph 固定在 R 4.4 兼容的 `2.1.4` build，以保留 SCENIC+ 层需要的 Python `lxml 5.3.1` / `libxml2 2.13` 栈。

hdWGCNA/plotting/network utilities 使用的 pinned R source dependencies：

```text
systemfonts                  1.3.2
tweenr                       2.0.3
WriteXLS                     6.8.0
rjson                        0.2.23
RhpcBLASctl                  0.23-42
graphlayouts                 1.2.3
tidygraph                    1.3.1 macOS/R 4.5; 1.3.0 Linux/R 4.4
ggforce                      0.5.0
proxy                        0.4-29
tester                       0.3.0
enrichR                      3.4
harmony                      2.0.4 macOS/R 4.5; 2.0.2 Linux/R 4.4
ggraph                       2.2.2
```

其他 pinned pip dependencies 见 `pip-constraints.txt`；安装器会 force-reinstall 该 pip 层，减少后续版本漂移。

`check_environment.sh` 会在加载 hdWGCNA 前设置 `options(enrichR.live = FALSE)`。这避免 Enrichr 网络服务慢、被阻断或不可用时造成误失败。它不会关闭 SCENIC+，只是防止包加载时的网络检查破坏环境验证。

## 排错

如果没有检测到 conda：

```bash
CONDA_ROOT=/path/to/miniforge3 bash install.sh
```

如果在临时目录解压并跳过了 relocation，可以之后手动放到 conda 下：

```bash
mkdir -p ~/miniforge3/share
cp -a scenicplus-grn-installer ~/miniforge3/share/
cd ~/miniforge3/share/scenicplus-grn-installer
bash install.sh
```

关闭 AutoZyme 安装但保留其余环境：

```bash
INSTALL_AUTOZYME=0 INSTALL_AUTOZYME_R=0 bash install.sh
```

保留 AutoZyme 安装，但在分析中关闭 runtime patch：

```bash
AUTOZYME_DISABLED=1 python your_analysis.py
AUTOZYME_DISABLED=1 Rscript your_analysis.R
```

只验证已有安装：

```bash
CONDA_ROOT=/path/to/conda
ENV_NAME=scenicplus-grn
bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/check_environment.sh" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME"
```
