# scenicplus-grn-workflow

用于 SCENIC+ gene regulatory network 分析的可复现 workflow toolkit，面向已经完成注释的 matched snRNA+snATAC / scMultiome 数据。它会创建独立的 `scenicplus-grn` conda 环境，安装固定版本的 Python、R、基因组学工具和 SCENIC+ 相关软件，并提供 `spgrn-*` 命令用于项目初始化、公共资源准备、pycisTopic、自定义 cisTarget database、SCENIC+ Snakemake 运行和结果整理。

推荐从 [GitHub Releases](https://github.com/JunyanKang/scenicplus-grn-workflow/releases) 下载 `scenicplus-grn-workflow-v*.tar.gz` 完整包安装；该包包含 `archives/vendor.tar.gz`，适合服务器和网络不稳定环境。GitHub 页面绿色 `Code` 下载的 ZIP 只包含源码，不包含离线 vendor archive，安装时需要能访问 GitHub 下载第三方源码。

## 适用场景

这个 workflow 适合需要在工作站或 Linux 服务器上复用同一套 SCENIC+ 分析环境，
并希望 workflow 脚本、版本记录和离线源码归档一起分发的项目。

适用于：

- 已注释的 scMultiome 对象和匹配的 ATAC fragments，
- metacell-based SCENIC+ 分析，
- 基于项目 consensus region universe 构建 custom cisTarget database，
- direct、orthology-mapped、audited generated 或 user-supplied species-specific motif2TF table 准备和审计，
- 可复现 rerun、稳定性检查和结果输出。

详细分析流程见：[SCENICPLUS_STEP_BY_STEP.md](https://github.com/JunyanKang/scenicplus-grn-workflow/blob/main/docs/SCENICPLUS_STEP_BY_STEP.md)。

## 支持平台

```text
Linux x86_64      glibc >= 2.17
macOS arm64       Apple Silicon
```

安装前需要已有 conda 风格目录，例如 Miniforge、Miniconda、Mambaforge、
Anaconda 或 `/opt/conda`。

## 快速安装

在目标机器解压 release archive：

```bash
tar -xzf scenicplus-grn-workflow.tar.gz
cd scenicplus-grn-workflow
```

运行安装：

```bash
CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

非交互安装：

```bash
ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

成功安装时日志结尾为：

```text
DONE: SCENIC+ environment is installed and checked.
```

## 后台安装

服务器安装耗时较长时，建议用 `tmux` 或 `nohup`。

使用 `tmux`：

```bash
tmux new -s scenicplus-install
cd scenicplus-grn-workflow
ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

按 `Ctrl-b` 然后 `d` 退出会话，之后重新进入：

```bash
tmux attach -t scenicplus-install
```

不用 `tmux`：

```bash
cd scenicplus-grn-workflow
mkdir -p logs
nohup env ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh \
  > logs/nohup_install_$(date +%Y%m%d_%H%M%S).out 2>&1 &
echo $! > logs/install.pid
```

追踪最新安装日志：

```bash
tail -f "$(ls -t logs/install_*.log | head -n 1)"
```

## 安装参数

| 变量 | 默认值 | 含义 |
|---|---:|---|
| `CONDA_ROOT` | 自动检测 | 使用的 conda 根目录。 |
| `ENV_NAME` | `scenicplus-grn` | 创建或更新的环境名。 |
| `ASSUME_YES` | `0` | 设为 `1` 后使用非交互模式。 |
| `FORCE` | `0` | 设为 `1` 后重建环境。 |
| `GITHUB_TRIES` | `3` | GitHub 失败多少次后使用本地归档。 |
| `INSTALL_MALLET` | `1` | 安装 pycisTopic MALLET topic-model backend。 |
| `AUTO_INSTALL_MAMBA` | `1` | 如果 base 里没有 mamba，自动安装 mamba。 |
| `PRECHECK_ONLY` | `0` | 只检查平台、路径和权限，不安装。 |
| `RELOCATE_INSTALLER` | `1` | 询问是否复制 workflow 包到 `$CONDA_ROOT/share/scenicplus-grn-workflow`。 |
| `LOG_DIR` | `logs/` | 安装日志目录。 |

示例：

```bash
ENV_NAME=scenicplus-grn-test CONDA_ROOT=/absolute/path/to/conda bash install.sh
GITHUB_TRIES=0 ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
FORCE=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
PRECHECK_ONLY=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

## 安装内容

- SCENIC+、pycisTopic、pycistarget 和 create_cisTarget_databases，
- Python 单细胞和基因组学依赖，
- R、Seurat、Signac、hdWGCNA 和 workflow 依赖，
- `samtools`、`tabix`、`bgzip`、`bedtools`、`macs2`、Snakemake 等命令行工具，
- MALLET 2.0.8 topic-model backend，
- `spgrn-*` workflow 命令，用于项目初始化、资源准备、pycisTopic、custom cisTarget、SCENIC+ 运行、QC 和结果输出。


## 验证安装

激活环境：

```bash
source /absolute/path/to/conda/bin/activate scenicplus-grn
```

运行检查：

```bash
spgrn-check
spgrn-check-workflow-installation
```

检查内容包括核心 Python imports、R packages、命令行工具、MALLET、workflow
脚本和已安装文档。

## 离线或网络受限安装

Git 仓库不跟踪 bundled source archive，因为这是体积较大的 release 归档文件。
`archives/vendor.tar.gz` 只随 GitHub Release package 分发，不会出现在
普通源码 checkout 或 GitHub "Code" 下载包里。

release package 包含：

```text
archives/vendor.tar.gz
```

`install.sh` 会在运行时解压该归档；GitHub 重试失败后使用其中的本地源码。
如需完全跳过 GitHub：

```bash
GITHUB_TRIES=0 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

## 许可证

本 workflow helper 和安装脚本使用 MIT License。release package 中的
`archives/vendor.tar.gz` 只是用于离线和可复现安装的第三方源码缓存；第三方组件仍遵循各自上游许可证或使用条款。尤其是 SCENIC+、pycisTopic 和 pycistarget 使用 academic non-commercial 许可，不是通用商业开源授权。
