# SCENIC+ GRN Conda Installer

This package installs a dedicated conda environment for SCENIC+/GRN analysis. It is designed to be copied to a workstation or Linux server as a small archive, unpacked anywhere, and then run with `bash install.sh`.

## Packaged Contents And Roles

```text
Installer entry points:
  install.sh                         Required. Main installer and bootstrap script.
  README.md                          User-facing installer guide.
  README.zh-CN.md                    Chinese installer guide.
  VERSION                            Installer package version.

Installer helper entry points:
  bin/check_environment.sh           Required. Environment and workflow self-check entry point.
  bin/initialize_scenicplus_project.sh
                                      User-facing one-step initializer after installation.
  bin/install_r.R                    Required R/hdWGCNA/Seurat/Signac installation layer for the metacell workflow.

Environment recipes and pinned Python layer:
  config/environment-linux-64.yml    Required on Linux x86_64, compatible with glibc >= 2.17.
  config/environment-macos-arm64.yml Required on Apple Silicon macOS.
  config/pip-constraints.txt         Required by the pinned pip supplement.
  config/scenicplus_config_template.yaml
                                      Required Snakemake config template used by workflow generators.
  config/locks/environment-linux-64.solved-lock.yml
                                      Linux dry-run conda lock record with build strings for auditing/debugging.

Offline/restricted-network source archives:
  archives/vendor.tar.gz             Required in release archives for robust installation when GitHub is slow or blocked.
                                      install.sh expands it to hidden runtime cache .vendor/.
  .vendor/github/                    Runtime-expanded pinned source archives used after GitHub retries fail.
  .vendor/mallet/                    Runtime-expanded MALLET 2.0.8 archive used by the default MALLET installation.

Installed workflow assets:
  scripts/                           Required executable workflow entry points after installation.
                                      Installed to $CONDA_PREFIX/share/scenicplus-grn/scripts/.
                                      Includes project parameter setup, raw-data sample-sheet generation,
                                      pycisTopic, cisTarget DB, SCENIC+ config, Snakemake and postprocessing wrappers.
  modules/                           Required internal helper modules imported by scripts; not user-facing commands.
                                      Installed to $CONDA_PREFIX/share/scenicplus-grn/modules/.

Documentation and audit records:
  docs/SCENICPLUS_STEP_BY_STEP.md    Strict matched snRNA+snATAC workflow guide.
  docs/SCENICPLUS_STEP_BY_STEP.zh-CN.md
                                      Chinese matched snRNA+snATAC workflow guide.
  docs/CHANGELOG.md                  Release history.
  docs/RELEASE_NOTES.md              Current release notes.
  docs/VERSION_LOCK.md               Human-readable record of pinned analysis software versions.
```

## Recommended Use

Copy the archive to the target machine, unpack it, and run the installer:

```bash
tar -xzf scenicplus-grn-installer.tar.gz
cd scenicplus-grn-installer
bash install.sh
```

The script first looks for a conda-style root. You can make this explicit:

```bash
CONDA_ROOT=/path/to/conda bash install.sh
```

If the unpacked installer is not already under the detected conda root, it will offer to copy itself to:

```text
$CONDA_ROOT/share/scenicplus-grn-installer
```

and continue from there. This keeps the reusable installer next to the conda installation rather than inside a project directory.

## What Happens When You Run It

`install.sh` performs these steps:

```text
1. Starts a timestamped log.
2. Detects a conda/miniforge/miniconda/mamba-style root.
3. Asks you to confirm the detected root.
4. If needed, offers to copy the installer to $CONDA_ROOT/share/scenicplus-grn-installer.
5. Checks write permission for $CONDA_ROOT, envs, pkgs, and share.
6. Bootstraps mamba into conda base if it is not already available.
7. Creates or updates a dedicated conda environment named scenicplus-grn.
8. Installs the conda/mamba-resolved Python, CLI, genomics, and R base layer.
9. Expands `archives/vendor.tar.gz` to `.vendor/` when needed, then installs a small pinned pip supplement plus the pinned SCENIC+/pycisTopic/pycistarget source layer. Each GitHub source is tried 3 times, then installed from the bundled `.vendor/github` archive if GitHub is unstable.
10. Installs MALLET 2.0.8 by default for the pycisTopic MALLET LDA backend, using a Python wrapper around the Java classes and a real `import-file` smoke test.
11. Installs the R/hdWGCNA layer required for the metacell workflow, using the bundled pinned source archive when GitHub is unavailable.
12. Runs check_environment.sh.
13. Copies the installer recipe into $CONDA_PREFIX/share/scenicplus-grn.
14. Runs the workflow asset checker from the installed environment copy.
```

By default, it creates or updates an independent environment named:

```text
scenicplus-grn
```

It does not install into your currently active environment unless you explicitly use `MODE=active`.

## Logs

A log is written for every run.

Default location:

```text
scenicplus-grn-installer/logs/install_YYYYMMDD_HHMMSS.log
```

If the unpacked directory is not writable, logs are written to:

```text
/tmp/scenicplus-grn-installer-logs
```

You can choose a log directory:

```bash
LOG_DIR=/path/to/logs bash install.sh
```

## Background Installation

The full installation can take a long time on a shared server. You do not need
to keep watching the terminal as long as the log is being written.

Recommended, if `tmux` is available:

```bash
tmux new -s scenicplus-install
cd scenicplus-grn-installer
ASSUME_YES=1 CONDA_ROOT=/path/to/miniconda3 bash install.sh
```

Detach from the session with `Ctrl-b` then `d`. Reattach later:

```bash
tmux attach -t scenicplus-install
```

Simple background mode without `tmux`:

```bash
cd scenicplus-grn-installer
mkdir -p logs
nohup env ASSUME_YES=1 CONDA_ROOT=/path/to/miniconda3 bash install.sh \
  > logs/nohup_install_$(date +%Y%m%d_%H%M%S).out 2>&1 &
echo $! > logs/install.pid
```

Check whether it is still running:

```bash
ps -p "$(cat logs/install.pid)"
```

Follow the newest installer log:

```bash
tail -f "$(ls -t logs/install_*.log | head -n 1)"
```

A successful run ends with:

```text
DONE: SCENIC+ environment is installed and checked.
```

If it fails, send the newest `logs/install_*.log` file for diagnosis.

## Main Commands

On Linux, `install.sh` automatically bootstraps `mamba` into conda `base` if
`$CONDA_ROOT/bin/mamba` is not already present. The bootstrap step uses the
existing conda solver for maximum compatibility with older Miniconda
installations; after mamba is available, the script uses mamba for the analysis
environment creation step.
This avoids the classic conda solver getting stuck in long conflict analysis.

To disable automatic mamba bootstrapping:

```bash
AUTO_INSTALL_MAMBA=0 bash install.sh
```

Create or update the dedicated environment:

```bash
bash install.sh
```

Use a specific conda root:

```bash
CONDA_ROOT=/path/to/miniforge3 bash install.sh
```

Create a differently named environment:

```bash
ENV_NAME=scenicplus-grn-test bash install.sh
```

Change how many times GitHub is tried before using bundled source archives:

```bash
GITHUB_TRIES=1 bash install.sh
```

Skip GitHub entirely and install the pinned source layer from bundled archives:

```bash
GITHUB_TRIES=0 bash install.sh
```

Rebuild the environment from scratch:

```bash
FORCE=1 bash install.sh
```

Skip MALLET installation:

```bash
INSTALL_MALLET=0 bash install.sh
```

MALLET is installed by default and is recommended for large pycisTopic topic
models. The installer verifies the wrapper with a small `import-file` run. If
MALLET is disabled, the workflow remains usable with `pycistopic.lda_backend=cgs`.

Run non-interactively with detected defaults:

```bash
ASSUME_YES=1 bash install.sh
```

Check detection, relocation, permissions, and platform recipe selection without installing:

```bash
PRECHECK_ONLY=1 bash install.sh
```

Do not copy the installer under the conda root:

```bash
RELOCATE_INSTALLER=0 bash install.sh
```

Install into the currently active conda environment instead:

```bash
source /path/to/miniforge3/bin/activate my-existing-env
MODE=active bash install.sh
```

The installer refuses to install into `base` unless you explicitly allow it:

```bash
ALLOW_BASE=1 MODE=active bash install.sh
```

## After Installation

Enter the environment and project parameters in the terminal:

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

`PROJECT_DIR` should be the dedicated SCENIC+ analysis root, not a larger study
repository unless all workflow-created `inputs/`, `work/`, `resources/`, `logs/`
and `results/` folders are intended to live there.

Run the one-step initializer. It checks the environment first, updates
`$PROJECT_DIR/scenicplus_project.env`, and initializes the project runtime files:

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
  bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/bin/initialize_scenicplus_project.sh"
```

After the initializer reports `PROJECT INITIALIZATION OK`, load the project
runtime variables:

```bash
source "$PROJECT_DIR/project_env.sh"
```

The project settings file is written to `$PROJECT_DIR/scenicplus_project.env`.
Workflow scripts detect available memory, current load and CPU count at launch
time to choose conservative parallel worker counts.

Expected core versions from the current recipe:

```text
pycisTopic 1.0.2
pycistarget 1.1
scenicplus 1.0a2
MACS2 2.2.9.1
Snakemake 7.32.4
MALLET 2.0.8
```

## Platform Detection

`install.sh` chooses the environment file automatically:

```text
macOS arm64      -> config/environment-macos-arm64.yml
Linux x86_64     -> config/environment-linux-64.yml
```

For Linux, the script exports `CONDA_SUBDIR=linux-64` if it is not already set. This recipe targets ordinary Linux x86_64 servers and uses glibc 2.17 as the minimum compatibility line. In practice, it should also work on newer glibc systems such as Rocky/Alma/CentOS Stream, Ubuntu, and Debian. Servers older than glibc 2.17 are not supported by many current conda-forge/bioconda packages. Linux ARM/aarch64 would need a separate recipe.

## What Gets Installed

The installer is version-locked to the successful local dry run on macOS. The
Linux recipe is conda-first: Python, CLI, genomics, SCENIC-adjacent packages and
R base are resolved by mamba from conda-forge/bioconda using a CentOS7/glibc
2.17-compatible lock. `config/environment-linux-64.yml` pins the top-level install
versions. `config/locks/environment-linux-64.solved-lock.yml` records the full dry-run
solve, including transitive dependency build strings, for auditing and exact
troubleshooting. Pip is used only for a small pinned supplement and exact GitHub
source commits. The SCENIC+/pycisTopic/pycistarget source commits remain
identical across platforms. Release archives also include exact source archives
as `archives/vendor.tar.gz`; `install.sh` expands this archive to `.vendor/`
before installation so that GitHub failures on restricted networks do not block
installation. For local bundled-source installation, the installer explicitly passes
the pinned package versions to setuptools-scm, so GitHub archive tarballs do not
need `.git` metadata to build reproducibly. See `docs/VERSION_LOCK.md` for the full
readable version summary.

Workflow scripts and helper modules are installed inside the independent conda
environment, not into a global system directory. After activation,
`$CONDA_PREFIX` points to the environment, so the scripts are available at:

```text
$CONDA_PREFIX/share/scenicplus-grn/scripts/
```

For convenience, the step-by-step guide defines:

```bash
export SCENICPLUS_HOME="$CONDA_PREFIX/share/scenicplus-grn"
```

and then calls scripts as `$SCENICPLUS_HOME/scripts/<script>.py`.

Core Python, command-line, and genomics layer:

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

AutoZyme is treated as a no-dependency overlay. It must not upgrade, downgrade,
or replace the pinned Scanpy, Seurat, SCENIC+, pycisTopic, or R package stack.
The installer defaults to `INSTALL_AUTOZYME=1` and `INSTALL_AUTOZYME_R=1`. Runtime
activation can be disabled for analysis scripts with `AUTOZYME_DISABLED=1`.
AutoZyme may print version-scope warnings when a patch was lifted against a
nearby upstream version; these warnings are logged and do not change package
versions.

Bundled GitHub source archives used after GitHub retries fail are stored inside
`archives/vendor.tar.gz` and expanded at runtime to:

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

Bundled archives are checksum-validated before use. If an uploaded archive is
corrupted, the installer stops with a checksum error instead of continuing with
a partial source tree.

R, Bioconductor, and hdWGCNA layer:

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

On Linux, `cmake`, `pkg-config`, `libuv`, `xz`, `liblzma-devel`, Seurat,
Signac, WGCNA, R igraph, and the heavy Seurat dependency stack are installed
through conda before the R source layer runs. This avoids failures on servers
whose system CMake is too old for CRAN packages such as `fs`, and avoids
source-compiling R igraph/Seurat on shared servers. R igraph is intentionally
held at the R 4.4-compatible `2.1.4` build to preserve the Python `lxml 5.3.1`
/ `libxml2 2.13` stack required by the SCENIC+ layer.

Pinned R source dependencies used by hdWGCNA/plotting/network utilities:

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

Additional pinned pip dependencies are listed in `config/pip-constraints.txt`; the installer force-reinstalls that pip layer to reduce future version drift.

`bin/check_environment.sh` sets `options(enrichR.live = FALSE)` before loading hdWGCNA. This avoids false failures when the Enrichr web service is slow, blocked, or unavailable on a server. It does not disable SCENIC+ itself; it only prevents a package attach-time network check from breaking the environment validation.

## Troubleshooting

If conda is not detected:

```bash
CONDA_ROOT=/path/to/miniforge3 bash install.sh
```

If you unpacked the archive somewhere temporary and skipped relocation, you can manually place it under conda later:

```bash
mkdir -p ~/miniforge3/share
cp -a scenicplus-grn-installer ~/miniforge3/share/
cd ~/miniforge3/share/scenicplus-grn-installer
bash install.sh
```

To disable AutoZyme installation while keeping the rest of the environment:

```bash
INSTALL_AUTOZYME=0 INSTALL_AUTOZYME_R=0 bash install.sh
```

To keep AutoZyme installed but disable runtime patch activation in an analysis:

```bash
AUTOZYME_DISABLED=1 python your_analysis.py
AUTOZYME_DISABLED=1 Rscript your_analysis.R
```

If you only need to validate an existing installation:

```bash
CONDA_ROOT=/path/to/conda
ENV_NAME=scenicplus-grn
bash "$CONDA_ROOT/envs/$ENV_NAME/share/scenicplus-grn/bin/check_environment.sh" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME"
```
