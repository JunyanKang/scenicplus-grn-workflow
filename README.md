# scenicplus-grn-installer

Reproducible conda installer and workflow launcher for SCENIC+ gene regulatory
network analysis from annotated matched snRNA+snATAC / scMultiome data.

The installer creates an isolated `scenicplus-grn` conda environment, installs
the pinned Python, R, genomics and SCENIC+ layers, and provides `spgrn-*`
commands for project initialization, resource preparation, pycisTopic,
custom cisTarget database construction, SCENIC+ Snakemake execution and
postprocessing.

## What This Is For

Use this package when you need a portable SCENIC+ analysis environment on a
workstation or Linux server and want the workflow scripts, version records and
offline source archives shipped together.

It is intended for:

- annotated scMultiome objects plus matched ATAC fragments,
- metacell-based SCENIC+ runs,
- custom cisTarget databases built from the project consensus region universe,
- reproducible reruns and output checks.

Detailed analysis instructions live in:

```text
docs/SCENICPLUS_STEP_BY_STEP.md
docs/SCENICPLUS_STEP_BY_STEP.zh-CN.md
```

## Supported Platforms

```text
Linux x86_64      glibc >= 2.17
macOS arm64       Apple Silicon
```

The installer expects an existing conda-style root such as Miniforge,
Miniconda, Mambaforge, Anaconda or `/opt/conda`.

## Quick Start

Unpack the release archive on the target machine:

```bash
tar -xzf scenicplus-grn-installer.tar.gz
cd scenicplus-grn-installer
```

Run the installer:

```bash
CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

For unattended installation:

```bash
ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

A successful run ends with:

```text
DONE: SCENIC+ environment is installed and checked.
```

## Background Installation

Long server installs should run in `tmux` or with `nohup`.

With `tmux`:

```bash
tmux new -s scenicplus-install
cd scenicplus-grn-installer
ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

Detach with `Ctrl-b` then `d`, and reattach later:

```bash
tmux attach -t scenicplus-install
```

Without `tmux`:

```bash
cd scenicplus-grn-installer
mkdir -p logs
nohup env ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh \
  > logs/nohup_install_$(date +%Y%m%d_%H%M%S).out 2>&1 &
echo $! > logs/install.pid
```

Follow the newest installer log:

```bash
tail -f "$(ls -t logs/install_*.log | head -n 1)"
```

## Installation Options

| Variable | Default | Meaning |
|---|---:|---|
| `CONDA_ROOT` | auto-detect | Conda root to use. |
| `ENV_NAME` | `scenicplus-grn` | Environment name to create or update. |
| `ASSUME_YES` | `0` | Set to `1` for non-interactive prompts. |
| `FORCE` | `0` | Set to `1` to recreate the environment. |
| `GITHUB_TRIES` | `3` | GitHub attempts before bundled source archives are used. |
| `INSTALL_MALLET` | `1` | Install MALLET backend for pycisTopic topic modeling. |
| `AUTO_INSTALL_MAMBA` | `1` | Bootstrap `mamba` into conda base if missing. |
| `PRECHECK_ONLY` | `0` | Check platform, paths and permissions without installing. |
| `RELOCATE_INSTALLER` | `1` | Offer to copy the installer to `$CONDA_ROOT/share/scenicplus-grn-installer`. |
| `LOG_DIR` | `logs/` | Installer log directory. |

Examples:

```bash
ENV_NAME=scenicplus-grn-test CONDA_ROOT=/absolute/path/to/conda bash install.sh
GITHUB_TRIES=0 ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
FORCE=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
PRECHECK_ONLY=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

## What Gets Installed

The environment contains:

- SCENIC+, pycisTopic, pycistarget and create_cisTarget_databases,
- Python single-cell and genomics libraries,
- R, Seurat, Signac, hdWGCNA and workflow dependencies,
- command-line tools including `samtools`, `tabix`, `bgzip`, `bedtools`,
  `macs2` and Snakemake,
- MALLET 2.0.8 by default,
- workflow scripts installed under:

```text
$CONDA_PREFIX/share/scenicplus-grn
```

Short commands are installed into:

```text
$CONDA_PREFIX/bin/spgrn-*
```

Version details are recorded in:

```text
docs/VERSION_LOCK.md
docs/RELEASE_NOTES.md
docs/CHANGELOG.md
```

## Verify The Installation

Activate the environment:

```bash
source /absolute/path/to/conda/bin/activate scenicplus-grn
```

Run the installed checks:

```bash
spgrn-check
spgrn-check-workflow-installation
```

The checks verify core Python imports, R packages, command-line tools, MALLET,
workflow scripts and installed documentation.

## Project Initialization

Create or update a project using the installed initializer:

```bash
# Replace these example values before running.
# CONDA_ROOT must be the real conda/miniforge/miniconda root.
export CONDA_ROOT=/absolute/path/to/conda
export ENV_NAME=scenicplus-grn
# PROJECT_DIR must be a dedicated SCENIC+ analysis directory.
export PROJECT_DIR=/absolute/path/to/grn_project
# ORGANISM must be one supported organism key.
export ORGANISM=mouse
export AUTOZYME=on
export ENSEMBL_RELEASE=115
# ANNOTATED_OBJECT must be the active annotated object for this analysis.
export ANNOTATED_OBJECT=/absolute/path/to/annotated_multiome_object.rds
# CELL_LABEL_COLUMN must be a real metadata column in ANNOTATED_OBJECT.
export CELL_LABEL_COLUMN=cell_annotation
# ATAC_INPUT_LAYOUT must match the real ATAC_DATA_ROOT structure.
export ATAC_INPUT_LAYOUT=split_ge_arc
export ATAC_DATA_ROOT=/absolute/path/to/raw_atac_data
spgrn-initialize
```

Then follow the step-by-step guide for resource preparation, pycisTopic,
custom cisTarget database construction, SCENIC+ inference and postprocessing:

```text
$CONDA_PREFIX/share/scenicplus-grn/docs/SCENICPLUS_STEP_BY_STEP.md
```

## Offline And Restricted-Network Installs

The Git repository does not track the bundled source archive because it is a
large binary release artifact. `archives/vendor.tar.gz` is distributed with the
GitHub Release package, not with a normal source checkout or GitHub "Code"
download.

Release packages include:

```text
archives/vendor.tar.gz
```

`install.sh` expands this archive at runtime and uses the bundled sources when
GitHub retries fail. To skip GitHub entirely:

```bash
GITHUB_TRIES=0 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

## Logs And Troubleshooting

Installer logs are written to:

```text
logs/install_YYYYMMDD_HHMMSS.log
```

If the unpacked directory is not writable, logs are written to:

```text
/tmp/scenicplus-grn-installer-logs
```

Useful checks:

```bash
tail -n 80 logs/install_*.log
CONDA_ROOT=/absolute/path/to/conda PRECHECK_ONLY=1 bash install.sh
source /absolute/path/to/conda/bin/activate scenicplus-grn
spgrn-check
```

For a failed installation, keep the newest installer log and the command that
was used to start the run.

## Repository Layout

Source checkout:

```text
install.sh                 Main installer.
bin/                       Installer bootstrap checks and R layer installer.
config/                    Conda recipes, pinned pip constraints and templates.
scripts/                   Installed workflow command implementations.
modules/                   Internal helper modules.
docs/                      Step-by-step guides, changelog and version records.
```

Release package only:

```text
archives/vendor.tar.gz     Bundled offline source archives.
```

For offline or restricted-network installation, download the GitHub Release
package or release asset rather than the source-only repository archive.
