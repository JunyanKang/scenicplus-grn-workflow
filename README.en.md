# scenicplus-grn-workflow

Reproducible workflow toolkit for SCENIC+ gene regulatory network analysis from annotated matched snRNA+snATAC / scMultiome data. It creates an isolated `scenicplus-grn` conda environment, installs the pinned Python, R, genomics and SCENIC+ layers, and provides `spgrn-*` commands for project initialization, resource preparation, pycisTopic, custom cisTarget database construction, SCENIC+ Snakemake execution and postprocessing.

Install from the complete `scenicplus-grn-workflow-v*.tar.gz` package on [GitHub Releases](https://github.com/JunyanKang/scenicplus-grn-workflow/releases). That package includes `archives/vendor.tar.gz` for servers and restricted-network environments. The green GitHub `Code` ZIP contains source only and does not include the offline vendor archive, so installation from that ZIP requires GitHub access for third-party source downloads.

## What This Is For

Use this package when you need a portable SCENIC+ analysis environment on a
workstation or Linux server and want the workflow scripts, version records and
offline source archives shipped together.

It is intended for:

- annotated scMultiome objects plus matched ATAC fragments,
- metacell-based SCENIC+ runs,
- custom cisTarget databases built from the project consensus region universe,
- direct, orthology-mapped, audited generated or user-supplied species-specific motif2TF table preparation and audit,
- reproducible reruns and output checks.

Detailed analysis instructions: [SCENICPLUS_STEP_BY_STEP.en.md](https://github.com/JunyanKang/scenicplus-grn-workflow/blob/main/docs/SCENICPLUS_STEP_BY_STEP.en.md).

The GitHub-default Chinese README and workflow guide are:

[README.md](https://github.com/JunyanKang/scenicplus-grn-workflow/blob/main/README.md) and [SCENICPLUS_STEP_BY_STEP.md](https://github.com/JunyanKang/scenicplus-grn-workflow/blob/main/docs/SCENICPLUS_STEP_BY_STEP.md).

## Supported Platforms

```text
Linux x86_64      glibc >= 2.17
macOS arm64       Apple Silicon
```

The workflow installer expects an existing conda-style root such as Miniforge,
Miniconda, Mambaforge, Anaconda or `/opt/conda`.

## Quick Start

Unpack the release archive on the target machine:

```bash
tar -xzf scenicplus-grn-workflow-v*.tar.gz
cd scenicplus-grn-workflow
```

Run the installation:

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
cd scenicplus-grn-workflow
ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

Detach with `Ctrl-b` then `d`, and reattach later:

```bash
tmux attach -t scenicplus-install
```

Without `tmux`:

```bash
cd scenicplus-grn-workflow
mkdir -p logs
nohup env ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh \
  > logs/nohup_install_$(date +%Y%m%d_%H%M%S).out 2>&1 &
echo $! > logs/install.pid
```

Follow the newest installation log:

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
| `RELOCATE_INSTALLER` | `1` | Offer to copy the workflow package to `$CONDA_ROOT/share/scenicplus-grn-workflow`. |
| `LOG_DIR` | `logs/` | Installer log directory. |

Examples:

```bash
ENV_NAME=scenicplus-grn-test CONDA_ROOT=/absolute/path/to/conda bash install.sh
GITHUB_TRIES=0 ASSUME_YES=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
FORCE=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
PRECHECK_ONLY=1 CONDA_ROOT=/absolute/path/to/conda bash install.sh
```

## What Gets Installed

- SCENIC+, pycisTopic, pycistarget and create_cisTarget_databases,
- Python single-cell and genomics libraries,
- R, Seurat, Signac, hdWGCNA and workflow dependencies,
- command-line tools including `samtools`, `tabix`, `bgzip`, `bedtools`,
  `macs2` and Snakemake,
- MALLET 2.0.8 topic-model backend,
- `spgrn-*` workflow commands for project initialization, resource preparation,
  pycisTopic, custom cisTarget, SCENIC+ execution, QC and output generation.


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

## License

The workflow helper and installation scripts are licensed under the MIT License. The
release package's `archives/vendor.tar.gz` is only an offline cache of
third-party source archives for reproducible installation; those components
remain governed by their upstream licenses or terms. In particular, SCENIC+,
pycisTopic and pycistarget use academic non-commercial terms rather than a
general commercial open-source grant.
