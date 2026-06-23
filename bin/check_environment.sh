#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-scenicplus-grn}"
CONDA_ROOT="${CONDA_ROOT:-}"
CONFIG_FILE=""
SKIP_R=0
SKIP_WORKFLOW_ASSETS=0

CACHE_ROOT="${SCENICPLUS_CHECK_CACHE_ROOT:-${TMPDIR:-/tmp}/scenicplus-grn-check-${ENV_NAME}}"
mkdir -p "$CACHE_ROOT/numba" "$CACHE_ROOT/matplotlib"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$CACHE_ROOT/numba}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$CACHE_ROOT/matplotlib}"

usage() {
  cat <<'EOF'
Usage:
  bash check_environment.sh --conda-root /path/to/miniforge3
  bash check_environment.sh --config scenicplus_project.env

Options:
  --config PATH              Settings file with CONDA_ROOT and ENV_NAME.
  --conda-root PATH          Conda/miniforge/miniconda/mamba/anaconda root.
  --env-name NAME            Environment name. Default: scenicplus-grn.
  --skip-r                   Skip R package smoke tests.
  --skip-workflow-assets     Skip installed workflow asset checks.
  -h, --help                 Show this message.
EOF
}

load_config() {
  local config="$1"
  if [[ ! -f "$config" ]]; then
    echo "ERROR: config file not found: $config" >&2
    exit 1
  fi
  set -a
  # shellcheck source=/dev/null
  source "$config"
  set +a
}

expand_path() {
  local path="$1"
  if [[ "$path" == "~" ]]; then
    printf "%s\n" "$HOME"
  elif [[ "$path" == ~/* ]]; then
    printf "%s/%s\n" "$HOME" "${path#~/}"
  else
    printf "%s\n" "$path"
  fi
}

infer_conda_root_from_prefix() {
  local prefix="${CONDA_PREFIX:-}"
  if [[ -n "$prefix" && -d "$prefix" ]]; then
    local parent grandparent
    parent="$(dirname "$prefix")"
    grandparent="$(dirname "$parent")"
    if [[ "$(basename "$parent")" == "envs" && -x "$grandparent/bin/conda" ]]; then
      printf "%s\n" "$grandparent"
      return 0
    fi
  fi
  return 1
}

infer_conda_root_from_script() {
  local script_path script_dir p
  script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
  script_dir="$(dirname "$script_path")"

  if [[ "$(basename "$script_dir")" == "bin" ]]; then
    p="$(dirname "$script_dir")"
    if [[ -x "$p/bin/conda" ]]; then
      printf "%s\n" "$p"
      return 0
    fi
  fi

  p="$script_dir"
  while [[ "$p" != "/" ]]; do
    if [[ "$(basename "$p")" == "envs" ]]; then
      p="$(dirname "$p")"
      if [[ -x "$p/bin/conda" ]]; then
        printf "%s\n" "$p"
        return 0
      fi
      break
    fi
    p="$(dirname "$p")"
  done
  return 1
}

detect_conda_root() {
  if [[ -n "$CONDA_ROOT" ]]; then
    expand_path "$CONDA_ROOT"
    return
  fi
  if infer_conda_root_from_prefix; then
    return
  fi
  if infer_conda_root_from_script; then
    return
  fi

  local candidates=(
    "$HOME/miniforge3"
    "$HOME/miniconda3"
    "$HOME/mambaforge"
    "$HOME/anaconda3"
    "/opt/conda"
  )
  local c
  for c in "${candidates[@]}"; do
    if [[ -x "$c/bin/conda" ]]; then
      printf "%s\n" "$c"
      return
    fi
  done
  echo "ERROR: could not find conda root. Pass --conda-root /path/to/miniforge3." >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --conda-root)
      CONDA_ROOT="$2"
      shift 2
      ;;
    --env-name)
      ENV_NAME="$2"
      shift 2
      ;;
    --skip-r)
      SKIP_R=1
      shift
      ;;
    --skip-workflow-assets)
      SKIP_WORKFLOW_ASSETS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -n "$CONFIG_FILE" ]]; then
  load_config "$CONFIG_FILE"
fi
if [[ "${CHECK_R:-1}" == "0" ]]; then
  SKIP_R=1
fi

CONDA_ROOT="$(detect_conda_root)"
CONDA_BIN="$CONDA_ROOT/bin/conda"
ENV_PREFIX="$CONDA_ROOT/envs/$ENV_NAME"

if [[ ! -x "$CONDA_BIN" ]]; then
  echo "ERROR: conda executable not found: $CONDA_BIN" >&2
  exit 1
fi
if [[ ! -d "$ENV_PREFIX" ]]; then
  echo "ERROR: environment not found: $ENV_PREFIX" >&2
  echo "Run the installer first, or pass --env-name if you used a different environment name." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$CONDA_ROOT/etc/profile.d/conda.sh"
set +u
conda activate "$ENV_PREFIX"
set -u

export SCENICPLUS_HOME="$CONDA_PREFIX/share/scenicplus-grn"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/scenicplus-grn-mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/scenicplus-grn-xdg-cache}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${TMPDIR:-/tmp}/scenicplus-grn-numba-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$NUMBA_CACHE_DIR"

echo "Conda root: $CONDA_ROOT"
echo "Environment: $CONDA_PREFIX"
echo "SCENICPLUS_HOME: $SCENICPLUS_HOME"
echo

python - <<'PY'
import os
import flatbuffers
import pycisTopic, pycistarget, scenicplus, scanpy, mudata, pandas
import pycisTopic.cistopic_class, pycisTopic.pseudobulk_peak_calling, pycisTopic.lda_models
import pycistarget.motif_enrichment_cistarget
import scenicplus.scenicplus_class, scenicplus.enhancer_to_gene, scenicplus.TF_to_gene
print("OK python core")
print("flatbuffers", getattr(flatbuffers, "__version__", "unknown"))
print("pycisTopic", pycisTopic.__version__)
print("pycistarget", pycistarget.__version__)
print("scenicplus", scenicplus.__version__)
print("scanpy", scanpy.__version__)
print("mudata", mudata.__version__)
print("pandas", pandas.__version__)
if os.environ.get("INSTALL_AUTOZYME", "1") == "1":
    import autozyme
    print("autozyme", getattr(autozyme, "__version__", "unknown"))
    autozyme.activate("scanpy")
    print("autozyme_scanpy", autozyme.status().get("scanpy"))
PY

macs2 --version
snakemake --version
scenicplus --help >/dev/null
pycistarget --help >/dev/null
python "$CONDA_PREFIX/opt/create_cisTarget_databases/create_cistarget_motif_databases.py" --help >/dev/null
cbust -h >/dev/null || true
echo "OK command-line tools"

if [[ "$SKIP_R" == "1" ]]; then
  echo "Skipping R package checks."
else
  Rscript - <<'RS'
options(enrichR.live = FALSE)
packages <- c(
  "Seurat", "Signac", "WGCNA", "hdWGCNA", "GenomicRanges",
  "GeneOverlap", "UCell", "impute", "preprocessCore"
)
if (identical(Sys.getenv("INSTALL_AUTOZYME_R", unset = "1"), "1")) {
  packages <- c(packages, "autozyme")
}
for (pkg in packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(pkg, as.character(packageVersion(pkg)), "\n")
}
cat(
  "MetacellsByGroups",
  exists("MetacellsByGroups", where = asNamespace("hdWGCNA"), inherits = FALSE),
  "\n"
)
if (identical(Sys.getenv("INSTALL_AUTOZYME_R", unset = "1"), "1")) {
  autozyme::activate("seurat")
  cat("autozyme_seurat", "checked", "\n")
}
RS
  echo "OK R packages"
fi

if [[ "$SKIP_WORKFLOW_ASSETS" == "1" ]]; then
  echo "Skipping workflow asset checks."
else
  if [[ ! -d "$SCENICPLUS_HOME" ]]; then
    echo "ERROR: installed workflow assets not found: $SCENICPLUS_HOME" >&2
    echo "Run install.sh again so scripts, docs, and config template are copied into the environment." >&2
    exit 1
  fi
  python "$SCENICPLUS_HOME/scripts/check_workflow_installation.py"
fi
