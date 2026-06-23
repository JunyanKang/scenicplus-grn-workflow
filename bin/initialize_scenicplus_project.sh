#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  CONDA_ROOT=/absolute/path/to/conda ENV_NAME=scenicplus-grn PROJECT_DIR=/absolute/path/to/project \
    ORGANISM=mouse AUTOZYME=on ENSEMBL_RELEASE=115 \
    ANNOTATED_OBJECT=/absolute/path/to/active_object.rds CELL_LABEL_COLUMN=cell_annotation \
    ATAC_INPUT_LAYOUT=split_ge_arc ATAC_DATA_ROOT=/absolute/path/to/atac_input_root \
    bash initialize_scenicplus_project.sh

This script validates the installed SCENIC+ environment, creates the project
directory layout, writes PROJECT_DIR/project_env.sh, and records the input
parameters in PROJECT_DIR/scenicplus_project.env.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

require_var() {
  local name="$1"
  local value="${!name:-}"
  if [[ -z "$value" ]]; then
    echo "ERROR: $name is required." >&2
    usage >&2
    exit 2
  fi
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

shell_quote() {
  local value="$1"
  printf "'%s'" "${value//\'/\'\"\'\"\'}"
}

update_config_setting() {
  local key="$1"
  local value="$2"
  local quoted
  quoted="$(shell_quote "$value")"
  if [[ -f "$CONFIG_FILE" ]]; then
    local tmp
    tmp="$(mktemp "$CONFIG_FILE.tmp.XXXXXX")"
    awk -v key="$key" -v value="$quoted" '
      BEGIN { replaced = 0 }
      $0 ~ "^" key "=" {
        print key "=" value
        replaced = 1
        next
      }
      { print }
      END {
        if (replaced == 0) {
          print key "=" value
        }
      }
    ' "$CONFIG_FILE" > "$tmp"
    mv "$tmp" "$CONFIG_FILE"
  else
    printf "%s=%s\n" "$key" "$quoted" >> "$CONFIG_FILE"
  fi
}

require_var CONDA_ROOT
require_var ENV_NAME
require_var PROJECT_DIR
require_var ORGANISM
require_var ANNOTATED_OBJECT
require_var CELL_LABEL_COLUMN
require_var ATAC_INPUT_LAYOUT
require_var ATAC_DATA_ROOT

AUTOZYME="${AUTOZYME:-on}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-115}"
SCENICPLUS_MAX_MEMORY_GB="${SCENICPLUS_MAX_MEMORY_GB:-auto}"

CONDA_ROOT="$(expand_path "$CONDA_ROOT")"
PROJECT_DIR="$(expand_path "$PROJECT_DIR")"
ANNOTATED_OBJECT="$(expand_path "$ANNOTATED_OBJECT")"
ATAC_DATA_ROOT="$(expand_path "$ATAC_DATA_ROOT")"

case "$ORGANISM" in
  human|mouse|cyno|rat|rabbit|chicken|zebrafish) ;;
  *)
    echo "ERROR: unsupported ORGANISM=$ORGANISM" >&2
    echo "Allowed: human mouse cyno rat rabbit chicken zebrafish" >&2
    exit 2
    ;;
esac
case "$AUTOZYME" in
  on|off) ;;
  *)
    echo "ERROR: AUTOZYME must be on or off." >&2
    exit 2
    ;;
esac
case "$ATAC_INPUT_LAYOUT" in
  split_ge_arc|cellranger_outs) ;;
  *)
    echo "ERROR: ATAC_INPUT_LAYOUT must be split_ge_arc or cellranger_outs." >&2
    exit 2
    ;;
esac
if [[ ! -f "$ANNOTATED_OBJECT" ]]; then
  echo "ERROR: ANNOTATED_OBJECT not found: $ANNOTATED_OBJECT" >&2
  exit 1
fi
if [[ ! -d "$ATAC_DATA_ROOT" ]]; then
  echo "ERROR: ATAC_DATA_ROOT not found: $ATAC_DATA_ROOT" >&2
  exit 1
fi

mkdir -p "$PROJECT_DIR"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd -P)"
ENV_PREFIX="$CONDA_ROOT/envs/$ENV_NAME"
SCENICPLUS_HOME="$ENV_PREFIX/share/scenicplus-grn"
CONFIG_FILE="$PROJECT_DIR/scenicplus_project.env"

if [[ ! -x "$CONDA_ROOT/bin/conda" ]]; then
  echo "ERROR: conda executable not found: $CONDA_ROOT/bin/conda" >&2
  exit 1
fi
if [[ ! -d "$ENV_PREFIX" ]]; then
  echo "ERROR: conda environment not found: $ENV_PREFIX" >&2
  exit 1
fi
if [[ ! -f "$SCENICPLUS_HOME/bin/check_environment.sh" ]]; then
  echo "ERROR: installed checker not found: $SCENICPLUS_HOME/bin/check_environment.sh" >&2
  exit 1
fi
if [[ ! -f "$SCENICPLUS_HOME/scripts/init_scenicplus_project.py" ]]; then
  echo "ERROR: project initializer not found: $SCENICPLUS_HOME/scripts/init_scenicplus_project.py" >&2
  exit 1
fi

mkdir -p "$PROJECT_DIR/logs"
echo "Running environment check..."
bash "$SCENICPLUS_HOME/bin/check_environment.sh" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  2>&1 | tee "$PROJECT_DIR/logs/check_environment.log"

touch "$CONFIG_FILE"
update_config_setting "CONDA_ROOT" "$CONDA_ROOT"
update_config_setting "ENV_NAME" "$ENV_NAME"
update_config_setting "PROJECT_DIR" "$PROJECT_DIR"
update_config_setting "ORGANISM" "$ORGANISM"
update_config_setting "AUTOZYME" "$AUTOZYME"
update_config_setting "ENSEMBL_RELEASE" "$ENSEMBL_RELEASE"
update_config_setting "SCENICPLUS_MAX_MEMORY_GB" "$SCENICPLUS_MAX_MEMORY_GB"
for optional_key in ANNOTATED_OBJECT CELL_LABEL_COLUMN ATAC_INPUT_LAYOUT ATAC_DATA_ROOT; do
  optional_value="${!optional_key:-}"
  if [[ -n "$optional_value" ]]; then
    if [[ "$optional_key" == "ANNOTATED_OBJECT" || "$optional_key" == "ATAC_DATA_ROOT" ]]; then
      optional_value="$(expand_path "$optional_value")"
    fi
    update_config_setting "$optional_key" "$optional_value"
  fi
done

echo "Initializing project parameters and directory layout..."
"$ENV_PREFIX/bin/python" "$SCENICPLUS_HOME/scripts/init_scenicplus_project.py" \
  --config "$CONFIG_FILE" \
  --project-dir "$PROJECT_DIR" \
  --organism "$ORGANISM" \
  --autozyme "$AUTOZYME" \
  --conda-root "$CONDA_ROOT" \
  --env-name "$ENV_NAME" \
  --ensembl-release "$ENSEMBL_RELEASE" \
  --max-memory-gb "$SCENICPLUS_MAX_MEMORY_GB" \
  2>&1 | tee "$PROJECT_DIR/logs/initialize_scenicplus_project.log"

echo "WROTE $CONFIG_FILE"

echo
echo "PROJECT INITIALIZATION OK"
echo "Next command:"
echo "  source $(shell_quote "$PROJECT_DIR/project_env.sh")"
