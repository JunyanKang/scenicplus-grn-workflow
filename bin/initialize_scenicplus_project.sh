#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  CONDA_ROOT=/absolute/path/to/conda ENV_NAME=scenicplus-grn PROJECT_DIR=/absolute/path/to/project \
    AUTOZYME=on bash initialize_scenicplus_project.sh

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

AUTOZYME="${AUTOZYME:-on}"
ORGANISM="${ORGANISM:-}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-}"
SCENICPLUS_MAX_MEMORY_GB="${SCENICPLUS_MAX_MEMORY_GB:-auto}"
MOTIF2TF_REFERENCE="${MOTIF2TF_REFERENCE:-}"
MOTIF2TF_TABLE="${MOTIF2TF_TABLE:-}"
ANNOTATED_OBJECT="${ANNOTATED_OBJECT:-}"
CELL_LABEL_COLUMN="${CELL_LABEL_COLUMN:-}"
ATAC_INPUT_LAYOUT="${ATAC_INPUT_LAYOUT:-}"
ATAC_DATA_ROOT="${ATAC_DATA_ROOT:-}"

CONDA_ROOT="$(expand_path "$CONDA_ROOT")"
PROJECT_DIR="$(expand_path "$PROJECT_DIR")"
if [[ -n "$ANNOTATED_OBJECT" ]]; then
  ANNOTATED_OBJECT="$(expand_path "$ANNOTATED_OBJECT")"
fi
if [[ -n "$ATAC_DATA_ROOT" ]]; then
  ATAC_DATA_ROOT="$(expand_path "$ATAC_DATA_ROOT")"
fi
if [[ -n "$MOTIF2TF_TABLE" ]]; then
  MOTIF2TF_TABLE="$(expand_path "$MOTIF2TF_TABLE")"
fi

if [[ -n "$MOTIF2TF_REFERENCE" ]]; then
  case "$MOTIF2TF_REFERENCE" in
    human|mouse|fly|chicken) ;;
    *)
      echo "ERROR: MOTIF2TF_REFERENCE must be human, mouse, fly, or chicken." >&2
      exit 2
      ;;
  esac
fi
case "$AUTOZYME" in
  on|off) ;;
  *)
    echo "ERROR: AUTOZYME must be on or off." >&2
    exit 2
    ;;
esac
if [[ -n "$ATAC_INPUT_LAYOUT" ]]; then
  case "$ATAC_INPUT_LAYOUT" in
    split_ge_arc|cellranger_outs) ;;
    *)
      echo "ERROR: ATAC_INPUT_LAYOUT must be split_ge_arc or cellranger_outs." >&2
      exit 2
      ;;
  esac
fi
if [[ -n "$ANNOTATED_OBJECT" && ! -f "$ANNOTATED_OBJECT" ]]; then
  echo "ERROR: ANNOTATED_OBJECT not found: $ANNOTATED_OBJECT" >&2
  exit 1
fi
if [[ -n "$ATAC_DATA_ROOT" && ! -d "$ATAC_DATA_ROOT" ]]; then
  echo "ERROR: ATAC_DATA_ROOT not found: $ATAC_DATA_ROOT" >&2
  exit 1
fi
if [[ -n "$MOTIF2TF_TABLE" && ! -f "$MOTIF2TF_TABLE" ]]; then
  echo "ERROR: MOTIF2TF_TABLE not found: $MOTIF2TF_TABLE" >&2
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
INIT_CHECK="${SPGRN_INITIALIZE_CHECK:-auto}"
if [[ "$INIT_CHECK" == "always" ]] || [[ "$INIT_CHECK" == "1" ]] || ! compgen -G "$PROJECT_DIR/logs/check_environment_*.log" >/dev/null; then
  echo "Running environment check..."
  "$ENV_PREFIX/bin/spgrn-check"
else
  echo "Skipping environment check; existing check log found. Set SPGRN_INITIALIZE_CHECK=always to rerun."
fi

touch "$CONFIG_FILE"
update_config_setting "CONDA_ROOT" "$CONDA_ROOT"
update_config_setting "ENV_NAME" "$ENV_NAME"
update_config_setting "PROJECT_DIR" "$PROJECT_DIR"
update_config_setting "AUTOZYME" "$AUTOZYME"
update_config_setting "SCENICPLUS_MAX_MEMORY_GB" "$SCENICPLUS_MAX_MEMORY_GB"
for optional_key in ORGANISM ENSEMBL_RELEASE MOTIF2TF_REFERENCE ANNOTATED_OBJECT CELL_LABEL_COLUMN ATAC_INPUT_LAYOUT ATAC_DATA_ROOT MOTIF2TF_TABLE; do
  optional_value="${!optional_key:-}"
  if [[ -n "$optional_value" ]]; then
    if [[ "$optional_key" == "ANNOTATED_OBJECT" || "$optional_key" == "ATAC_DATA_ROOT" || "$optional_key" == "MOTIF2TF_TABLE" ]]; then
      optional_value="$(expand_path "$optional_value")"
    fi
    update_config_setting "$optional_key" "$optional_value"
  fi
done

echo "Initializing project parameters and directory layout..."
init_cmd=(
  "$ENV_PREFIX/bin/python" "$SCENICPLUS_HOME/scripts/init_scenicplus_project.py"
  --config "$CONFIG_FILE"
  --project-dir "$PROJECT_DIR"
  --autozyme "$AUTOZYME"
  --conda-root "$CONDA_ROOT"
  --env-name "$ENV_NAME"
  --max-memory-gb "$SCENICPLUS_MAX_MEMORY_GB"
)
if [[ -n "$ORGANISM" ]]; then
  init_cmd+=(--organism "$ORGANISM")
fi
if [[ -n "$ENSEMBL_RELEASE" ]]; then
  init_cmd+=(--ensembl-release "$ENSEMBL_RELEASE")
fi
if [[ -n "$MOTIF2TF_REFERENCE" ]]; then
  init_cmd+=(--motif2tf-reference "$MOTIF2TF_REFERENCE")
fi
"${init_cmd[@]}" 2>&1 | tee "$PROJECT_DIR/logs/initialize_scenicplus_project.log"

echo "WROTE $CONFIG_FILE"

echo
echo "PROJECT INITIALIZATION OK"
echo "Next command:"
echo "  source $(shell_quote "$PROJECT_DIR/project_env.sh")"
