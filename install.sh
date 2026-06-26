#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${BIN_DIR:-$SCRIPT_DIR/bin}"
CONFIG_DIR="${CONFIG_DIR:-$SCRIPT_DIR/config}"
DOCS_DIR="${DOCS_DIR:-$SCRIPT_DIR/docs}"
LOCKS_DIR="${LOCKS_DIR:-$CONFIG_DIR/locks}"
ENV_NAME="${ENV_NAME:-scenicplus-grn}"
MODE="${MODE:-new}"
INSTALL_R="${INSTALL_R:-1}"
FORCE="${FORCE:-0}"
ALLOW_BASE="${ALLOW_BASE:-0}"
ASSUME_YES="${ASSUME_YES:-0}"
AUTO_INSTALL_MAMBA="${AUTO_INSTALL_MAMBA:-1}"
RELOCATE_INSTALLER="${RELOCATE_INSTALLER:-1}"
PRECHECK_ONLY="${PRECHECK_ONLY:-0}"
LOG_DIR="${LOG_DIR:-$SCRIPT_DIR/logs}"
GITHUB_TRIES="${GITHUB_TRIES:-3}"
ARCHIVES_DIR="${ARCHIVES_DIR:-$SCRIPT_DIR/archives}"
VENDOR_DIR="${VENDOR_DIR:-$SCRIPT_DIR/.vendor}"
VENDOR_ARCHIVE="${VENDOR_ARCHIVE:-$ARCHIVES_DIR/vendor.tar.gz}"
VENDOR_GITHUB_DIR="${VENDOR_GITHUB_DIR:-$VENDOR_DIR/github}"
R_INSTALL_NCPUS="${R_INSTALL_NCPUS:-2}"
INSTALL_AUTOZYME="${INSTALL_AUTOZYME:-1}"
INSTALL_AUTOZYME_R="${INSTALL_AUTOZYME_R:-1}"
INSTALL_MALLET="${INSTALL_MALLET:-1}"

PYCISTOPIC_COMMIT="219225df56b32738d82cd14532b187a1483de04f"
PYCISTARGET_COMMIT="5aa517604e4842539a7531c16905825dc7cb80fb"
SCENICPLUS_COMMIT="e82b82f14b76618b850dfe442efc2421bb34f3b4"
CISTARGET_DB_COMMIT="304d5dc1b15e5c923908a50a1ec291c3faaccf9c"
CBUST_COMMIT="5911cd6201b767a43316ce613afc6c9255dc3511"
LOOMXPY_COMMIT="61995ff10940968eac2cee8fe48300ab477a15d0"
AUTOZYME_COMMIT="35f91f2229eb44d82710470803865d3c15102716"

PYCISTOPIC_VERSION="1.0.2"
PYCISTARGET_VERSION="1.1"
SCENICPLUS_VERSION="1.0a2"
LOOMXPY_VERSION="0.4.2"
AUTOZYME_VERSION="0.3.1"
MALLET_VERSION="2.0.8"
MALLET_URL="https://mallet.cs.umass.edu/dist/mallet-${MALLET_VERSION}.zip"

usage() {
  cat <<'EOF'
Usage:
  bash install.sh

Recommended workflow:
  tar -xzf scenicplus-grn-workflow.tar.gz
  cd scenicplus-grn-workflow
  bash install.sh

The script detects a conda/miniforge/miniconda/mamba-style root, asks for
confirmation, checks write permissions, then creates or updates a dedicated
scenicplus-grn conda environment. If the workflow package was unpacked in a random
directory, it offers to copy itself to:
  $CONDA_ROOT/share/scenicplus-grn-workflow

Modes:
  MODE=new      Create/update a dedicated scenicplus-grn conda environment.
  MODE=active   Install/update the currently activated conda environment.

Options:
  ENV_NAME=scenicplus-grn       Target env name for MODE=new.
  CONDA_ROOT=/absolute/path/to/conda     Conda/miniforge/miniconda/mamba root.
  FORCE=1                       Remove and recreate existing ENV_NAME.
  INSTALL_R=0                   Skip optional R/hdWGCNA layer.
  GITHUB_TRIES=3                GitHub attempts before bundled source archives are used.
  R_INSTALL_NCPUS=2             R package install parallelism.
  INSTALL_AUTOZYME=1            Install Python AutoZyme from pinned local archive with --no-deps.
  INSTALL_AUTOZYME_R=1          Install R AutoZyme only if its required R deps are already present.
  INSTALL_MALLET=1              Install MALLET 2.0.8 for pycisTopic's MALLET LDA backend.
  ALLOW_BASE=1                  Allow MODE=active installation into base.
  ASSUME_YES=1                  Accept detected defaults without prompts.
  AUTO_INSTALL_MAMBA=0          Do not bootstrap mamba into conda base.
  RELOCATE_INSTALLER=0          Do not offer to copy workflow package under conda root.
  PRECHECK_ONLY=1               Check detection/permissions and exit before install.
  LOG_DIR=/path/to/logs         Override installation log directory.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

setup_logging() {
  local stamp
  stamp="$(date +%Y%m%d_%H%M%S)"
  if ! mkdir -p "$LOG_DIR" 2>/dev/null; then
    LOG_DIR="${TMPDIR:-/tmp}/scenicplus-grn-workflow-logs"
    mkdir -p "$LOG_DIR"
  fi
  LOG_FILE="$LOG_DIR/install_${stamp}_$$.log"
  touch "$LOG_FILE"
  exec > >(tee -a "$LOG_FILE") 2>&1
  echo "SCENIC+ GRN workflow install log: $LOG_FILE"
  echo "Started: $(date)"
  echo "Installer directory: $SCRIPT_DIR"
}

die() {
  echo "ERROR: $*" >&2
  if [[ -n "${LOG_FILE:-}" ]]; then
    echo "See log: $LOG_FILE" >&2
  fi
  exit 1
}

validate_settings() {
  if [[ ! "$GITHUB_TRIES" =~ ^[0-9]+$ ]]; then
    die "GITHUB_TRIES must be a non-negative integer, got: $GITHUB_TRIES"
  fi
  if [[ ! "$R_INSTALL_NCPUS" =~ ^[0-9]+$ || "$R_INSTALL_NCPUS" == "0" ]]; then
    die "R_INSTALL_NCPUS must be a positive integer, got: $R_INSTALL_NCPUS"
  fi
  export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
  export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
  export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
  export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
  export R_INSTALL_NCPUS
  export GITHUB_TRIES
  export INSTALL_AUTOZYME
  export INSTALL_AUTOZYME_R
  export INSTALL_MALLET
  export R_REMOTES_NO_ERRORS_FROM_WARNINGS="${R_REMOTES_NO_ERRORS_FROM_WARNINGS:-true}"
  export MAKEFLAGS="${MAKEFLAGS:--j1}"
}

confirm() {
  local prompt="$1"
  local default="${2:-y}"
  local answer
  if [[ "$ASSUME_YES" == "1" || ! -t 0 ]]; then
    [[ "$default" =~ ^[Yy]$ ]]
    return
  fi
  if [[ "$default" =~ ^[Yy]$ ]]; then
    read -r -p "$prompt [Y/n] " answer
    [[ -z "$answer" || "$answer" =~ ^[Yy]$ ]]
  else
    read -r -p "$prompt [y/N] " answer
    [[ "$answer" =~ ^[Yy]$ ]]
  fi
}

prompt_value() {
  local prompt="$1"
  local value
  if [[ ! -t 0 ]]; then
    return 1
  fi
  read -r -p "$prompt" value
  [[ -n "$value" ]] || return 1
  printf "%s\n" "$value"
}

resolve_path() {
  local path="$1"
  if [[ -d "$path" ]]; then
    (cd "$path" && pwd)
  else
    return 1
  fi
}

is_conda_root() {
  local root="$1"
  [[ -x "$root/bin/conda" ]]
}

root_from_conda_exe() {
  local exe="$1"
  local bindir
  bindir="$(cd "$(dirname "$exe")" && pwd)"
  (cd "$bindir/.." && pwd)
}

candidate_roots() {
  {
    if [[ -n "${CONDA_ROOT:-}" ]]; then
      printf "%s\n" "$CONDA_ROOT"
    fi
    if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE:-}" ]]; then
      root_from_conda_exe "$CONDA_EXE"
    fi
    if command -v conda >/dev/null 2>&1; then
      root_from_conda_exe "$(command -v conda)"
    fi
    printf "%s\n" \
      "$HOME/miniforge3" \
      "$HOME/miniconda3" \
      "$HOME/mambaforge" \
      "$HOME/anaconda3" \
      "$HOME/conda" \
      "/opt/conda" \
      "/usr/local/miniforge3" \
      "/usr/local/miniconda3"
    find "$HOME" -maxdepth 3 -type f -path "*/bin/conda" 2>/dev/null | while read -r exe; do
      root_from_conda_exe "$exe"
    done
  } | awk 'NF && !seen[$0]++'
}

detect_conda_root() {
  local candidate
  while IFS= read -r candidate; do
    [[ -n "$candidate" ]] || continue
    candidate="${candidate/#\~/$HOME}"
    if candidate="$(resolve_path "$candidate" 2>/dev/null)" && is_conda_root "$candidate"; then
      printf "%s\n" "$candidate"
      return 0
    fi
  done < <(candidate_roots)
  return 1
}

choose_conda_root() {
  local detected="${1:-}"
  local manual
  if [[ -n "$detected" ]]; then
    echo "Detected conda root: $detected"
    if confirm "Use this conda root?" "y"; then
      CONDA_ROOT="$detected"
      return
    fi
  fi

  manual="$(prompt_value "Enter conda/miniforge/miniconda/mamba root path, or Ctrl-C to stop: ")" || \
    die "Could not determine conda root. Set CONDA_ROOT=/path/to/miniforge3 and rerun."
  manual="${manual/#\~/$HOME}"
  manual="$(resolve_path "$manual")" || die "Conda root does not exist: $manual"
  is_conda_root "$manual" || die "No executable bin/conda found under: $manual"
  CONDA_ROOT="$manual"
}

path_is_inside() {
  local child="$1"
  local parent="$2"
  [[ "$child" == "$parent" || "$child" == "$parent/"* ]]
}

copy_distribution_assets() {
  local target_dir="$1"
  rm -rf "$target_dir/bin" "$target_dir/config" "$target_dir/docs"
  mkdir -p "$target_dir/bin"
  mkdir -p "$target_dir/config/locks"
  mkdir -p "$target_dir/docs"
  mkdir -p "$target_dir/archives"
  cp -f \
    "$SCRIPT_DIR/LICENSE" \
    "$SCRIPT_DIR/THIRD_PARTY_NOTICES.md" \
    "$SCRIPT_DIR/VERSION" \
    "$SCRIPT_DIR/README.md" \
    "$SCRIPT_DIR/README.en.md" \
    "$SCRIPT_DIR/install.sh" \
    "$target_dir/"
  cp -f \
    "$BIN_DIR/check_environment.sh" \
    "$BIN_DIR/initialize_scenicplus_project.sh" \
    "$BIN_DIR/install_r.R" \
    "$target_dir/bin/"
  cp -f \
    "$CONFIG_DIR/pip-constraints.txt" \
    "$CONFIG_DIR/environment-linux-64.yml" \
    "$CONFIG_DIR/environment-macos-arm64.yml" \
    "$CONFIG_DIR/scenicplus_config_template.yaml" \
    "$target_dir/config/"
  cp -f "$LOCKS_DIR/environment-linux-64.solved-lock.yml" "$target_dir/config/locks/"
  cp -f \
    "$DOCS_DIR/VERSION_LOCK.md" \
    "$DOCS_DIR/CHANGELOG.md" \
    "$DOCS_DIR/RELEASE_NOTES.md" \
    "$DOCS_DIR/SCENICPLUS_STEP_BY_STEP.md" \
    "$DOCS_DIR/SCENICPLUS_STEP_BY_STEP.en.md" \
    "$target_dir/docs/"
  if [[ -s "$VENDOR_ARCHIVE" ]]; then
    cp -f "$VENDOR_ARCHIVE" "$target_dir/archives/vendor.tar.gz"
  elif [[ -d "$VENDOR_DIR" ]]; then
    rm -rf "$target_dir/.vendor"
    cp -a "$VENDOR_DIR" "$target_dir/.vendor"
  fi
  if [[ -d "$SCRIPT_DIR/modules" ]]; then
    rm -rf "$target_dir/modules"
    cp -a "$SCRIPT_DIR/modules" "$target_dir/modules"
  fi
  if [[ -d "$SCRIPT_DIR/scripts" ]]; then
    rm -rf "$target_dir/scripts"
    cp -a "$SCRIPT_DIR/scripts" "$target_dir/scripts"
  fi
}

copy_to_conda_share_if_needed() {
  local target_dir="$CONDA_ROOT/share/scenicplus-grn-workflow"
  if [[ "$RELOCATE_INSTALLER" != "1" || "${SCENICPLUS_INSTALLER_RELOCATED:-0}" == "1" ]]; then
    return
  fi
  if path_is_inside "$SCRIPT_DIR" "$CONDA_ROOT"; then
    echo "Workflow package is already under conda root."
    return
  fi

  echo "Workflow package is not under the detected conda root."
  echo "Current workflow package directory: $SCRIPT_DIR"
  echo "Recommended workflow package directory: $target_dir"
  if ! confirm "Copy workflow package to the recommended directory and continue there?" "y"; then
    echo "Continuing from current directory."
    return
  fi

  check_writable_dir "$CONDA_ROOT"
  mkdir -p "$CONDA_ROOT/share"
  mkdir -p "$target_dir"
  copy_distribution_assets "$target_dir"
  chmod +x "$target_dir/install.sh" "$target_dir/bin/check_environment.sh" "$target_dir/bin/initialize_scenicplus_project.sh" "$target_dir/bin/install_r.R" "$target_dir/bin/run_python_entrypoint.py"
  echo "Copied workflow package to: $target_dir"
  echo "Re-running from recommended directory."
  export SCENICPLUS_INSTALLER_RELOCATED=1
  exec bash "$target_dir/install.sh" "$@"
}

install_env_wrappers() {
  local env_bin="$CONDA_PREFIX/bin"
  local share_dir="$CONDA_PREFIX/share/scenicplus-grn"
  local script_dir="$share_dir/scripts"
  local script base cmd wrapper
  mkdir -p "$env_bin"
  rm -f "$env_bin"/scenicplus-grn-* "$env_bin"/spgrn-*

  cat > "$env_bin/spgrn-check" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
exec "$(dirname "${BASH_SOURCE[0]}")/../share/scenicplus-grn/bin/check_environment.sh" "$@"
EOF
  cat > "$env_bin/spgrn-initialize" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
exec "$(dirname "${BASH_SOURCE[0]}")/../share/scenicplus-grn/bin/initialize_scenicplus_project.sh" "$@"
EOF

  if [[ -d "$script_dir" ]]; then
    for script in "$script_dir"/*.py; do
      [[ -f "$script" ]] || continue
      base="$(basename "$script" .py)"
      cmd="spgrn-${base//_/-}"
      wrapper="$env_bin/$cmd"
      cat > "$wrapper" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "\$(dirname "\${BASH_SOURCE[0]}")/python" "\$(dirname "\${BASH_SOURCE[0]}")/../share/scenicplus-grn/bin/run_python_entrypoint.py" "\$(dirname "\${BASH_SOURCE[0]}")/../share/scenicplus-grn/scripts/$base.py" "\$@"
EOF
      chmod +x "$wrapper"
    done
    for script in "$script_dir"/*.R; do
      [[ -f "$script" ]] || continue
      base="$(basename "$script" .R)"
      cmd="spgrn-${base//_/-}"
      wrapper="$env_bin/$cmd"
      cat > "$wrapper" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "\$(dirname "\${BASH_SOURCE[0]}")/Rscript" "\$(dirname "\${BASH_SOURCE[0]}")/../share/scenicplus-grn/scripts/$base.R" "\$@"
EOF
      chmod +x "$wrapper"
    done
  fi
  chmod +x "$env_bin/spgrn-check" "$env_bin/spgrn-initialize"
}

check_writable_dir() {
  local dir="$1"
  local test_file
  mkdir -p "$dir" || die "Cannot create directory: $dir"
  test_file="$dir/.scenicplus_write_test_$$"
  touch "$test_file" 2>/dev/null || die "No write permission for: $dir"
  rm -f "$test_file"
}

check_permissions() {
  echo "Checking permissions under conda root: $CONDA_ROOT"
  check_writable_dir "$CONDA_ROOT"
  check_writable_dir "$CONDA_ROOT/envs"
  check_writable_dir "$CONDA_ROOT/pkgs"
  check_writable_dir "$CONDA_ROOT/share"
}

select_env_file() {
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64)
      echo "$CONFIG_DIR/environment-macos-arm64.yml"
      ;;
    Linux-*)
      export CONDA_SUBDIR="${CONDA_SUBDIR:-linux-64}"
      echo "$CONFIG_DIR/environment-linux-64.yml"
      ;;
    *)
      die "Unsupported platform: $(uname -s)-$(uname -m)"
      ;;
  esac
}

conda_env_exists() {
  "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"
}

conda_env_supports_libmamba_solver() {
  "$CONDA_BIN" env create --help 2>/dev/null | grep -q "libmamba"
}

setup_solver_command() {
  MAMBA_BIN="$CONDA_ROOT/bin/mamba"
  SOLVER_KIND="mamba"
  if [[ -x "$MAMBA_BIN" ]]; then
    SOLVER_LABEL="$MAMBA_BIN"
    return
  fi

  if [[ "$AUTO_INSTALL_MAMBA" == "1" ]]; then
    echo "No mamba executable found at: $MAMBA_BIN"
    echo "Bootstrapping mamba into conda base from conda-forge using the existing conda solver."
    "$CONDA_BIN" install -n base -c conda-forge mamba -y
    MAMBA_BIN="$CONDA_ROOT/bin/mamba"
    [[ -x "$MAMBA_BIN" ]] || die "mamba installation did not create: $MAMBA_BIN"
    SOLVER_KIND="mamba"
    SOLVER_LABEL="$MAMBA_BIN"
    return
  fi

  if conda_env_supports_libmamba_solver; then
    echo "AUTO_INSTALL_MAMBA=0 and no mamba executable found; falling back to conda --solver libmamba."
    SOLVER_KIND="conda-libmamba"
    SOLVER_LABEL="$CONDA_BIN --solver libmamba"
    return
  fi

  die "No mamba executable found and this conda does not expose --solver libmamba. Rerun with AUTO_INSTALL_MAMBA=1 or install mamba into base."
}

solver_env_create() {
  if [[ "$SOLVER_KIND" == "mamba" ]]; then
    "$MAMBA_BIN" env create -n "$ENV_NAME" -f "$ENV_FILE"
  else
    "$CONDA_BIN" env create -n "$ENV_NAME" -f "$ENV_FILE" --solver libmamba
  fi
}

solver_env_update_name() {
  if [[ "$SOLVER_KIND" == "mamba" ]]; then
    "$MAMBA_BIN" env update -n "$ENV_NAME" -f "$ENV_FILE"
  else
    "$CONDA_BIN" env update -n "$ENV_NAME" -f "$ENV_FILE" --solver libmamba
  fi
}

solver_env_update_prefix() {
  if [[ "$SOLVER_KIND" == "mamba" ]]; then
    "$MAMBA_BIN" env update -p "$CONDA_PREFIX" -f "$ENV_FILE"
  else
    "$CONDA_BIN" env update -p "$CONDA_PREFIX" -f "$ENV_FILE" --solver libmamba
  fi
}

stage_if_running_from_target_prefix() {
  local target_prefix="$CONDA_ROOT/envs/$ENV_NAME"
  if [[ "${STAGED_INSTALLER:-0}" == "1" ]]; then
    return
  fi
  if [[ "$SCRIPT_DIR" == "$target_prefix/"* ]] && [[ "$FORCE" == "1" || ! -d "$target_prefix/conda-meta" ]]; then
    local staged_parent
    staged_parent="$(mktemp -d)"
    local staged="$staged_parent/scenicplus-grn-workflow"
    cp -a "$SCRIPT_DIR" "$staged"
    if [[ ! -d "$target_prefix/conda-meta" ]]; then
      rm -rf "$target_prefix"
    fi
    export STAGED_INSTALLER=1
    exec bash "$staged/install.sh" "$@"
  fi
}

extract_archive_to_dest() {
  local archive="$1"
  local dest="$2"
  local tmpdir
  local topdir

  [[ -s "$archive" ]] || return 1
  tmpdir="$(mktemp -d)"
  tar -xzf "$archive" -C "$tmpdir"
  topdir="$(find "$tmpdir" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
  if [[ ! -d "$topdir" ]]; then
    rm -rf "$tmpdir"
    die "Archive did not contain a top-level source directory: $archive"
  fi
  rm -rf "$dest"
  mkdir -p "$(dirname "$dest")"
  mv "$topdir" "$dest"
  rm -rf "$tmpdir"
}

ensure_vendor_available() {
  if [[ -d "$VENDOR_GITHUB_DIR" && -d "$VENDOR_DIR/mallet" ]]; then
    return
  fi
  if [[ ! -d "$VENDOR_DIR" && -d "$SCRIPT_DIR/vendor" ]]; then
    echo "Using legacy uncompressed vendor directory: $SCRIPT_DIR/vendor"
    VENDOR_DIR="$SCRIPT_DIR/vendor"
    VENDOR_GITHUB_DIR="$VENDOR_DIR/github"
    return
  fi
  if [[ -s "$VENDOR_ARCHIVE" ]]; then
    local tmpdir
    echo "Expanding bundled vendor archive: $VENDOR_ARCHIVE"
    rm -rf "$VENDOR_DIR"
    tmpdir="$(mktemp -d)"
    tar -xzf "$VENDOR_ARCHIVE" -C "$tmpdir"
    if [[ ! -d "$tmpdir/vendor" ]]; then
      rm -rf "$tmpdir"
      die "Vendor archive did not contain a top-level vendor directory: $VENDOR_ARCHIVE"
    fi
    mv "$tmpdir/vendor" "$VENDOR_DIR"
    rm -rf "$tmpdir"
  fi
  if [[ ! -d "$VENDOR_GITHUB_DIR" ]]; then
    echo "No bundled vendor/github directory found. GitHub fallback archives will be unavailable."
  fi
  if [[ ! -d "$VENDOR_DIR/mallet" ]]; then
    echo "No bundled vendor/mallet directory found. MALLET will be downloaded if INSTALL_MALLET=1."
  fi
}

sha256_file() {
  local file="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file" | awk '{print $1}'
  else
    shasum -a 256 "$file" | awk '{print $1}'
  fi
}

verify_vendor_archive() {
  local archive="$1"
  local manifest="$VENDOR_GITHUB_DIR/SHA256SUMS"
  local name
  local expected
  local observed

  [[ -s "$archive" ]] || return 1
  [[ -s "$manifest" ]] || return 0

  name="$(basename "$archive")"
  expected="$(awk -v n="$name" '$2 == n {print $1}' "$manifest")"
  [[ -n "$expected" ]] || die "No checksum entry for bundled source archive: $name"
  observed="$(sha256_file "$archive")"
  [[ "$observed" == "$expected" ]] || die "Checksum mismatch for bundled source archive: $name"
}

git_clone_checkout() {
  local repo_url="$1"
  local dest="$2"
  local ref="$3"
  local archive="${4:-}"
  local branch="${5:-}"
  local attempt
  local -a clone_cmd

  for ((attempt = 1; attempt <= GITHUB_TRIES; attempt++)); do
    rm -rf "$dest"
    echo "Cloning $repo_url at $ref (attempt $attempt/$GITHUB_TRIES)."
    if [[ -n "$branch" ]]; then
      clone_cmd=(git -c http.lowSpeedLimit=1000 -c http.lowSpeedTime=60 clone --filter=blob:none -b "$branch" "$repo_url" "$dest")
    else
      clone_cmd=(git -c http.lowSpeedLimit=1000 -c http.lowSpeedTime=60 clone --filter=blob:none "$repo_url" "$dest")
    fi

    if "${clone_cmd[@]}" && git -C "$dest" checkout "$ref"; then
      return 0
    fi

    echo "Clone or checkout failed for $repo_url; retrying after ${attempt}s."
    sleep "$attempt"
  done

  if [[ -n "$archive" && -s "$archive" ]]; then
    echo "GitHub clone failed after $GITHUB_TRIES attempts; using bundled source archive: $archive"
    verify_vendor_archive "$archive"
    extract_archive_to_dest "$archive" "$dest"
    return 0
  fi

  die "Failed to clone/checkout $repo_url at $ref after $GITHUB_TRIES attempts, and no bundled source archive was available."
}

pip_install_github_or_local() {
  local requirement="$1"
  local archive="$2"
  local label="$3"
  local scm_env_name="${4:-}"
  local scm_version="${5:-}"
  local attempt
  local -a pip_env

  pip_env=()
  if [[ -n "$scm_env_name" && -n "$scm_version" ]]; then
    pip_env=(env "${scm_env_name}=${scm_version}")
  fi

  for ((attempt = 1; attempt <= GITHUB_TRIES; attempt++)); do
    echo "Installing $label from GitHub (attempt $attempt/$GITHUB_TRIES)."
    if "${pip_env[@]}" python -m pip install --retries 3 --timeout 120 --no-deps --force-reinstall -c "$CONFIG_DIR/pip-constraints.txt" "$requirement"; then
      return 0
    fi
    echo "GitHub pip install failed for $label; retrying after ${attempt}s."
    sleep "$attempt"
  done

  if [[ -s "$archive" ]]; then
    echo "GitHub pip install failed after $GITHUB_TRIES attempts; installing $label from bundled source archive: $archive"
    verify_vendor_archive "$archive"
    "${pip_env[@]}" python -m pip install --no-deps --force-reinstall -c "$CONFIG_DIR/pip-constraints.txt" "$archive"
    return 0
  fi

  die "Failed to install $label from GitHub after $GITHUB_TRIES attempts, and no bundled source archive was available."
}


install_autozyme_python() {
  if [[ "$INSTALL_AUTOZYME" != "1" ]]; then
    echo "INSTALL_AUTOZYME=0: skipping Python AutoZyme."
    return
  fi
  local workdir="$1"
  local archive="$VENDOR_GITHUB_DIR/autozyme-${AUTOZYME_COMMIT}.tar.gz"
  local src="$workdir/autozyme"
  [[ -s "$archive" ]] || die "Missing bundled AutoZyme archive: $archive"
  echo "Installing Python AutoZyme $AUTOZYME_VERSION from pinned bundled archive as editable no-deps overlay."
  verify_vendor_archive "$archive"
  extract_archive_to_dest "$archive" "$src"
  rm -rf "$CONDA_PREFIX/opt/autozyme"
  cp -a "$src" "$CONDA_PREFIX/opt/autozyme"
  python -m pip install --no-deps --no-build-isolation --force-reinstall -e "$CONDA_PREFIX/opt/autozyme/autozyme_py"
}

install_mallet() {
  if [[ "$INSTALL_MALLET" != "1" ]]; then
    echo "INSTALL_MALLET=0: skipping MALLET."
    return
  fi
  local archive="$VENDOR_DIR/mallet/mallet-${MALLET_VERSION}.zip"
  local workdir="$1"
  local zip_path="$archive"
  if [[ ! -s "$zip_path" ]]; then
    zip_path="$workdir/mallet-${MALLET_VERSION}.zip"
    echo "Bundled MALLET archive not found; downloading $MALLET_URL"
    curl -L --retry 3 --retry-delay 5 -o "$zip_path" "$MALLET_URL"
  else
    echo "Installing MALLET $MALLET_VERSION from bundled archive: $archive"
  fi
  command -v unzip >/dev/null 2>&1 || die "unzip is required to install MALLET."
  rm -rf "$workdir/mallet-${MALLET_VERSION}" "$CONDA_PREFIX/opt/mallet-${MALLET_VERSION}"
  unzip -q "$zip_path" -d "$workdir"
  [[ -x "$workdir/mallet-${MALLET_VERSION}/bin/mallet" ]] || chmod +x "$workdir/mallet-${MALLET_VERSION}/bin/mallet"
  cp -a "$workdir/mallet-${MALLET_VERSION}" "$CONDA_PREFIX/opt/mallet-${MALLET_VERSION}"
  chmod +x "$CONDA_PREFIX/opt/mallet-${MALLET_VERSION}/bin/mallet"
  cat > "$CONDA_PREFIX/bin/mallet" <<EOF
#!/usr/bin/env python3
import os
import subprocess
import sys

MALLET_ROOT = "$CONDA_PREFIX/opt/mallet-${MALLET_VERSION}"
CLASS_MAP = {
    "import-dir": "cc.mallet.classify.tui.Text2Vectors",
    "import-file": "cc.mallet.classify.tui.Csv2Vectors",
    "import-svmlight": "cc.mallet.classify.tui.SvmLight2Vectors",
    "info": "cc.mallet.classify.tui.Vectors2Info",
    "train-classifier": "cc.mallet.classify.tui.Vectors2Classify",
    "classify-dir": "cc.mallet.classify.tui.Text2Classify",
    "classify-file": "cc.mallet.classify.tui.Csv2Classify",
    "classify-svmlight": "cc.mallet.classify.tui.SvmLight2Classify",
    "train-topics": "cc.mallet.topics.tui.TopicTrainer",
    "infer-topics": "cc.mallet.topics.tui.InferTopics",
    "evaluate-topics": "cc.mallet.topics.tui.EvaluateTopics",
    "prune": "cc.mallet.classify.tui.Vectors2Vectors",
    "split": "cc.mallet.classify.tui.Vectors2Vectors",
    "bulk-load": "cc.mallet.util.BulkLoader",
}


def print_help() -> int:
    sys.stderr.write(
        "Mallet 2.0 commands:\\n\\n"
        "  import-dir\\n  import-file\\n  import-svmlight\\n  info\\n"
        "  train-classifier\\n  classify-dir\\n  classify-file\\n"
        "  classify-svmlight\\n  train-topics\\n  infer-topics\\n"
        "  evaluate-topics\\n  prune\\n  split\\n  bulk-load\\n\\n"
        "Include --help TRUE with any option for more information.\\n"
    )
    return 1


def main() -> int:
    args = sys.argv[1:]
    if not args:
        return print_help()
    command = args.pop(0)
    if command == "run":
        if not args:
            sys.stderr.write("mallet run requires a Java class name.\\n")
            return 1
        klass = args.pop(0)
    else:
        klass = CLASS_MAP.get(command)
        if klass is None:
            sys.stderr.write(f"Unrecognized command: {command}\\n")
            return print_help()
    classpath = os.pathsep.join(
        [
            os.path.join(MALLET_ROOT, "class"),
            os.path.join(MALLET_ROOT, "lib", "mallet-deps.jar"),
            os.environ.get("CLASSPATH", ""),
        ]
    )
    memory = os.environ.get("MALLET_MEMORY", "1g")
    java = os.environ.get("JAVA", "java")
    java_cmd = [
        java,
        f"-Xmx{memory}",
        "-ea",
        "-Djava.awt.headless=true",
        "-Dfile.encoding=UTF-8",
        "-server",
        "-classpath",
        classpath,
        klass,
        *args,
    ]
    return subprocess.call(java_cmd)


if __name__ == "__main__":
    raise SystemExit(main())
EOF
  chmod +x "$CONDA_PREFIX/bin/mallet"
  [[ -s "$CONDA_PREFIX/opt/mallet-${MALLET_VERSION}/lib/mallet-deps.jar" ]] || \
    die "MALLET installed but lib/mallet-deps.jar is missing."
  [[ -s "$CONDA_PREFIX/opt/mallet-${MALLET_VERSION}/class/cc/mallet/topics/tui/TopicTrainer.class" ]] || \
    die "MALLET installed but TopicTrainer.class is missing."
  python - <<PY
import pathlib
import shutil
import subprocess
import sys
import tempfile

mallet = pathlib.Path("$CONDA_PREFIX/bin/mallet")
tmp = pathlib.Path(tempfile.mkdtemp())
try:
    (tmp / "toy.txt").write_text("doc1\\talpha beta gamma\\ndoc2\\tbeta delta\\n")
    out = tmp / "toy.mallet"
    proc = subprocess.run(
        [
            str(mallet),
            "import-file",
            "--input",
            str(tmp / "toy.txt"),
            "--output",
            str(out),
            "--keep-sequence",
            "TRUE",
        ],
        text=True,
        capture_output=True,
        timeout=30,
    )
    if proc.returncode != 0 or not out.exists() or out.stat().st_size == 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit("MALLET import-file smoke test failed.")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
PY
  echo "MALLET $MALLET_VERSION installed: $CONDA_PREFIX/bin/mallet"
}

install_source_layer() {
  if [[ -z "${CONDA_PREFIX:-}" ]]; then
    die "CONDA_PREFIX is not set. Activate the target environment first."
  fi

  mkdir -p "$CONDA_PREFIX/opt"
  local workdir
  workdir="$(mktemp -d)"
  trap 'rm -rf "${workdir:-}"' RETURN

  python -m pip install --retries 5 --timeout 120 --upgrade "setuptools==80.10.2" "wheel==0.47.0"

  python -m pip install --retries 5 --timeout 120 --no-deps --force-reinstall -c "$CONFIG_DIR/pip-constraints.txt" \
    -r "$CONFIG_DIR/pip-constraints.txt"

  pip_install_github_or_local \
    "loomxpy @ git+https://github.com/aertslab/LoomXpy@${LOOMXPY_COMMIT}" \
    "$VENDOR_GITHUB_DIR/LoomXpy-${LOOMXPY_COMMIT}.tar.gz" \
    "LoomXpy" \
    "SETUPTOOLS_SCM_PRETEND_VERSION_FOR_LOOMXPY" \
    "$LOOMXPY_VERSION"

  git_clone_checkout https://github.com/aertslab/pycisTopic.git "$workdir/pycisTopic" "$PYCISTOPIC_COMMIT" "$VENDOR_GITHUB_DIR/pycisTopic-${PYCISTOPIC_COMMIT}.tar.gz"
  git_clone_checkout https://github.com/aertslab/pycistarget.git "$workdir/pycistarget" "$PYCISTARGET_COMMIT" "$VENDOR_GITHUB_DIR/pycistarget-${PYCISTARGET_COMMIT}.tar.gz"
  git_clone_checkout https://github.com/aertslab/scenicplus.git "$workdir/scenicplus" "$SCENICPLUS_COMMIT" "$VENDOR_GITHUB_DIR/scenicplus-${SCENICPLUS_COMMIT}.tar.gz"

  SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PYCISTOPIC="$PYCISTOPIC_VERSION" \
  SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PYCISTARGET="$PYCISTARGET_VERSION" \
  SETUPTOOLS_SCM_PRETEND_VERSION_FOR_SCENICPLUS="$SCENICPLUS_VERSION" \
  python -m pip install --no-deps --force-reinstall \
    "$workdir/pycisTopic" \
    "$workdir/pycistarget" \
    "$workdir/scenicplus"

  rm -rf "$CONDA_PREFIX/opt/pycisTopic" "$CONDA_PREFIX/opt/pycistarget" "$CONDA_PREFIX/opt/scenicplus"
  cp -a "$workdir/pycisTopic" "$CONDA_PREFIX/opt/pycisTopic"
  cp -a "$workdir/pycistarget" "$CONDA_PREFIX/opt/pycistarget"
  cp -a "$workdir/scenicplus" "$CONDA_PREFIX/opt/scenicplus"

  rm -rf "$CONDA_PREFIX/opt/create_cisTarget_databases"
  git_clone_checkout https://github.com/aertslab/create_cisTarget_databases.git "$CONDA_PREFIX/opt/create_cisTarget_databases" "$CISTARGET_DB_COMMIT" "$VENDOR_GITHUB_DIR/create_cisTarget_databases-${CISTARGET_DB_COMMIT}.tar.gz"

  git_clone_checkout https://github.com/ghuls/cluster-buster.git "$workdir/cluster-buster" "$CBUST_COMMIT" "$VENDOR_GITHUB_DIR/cluster-buster-${CBUST_COMMIT}.tar.gz" "change_f4_output"
  make -C "$workdir/cluster-buster" clean || true
  make -C "$workdir/cluster-buster" cbust CXXFLAGS="-Wall -std=c++11 -O3 -D NDEBUG"
  cp -f "$workdir/cluster-buster/cbust" "$CONDA_PREFIX/bin/cbust"
  chmod +x "$CONDA_PREFIX/bin/cbust"
  rm -rf "$CONDA_PREFIX/opt/cluster-buster"
  cp -a "$workdir/cluster-buster" "$CONDA_PREFIX/opt/cluster-buster"

  install_autozyme_python "$workdir"
  install_mallet "$workdir"

  rm -rf "$workdir"
  trap - RETURN
}

copy_recipe_into_env() {
  mkdir -p "$CONDA_PREFIX/share/scenicplus-grn"
  copy_distribution_assets "$CONDA_PREFIX/share/scenicplus-grn"
  if [[ -d "$CONDA_PREFIX/share/scenicplus-grn/scripts" ]]; then
    find "$CONDA_PREFIX/share/scenicplus-grn/scripts" -type f -name "*.py" -exec chmod +x {} \;
  fi
  chmod +x "$CONDA_PREFIX/share/scenicplus-grn/install.sh" "$CONDA_PREFIX/share/scenicplus-grn/bin/check_environment.sh" "$CONDA_PREFIX/share/scenicplus-grn/bin/initialize_scenicplus_project.sh" "$CONDA_PREFIX/share/scenicplus-grn/bin/install_r.R"
  install_env_wrappers
}

run_checks() {
  local check_args=(--conda-root "$CONDA_ROOT" --env-name "$ENV_NAME" --skip-workflow-assets)
  if [[ "$INSTALL_R" == "1" ]]; then
    bash "$BIN_DIR/check_environment.sh" "${check_args[@]}"
  else
    bash "$BIN_DIR/check_environment.sh" "${check_args[@]}" --skip-r
  fi
}

run_workflow_asset_check() {
  export SCENICPLUS_HOME="$CONDA_PREFIX/share/scenicplus-grn"
  if [[ "$INSTALL_MALLET" == "1" ]]; then
    SCENICPLUS_REQUIRE_MALLET=1 python "$SCENICPLUS_HOME/scripts/check_workflow_installation.py"
  else
    SCENICPLUS_REQUIRE_MALLET=0 python "$SCENICPLUS_HOME/scripts/check_workflow_installation.py"
  fi
}

main() {
  setup_logging
  validate_settings

  local detected_root
  detected_root="$(detect_conda_root || true)"
  choose_conda_root "$detected_root"
  CONDA_BIN="$CONDA_ROOT/bin/conda"

  copy_to_conda_share_if_needed "$@"
  ensure_vendor_available
  check_permissions

  ENV_FILE="$(select_env_file)"
  echo "Using environment file: $ENV_FILE"
  echo "Using conda command: $CONDA_BIN"
  echo "Target environment: $ENV_NAME"

  if [[ "$PRECHECK_ONLY" == "1" ]]; then
    if [[ -x "$CONDA_ROOT/bin/mamba" ]]; then
      echo "PRECHECK_ONLY=1: mamba is already available at $CONDA_ROOT/bin/mamba."
    elif [[ "$AUTO_INSTALL_MAMBA" == "1" ]]; then
      echo "PRECHECK_ONLY=1: mamba is not present and would be bootstrapped into conda base during installation."
    elif conda_env_supports_libmamba_solver; then
      echo "PRECHECK_ONLY=1: mamba is not present; installer would use conda --solver libmamba."
    else
      echo "PRECHECK_ONLY=1: mamba is not present and AUTO_INSTALL_MAMBA=0 would block installation."
    fi
    echo "PRECHECK_ONLY=1: detection, relocation, permissions, and platform recipe checks passed."
    echo "Log saved to: $LOG_FILE"
    exit 0
  fi

  setup_solver_command
  echo "Using solver command: $SOLVER_LABEL"

  case "$MODE" in
    new)
      stage_if_running_from_target_prefix "$@"

      if conda_env_exists; then
        if [[ "$FORCE" == "1" ]]; then
          "$CONDA_BIN" env remove -n "$ENV_NAME" -y
        else
          echo "Environment $ENV_NAME already exists; updating it."
        fi
      fi

      if conda_env_exists; then
        solver_env_update_name
      else
        solver_env_create
      fi

      eval "$("$CONDA_BIN" shell.bash hook)"
      set +u
      conda activate "$ENV_NAME"
      set -u
      ;;
    active)
      if [[ -z "${CONDA_PREFIX:-}" ]]; then
        die "No active conda environment detected."
      fi
      if [[ "${CONDA_DEFAULT_ENV:-}" == "base" && "$ALLOW_BASE" != "1" ]]; then
        die "Refusing to install into conda base. Set ALLOW_BASE=1 to override."
      fi
      solver_env_update_prefix
      ;;
    *)
      usage
      die "Unknown MODE=$MODE"
      ;;
  esac

  install_source_layer
  if [[ "$INSTALL_R" == "1" ]]; then
    export SCENICPLUS_VENDOR_GITHUB_DIR="$VENDOR_GITHUB_DIR"
    export SCENICPLUS_AUTOZYME_ARCHIVE="$VENDOR_GITHUB_DIR/autozyme-${AUTOZYME_COMMIT}.tar.gz"
    Rscript "$BIN_DIR/install_r.R"
  fi
  run_checks
  copy_recipe_into_env
  run_workflow_asset_check

  echo "DONE: SCENIC+ environment is installed and checked."
  echo "Log saved to: $LOG_FILE"
}

main "$@"
