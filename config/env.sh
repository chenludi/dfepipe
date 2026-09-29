#!/bin/bash
# Software environment. EDIT TOOL VERSIONS HERE ONLY.
#
# Sourced by every stage script after config/pipeline.conf. Exposes
# dfepipe_init_modules, load_modules <mod>... and activate_conda.
BEDTOOLS_MOD="${BEDTOOLS_MOD:-bedtools}"
R_MOD="${R_MOD:-R}"
CONDA_MOD="${CONDA_MOD:-conda}"
# Not asserted here: sourcing env.sh must not require conda, because stages
# like prep need only bedtools. activate_conda does the assertion instead.
CONDA_ENV="${CONDA_ENV:-}"

# SITE-SPECIFIC: the one line below is the lmod init script's location on UCLA
# Hoffman2. On another cluster, change MODULE_INIT (or export MODULESHOME so the
# init is skipped entirely). It is a software path, not a data path, so it is
# allowed to be hardcoded; data paths live in pipeline.conf and profiles/.
MODULE_INIT="${MODULE_INIT:-/u/local/Modules/default/init/bash}"

# lmod and conda init scripts misbehave under `set -eu`, so strictness is
# toggled off around them and restored exactly.
#
# `$(set +o)` runs in a command-substitution subshell, and bash reports errexit
# as OFF there however the caller had it -- so saving state from that alone
# silently drops `set -e` for the rest of the script. errexit is therefore read
# from $- in THIS shell and restored separately.
_dfepipe_strict_off() {
  _dfe_saved="$(set +o)"
  case "$-" in
    *e*) _dfe_saved_e="set -o errexit" ;;
    *)   _dfe_saved_e="set +o errexit" ;;
  esac
  set +eu
}
_dfepipe_strict_restore() { eval "$_dfe_saved"; eval "$_dfe_saved_e"; }

dfepipe_init_modules() {
  _dfepipe_strict_off
  [[ -n "${MODULESHOME:-}" ]] || source "$MODULE_INIT"
  _dfepipe_strict_restore
}

load_modules() {
  dfepipe_init_modules
  _dfepipe_strict_off
  local m; for m in "$@"; do module load "$m"; done
  _dfepipe_strict_restore
}

# The scheduler is not on PATH in a non-interactive shell (ssh host 'cmd', and
# anything a hook or automation runs), so a launcher must call this before qsub.
# SGE_ROOT and the arch subdirectory are site-specific; override them in a profile.
SGE_ROOT_DIR="${SGE_ROOT_DIR:-${SGE_ROOT:-/u/systems/UGE8.6.4}}"
SGE_ARCH="${SGE_ARCH:-lx-amd64}"
# SGE_CELL matters: qsub refuses to start without its cell directory, so putting
# the binary on PATH is not enough. settings.sh sets the whole environment.
SGE_CELL_NAME="${SGE_CELL_NAME:-hoffman2}"
ensure_scheduler() {
  command -v qsub >/dev/null 2>&1 && [[ -n "${SGE_CELL:-}" ]] && return 0
  local settings="${SGE_ROOT_DIR}/${SGE_CELL_NAME}/common/settings.sh"
  if [[ -r "$settings" ]]; then
    _dfepipe_strict_off
    # shellcheck disable=SC1090
    source "$settings"
    _dfepipe_strict_restore
  fi
  if ! command -v qsub >/dev/null 2>&1; then
    local d="${SGE_ROOT_DIR}/bin/${SGE_ARCH}"
    [[ -x "${d}/qsub" ]] || { echo "qsub not found (tried ${settings} and ${d})" >&2; return 1; }
    export SGE_ROOT="${SGE_ROOT_DIR}" SGE_CELL="${SGE_CELL_NAME}"
    export PATH="${d}:${PATH}"
  fi
}

activate_conda() {
  : "${CONDA_ENV:?CONDA_ENV not set: path to the conda env (see config/profiles/)}"
  dfepipe_init_modules      # `module` may not exist yet in a non-interactive shell
  _dfepipe_strict_off
  module load "$CONDA_MOD"
  conda activate "$CONDA_ENV" || source activate "$CONDA_ENV"
  _dfepipe_strict_restore
}
