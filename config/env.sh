#!/bin/bash
# Software environment. EDIT TOOL VERSIONS HERE ONLY.
#
# Sourced by every stage script after config/pipeline.conf. Exposes
# dfepipe_init_modules, load_modules <mod>... and activate_conda.
BEDTOOLS_MOD="${BEDTOOLS_MOD:-bedtools}"
R_MOD="${R_MOD:-R}"
CONDA_MOD="${CONDA_MOD:-conda}"
CONDA_ENV="${CONDA_ENV:?CONDA_ENV not set: path to the dfepipe conda env}"

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

activate_conda() {
  dfepipe_init_modules      # `module` may not exist yet in a non-interactive shell
  _dfepipe_strict_off
  module load "$CONDA_MOD"
  conda activate "$CONDA_ENV" || source activate "$CONDA_ENV"
  _dfepipe_strict_restore
}
