#!/bin/bash
# Create the dfepipe conda env and apply the dadi / dadi-cli patches.
#   bash setup/setup.sh [env_prefix]
# Idempotent: re-running with the env present only re-verifies the patches.
#
# The patches add the neugamma / gammaneud1 probability density
# functions. They were cut from the dadi0116 env on hoffman2 against pristine
# dadi 2.3.3 and dadi-cli 0.9.5; see records/260929_extract_patches.sh.
#
# APPLY_CDF_INTEGRATION=1 additionally applies the optional CDF cache-integration
# patch. The default (stock `quad` integration on the exterior tails) is what every
# reference result in this package was produced with; the evidence for the variant
# is in setup/patches/optional/reference_integration/.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prefix="${1:-${HOME}/condaenv/dfepipe}"

# conda first: it is what both our clusters have and what built the environment
# that produced the published results. micromamba is the fallback for a cluster
# with no conda at all - it is a single static binary needing no base install.
if command -v conda >/dev/null 2>&1; then
  CREATE=(conda env create -p "$prefix" -f "${here}/environment.yml")
elif command -v micromamba >/dev/null 2>&1; then
  CREATE=(micromamba create -y -p "$prefix" -f "${here}/environment.yml")
else
  echo "neither conda nor micromamba found on PATH" >&2
  echo "  (on a cluster: try 'module load conda' first)" >&2
  exit 1
fi

if [[ ! -d "$prefix" ]]; then
  "${CREATE[@]}"
fi

# dadi and dadi-cli come from the vendored sdists, never from PyPI, so the exact
# version that produced the published results is what gets built here. They are
# sdists because dadi has C/Cython extensions that must be compiled.
vendor="${here}/vendor"
for pkg in "${vendor}"/dadi-*.tar.gz "${vendor}"/dadi_cli-*.tar.gz; do
  [[ -e "$pkg" ]] || { echo "missing vendored sdist: $pkg" >&2; exit 1; }
done
if ! "$prefix/bin/python" -c "import dadi, dadi_cli" >/dev/null 2>&1; then
  "$prefix/bin/pip" install --no-index --no-build-isolation \
      "${vendor}"/dadi-*.tar.gz "${vendor}"/dadi_cli-*.tar.gz
fi

sp="$(ls -d "$prefix"/lib/python3.*/site-packages 2>/dev/null | head -1)"
[[ -n "$sp" && -d "$sp" ]] || { echo "no site-packages under $prefix" >&2; exit 1; }

# Apply one patch, or confirm it is already applied. Anything else is fatal:
# a half-patched env silently produces wrong DFEs.
#
# Every call passes -N (--forward), and the already-applied case is tested first.
# Without -N, `patch` notices that an already-applied patch reverses cleanly and
# applies it in REVERSE, so a second run of this script would quietly strip the
# patches back out. -N makes patch refuse instead.
#
# The three-way ladder below was checked against GNU patch 2.7.1 (hoffman2) and
# Apple patch 2.0-12u11 (macOS), which disagree on `-R -N` for an already-patched
# file. On both, pristine -> apply, patched -> already applied, and anything else
# (upstream moved, wrong file) -> all three fail and the run stops.
apply_patch() {
  local p="$1" f="$2"
  [[ -f "$f" ]] || { echo "patch target missing: $f" >&2; exit 1; }
  if patch -R -N --dry-run -s -p0 "$f" < "$p" >/dev/null 2>&1; then
    echo "already applied: $(basename "$p")"
  elif patch -N --dry-run -s -p0 "$f" < "$p" >/dev/null 2>&1; then
    patch -N -p0 "$f" < "$p"
    echo "applied $(basename "$p")"
  elif patch -R --dry-run -s -p0 "$f" < "$p" >/dev/null 2>&1; then
    echo "already applied: $(basename "$p")"
  else
    echo "PATCH FAILED for $f -- upstream has moved; regenerate the patch" >&2
    exit 1
  fi
}

for p in "${here}"/patches/*.patch; do
  case "$(basename "$p" .patch)" in
    dadi_cli_Pdfs)   apply_patch "$p" "$sp/dadi_cli/Pdfs.py" ;;
    dadi_DFE_PDFs)   apply_patch "$p" "$sp/dadi/DFE/PDFs.py" ;;
    *) echo "unknown patch $p" >&2; exit 1 ;;
  esac
done

if [[ "${APPLY_CDF_INTEGRATION:-0}" == "1" ]]; then
  apply_patch "${here}/patches/optional/cache1d_cdf_integration.patch" \
              "$sp/dadi/DFE/Cache1D_mod.py"
else
  echo "skipping optional CDF cache integration (APPLY_CDF_INTEGRATION=1 to apply)"
fi

# Goal tests: not "is dadi 2.3.3 installed" but "does this environment do what
# the patch is for". These are also what to run against a future dadi to decide
# whether the patch is still needed.
"$prefix/bin/python" "$(dirname "$here")/scripts/checks/test_env.py" --prefix "$prefix"

# Capture the exact build so it can be recreated rather than merely approximated.
mkdir -p "${here}/lock"
lockname="${LOCK_NAME:-$(hostname -s)}"
if command -v conda >/dev/null 2>&1; then
  conda list --explicit -p "$prefix" > "${here}/lock/${lockname}.txt" 2>/dev/null || true
fi

echo "env ready: $prefix"
echo
echo "Put this in your site profile (config/profiles/<cluster>.conf):"
echo "  export CONDA_ENV=\"$prefix\""
