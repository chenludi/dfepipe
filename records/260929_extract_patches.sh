#!/bin/bash
# 2026/09/29  Extract the neugamma / neugammazero / gammaneud1 additions to dadi and
# dadi-cli out of the dadi0116 conda env on hoffman2 and write them as patch files.
#
#   bash records/260929_extract_patches.sh        # run from the package root
#
# Produces:
#   setup/patches/dadi_cli_Pdfs.patch   -> dadi_cli/Pdfs.py    (dadi-cli 0.9.5)
#   setup/patches/dadi_DFE_PDFs.patch   -> dadi/DFE/PDFs.py    (dadi 2.3.3)
#
# The dadi0116 env keeps `.original` backups next to both customised files, so the
# patches are a plain `diff -u original patched`. But a backup is only a valid patch
# base if it really is the pristine upstream file, and nothing in the env guarantees
# that, so step 2 downloads the two sdists from PyPI and requires byte equality
# before the diffs are taken. If that check ever fails, regenerate the patches
# against the sdists rather than against the backups.
#
# Requires: passwordless ssh to the host alias `hoffman2`, network access to PyPI.
set -euo pipefail

# Path relative to the remote home directory, so no cluster absolute path is
# baked into the package (see tests/test_config.py).
b=project-klohmuel/noncdoingdfe/condaenv/dadi0116/lib/python3.12/site-packages

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# -- 1. pull the customised files and the env's own .original backups -----------
# NB the local names must not differ only in case: dadi_cli/Pdfs.py and
# dadi/DFE/PDFs.py would silently overwrite each other on a case-insensitive
# filesystem (macOS), so they are prefixed cli_ and dfe_ here.
ssh hoffman2 "cat $b/dadi_cli/Pdfs.py"          > "$work/cli_Pdfs.py"
ssh hoffman2 "cat $b/dadi_cli/Pdfs.py.original" > "$work/cli_Pdfs.py.original"
ssh hoffman2 "cat $b/dadi/DFE/PDFs.py"          > "$work/dfe_PDFs.py"
ssh hoffman2 "cat $b/dadi/DFE/PDFs.py.original" > "$work/dfe_PDFs.py.original"
for f in cli_Pdfs.py cli_Pdfs.py.original dfe_PDFs.py dfe_PDFs.py.original; do
  [[ -s "$work/$f" ]] || { echo "fetch failed or empty: $f" >&2; exit 1; }
done
# The env's dist-info, for the record:
ssh hoffman2 "ls -d $b/dadi-*.dist-info $b/dadi_cli-*.dist-info" | xargs -n1 basename

# -- 2. prove the .original backups are pristine dadi 2.3.3 / dadi-cli 0.9.5 ---
(
  cd "$work"
  pip3 download dadi-cli==0.9.5 --no-deps --no-binary :all: -d . >/dev/null
  sdist=$(curl -fsS https://pypi.org/pypi/dadi/2.3.3/json |
    python3 -c 'import json,sys; print([u["url"] for u in json.load(sys.stdin)["urls"] if u["packagetype"]=="sdist"][0])')
  curl -fsSL "$sdist" -o dadi-2.3.3.tar.gz
  tar xzf dadi_cli-0.9.5.tar.gz
  tar xzf dadi-2.3.3.tar.gz
)
cmp "$work/dadi_cli-0.9.5/dadi_cli/Pdfs.py" "$work/cli_Pdfs.py.original"
cmp "$work/dadi-2.3.3/dadi/DFE/PDFs.py"     "$work/dfe_PDFs.py.original"
echo "verified: both .original backups are byte-identical to pristine upstream"

# -- 3. write the patches ------------------------------------------------------
# -L gives the patch stable, package-relative headers instead of mktemp paths.
# diff exits 1 when the files differ, which is the expected outcome here.
mkdir -p setup/patches
diff -u -L dadi_cli/Pdfs.py.original -L dadi_cli/Pdfs.py \
  "$work/cli_Pdfs.py.original" "$work/cli_Pdfs.py" > setup/patches/dadi_cli_Pdfs.patch || true
diff -u -L dadi/DFE/PDFs.py.original -L dadi/DFE/PDFs.py \
  "$work/dfe_PDFs.py.original" "$work/dfe_PDFs.py" > setup/patches/dadi_DFE_PDFs.patch || true

# -- 4. verify the patches -----------------------------------------------------
for p in setup/patches/dadi_cli_Pdfs.patch setup/patches/dadi_DFE_PDFs.patch; do
  [[ -s "$p" ]] || { echo "empty patch: $p" >&2; exit 1; }
  grep -q neugamma "$p" || { echo "no neugamma in $p" >&2; exit 1; }
done
echo "patches mentioning neugammazero:"
grep -l neugammazero setup/patches/*.patch

# Apply to a throwaway copy of pristine upstream: this is exactly what setup.sh
# will do to a freshly created env.
patch --dry-run -s -p0 "$work/dadi_cli-0.9.5/dadi_cli/Pdfs.py" < setup/patches/dadi_cli_Pdfs.patch
patch --dry-run -s -p0 "$work/dadi-2.3.3/dadi/DFE/PDFs.py"     < setup/patches/dadi_DFE_PDFs.patch
echo "both patches apply cleanly to pristine upstream"

wc -l setup/patches/*.patch
