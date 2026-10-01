#!/bin/bash
# Copy the package to a cluster. Nothing else; it runs no pipeline stage.
#
#   bash sync/push.sh <host>:<path>        push the package
#   DRY_RUN=1 bash sync/push.sh <host>:<path>   show what would change
#   CHECK=1   bash sync/push.sh <host>:<path>   compare only, change nothing
#
#   bash sync/push.sh hoffman2:project-klohmuel/noncdoingdfe/test_0118/dfepipe
#
# This exists because the package is meant to run on a cluster and every user
# was otherwise inventing their own rsync, with their own guess at what to
# exclude. Getting that wrong is not harmless: a stale copy means the behaviour
# you reasoned about locally is not the behaviour that runs, and the one time it
# happened here a run of a half-updated package wasted an array.
#
# Checksums, not timestamps (-c): the mtimes of a fresh checkout mean nothing,
# and a stage that silently did not update is the failure this guards against.
# Permissions and times are left alone, because cluster filesystems often refuse
# them and the package does not care.
#
# Output is NOT synced. Results live under $OUT on the cluster and are never
# copied back by this script - they are large, and overwriting a result with a
# local stale copy would be the worst thing this could do. For the same reason
# --delete is safe here: it only ever touches the package directory.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(dirname "$here")"

if [[ $# -ne 1 ]]; then
  echo "usage: push.sh <host>:<path>" >&2
  echo "       e.g. push.sh hoffman2:project-klohmuel/noncdoingdfe/test_0118/dfepipe" >&2
  exit 2
fi
dest="$1"
[[ "$dest" == *:* ]] || { echo "push.sh: destination needs <host>:<path>" >&2; exit 2; }

# Everything the package needs to run, and nothing that is generated.
EXCLUDES=(
  --exclude '.git/'            # history is not needed to run
  --exclude '__pycache__/'
  --exclude '.pytest_cache/'
  --exclude '*.pyc'
  --exclude '.DS_Store'
  --exclude 'output/'          # should never exist here, but be sure
)

# --delete matters as much as the copy: a script deleted here must not linger
# there, or it can still be sourced. The destination is the package's own
# directory, and output lives under $OUT elsewhere, so nothing else is at risk.
# --delete-excluded is deliberately NOT set: the excluded paths (.git, caches)
# are left alone if the remote happens to have them.
flags=(-rlc --delete --no-perms --no-times --itemize-changes)
[[ "${DRY_RUN:-0}" == "1" || "${CHECK:-0}" == "1" ]] && flags+=(--dry-run)

echo "push.sh: ${PKG} -> ${dest}"
[[ "${DRY_RUN:-0}" == "1" ]] && echo "push.sh: DRY_RUN=1, nothing will be written"
[[ "${CHECK:-0}" == "1" ]] && echo "push.sh: CHECK=1, comparing only"

# The trailing slash matters: contents of PKG into dest, not PKG as a subdir.
out="$(rsync "${flags[@]}" "${EXCLUDES[@]}" "${PKG}/" "${dest}/")"
printf '%s\n' "$out" | sed 's/^/  /'

n="$(printf '%s\n' "$out" | grep -c '^[<>ch.]' || true)"
if [[ "${CHECK:-0}" == "1" ]]; then
  if [[ "$n" -eq 0 ]]; then
    echo "push.sh: the remote copy matches."
  else
    echo "push.sh: ${n} file(s) differ. Push before you submit anything." >&2
    exit 1
  fi
  exit 0
fi
echo "push.sh: ${n} file(s) transferred."

# The one thing worth verifying after a push: that the entry point and the config
# really are what is here, since everything else is reasoned about from them.
if [[ "${DRY_RUN:-0}" != "1" ]]; then
  host="${dest%%:*}"; path="${dest#*:}"
  for f in config/pipeline.conf scripts/chains/run_all.sh; do
    l="$(shasum -a256 "${PKG}/${f}" 2>/dev/null | cut -c1-16)"
    r="$(ssh -o BatchMode=yes "$host" "sha256sum '${path}/${f}' 2>/dev/null | cut -c1-16" 2>/dev/null | tr -d '\r')"
    if [[ -n "$r" && "$l" == "$r" ]]; then
      printf '  verified  %s\n' "$f"
    else
      printf '  MISMATCH  %s (local %s, remote %s)\n' "$f" "$l" "${r:-absent}" >&2
      exit 1
    fi
  done
fi
