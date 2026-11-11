#!/bin/bash
# Copy the package to a cluster. Nothing else; it runs no pipeline stage.
#
#   bash sync/push.sh <host>:<path>        push the package
#   DRY_RUN=1 bash sync/push.sh <host>:<path>   show what would change
#   CHECK=1   bash sync/push.sh <host>:<path>   compare only, change nothing
#   ALLOW_RUNNING=1 bash sync/push.sh ...       push even with jobs in the queue
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
# A push is refused while you have jobs in the queue on that host, because a
# running stage reads its scripts from this directory as it goes. That is not
# hypothetical: pushing run_cache.sh twelve minutes into a cache job left the job
# running the pre-push copy, so its checkpoint never ran, no check file was
# written, and the gate then stopped bestdfe, summary and bestpdf. The cache
# itself was fine; three stages stopped because of the push. ALLOW_RUNNING=1
# overrides, for when the queued jobs are unrelated to the package.
#
# Output is not synced. Results live under $OUT on the cluster and are never
# copied back by this script - they are large, and overwriting a result with a
# local stale copy would be the worst thing this could do. For the same reason
# --delete is safe here: it only ever touches the package directory.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(dirname "$here")"

# DFEPIPE_REMOTE is the default destination, so a site profile can record it
# once instead of everyone retyping it. An argument overrides it.
dest="${1:-${DFEPIPE_REMOTE:-}}"
if [[ $# -gt 1 || -z "$dest" ]]; then
  echo "usage: push.sh [<host>:<path>]" >&2
  echo "       with no argument, uses \$DFEPIPE_REMOTE from your site profile" >&2
  echo "       e.g. push.sh hoffman2:project-klohmuel/noncdoingdfe/test_0118/dfepipe" >&2
  exit 2
fi
[[ "$dest" == *:* ]] || { echo "push.sh: destination needs <host>:<path>" >&2; exit 2; }

# ---------------------------------------------------- refuse to push mid-run --
# Only for a real push: CHECK and DRY_RUN write nothing, so they are always safe.
if [[ "${DRY_RUN:-0}" != "1" && "${CHECK:-0}" != "1" \
      && "${ALLOW_RUNNING:-0}" != "1" ]]; then
  _host="${dest%%:*}"
  # Found the same way config/env.sh:52-71 finds it, because a non-interactive
  # ssh does not get the login shell's PATH and qstat is usually not on it.
  # The remote snippet writes no temp file on purpose: $TMPDIR and /tmp are not
  # writable on every login node (they are not on hoffman2), and a probe that
  # fails for that reason would be indistinguishable from an idle queue.
  _q="$(ssh -o BatchMode=yes -o ConnectTimeout=30 "$_host" '
      export SGE_ROOT="${SGE_ROOT:-/u/systems/UGE8.6.4}"
      export SGE_CELL="${SGE_CELL:-hoffman2}"
      for c in "$(command -v qstat 2>/dev/null)" \
               "${SGE_ROOT}/bin/lx-amd64/qstat" /u/local/bin/qstat; do
        [ -n "$c" ] && [ -x "$c" ] || continue
        out=$("$c" -u "$USER" 2>/dev/null); rc=$?
        # A blank qstat means an empty queue only when it also succeeded. A
        # failed qstat prints nothing either, and reading that as "idle" is the
        # mistake this whole block exists to prevent.
        [ $rc -eq 0 ] || continue
        printf "%s" "$out" | awk "NR>2" | grep -c . || true
        exit 0
      done
      echo unknown' 2>/dev/null | tail -1 | tr -d "[:space:]" || true)"
  # `|| true` so an ssh that cannot connect at all leaves _q empty and falls into
  # the message below. Without it, set -e aborts here with ssh's own status and
  # no explanation, which is the least useful thing this script could do.
  if [[ "$_q" == "unknown" || -z "$_q" ]]; then
    echo "push.sh: could not ask ${_host} whether you have jobs running." >&2
    echo "         qstat was not runnable there, and a push during a run can" >&2
    echo "         make a job read half-updated scripts. Check by hand, then" >&2
    echo "         ALLOW_RUNNING=1 bash sync/push.sh ${dest}" >&2
    exit 1
  fi
  if [[ "$_q" != "0" ]]; then
    echo "push.sh: ${_q} job(s) of yours are in the queue on ${_host}." >&2
    echo "         A running stage reads its scripts from the destination as it" >&2
    echo "         goes, so pushing now can leave it running a mixture. Wait for" >&2
    echo "         them, or ALLOW_RUNNING=1 bash sync/push.sh ${dest}" >&2
    exit 1
  fi
fi

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
# --delete-excluded is deliberately not set: the excluded paths (.git, caches)
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
  # A leading ~ has to reach the remote shell unquoted to expand, but the path
  # still needs quoting against spaces. Hand the remote $HOME instead and let it
  # expand that inside double quotes. Single-quoting the whole path, as this did
  # until 2026-10-25, made every push to a ~ destination report "remote absent" -
  # a verification that always fails is one you learn to ignore.
  rpath="${path}"
  case "$rpath" in
    "~/"*) rpath="\$HOME/${rpath#\~/}" ;;
    "~")   rpath="\$HOME" ;;
  esac
  bad=0
  for f in config/pipeline.conf scripts/chains/run_all.sh; do
    l="$(shasum -a256 "${PKG}/${f}" 2>/dev/null | cut -c1-16)"
    r="$(ssh -o BatchMode=yes "$host" "sha256sum \"${rpath}/${f}\" 2>/dev/null | cut -c1-16" 2>/dev/null | tr -d '\r')"
    if [[ -n "$r" && "$l" == "$r" ]]; then
      printf '  verified  %s\n' "$f"
    else
      printf '  MISMATCH  %s (local %s, remote %s)\n' "$f" "$l" "${r:-absent}" >&2
      bad=1
    fi
  done
  [[ "$bad" -eq 0 ]] || exit 1
fi
