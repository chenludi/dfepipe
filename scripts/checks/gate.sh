#!/bin/bash
# gate.sh <upstream_check_file> <upstream_stage_name>
# Warns and ALWAYS exits 0. Re-running one stage against older output is a normal
# workflow here, so an upstream failure is reported, never enforced.
set -uo pipefail
f="${1:?usage: gate.sh <check_file> <stage>}"
stage="${2:?usage: gate.sh <check_file> <stage>}"

if [[ ! -s "$f" ]]; then
  echo "WARNING: no check file for upstream stage '${stage}' (${f}); proceeding unverified" >&2
  exit 0
fi
bad=$(awk -F'\t' 'NR>1 && $2=="FAIL" {printf "%s ", $1}' "$f")
if [[ -n "$bad" ]]; then
  echo "WARNING: upstream stage '${stage}' FAILED for: ${bad}" >&2
  echo "         see ${f}; proceeding anyway" >&2
fi
exit 0
