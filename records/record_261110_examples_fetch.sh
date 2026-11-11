#!/bin/bash
# 2026-11-10  Fetch the small "expected output" files for examples/.
#
# Read-only on the cluster. It copies down only the files small enough to live in
# the repository - check files, the best-fit parameters, the discrete DFE and its
# csv, and two spectra plots. Beds, VCFs, caches and the per-fit InferDFE files
# stay where they are: the ELS bed alone is ~25 MB and a cache is larger.
#
# The run keys are part of the paths, and each stage's directory name carries its
# upstream key as well as its own, so the paths below are the provenance:
#
#   ELS    prep_261025 / fs_261025 / demo_261025 / cache_261025 / dfe_261025 / sum_261025
#   probe  (no prep)   / fs_pp261026 / demo_pp261026 / cache_pp261026 / dfe_pp261026 / sum_pp261026
#
# The probe ran FROM=sfs, so it has no check_prep. Its shared check files are
# rewritten by any later run that reuses fs_pp261026 - the 39-bed run submitted
# on 2026-11-10 does - so what lands here is a snapshot of the probe-only state,
# taken while that run was still queued.
#
# Only encode_els is published under examples/. The conservation-bin block is
# kept because it is the one that was fetched and checked on 2026-11-10, and it
# is the shape any later example takes:
#
#   bash records/record_261110_examples_fetch.sh
#   EXAMPLES="encode_els conservation_bin" bash records/...
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

HOST="${HOST:-hoffman2}"
EXAMPLES="${EXAMPLES:-encode_els}"
want () { printf '%s\n' $EXAMPLES | grep -qx "$1"; }
O="project-klohmuel/noncdoingdfe/test_0118/output_dfepipe"

get () {   # get <remote-relative-path> <local-path>
  local src="$1" dst="$2"
  mkdir -p "$(dirname "$dst")"
  if scp -q "${HOST}:${O}/${src}" "$dst" 2>/dev/null; then
    echo "  $(printf '%7s' "$(wc -c < "$dst" | tr -d ' ')")  $dst"
  else
    echo "  MISSING  ${src}" >&2
    return 1
  fi
}

# ------------------------------------------------------------------- ELS ------
if want encode_els; then
b=GRCh38-ELS
d=examples/encode_els/expected
echo "examples/encode_els  <- ${b}, run keys *_261025"
get "prep/check_prep.prep_261025.tsv"                              "$d/prep/check_prep.prep_261025.tsv"
get "sfs/check_sfs.fs_261025.tsv"                                  "$d/sfs/check_sfs.fs_261025.tsv"
get "demo/check_demo.fs_261025-three_epoch-demo_261025.tsv"        "$d/demo/check_demo.fs_261025-three_epoch-demo_261025.tsv"
get "demo/${b}/fs_261025-three_epoch-demo_261025/best_inference.txt" "$d/demo/best_inference.txt"
get "demo/${b}/fs_261025-three_epoch-demo_261025/sfs_bestfit.pdf"  "$d/demo/sfs_bestfit.pdf"
get "cache/${b}/demo_261025-1-cache_261025/check_cache.tsv"        "$d/cache/check_cache.tsv"
get "dfe/check_dfe.cache_261025-1-gamma-dfe_261025.tsv"            "$d/dfe/check_dfe.cache_261025-1-gamma-dfe_261025.tsv"
g="dfe/${b}/cache_261025-1-gamma-dfe_261025"
n="dfe_exHET_noq.GRCh38-ELS_qui2.chrs_1t22_gamma"
get "${g}/${n}.pdf"            "$d/dfe/${n}.pdf"
get "${g}/${n}.csv"            "$d/dfe/${n}.csv"
get "${g}/sfs_rank1_poisson.pdf" "$d/dfe/sfs_rank1_poisson.pdf"
get "${g}/spec.bestdfe.txt"      "$d/dfe/spec.bestdfe.txt"
get "summary/check_summary.sum_261025.tsv"                         "$d/summary/check_summary.sum_261025.tsv"
get "summary/sum_261025/bestpdf_sumdfe.tsv"                        "$d/summary/bestpdf_sumdfe.tsv"
fi

# --------------------------------------------------- the conservation bin -----
if want conservation_bin; then
b=phyloP243way_0-5_phyloP447way_40-100
d=examples/conservation_bin/expected
echo "examples/conservation_bin  <- ${b}, run keys *_pp261026"
get "sfs/check_sfs.fs_pp261026.tsv"                                     "$d/sfs/check_sfs.fs_pp261026.tsv"
get "demo/check_demo.fs_pp261026-three_epoch-demo_pp261026.tsv"         "$d/demo/check_demo.fs_pp261026-three_epoch-demo_pp261026.tsv"
get "demo/${b}/fs_pp261026-three_epoch-demo_pp261026/best_inference.txt" "$d/demo/best_inference.txt"
get "demo/${b}/fs_pp261026-three_epoch-demo_pp261026/sfs_bestfit.pdf"   "$d/demo/sfs_bestfit.pdf"
get "cache/${b}/demo_pp261026-1-cache_pp261026/check_cache.tsv"         "$d/cache/check_cache.tsv"
get "dfe/check_dfe.cache_pp261026-1-gamma-dfe_pp261026.tsv"             "$d/dfe/check_dfe.cache_pp261026-1-gamma-dfe_pp261026.tsv"
g="dfe/${b}/cache_pp261026-1-gamma-dfe_pp261026"
n="dfe_exHET_noq.${b}_qui2.chrs_1t22_gamma"
get "${g}/${n}.pdf"              "$d/dfe/${n}.pdf"
get "${g}/${n}.csv"              "$d/dfe/${n}.csv"
get "${g}/sfs_rank1_poisson.pdf" "$d/dfe/sfs_rank1_poisson.pdf"
get "${g}/spec.bestdfe.txt"      "$d/dfe/spec.bestdfe.txt"
get "summary/check_summary.sum_pp261026.tsv"                            "$d/summary/check_summary.sum_pp261026.tsv"
get "summary/sum_pp261026/bestpdf_sumdfe.tsv"                           "$d/summary/bestpdf_sumdfe.tsv"
fi

echo
echo "total: $(find examples -type f | wc -l | tr -d ' ') files, $(du -sh examples | cut -f1)"
