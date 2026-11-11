#!/bin/bash
# 2026-10-26  ENCODE v4: PLS, ELS, ELS ∩ Prom, ELS − Prom.
#
# Run this on the cluster, from the package root. It stages the beds, checks
# them, and dry-runs the chain. It submits nothing unless SUBMIT=1.
#
#   bash records/record_261026_els_upd_run.sh              stage, check, dry-run
#   SUBMIT=1 bash records/record_261026_els_upd_run.sh     the same, then submit
#
# Everything it does is echoed and appended to the log named at the end, so what
# was run on the cluster is recoverable afterwards.
#
# The four beds, all ENCODE v4. Both source files were verified byte-identical
# to Chenlu's own downloads on 2026-10-25:
#   GRCh38-cCREs.PLS.bed  sha256 65443e5ace6db2f1... 47,532 int    13,461,873 bp
#   GRCh38-cCREs.ELS.bed  sha256 e09eb2709fd686f3... 1,718,669 int 458,253,930 bp
# The ELS file pools dELS and pELS (column 6), as the old GRCh38-ELS.bed did.
#
# Prom = the union of the ChromHMM groups PromF, TSS and BivProm, which is how
# the supplement defines promoters. The two ELS parts were built by
# ../../record_260928_ELS-ChromHMMProm_length.sh on 2026-09-28 and are not
# rebuilt here; the script re-checks that they still partition ELS exactly.
#
# The staged beds are copies. Each source's sha256 is recorded in the log at
# staging time, so which bytes a run used stays recoverable either way.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source config/profiles/hoffman_replicate.conf
source config/pipeline.conf

export PREP_RUN=prep_261026 FS_RUN=fs_261026 DEMO_RUN=demo_261026 \
       CACHE_RUN=cache_261026 DFE_RUN=dfe_261026 SUMMARY_RUN=sum_261026

enc="$HOME/project-klohmuel/noncdoingdfe/data/encode/ENCODEupdate_202601"
drv="$REFS_DIR/result_small/mv_260928_ELS-ChromHMMProm"
stage="$REFS_DIR/data_bed/els_prom_upd_261026"
log="$REFS_DIR/submit_261026.log"

say () { echo "$@" | tee -a "$log"; }

say "=============================================================="
say " ENCODE v4 run   $(date '+%Y-%m-%d %H:%M:%S')   host=$(hostname)"
say " run keys: ${PREP_RUN} ${FS_RUN} ${DEMO_RUN} ${CACHE_RUN} ${DFE_RUN} ${SUMMARY_RUN}"
say "=============================================================="

# ------------------------------------------------------------------ 1. stage --
# Copied, not symlinked, so this run's inputs are a snapshot. The two derived
# beds live in an output directory and re-running
# record_260928_ELS-ChromHMMProm_length.sh overwrites them in place; a symlink
# would follow that silently, with spec.txt still recording the same path and
# nothing to notice the change by. The source path, size, mtime and sha256 of
# each go in the log, which is what a symlink would have shown in `ls -l`.
mkdir -p "$stage"
copy () {  # copy <src> <dest-name>
  [[ -s "$1" ]] || { say "MISSING source: $1"; exit 1; }
  cp -f "$1" "${stage}/$2"
  say "  $2"
  say "      <- $1"
  say "      $(stat -c '%s bytes, mtime %y' "$1" | cut -d. -f1)"
  say "      sha256 $(sha256sum "$1" | cut -c1-32)..."
}
say
say "1. staging into ${stage}  (copies)"
copy "${enc}/GRCh38-cCREs.PLS.bed"      GRCh38-PLSv4.bed
copy "${enc}/GRCh38-cCREs.ELS.bed"      GRCh38-ELSv4.bed
copy "${drv}/upd.raw.ELS_x_Prom.bed"    GRCh38-ELSv4_x_Prom.bed
copy "${drv}/upd.raw.ELS_not_Prom.bed"  GRCh38-ELSv4_not_Prom.bed

# --------------------------------------------------- 2. the parts must add up --
# If the two parts no longer partition ELS, the derived beds are stale relative
# to the ELS file beside them and the contrast would not mean what it says.
say
say "2. do the two parts still partition ELS?"
bp () { awk '!/^(#|track|browser)/ && NF>=3 {s+=$3-$2} END{printf "%d", s+0}' "$1"; }
a=$(bp "${stage}/GRCh38-ELSv4_x_Prom.bed")
b=$(bp "${stage}/GRCh38-ELSv4_not_Prom.bed")
t=$(bp "${stage}/GRCh38-ELSv4.bed")
say "   ELS_x_Prom    ${a}"
say "   ELS_not_Prom  ${b}"
say "   sum           $((a + b))"
say "   ELS           ${t}"
if [[ $((a + b)) -ne "$t" ]]; then
  say "   MISMATCH by $((a + b - t)) bp - the parts do not partition ELS. Stopping."
  say "   Rebuild them with record_260928_ELS-ChromHMMProm_length.sh first."
  exit 1
fi
say "   difference    0   ok"

beds=("${stage}/GRCh38-PLSv4.bed"
      "${stage}/GRCh38-ELSv4.bed"
      "${stage}/GRCh38-ELSv4_x_Prom.bed"
      "${stage}/GRCh38-ELSv4_not_Prom.bed")

# ---------------------------------------------------------------- 3. preflight --
say
say "3. check_bed on all four"
python3 scripts/checks/check_bed.py "${beds[@]}" 2>&1 | tee -a "$log"

# ----------------------------------------------------------------- 4. dry run --
say
say "4. dry run"
DRY_RUN=1 bash scripts/chains/run_all.sh "${beds[@]}" 2>&1 | tail -20 | tee -a "$log"

# ------------------------------------------------------------------ 5. submit --
if [[ "${SUBMIT:-0}" != "1" ]]; then
  say
  say "Nothing was submitted. To submit:"
  say "   SUBMIT=1 bash records/record_261026_els_upd_run.sh"
  exit 0
fi

say
say "5. submitting"
bash scripts/chains/run_all.sh "${beds[@]}" 2>&1 | tee -a "$log"

say
say "job ids as the queue sees them:"
qstat -u "$USER" 2>&1 | tee -a "$log"
say
say "logged to ${log}"
say "the table will be \$OUT/summary/${SUMMARY_RUN}/bestpdf_sumdfe.tsv"
