#!/bin/bash
# 2026/09/29  Every cluster action taken while building and validating the dfe
# stage, in the order it was taken. Not meant to be run end to end: it is the
# record the plan's Global Constraints ask for. Run from the package root on a
# machine with ssh access to hoffman2.
#
# Validation bed: phastCons17way_0-5_nomiss, against the published results in
#   $REFS_DIR/result_dadi/has_qui.exHET_noq.phastCons17way_0-5_nomiss_qui2.chrs_1t22_{neu,gamma,neugamma}
# Nothing under $REFS_DIR was written except the package directory itself and
# $OUT; the reference result directories were only read from.
set -euo pipefail

pkg=project-klohmuel/noncdoingdfe/test_0118/dfepipe

# --- 1. push the new stage files -------------------------------------------
rsync -avc --no-perms --no-times \
  scripts/dfe/run_fitneu.sh scripts/dfe/run_dfe.sh scripts/dfe/run_bestdfe.sh \
  "hoffman2:${pkg}/scripts/dfe/"
rsync -avc --no-perms --no-times \
  scripts/sge/submit_dfe.sh scripts/sge/worker_dfe.sh "hoffman2:${pkg}/scripts/sge/"
rsync -avc --no-perms --no-times scripts/checks/check_dfe.py "hoffman2:${pkg}/scripts/checks/"
rsync -avc --no-perms --no-times tests/test_check_dfe.py "hoffman2:${pkg}/tests/"

# --- the run keys every step below used ------------------------------------
read -r -d '' env_block <<'ENV' || true
source config/profiles/hoffman_replicate.conf
export PREP_RUN=prep_260929 FS_RUN=fs_260929 DEMO_RUN=demo_260929
export CACHE_RUN=cache_260929 DFE_RUN=dfe_260929 SUMMARY_RUN=sum_260929
# The cache is the user's existing one, not one this pipeline built.
export DFE_CACHE_FILE=$REFS_DIR/result_dadi/quiecent_ctrl.exHET_noq.phastCons17way_0-5_nomiss_qui2.chrs_1t22-demo_three_epoch-1/cache_1.spectra.bpkl
ENV

# --- 2. the neutral null ----------------------------------------------------
# Result: ll = -28858.717221247774, mutratio = 0.10908421521059425
# Published: ll = -28858.706648022526, same mutratio. The 0.011 ll difference is
# our own demographic fit (theta 503196.60542 vs 503196.65647), not the port.
ssh hoffman2 "cd ${pkg} && ${env_block}
  bash scripts/dfe/run_fitneu.sh phastCons17way_0-5_nomiss.bed"

# --- 3. the command packs, and a diff against the original's -----------------
# 190 gamma starting points in 4 packs, 210 neugamma in 5: the same counts the
# published run produced (record_251006_repeat_phyloP243_0-5.sh submitted
# packs 0/50/100/150 for gamma and 0/50/100/150/200 for neugamma).
ssh hoffman2 "cd ${pkg} && ${env_block}
  DFE_GENERATE_ONLY=1 bash scripts/dfe/run_dfe.sh gamma phastCons17way_0-5_nomiss.bed
  sed -n 29p \$OUT/dfe/phastCons17way_0-5_nomiss/cache_260929-gamma-dfe_260929/dadi-cli_dfe.150.sh"
# Compare with the original's line for the same starting point:
#   sed -n 29p \$REFS_DIR/result_dadi/has_qui.exHET_noq.phastCons17way_0-5_nomiss_qui2.chrs_1t22_gamma/dadi-cli_dfe.150.sh
# Identical but for the paths, the demographic .bestfits key (109 vs 29, our
# demo run's winner) and `--cpus` spelled in full instead of `--cpu`.

# --- 4. one fit by hand, the starting point that won the published run --------
# Result: the converged block's last row, ll -610.1732410661816,
# shape 0.16902608292081361, scale 65.56182508865463, theta 54890.806799339436.
ssh hoffman2 "cd ${pkg} && ${env_block}
  source /u/local/Modules/default/init/bash; module load conda; conda activate \$CONDA_ENV
  o=\$OUT/dfe/phastCons17way_0-5_nomiss/cache_260929-gamma-dfe_260929
  bash -c \"\$(sed -n 29p \$o/dadi-cli_dfe.150.sh)\""

# --- 5. the full arrays -----------------------------------------------------
# 4 gamma tasks + 5 neugamma tasks, one per pack. SGE_DFE_RT cut to 2h while
# developing; the default is 23:50:00.
ssh hoffman2 "cd ${pkg} && ${env_block}
  export SGE_DFE_RT=2:00:00
  FORCE=1 bash scripts/sge/submit_dfe.sh gamma phastCons17way_0-5_nomiss.bed
  bash scripts/sge/submit_dfe.sh neugamma phastCons17way_0-5_nomiss.bed"

# --- 6. run_bestdfe.sh against the PUBLISHED fits ----------------------------
# The strongest check of the collection / best-pick / discretisation chain: the
# reference .InferDFE.bestfits are copied into a run of their own (DFE_RUN=
# refcopy) and run_bestdfe.sh must rebuild the published large_ll.txt and DFE.
# Result: identical filename / Log(likelihood) / shape / scale / theta columns
# and the same row counts for both pdfs; scale_div2Na, sdadi, Ne and the DFE
# proportions differ only in the 8th significant figure, from our Na
# (18460.05516977379) against theirs (18460.05704253424).
ssh hoffman2 "cd ${pkg} && source config/profiles/hoffman_replicate.conf
  for pdf in gamma neugamma; do
    src=\$REFS_DIR/result_dadi/has_qui.exHET_noq.phastCons17way_0-5_nomiss_qui2.chrs_1t22_\$pdf
    dst=\$OUT/dfe/phastCons17way_0-5_nomiss/cache_260929-\$pdf-refcopy
    mkdir -p \$dst && cp \$src/*.InferDFE.bestfits \$dst/
  done"
ssh hoffman2 "cd ${pkg} && ${env_block}
  export DFE_RUN=refcopy
  for pdf in gamma neugamma; do
    bash scripts/dfe/run_bestdfe.sh \$pdf phastCons17way_0-5_nomiss.bed
  done"
