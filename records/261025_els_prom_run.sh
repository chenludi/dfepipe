#!/bin/bash
# 2026-10-25  ELS, ELS ∩ Prom, ELS − Prom. Not a script to run: the record of a run.
#
# Three beds in one chain, so each shares nothing but the reference data: every
# bed gets its own control set, demography, cache and DFE.
#
#   Prom = TSS + PromF + BivProm (ChromHMM), as the supplement defines it.
#   Inputs are the `old` ENCODE release (the manuscript's) at the `raw` level,
#   so prep applies HET/Quies/CDS/mask identically to all three. The two parts
#   partition ELS exactly: 32,093,548 + 230,309,010 = 262,402,558 bp, diff 0.
#   Built by ../../record_260928_ELS-ChromHMMProm_length.sh, not rebuilt here.
#
# check_bed: WARN on all three, both reasons benign - chrX/chrY present (2.5% of
# bp, dropped at concat since CHROMS is chr1-22, same as the published ELS run),
# and 66 unsorted intervals in the raw ENCODE bed, which prep sorts anyway.
# No interval exceeds 2 x CTRL_FLANK (max 917 bp), so the long-interval control
# window question (todo/tocheck.md item 1) does not arise for these beds.

stage=$REFS_DIR/data_bed/els_prom_261025     # GRCh38-ELS{,_x_Prom,_not_Prom}.bed

: '
cd ~/project-klohmuel/noncdoingdfe/test_0118/dfepipe
source config/profiles/hoffman_replicate.conf
export PREP_RUN=prep_261025 FS_RUN=fs_261025 DEMO_RUN=demo_261025 \
       CACHE_RUN=cache_261025 DFE_RUN=dfe_261025 SUMMARY_RUN=sum_261025
export SGE_SFS_RT=8:00:00
d=$REFS_DIR/data_bed/els_prom_261025
bash scripts/chains/run_all.sh $d/GRCh38-ELS.bed $d/GRCh38-ELS_x_Prom.bed $d/GRCh38-ELS_not_Prom.bed
'

# Submitted twice. The first attempt, 22:12, was abandoned at 22:30 - not
# because anything failed, but because Chenlu pointed out that running the beds
# one after another inside a task is wasted wall clock. It was:
#
#   prep 15278009, sfs 15278010 (1-22), concat ...011, demo ...012 (1-3),
#   bestdemo ...013, launch:cache ...014
#
# sfs needed SGE_SFS_RT raised to 8h by hand and demo needed
# `qalter -l h_data=3G,h_rt=8:00:00 15278012` before it ran, both because one
# task did all three beds. qalter also reported "modified environment", so the
# env was checked afterwards: DEMO_P0_FILES, DEMO_P0_DIR, REPO_DIR, the run
# keys, GATE_FILES and CONDA_ENV all survived, DEMO_P0_FILES with its spaces
# intact - which settles an open question, qsub -v does preserve spaces here.
#
# The arrays were then split over the product, outputs deleted, and the whole
# thing rerun from prep. The defect: sfs and demo both looped over the beds
# inside one task, so a task's cost
# scaled with the bed count while its h_rt stayed one number for the stage.
# Splitting needed no change to the output layout, because every file was
# already keyed by bed and by chromosome or marker.
#
# phase 1, submitted 22:35, on the split arrays
#   prep           15278136   one job, all 3 beds (not an array: it writes the
#                             stage's check file itself, so tasks would race)
#   sfs            15278137   1-66   = 22 chromosomes x 3 beds
#   concat         15278138
#   demo           15278139   1-9    = 3 p0 files x 3 beds
#   bestdemo       15278140
#   launch:cache   15278141
#
# Both arrays run on the shipped defaults now - sfs h_rt 3:40, demo 3:10 - with
# no override, because a task is one bed again. That is the test that the split
# did what it was meant to.
#
# Expected: sfs tasks about as long as run 1's 14 min, demo about 1h44 rather
# than 5h12, so phase 1 lands in roughly two hours instead of six.
#
# 12 .sfs.chr* scratch directories were left behind when the first attempt was
# qdel'd, confirming that the EXIT trap in run_sfs_chrom.sh does not fire on
# SIGKILL. They went with the deleted output tree.
#
# Difference from the published ELS, known in advance: run_prep.sh:112 merges,
# and the published ELS predates the 2024-11-04 merge. Measured on ELS at -26 bp
# in 202.7 Mb, 1.3e-7. todo/tocheck.md item 1.
