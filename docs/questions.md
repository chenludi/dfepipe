# Findings from the port

What was found while porting the original scripts, kept because it explains why
this package behaves as it does. Nothing here was acted on in the originals: no
original script and no manuscript file was changed.

Item numbers are stable - `docs/PIPELINE.md` and `scripts/sfs/run_sfs_chrom.sh`
cite them - so gaps in the numbering are items that were resolved and removed.

## The one deliberate numerical difference

**The original's reported "best fit" is the worst converged row of the winning
file.** `plot_converge_dfe_0430.py:282-321` does `data_dict.update({filename:
...})` for every row of the `# Converged results` block without breaking, so the
last row survives, and dadi-cli writes that block best-first. This package takes
the best row instead.

The cost in the case measured was 7e-9 in log-likelihood and the sixth decimal
of shape, but it is systematic across every published fit and unbounded for a
poorly converged file. How far it moved published results is open:
`todo/tocheck.md` item 2.

Two other scripts in the same pipeline already take the first row
(`dadi_cal_ll_neu_0723.py:112-127`, `make_p0_dfe_0430.py:206-221`).

## Errors in the original scripts

Recorded, not ported. Still true unless marked otherwise.

6. `dadi-cli_demog_fsfile.sh`: the `scripts/` copy builds `out_demog` from
   variables assigned later (line 32 vs 38/43), so its files would lack the
   marker `max_ll_key` needs; the cluster copy has the `dadi-cli InferDM` call
   deleted. Both edits postdate the published results.
7. `--output` / `--cpu` / `--sample-size` are argparse abbreviations, not real
   options (`--output-prefix`, `--cpus`, `--sample-sizes`).
8. `get_fs_mut_bychr.sh` writes temp files into the mutation-rate reference
   directory, so 22 array tasks can collide.
10. `piperef_stddfe_4_plotdemo.sh:52` `dadi-cli Plot` evidently never ran - no
    `*.bestfits.pdf` exists beside any published result.
11. Destructive paths: `rm -r $dirout` in steps 5 and 7; the second `.bak`
    destruction in `make_p0_dfe_0430.py:284-295`; the unconditional `.bak`
    clobber in `plot_DFE_single.py:132-156`; the shared `plots_dfe_mutload/`
    copy keyed only by bed and pdf.
12. Four "Error: ... does not exist" messages whose `exit 1` is commented out,
    so a missing input surfaces later as a traceback.
14. `plot_converge_dfe_0430.py:370` IndexErrors when two or fewer starting
    points converged; `:404` intends a `large_ll.txt` backup that never happens;
    `:494` rounds figure names to 3 decimals, so near-identical fits overwrite
    each other. Checked: nothing published was affected - the two thinnest published
    runs converged 3 and 4 starting points and both completed.
15. `piperef_stddfe_8_plotdfe.sh:69` hardcodes `cachenumber=1`, so step 8 can
    only read `cache_1.spectra.bpkl` whatever step 6 built.
16. The published beds came through `filt_nc-unmask.sh`, not through the chain
    the step-1 script suggests: `piperef_stddfe_1_rmHETqui.sh` stops after HET
    and quiescent and only counts lengths, while `submit_exHETquiencmask.sh`
    repeats those two subtracts and then subtracts CDS and intersects the mask.
    The `s.nc_` prefix on every published file is that script's output. This
    package reproduces those beds exactly, which is the evidence that the
    reading is right.

    The originals' HET and quiescent paths are stale - they name
    `chromHMM_100_states/temp_HET.bed`, which no longer exists - so a profile
    must point at the files that do. The four filters are set operations, so
    their order changes only the intermediate counts in `lengths.tsv`.
17. GNU-only: `stat -c %s` (`piperef_stddfe_2_concatefs.sh`), `sed -i`
    (`piperef_stddfe_9_sumdfe.sh`).
18. Checked and **not** a bug: `dadi_cal_ll_neu_0723.py:206-214` passes an
    unfolded model against folded data, but `dadi.Inference.ll_per_bin` folds
    the model itself. Verified identical to 1e-11.

## Manuscript items

Recorded, not changed.

19. **Control window.** The supplement says "100 kb window centered on the
    middle of the test sites"; every run used `mid +/- 25000`, i.e. 50 kb.
20. **Two published csv files for one bed came from different script versions.**
    `dfe_..._gamma.csv` has `0.0` in the `>=0.5` bin, the current CDF form;
    `dfe_..._neugamma.csv` has `9.61e-126`, exactly the commented-out `quad`
    integration, so it predates the 2025-03-25 switch. The summary table was
    assembled from pre-switch csv files. This package uses the CDF form, which
    matches the original's current code.
21. **Mixed provenance** in the published results, split by date:
    - windowing: unmerged before 2024-11-04 (ELS, EnhA, PromF, TSS, BivProm,
      PLS, the phastCons17 bins), merged after (PromF_TSS_BivProm, 2025 bins).
    - dadi: demography and cache under 2.3.0, DFE fits under 2.3.3.

    Both are about 1e-6 on everything measured.

## How much a DFE actually moves

The scale to judge any difference against. Spread across the 15 EnhA bootstrap
replicates (Data S2a, Fig. S14): neutral bin 0.0161, other bins 0.0076-0.0089.
The two validation runs reproduce the published bins to 8.2e-08 and 2.5e-05
relative, both far inside that spread.
