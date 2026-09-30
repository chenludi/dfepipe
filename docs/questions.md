# Questions for Chenlu

Things found while porting that need a decision or a check. Nothing here has
been acted on, and no original script or manuscript file has been changed.

## Decisions needed

1. **Submit the cache job?** ~10 h on 10 cores, `highp`. The DFE work does not
   need it — it uses your existing `cache_1.spectra.bpkl`. Running ours would
   validate the cache stage by comparison.
2. **`CACHE_MIN_SIZE_FRAC = 0.5`** — fail a cache below half its expected size
   (`gamma_pts x (projection+1) x 8`, ~40% pickle overhead). Is half the right
   margin? Measured: 300x161 -> 560 kB, 300x101 -> 335 kB, 1000x161 -> 1.8 MB.
3. **`HASQUI_MIN_FRAC = 0.90`**, not the 0.95 you suggested: ELS is 0.9105 and
   the Methods say "at least 90%", so 0.95 fails a published analysis.
4. **Per-chromosome `hasqui_ratio`** — the original divides a per-chromosome
   numerator by the genome-wide bed length, so every `.track` value is ~1/22 of
   the truth. Nothing reads it. We report genome-over-genome in the summary
   instead. Want a correct per-chromosome column too?

## Read first — affects every published DFE

5. **The reported "best fit" is the WORST converged row of the winning file.**
   `plot_converge_dfe_0430.py:282-321` does `data_dict.update({filename: ...})`
   for every row of the `# Converged results` block without breaking, so the
   last row survives. dadi-cli writes that block best-first. Verified on
   `178.InferDFE.bestfits` for phastCons17way_0-5_nomiss: rows are
   -610.1733440069229 / -610.1733440099965 / **-610.1733440139265**, and the
   published value is the third. Cost here is 7e-9 in ll and the 6th decimal of
   shape, but it is systematic across every published fit and is not bounded
   for a poorly converged file. Ported as-is. **Do you want it kept, or the
   first (best) row taken?**
   Note two other scripts in the same pipeline take the FIRST row
   (`dadi_cal_ll_neu_0723.py:112-127`, `make_p0_dfe_0430.py:206-221`), so the
   conventions already disagree internally.

6. **A DFE pack can need 50 h inside a 23:50 h job.** `PACKNUMBER=50` commands
   x `timeout 60m` each vs `h_rt=23:50:00`
   (`submit_dadi-cli_dfebash.sh:8,66`). Not hypothetical: the reference gamma
   directory has 190 starting points but only 144 `.bestfits` - 46 died
   unfinished; neugamma 210 vs 175. **Cut PACKNUMBER, cut the timeout, or raise
   h_rt?**

7. **The two published csv files for this bed came from different script
   versions.** `dfe_..._gamma.csv` has `0.0` in the `>=0.5` bin (the current CDF
   form); `dfe_..._neugamma.csv` has `9.61e-126`, which is exactly the
   commented-out `quad` integration. So the neugamma csv predates the
   2025/03/25 switch and the gamma csv does not. We use the CDF form.

## Errors found in the original scripts (not changed)

5. **`dadi-cli_demog_fsfile.sh` is broken, in two different ways.** The copy in
   `scripts/` builds `out_demog` from `$model_demog` and `$p0marker` at line 32,
   before line 38 sets `p0marker` and line 43 sets `model_demog` — so every file
   it writes would lack the marker that `max_ll_key` depends on. The cluster copy
   has the assignment in the right place but has had the `dadi-cli InferDM` call
   itself deleted. Both edits postdate your results; neither copy runs correctly
   today.
6. **`--output`/`--cpu`/`--sample-size` are argparse abbreviations**, not real
   dadi-cli options (`--output-prefix`, `--cpus`, `--sample-sizes`). They work
   today and would break on any upstream change.
7. **`get_fs_mut_bychr.sh` writes its temp files into the mutation-rate
   reference directory** (`$dirmut`), so 22 concurrent array tasks can collide.
8. **`windctrl_splitchr.sh:122` `hasqui_ratio`** — see 4.
9. **`plot_bestll_demo_args_0421.py:119-122`** lists `bottlegrowth_1d` with two
   parameters; it has three, so that model would abort at :163. We read the
   parameter names from the `.bestfits` header instead.
10. **`pipeline.conf` said `bottlegrowth`**; dadi-cli calls it `bottlegrowth_1d`.
    Fixed in our config only.
11. **`dadi-cli Plot` in `piperef_stddfe_4_plotdemo.sh:52` evidently never ran** —
    no `*.bestfits.pdf` exists beside any published result. Ours writes it.
12. **`piperef_stddfe_5_fitneu.sh` and `_7_fitgamneug.sh` begin with
    `rm -r $dirout`.** Not ported: runs are keyed, so nothing is deleted
    implicitly.
13. **`piperef_stddfe_8_plotdfe.sh` has no `neugammazero` branch** and exits 1.
14. **`piperef_stddfe_2_concatefs.sh` uses GNU `stat -c %s`** and
    `piperef_stddfe_9_sumdfe.sh` uses GNU `sed -i`; both are Linux-only.
15. **`piperef_stddfe_1_rmHETqui.sh` cannot produce the CDS/mask filters** — it
    stops after HET and quiescent. The chain that does all four is
    `submit_exHETquiencmask.sh` then `filt_nc-unmask.sh`.

18. Destructive paths not ported: `rm -r $dirout` in steps 5 and 7; the second
    `.bak` destruction in `make_p0_dfe_0430.py:284-295`; the unconditional
    `.bak` clobber in `plot_DFE_single.py:132-156`; the shared
    `plots_dfe_mutload/` copy keyed only by bed and pdf.
19. Four "Error: ... does not exist" messages whose `exit 1` is commented out,
    so a missing input surfaces later as a traceback.
20. `plot_converge_dfe_0430.py:370` IndexErrors when two or fewer starting
    points converged; `:404` intended a `large_ll.txt` backup that never
    happens; `:494` rounds figure names to 3 decimals so near-identical fits
    overwrite each other (4 pdfs for 10 plotted fits in the reference dir).
21. `piperef_stddfe_8_plotdfe.sh:69` hardcodes `cachenumber=1`, so step 8 can
    only read `cache_1.spectra.bpkl` whatever step 6 built.
22. Checked and NOT a bug: `dadi_cal_ll_neu_0723.py:206-214` passes an unfolded
    model against folded data, but `dadi.Inference.ll_per_bin` folds the model
    itself. Verified: identical to 1e-11.

## Manuscript items (recorded, not changed)

16. **Control window**: the supplement says "100 kb window centered on the middle
    of the test sites"; every run used `mid ± 25000`, i.e. 50 kb. Added as item
    28 of `ToCorrect_260922.docx`.
17. **Mixed provenance.** Results were built with two different control-pairing
    rules and two different dadi versions, split by date:
    - windowing: unmerged before 2024-11-04 (ELS, EnhA, PromF, TSS, BivProm,
      PLS, phastCons17 bins), merged after (PromF_TSS_BivProm, the 2025 phyloP
      and phastCons bins).
    - dadi: demography and cache under 2.3.0, DFE fits under 2.3.3.
    Both differences are ~1e-6 on the numbers we have measured.
