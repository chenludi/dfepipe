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
