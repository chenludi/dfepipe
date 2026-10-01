# Questions for Chenlu

Found while porting. Nothing acted on; no original script or manuscript file
changed.

## Read first — affects every published DFE

**A. The reported "best fit" is the WORST converged row of the winning file.**
`plot_converge_dfe_0430.py:282-321` does `data_dict.update({filename: ...})` for
every row of the `# Converged results` block without breaking, so the last row
survives. dadi-cli writes that block best-first. Verified on
`178.InferDFE.bestfits` (phastCons17way_0-5_nomiss): the rows are
-610.1733440069229 / -610.1733440099965 / **-610.1733440139265**, and the
published value is the third. Cost here is 7e-9 in ll and the 6th decimal of
shape, but it is systematic across every published fit and unbounded for a
poorly converged file. **RESOLVED 2026-09-30: take the best row.** This is now
the only deliberate numerical difference between the old and new pipelines;
the impact on published values is a todo in docs/tocheck.md.
Two other scripts in the same pipeline already take the FIRST row
(`dadi_cal_ll_neu_0723.py:112-127`, `make_p0_dfe_0430.py:206-221`).

**B. Only 4 starting points were lost to the clock, not 46.** Diagnosed
2026-09-30 for phastCons17way 0-5 gamma (190 starting points):

| outcome | count |
|---|---|
| converged, wrote `.bestfits` | 144 |
| ran all 100 optimizations, no valid likelihood | 42 |
| never started (h_rt kill, contiguous 46-49) | 4 |

The 42 are not a timing problem: every row of their `.opts` files has `--` for
the likelihood and `scale` pinned at the 1e-10 lower bound, i.e. the optimiser
degenerated. They start from the smallest shapes in the grid
(`make_p0_dfe_0430.py` log-spaces shape 1e-10 to 2); e.g. point 62 starts at
shape 4.8e-05.

Still true that `PACKNUMBER=50` x `timeout 60m` = 50 h exceeds `h_rt=23:50:00`
(`submit_dadi-cli_dfebash.sh:8,66`), and that is what cost the 4. **Raise h_rt,
or cut PACKNUMBER?** Separately: **should the p0 grid still spend ~22% of its
starting points on shapes that cannot produce a likelihood?**

**C. dELS and pELS rest on 3 and 4 converged fits.** Measured 2026-09-30 over
60 published gamma runs. Nearly every starting point for those two produced no
output at all - the packs exist, so the jobs were submitted, but almost nothing
ran.

| run | converged | produced any output | starting points |
|---|---|---|---|
| ELS | 132 | 144 | 190 |
| phastCons17way 0-5 | 144 | 186 | 190 |
| CDS | 135 | 143 | 190 |
| PLS | 96 | 138 | 190 |
| **dELS** | **3** | **3** | 158 |
| **pELS** | **4** | **4** | 180 |
| Exc/Inh.PeakCalls | 3 | 3 | ~150-180 |

Not wrong - a converged fit is a converged fit - but far less explored than the
rest, and both appear in Data S2c and Fig. S2. **Worth re-running?**

## Decisions needed

0. **Cache runtime.** Our run took **20 min**, not the 10:10:00 your script
   requested. Cost is data-independent (grid size ~quadratic, gamma_pts linear,
   sample size minor, sequence length none), so one recommendation covers every
   bed. `SGE_CACHE_RT` is now 2:00:00, and `submit_cache.sh` estimates from
   40 core-seconds per spectrum and raises the request itself for a bigger
   marker (`large`: ~12 h estimated, 36 h requested). A killed build now FAILs
   with "raise SGE_CACHE_RT and re-run - a partial cache cannot be resumed".
   **Still on `highp`: keep it now the job is 20 min rather than 10 h?**
1. ~~`CACHE_MIN_SIZE_FRAC`~~ **RESOLVED: 0.5 kept.**
2. ~~`HASQUI_MIN_FRAC`~~ **RESOLVED: 0.90 kept.**
3. ~~Per-chromosome `hasqui_ratio`~~ **RESOLVED: added**, against that
   chromosome's own prep bp, so one bad chromosome is visible.
4. ~~DFE run directory naming~~ **RESOLVED: marker added**, now
   `<CACHE_RUN>-<CACHE_MARKER>-<pdf>-<DFE_RUN>`.
5. ~~Seven config keys the stages needed and `config/` lacks.~~ **Resolved
   2026-09-30:** every tunable is now in `config/pipeline.conf` and the stages
   assert rather than default. `neugammazero` sits in its own section 6, "not
   used by this pipeline". `scripts/checks/conf.py` reads the config so a
   check's argparse default is the configured value rather than a second copy.

### Outcome vocabulary and verdicts — resolved 2026-09-30

- `bestfit` (wrote a `.bestfits`) and `converged` (that file carries a converged
  block) are now separate. The convergence warning uses `converged`, the
  stricter number; keying it on `bestfit` was about 1.6x too generous.
- `DFE_MIN_CONVERGED_N = 10` **stays**. Few converged fits is worth flagging
  whatever the fraction says: it usually means poor starting points, a hard
  surface, or a grid the optimiser never explored.
- A fit that never started optimising (no `.opts`) is `error`, and every
  starting point erroring is a FAIL: that is systematic — missing cache,
  unreadable `--demo-popt`, wrong environment — not a hard likelihood surface.
  A `timeout` is `fail_optimize`, not worth retrying, and `RESUME=1` covers only
  `not_run`.
- `run_concat.sh` refuses when any chromosome has no usable `.fs`/`.tab` rather
  than quietly leaving it out; `CHROMS` is the deliberate way to run a subset.

## Errors in the original scripts (not changed, not ported)

6. `dadi-cli_demog_fsfile.sh`: the `scripts/` copy builds `out_demog` from
   variables assigned later (line 32 vs 38/43), so its files would lack the
   marker `max_ll_key` needs; the cluster copy has the `dadi-cli InferDM` call
   deleted. Both edits postdate your results.
7. `--output` / `--cpu` / `--sample-size` are argparse abbreviations, not real
   options (`--output-prefix`, `--cpus`, `--sample-sizes`).
8. `get_fs_mut_bychr.sh` writes temp files into the mutation-rate reference
   directory, so 22 array tasks can collide.
9. `plot_bestll_demo_args_0421.py:119-122` gives `bottlegrowth_1d` two
   parameters; it has three, so that model aborts at :163. (Our config also said
   `bottlegrowth`; dadi-cli calls it `bottlegrowth_1d`.)
10. `piperef_stddfe_4_plotdemo.sh:52` `dadi-cli Plot` evidently never ran — no
    `*.bestfits.pdf` exists beside any published result.
11. Destructive paths: `rm -r $dirout` in steps 5 and 7; the second `.bak`
    destruction in `make_p0_dfe_0430.py:284-295`; the unconditional `.bak`
    clobber in `plot_DFE_single.py:132-156`; the shared `plots_dfe_mutload/`
    copy keyed only by bed and pdf.
12. Four "Error: ... does not exist" messages whose `exit 1` is commented out,
    so a missing input surfaces later as a traceback.
13. `piperef_stddfe_8_plotdfe.sh` has no `neugammazero` branch and exits 1 —
    the pdf step 7 was extended to fit. Moot for us: neugammazero is out of the
    pipeline (`todo/toimplement_pdf_neugammazero.md`).
14. `plot_converge_dfe_0430.py`: `:370` IndexErrors when <=2 starting points
    converged; `:404` intends a `large_ll.txt` backup that never happens; `:494`
    rounds figure names to 3 decimals so near-identical fits overwrite each
    other (4 pdfs for 10 plotted fits in the reference dir).
15. `piperef_stddfe_8_plotdfe.sh:69` hardcodes `cachenumber=1`, so step 8 can
    only read `cache_1.spectra.bpkl` whatever step 6 built.
16. `piperef_stddfe_1_rmHETqui.sh` cannot produce the CDS/mask filters — it
    stops after HET and quiescent. The full chain is `submit_exHETquiencmask.sh`
    then `filt_nc-unmask.sh`.
17. GNU-only: `stat -c %s` (`piperef_stddfe_2_concatefs.sh`), `sed -i`
    (`piperef_stddfe_9_sumdfe.sh`).
18. Checked and NOT a bug: `dadi_cal_ll_neu_0723.py:206-214` passes an unfolded
    model against folded data, but `dadi.Inference.ll_per_bin` folds the model
    itself. Verified identical to 1e-11.

## Manuscript items (recorded, not changed)

19. **Control window**: the supplement says "100 kb window centered on the middle
    of the test sites"; every run used `mid +/- 25000`, i.e. 50 kb. Added as
    item 28 of `ToCorrect_260922.docx`.
20. **Two published csv files for one bed came from different script versions.**
    `dfe_..._gamma.csv` has `0.0` in the `>=0.5` bin (current CDF form);
    `dfe_..._neugamma.csv` has `9.61e-126`, exactly the commented-out `quad`
    integration, so it predates the 2025/03/25 switch. We use the CDF form.
21. **Mixed provenance**, split by date:
    - windowing: unmerged before 2024-11-04 (ELS, EnhA, PromF, TSS, BivProm,
      PLS, phastCons17 bins), merged after (PromF_TSS_BivProm, 2025 bins).
    - dadi: demography and cache under 2.3.0, DFE fits under 2.3.3.
    Both ~1e-6 on everything measured.

22. dadi warned during our cache build: "Most negative value is -178.88 ... If
    negative values are very negative (<-1), rerun with larger values for
    --grids". Your cache has the same 11 unusable spectra, so this is
    reproduced, not introduced - but dadi is suggesting grids larger than
    1500 1520 1540.

## Reference: how much a DFE actually moves

Spread across your 15 EnhA bootstrap replicates (Data S2a, Fig. S14):
neutral bin 0.0161, other bins 0.0076-0.0089. Our port reproduces the published
bins to 8 significant figures, six orders of magnitude inside that.
