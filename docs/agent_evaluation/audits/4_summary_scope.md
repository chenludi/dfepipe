# Audit 4: summary stage and scope creep

**Scope.** (A) `scripts/summary/run_summary.sh` and `run_bestpdf.sh`, checked against `piperef_stddfe_9_sumdfe.sh`, `summary_getul.py` and `summary_getbestpdf.py`. (B) Everything in the package that has no counterpart in the original scripts. Read-only.

The core arithmetic is faithful: u, l and ul, the model-selection comparison, and `num_nonneutral` all reproduce a published row exactly, and `config/cutvalue_modelsel.csv` is byte-identical to the original table.

**Main findings.**

- **`HASQUI_MIN_FRAC=0.90` FAILed and blocked the chain**, though its docstring said it "never blocks". 8 of 280 published beds (Tx, TxEnh, TxEx, TSS and others) fall below 0.90, so reproducing them would have needed an override. **Blocked known-good results.** (Since fixed; see Audit 5.)

- **`MODELSEL_SMALL_L` originally shipped at 1,500,000.** At that setting, beds shorter than 1.5 Mb were judged by `chi_5%` instead of `delta_50`, which can flip a reported pdf from neutral to gamma. It was set to 0 during the audit, with a note that the figure was invented while porting. **Changes the reported model when enabled.**

- **Output column names break the original figure scripts.** `pdf_plotting` and `best_bin_<label>` replace `pdf_plotting_delta_50` and `delta_50_<label>`, which the manuscript plotting scripts read. Three published columns (`pdf`, the DFE row name and `File`) are dropped, and `ul` moves. None of this is in the Differences list.

- **A crash in an invented check blocks the final table.** The stages treat their own check as advisory, but `gate.py` blocks on a missing check file. If `check_summary.py` crashes, `bestpdf_sumdfe.tsv` never appears even though `sumdfe.tsv` is correct. This has happened once already.

- **Deadline guard** (`run_dfe.sh`): up to 20% of starting points can be skipped based on queue timing while the run still reports PASS. **Could change the best fit.**

- **Interior-bin non-finite residual → FAIL** (`residuals.py`). This can be a correct outcome for a strongly selected target rather than a sign of breakage, yet it is the only residual-derived FAIL.

- **Undocumented summary deviations**: the "exactly one log file" check is gone (the first match is taken silently); the log glob pattern changed; `.bed` stripping is corrected rather than ported; `pdf_best_allconfidence` gained an extra condition; a missing `.etal` marker aborts every bed where the original failed one.

- **Wrong citations**: several line numbers in `run_summary.sh` and `run_bestpdf.sh` are off by up to ~100 lines; `.etal` is attributed to `concatefs.R` when it comes from `concate_mut.R`; the `ROULETTE_SCALE` citations in `pipeline.conf` are off. Two claimed hazards (identical rows from two beds collapsing; a "shift by one") do not happen as described.

- **Invented machinery judged over-engineered**: `conf.py`, a second source of truth for thresholds whose regex fails on some keys; `FORCE=1` with several unrelated meanings; the `DFE_TARGET_TASK_HOURS` sizing keys; the experiment write-up (including a .docx) under `setup/patches/optional/`. Judged justified: `check_bed.py`, `gate.py`, the run-key directories, the `RESUME` guards, the profiles and `test_env.py`.

- **Invented thresholds**: 16 warn/fail thresholds in `pipeline.conf` have no original. `MIN_SNPS=500` turns into a FAIL what the original only warned about. `CACHE_MIN_USABLE_FRAC=0.90` sits only about 7 points below the expected 96.3%.

**Not checked.**

- Nothing was executed. The calibration data behind `RESID_*`, `G2_*` and `DFE_SECONDS_PER_FIT` is not in the package.

- Scheduler behaviour, and whether the cluster copy matches the local one.

- Because two of three line citations checked in `run_summary.sh` were wrong, the package's line citations in general should be treated as unverified.
