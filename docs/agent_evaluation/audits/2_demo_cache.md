# Audit 2: demography and cache stages

**Scope.** `scripts/demo/run_demo.sh`, `scripts/demo/run_bestdemo.sh` and `scripts/cache/run_cache.sh`, checked against the originals (`dadi-cli_demog_fsfile.sh`, `plot_bestll_demo_args_0421.py`, `piperef_stddfe_4_plotdemo.sh`, `dadi-cli_cache_simp.sh`, `submit_dadi-cli_cache_simp.sh`) and against published output directories. Read-only.

The dadi-cli command lines for both InferDM and GenerateCache match the originals argument for argument, and the Na formula reproduces a published value to 16 digits.

**Main findings.**

- **False claim that the model spectrum matches the figure** (`run_bestdemo.sh`). The header says the saved model spectrum is the one `dadi-cli Plot` draws. The port folds the model before scaling; dadi-cli Plot does not. So `bestfit_model.fs`, and the residual and G2 numbers computed from it, are not the spectrum shown in `sfs_bestfit.pdf`. **Changes reported residuals and WARN/FAIL flags**, but not the demographic parameters, Na, theta or the chosen fit. The folded version is arguably the more correct one; only the claim of identity is wrong.

- **False justification for a renamed figure** (`run_bestdemo.sh`). The header says the original's per-file PDF was never produced. 667 such PDFs exist beside published results. The rename itself is fine.

- **Dropped model guard** (`run_bestdemo.sh`). The original exited if the parameter count did not match the model. The port reads parameter names from each file, so a run directory with fits from two models could silently pick a winner across models. **Could change a result**, but only in a mixed directory.

- **Contradictory runtime claim** (`run_cache.sh`). The header says "~10 hours on 10 cores", but that is the original's h_rt request; the package elsewhere measures 20 min and requests 2 h. The `highp` queue request was also dropped.

- **Uncited provenance**: two of the three `config/cache_markers.tsv` rows come from `pipefs_exinitials_6_cache.sh`, which is not cited. The values were verified to match.

- **Latent spurious WARN**: for the absolute-range cache markers, `run_cache.sh` still passes `GAMMA_S_HIGH` to `check_cache.py`, which compares it with an unrelated bound.

- **Shared knobs**: `MAXEVAL` and `NUM_CONVERGENCE` are single keys used by both demo and DFE. Their values match today, but changing one for the DFE silently changes the demography. `SGE_DEMO_PE`/`SGE_CACHE_PE` are labelled "no effect on results" but feed dadi-cli `--cpus`.

- **Undocumented deviations**: the Na fallback (`Na == "NA"`) is dropped; the original's 600000-byte cache-size warning becomes a fraction-based FAIL; marker rows can override `GAMMA_PTS`/`GRIDS`; the all-fits `ll_Na.pdf` is overwritten by the post-cutoff one; the scan for the `etal` row covers the whole table; malformed p0 lines are rejected. All cosmetic or defensive.

- **New thresholds with no original**: `CACHE_MIN_SIZE_FRAC=0.5` (no stated basis), `CACHE_MIN_USABLE_FRAC=0.9` (the expected value is 0.963), and the `RESID_*`/`G2_PER_SNP_*` thresholds, calibrated on 137 earlier fits with a single known bad fit. These change flags, not values.

- Minor citation slips: the GenerateCache call sits outside its cited line range; one citation is off by a line.

**Not checked.**

- Nothing was executed. The size of the folded-vs-unfolded residual difference was not computed.

- The p0 grid files (they live on the cluster) and an end-to-end reproduction of `best_inference.txt`.

- Which of the two cache-runtime statements is correct.
