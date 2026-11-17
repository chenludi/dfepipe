# Audit 3: DFE stage

**Scope.** `scripts/dfe/run_fitneu.sh`, `run_dfe.sh` and `run_bestdfe.sh`, checked against `piperef_stddfe_5_fitneu.sh`, `dadi_cal_ll_neu_0723.py`, `piperef_stddfe_7_fitgamneug.sh`, `make_p0_dfe_0430.py`, `submit_dadi-cli_dfebash.sh`, `piperef_stddfe_8_plotdfe.sh`, `plot_converge_dfe_0430.py` and `plot_DFE_single.py`. Read-only.

Apart from one exception, the numerics are a faithful port. The starting-point grid, the InferDFE argument list (flag for flag), the neutral likelihood, the likelihood cutoff and the CDF-based DFE discretisation all match the originals. Every line citation checked resolved correctly.

**Main findings.**

- **Halved timeout, not listed as a difference.** `DFE_TIMEOUT` is 30m; the original used 60m. `run_dfe.sh` cites the exact original line holding `60m` and says it runs fits "as" that script does. **Could change a result:** a fit that would converge between 30 and 60 min is lost. The timing evidence comes only from fits that succeeded, so it cannot rule this out.

- **Timed-out starting points are never retried** (`run_dfe.sh`), even under `RESUME=1`, unless `FORCE=1`. Together with the 30m timeout, this permanently excludes a starting point at half the original budget. **Could change a result.**

- **Deadline guard skips fits** that cannot finish before the task's h_rt and logs them `not_run`. Which starting points get optimised now depends on queue timing and earlier runs, so a resume can pick a different winning fit than a fresh run. **Could change a result.** Not in the Differences list.

- **False justification for the neutral grid** (`run_fitneu.sh`). The docstring says the neutral model runs on the same grid as the selection cache. It does not (354/388/422 against 1500/1520/1540), here or in the original. No number changes, but a reader comparing the two likelihoods would be misled.

- **Best fit instead of last converged fit** (`run_bestdfe.sh`). This is the one deliberate numerical difference, and it is disclosed accurately: about 7e-9 in ll on the reference bed, but unbounded in general.

- **False reason for keeping a factor of 2** (`run_dfe.sh`). Multiplying by 2 is exact in floating point, so keeping it does not change bit-for-bit results as claimed. Harmless.

- **Dropped knobs and models**: `mutratioadj`, `muadj_ctrl` and `muadj_test` were removed (all default to no adjustment); `snm`/`snm_1d` are rejected by the bash gate but still listed in the Python. No current effect.

- **Undocumented output changes**: output-file rotation is dropped (files are overwritten); the neu log body has extra fields; `large_ll.txt` uses LF instead of CRLF; on a malformed bestfits file the parser stops reading where the original kept going; there are new `spec.txt` files. None changes a number.

- **Miscited provenance** for `DFE_SECONDS_PER_FIT=330`: the cited line is the timeout. This setting affects scheduling only.

**Not checked.**

- Anything that needs dadi or dadi-cli to be installed: the auto-fold behaviour, the `cache1d.integrate` signature, the patched `neugamma` pdf, and `--cpu` abbreviation resolution.

- The quoted measurements (1e-11 fold equivalence; the ordering of rows near −610.17), which rest on project notes rather than primary files.

- Whether the generated command packs are byte-identical to the published ones (none were available to diff against).
