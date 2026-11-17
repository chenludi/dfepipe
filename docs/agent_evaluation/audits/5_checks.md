# Audit 5: the checks layer

**Scope.** The 11 files in `scripts/checks/` (`checkfile.py`, `conf.py`, `gate.py`, `residuals.py`, `check_bed.py`, and the six stage checks). Each was checked for docstring-versus-code agreement, threshold boundaries, crash paths, and whether it can alter results. The functions were called on synthetic input; the test suite passed (360 tests).

**Main findings.**

- **No check can alter a pipeline result.** The only write path in the whole layer is `checkfile.py`'s writer for the check `.tsv`, and there are no deletes, renames or subprocess calls. No PASS/WARN/FAIL inversion remains, and residual and deviance thresholds only ever WARN, even at 1e12 or inf.

- **NaN passes every threshold** (5 files). A NaN SNP count, `Na`, log-likelihood or DFE bin reads as PASS. A NaN `Na` would break the cache stage downstream, and negative counts also get through. **Could let a bad result pass unflagged.**

- **Crash on a ragged summary row** (`check_summary.py`). The FAIL row and the per-bed row share an item name, so `write_check` refuses both and no check file is written. The one structural FAIL this module advertises crashes in its most likely form, which then blocks `bestpdf`.

- **Other crashes that leave no check file (and so block the chain)**: a NUL byte in any CSV-read file (5 readers, including `gate.py` itself, where `FORCE=1` cannot override); a missing `--out` directory in `check_prep` and `check_summary`; duplicate items from two beds with the same basename, `--chroms "1 chr1"`, or beds that shorten to the same name.

- **Silent non-blocking cases in `gate.py`**: a check file with only a header, or the wrong header, parses to nothing and passes. An unrecognised status loses to a duplicate PASS row, despite the "worst status wins" docstring.

- **Weak cache checks** (`check_cache.py`). The gamma-range test warns only below 10% of the requested bound (the comment implies 90%), and the size floor passes a cache built with half the requested gamma points. `n_gamma` is never compared with `--gamma-pts`.

- **`check_demo` never flags poor convergence**: with only one converged fit out of 20, it returns PASS. `check_dfe` does have this check.

- **Threshold drift**: the `HASQUI_MIN_FRAC` fallback is 0.90 but the configured value is 0.75, so a check run by hand and a pipeline-run check disagree. `CONVERGENCE_LL_SPREAD`, `RANGE_TOLERANCE`, `PICKLE_OVERHEAD` and `OUT_OF_RANGE_SLACK` are hardcoded, which contradicts `conf.py`'s "one place" claim. A cache-size measurement is quoted with two different values.

- **Docstring overclaims**: `check_demo`, `check_cache` and `gate.py` say WARNs are "compared with past runs", which is false for thresholds that are plain judgements. `gate.py` says reasons are quoted, but they are actually stripped. Several FAIL and WARN lists are incomplete.

- **Inconsistencies**: `check_bed` names items with `.bed` while other checks do not, so their rows cannot be joined; `_bed_base` is copied into six files; there are two `.bestfits` parsers; `check_sfs` drops every line beginning with `t`; negative coordinates are reported as "out of order"; `--ctrl-flank 0` silently disables the long-interval test.

**Not checked.**

- The dadi branch of the residual code (dadi was not installed). The fallback was checked against the vendored dadi source and matches to machine precision.

- Real pipeline outputs: all inputs were synthetic.

- How the gates behave across the whole chain, and concurrent writes to the same check file.
