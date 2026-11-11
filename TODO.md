# Known-unfinished work

Not hidden, and worth reading before telling anyone this package is finished.
Each entry records a dated incident, and they are the reason most of the guards
in `AGENTS.md` exist.

Detail lives in `todo/`, one file per topic. Finished work is removed from
there rather than archived, so `git log` is where it went. Disk policy is
`docs/tidyup.md`.

| file | what |
|---|---|
| `todo/tocheck.md` | open questions about the pipeline's own behaviour |
| `todo/vendor_and_docs.md` | what is still outside the package, and the flag/naming wrinkles |
| `todo/toimplement_pdf_neugammazero.md` | the neugammazero pdf, removed 2026-09-30, and what putting it back needs |

## Open questions

The numbering matches the sections of `todo/tocheck.md`, which `AGENTS.md`
cites by number.

1. **Does `has_qui` change when the test bed is merged before windowing?**
   `MERGE_BEFORE_WINDOW` currently has no effect, because `run_prep.sh` already
   merges, so the bed reaching sfs is merged either way. Verified on chr22:
   re-merging changes neither the interval count nor the bp. It does matter for
   a run that starts at sfs and skips prep.
2. **How much did the worst-converged-row bug move the published results?**
   This pipeline takes the *best* converged fit; the original took the last,
   which was the worst. That is the one deliberate numerical difference, about
   1e-7 in practice but unbounded for a poorly converged file.
3. Use measured fit times instead of `DFE_SECONDS_PER_FIT`.

## One finding from the acceptance runs

The runs of 2026-10-25 and 2026-10-26 were the clean end-to-end test: one dated
key set each, no `DFE_CACHE_FILE`, no `RESUME`, no `FORCE`, driven by the chain
scripts, every check PASS or an understood WARN, and ELS matching the published
Data S2 row to 1e-5. One thing did not come out clean, and is still open:

- A neugamma fit for `GRCh38-ELSv4_not_Prom` hit `DFE_TIMEOUT` and the deadline
  guard never fired, so it was killed rather than left as a `not_run` row for a
  resume to pick up. 207-210 of 210 fits finished; the guard is meant to make
  that impossible.

## Still outside the package

`todo/vendor_and_docs.md` records what is not yet vendored, and the `FORCE=1`
overloading: one flag means six things, which is why `SKIP_GATE=1` exists as the
narrow override.
