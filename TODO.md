# Known-unfinished work

Not hidden, and worth reading before telling anyone this package is finished.
Each entry records a dated incident, and they are the reason most of the guards
in `AGENTS.md` exist.

Detail lives in `todo/`, one file per topic. Resolved questions are in
`docs/questions.md`; disk policy is `docs/tidyup.md`.

| file | what |
|---|---|
| `todo/tocheck.md` | open questions about the pipeline's own behaviour |
| `todo/final_clean_rerun.md` | the acceptance test: delete the development output and run the whole chain once, untouched |
| `todo/vendor_and_docs.md` | what is still outside the package, and the three flag/naming wrinkles |
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
4. Two caches got mixed inside one DFE run directory.
5. `RESUME=1` started a fresh run instead of resuming.
6. The demo stage's `h_rt` does not scale with the number of beds.

## Still outside the package

`todo/vendor_and_docs.md` records what is not yet vendored, and the `FORCE=1`
overloading: one flag means six things, which is why `SKIP_GATE=1` exists as the
narrow override.
