# Known-unfinished work

Not hidden, and worth reading before telling anyone this package is finished.
Each entry records a dated incident, and they are the reason most of the guards
in `AGENTS.md` exist.

This is the index. The detail that is too long to sit here, and that shipped
code cites by path, is in two files beside it:

| file | what |
|---|---|
| `tocheck.md` | the open questions in full, numbered; `config/pipeline.conf`, `check_bed.py`, `run_sfs_chrom.sh` and `run_bestdfe.sh` cite its numbers |
| `toimplement_pdf_neugammazero.md` | the neugammazero pdf, removed 2026-09-30, and what putting it back needs; `run_dfe.sh` cites it |

Finished work is removed rather than archived, so `git log` is where it went.
Disk policy is `docs/tidyup.md`.

## Open questions

Developer notes. `todo/tocheck.md` holds the detail and keeps its own numbering,
which shipped code cites by number.

1. **The control window is centred on the interval's midpoint, and its width is
   fixed.** So an interval longer than `2 x CTRL_FLANK` encloses its own window
   entirely, finds no neutral control, and is dropped from the analysis without
   appearing in the result. That is the caveat on the current scripts and it
   needs fixing: windowing from the interval edges would give every interval a
   window of its own length. The same root decides what the midpoint is in the
   first place, which is why merging the bed before windowing moves `has_qui` on
   two chromosomes. `todo/tocheck.md` item 1.
2. Use measured fit times instead of `DFE_SECONDS_PER_FIT`. `todo/tocheck.md`
   item 3.

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

- [ ] `DEMO_P0_DIR`, the demographic starting-point grids, live under
      `$REFS_DIR/dadi_initialparams/<model>_<set>/`. Only `three_epoch_1` and
      `two_epoch_1..6` exist, which is *why* `pipeline.conf` accepts only those
      two models. They are small text files - vendoring `three_epoch_1` and
      `two_epoch_1` would make the supported model set a property of the package
      rather than of one cluster's filesystem.
- [ ] The reference data proper (VCFs, mask, CDS, ChromHMM, Roulette rates, the
      control beds) cannot be vendored and should not be. What is missing is a
      *manifest*: a list of what must exist, with the format of each, so a new
      site knows what to assemble. `config/profiles/template.conf` names them;
      it does not describe their formats.

## Still unresolved, and worth knowing

- [ ] `FORCE=1` is still overloaded, though less than this note first said.
      It was written as three meanings; it turned out to be six, and the
      worst of them is now fixed. `SKIP_GATE=1` (2026-10) bypasses an
      upstream gate and does nothing else, so clearing a gate no longer
      rewrites finished work - which it did, because the gate's own message
      used to recommend `FORCE=1`. What remains conflated is redo every DFE
      fit (`run_dfe.sh`) with overwrite a finished table (`run_summary.sh`,
      `run_bestpdf.sh`), so `FORCE=1 RESUME=1` still discards the fits
      `RESUME` was meant to keep. Separate flags would finish the job.
- [ ] `RESUME=1` is honoured only by `run_dfe.sh`. So the natural recovery
      command, `RESUME=1 FROM=dfe bash scripts/chains/run_all.sh <bed>`, is
      *expected* to log a `fitneu` failure and usually a `summary` failure,
      because those refuse existing output and ignore `RESUME`. Harmless, and
      confusing until someone tells you.
- [ ] `check_cache.tsv` lives inside each bed's cache directory while every
      other check file is `<OUT>/<stage>/check_<stage>.<run>.tsv`.
