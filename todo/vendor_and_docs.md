# To do: the last self-containment gaps

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
