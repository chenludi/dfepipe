# To do: the last self-containment gaps

Closed on 2026-10-01, recorded here because each was found the same way - by
giving the package to an agent with no other context and watching where it had
to guess:

- the two R scripts `run_concat.sh` shelled out to are now vendored in
  `scripts/sfs/R/`, byte-identical, with provenance
- the input bed is now preflighted (`scripts/checks/check_bed.py`)
- `sync/push.sh` ships, so nobody has to invent an rsync. It uses checksums and
  `--delete`, and `CHECK=1` compares without writing
- `config/profiles/hoffman_klohmuel.conf`, which could not run and called itself
  "the worked example", is replaced by `config/profiles/template.conf`
- `RESID_REFERENCE` is gone: nothing read it and the file it named did not exist
- `MUTRATE_PATTERN` is overridable now; the literal `}` in `{chrnum}` had been
  closing the `${VAR:-default}` expansion early

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

- [ ] `FORCE=1` means three different things: bypass an upstream gate
      (`gate.py`), redo every DFE fit (`run_dfe.sh`), and overwrite a finished
      table (`run_summary.sh`, `run_bestpdf.sh`). There is no way to ask for
      only one, so `FORCE=1 RESUME=1` discards the fits `RESUME` was meant to
      keep. Separate flags would fix it.
- [ ] `RESUME=1` is honoured only by `run_dfe.sh`. So the natural recovery
      command, `RESUME=1 FROM=dfe bash scripts/chains/run_all.sh <bed>`, is
      *expected* to log a `fitneu` failure and usually a `summary` failure,
      because those refuse existing output and ignore `RESUME`. Harmless, and
      confusing until someone tells you.
- [ ] `check_cache.tsv` lives inside each bed's cache directory while every
      other check file is `<OUT>/<stage>/check_<stage>.<run>.tsv`.
