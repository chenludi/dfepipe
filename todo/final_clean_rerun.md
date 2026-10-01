# Final clean rerun: the acceptance test for the package

Development has run stages one at a time, resumed part-finished arrays and
swapped a cache mid-run, so the current output tree carries the history of
building it rather than the behaviour of the finished pipeline. Before the
package is called done, delete the development output and run the whole chain
once, untouched, from the bed files.

That run - not the development one - is the evidence the package works.

## Why it is needed, concretely

- `dfe_260929` mixes two caches: gamma 171 reference / 59 pipeline, neugamma
  99 / 52 (`tocheck.md` item 4). Its `large_ll.txt` ranks both together.
- Several stages were run with scripts older than the current ones: the DFE
  arrays of 2026-09-30 predate the outcome log, the severity column and the
  h_rt deadline guard, so none of that machinery has yet run on real data.
- `fit_outcomes.tsv` has never been written by a completed run.
- Run keys are a patchwork (`prep_260929`, `fs_260929`, `fs_260929nm`, ...).

## The run

- [ ] One set of run keys for every stage, dated the day of the run.
- [ ] Delete the development output root, or point `OUT` somewhere new. Keep
      `$REFS_DIR/result_dadi` untouched - it holds the published results the run
      is compared against.
- [ ] No `DFE_CACHE_FILE`: build the cache from scratch so the chain is tested,
      not short-circuited.
- [ ] No `RESUME`, no `FORCE`. If a stage needs either, that is a finding about
      the pipeline and belongs in the write-up.
- [ ] Run through the chain scripts, not stage by stage by hand.
- [ ] Record the submission commands and job IDs in `records/`.

## What to check afterwards

- [ ] Every stage's check file is PASS or an understood WARN, with no FAIL.
- [ ] The DFE bins match the published values inside the bootstrap spread
      (neutral bin 0.016, other bins ~0.008). Anything larger is a real finding.
- [ ] `fit_outcomes.tsv` exists for every DFE run, accounts for every starting
      point, and its severity counts are sane.
- [ ] No task was killed by `h_rt`: the deadline guard should have left
      `not_run` rows instead, and if it did, note how many and why.
- [ ] One cache per DFE run directory - `spec.txt`'s `cachefile=` matches what
      every `.opts` in that directory used.
