# Round 2, scenario 4: failure recovery

**Request.** The scheduler killed some DFE array tasks partway through a run (run keys `*_a`). The user asked how to finish without redoing work, and what to do about one check file that said WARN and another that said FAIL. The read-only agent had to give the commands, every fit outcome and whether it is worth re-running, WARN versus FAIL, the risks of resuming, confusions, and confidence.

**How it found its way in.** It started from `AGENTS.md`, whose §5 and §6 answered the question directly, then checked those answers against the scripts.

## What the agent got right

- Its procedure: read the check files first, count the outcomes in `fit_outcomes.tsv`, run `CHECK=1` on the sync, run `DRY_RUN=1`, then `RESUME=1 FROM=dfe bash scripts/chains/run_all.sh <bed>` with a new `SUMMARY_RUN`. No `FORCE=1`, and nothing else changed.
- It explained all six outcomes correctly by severity. Only `not_run` (severity `rerun`) is worth re-running, `timeout` is deliberately skipped, and `error` needs a fix first.
- It explained WARN versus FAIL. A WARN never blocks. A FAIL, or a missing check file, stops only the next job in the chain, through `gate.py`, not through `-hold_jid`.
- It listed the existing guards: refusing to resume with a different cache, refusing `RESUME=1` when there is nothing to resume, and the append-only `fit_outcomes.tsv` with `flock`.

## Problems it found

- The docs and the code disagreed. `AGENTS.md` and `docs/PIPELINE.md` said `RESUME=1` re-runs only `not_run`. `run_dfe.sh` actually skips only points with a `.bestfits` or a logged `timeout`, so the `no_bestfit` and `error` points run again.
- A live bug: `check_dfe.py` (`count_outcomes`) counted points that the deadline guard logged as `not_run` as `error`/`broken`. That can raise a false FAIL ("all N starting points exited..."), and no test covers it.
- The cache-mismatch guard compares only the cache path, so a cache rebuilt at the same path gets past it.
- Marker identity is only an index into the starting-point grid. Changing the bounds or bins, or re-running demo, between the run and the resume silently pairs old fits with new starting points.
- Packs are regenerated on every submission and stale ones are not cleared, so a changed `PACKNUMBER` can run some points twice. `spec.txt` is rewritten each time and stops describing the original run.
- `FORCE=1` has three meanings (bypass a gate, redo every fit, overwrite a table), so `FORCE=1 RESUME=1` destroys what the resume was keeping.
- The recovery command in `AGENTS.md` left the required new `SUMMARY_RUN` out of the command line. A FAIL for one pdf blocks the summary for the other, and nowhere says whether that is intended.
- The hand-run override path was under-documented: run on the login node, or wrap in `qsub`?

**Outcome.** A correct resume procedure, plus one real doc-versus-code disagreement and one classification bug, found by reading the code.
