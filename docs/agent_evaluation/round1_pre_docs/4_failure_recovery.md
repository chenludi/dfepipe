# Round 1, scenario 4: failure recovery

**Request.** A user whose DFE array tasks had been killed by the scheduler
wanted to finish the run without redoing work, and asked what to do about one
WARN and one FAIL check file. The report had to give resume commands, the
meaning of each per-starting-point outcome, WARN vs FAIL semantics, and risks.

**What the agent got right.**
- The resume: `RESUME=1 FROM=dfe` with a new `SUMMARY_RUN`, because the
  summary stages refuse to overwrite an existing table.
- The full outcome table in `fit_outcomes.tsv`, and that only `not_run` is
  worth re-running while `timeout` is not.
- WARN never blocks; FAIL stops the next job only inside a chain, through the
  gate on the upstream check file.
- The main trap: `FORCE=1 RESUME=1` discards the fits the resume was meant to keep.

**What it could not find or had to guess.**
- No documentation of the resume procedure; `docs/PIPELINE.md` was missing and
  the recipe was reconstructed from script headers and `todo/` post-mortems.
- `RESUME` was honoured only by `run_dfe.sh`; `fitneu` and the summary stages
  would log failures on a resume, and nothing said this was expected.
- `FORCE=1` meant three different things (bypass a gate, redo every fit,
  overwrite a table) with no narrow gate-only override.
- No `submit_bestdfe.sh`, so re-picking the best fit alone meant hand-building
  a `qsub` line.
- `CACHE_MARKER` was not exported by the profiles, so hand-built paths could
  silently omit it; `check_dfe.py` and `records/` used the old naming.
- `docs/tocheck.md` was cited but the file lived at `todo/tocheck.md`.
- Only the cache path was reconciled on resume; changed bounds, `PACKNUMBER`
  or demography could silently mix incomparable fits.
- The resume machinery had never run on real data, per `todo/final_clean_rerun.md`.

**Outcome.** Correct and careful answer; the agent was fairly confident, less
so on the exact commands, since no documented recipe existed to check against.
