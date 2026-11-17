# Audit 6: the SGE and chain layer

**Scope.** `scripts/chains/run_all.sh`, `scripts/sge/submit_*.sh`, `worker_*.sh`, `gate_upstream.sh`, `config/env.sh` and the site profiles: job submission, gating, resources, and `FORCE`/`RESUME`/`DRY_RUN` behaviour on a UGE/SGE cluster. Read-only; no qsub was run.

**Main findings.**

- **The DFE array never checks its upstream gate** (`worker_dfe.sh`). The gate call sits inside the `neu` branch, and the chain never sends `neu` through that script. A failed or missing cache check therefore does not stop the fits; only `fitneu` and `bestdfe` stop. This contradicts `run_all.sh`'s "every job gates" and explains the shape of a past incident.

- **Stale command packs after a `PACKNUMBER` change** (`run_dfe.sh`, `submit_dfe.sh`). Old packs are not removed, so with `FORCE=1` two array tasks run overlapping starting points against the same output prefix and race on the same `.bestfits` files. The package itself recommends both triggers: "Re-run with FORCE=1 PACKNUMBER=…", and setting `DFE_TARGET_TASK_HOURS` on a resume. **Can overwrite finished fits.**

- **`FORCE=1` has six meanings, not three.** It overrides a gate, overrides the bed preflight, overwrites finished outputs, redoes finished fits, rebuilds an existing cache, and re-attempts timed-out points. It applies to every later stage and phase, so using it to get past one gate (which the gate message suggests) rewrites finished work downstream. `FORCE=1 RESUME=1` discards the fits that `RESUME` was meant to keep.

- **Possible lost phase 3** (`run_all.sh` `launch_next`). `-hold_jid` is passed even when empty, which happens when the cache stage has nothing to do. Whether qsub rejects an empty value, which would stop the chain silently, is suspected but not tested.

- **Failure on the second pdf orphans the first.** If `submit_dfe.sh` fails for the second pdf, the first pdf's array runs with nothing scheduled to collect it, and `summary`/`bestpdf` are never submitted.

- **Silent abort** in `submit_dfe.sh` when a pack has no fit lines: `grep -c` returns 1 under `set -e`, with no message.

- **Spaces in paths shift the task-to-bed mapping**: lists are joined with spaces in the submit scripts and split on whitespace in the workers, so a bed path containing a space maps tasks to the wrong beds. Suspected: SGE `-v` may not keep spaces inside values such as `PDFS`, which could mean no DFE array is submitted at all.

- **Gate cascades**: one crashed check blocks every stage below it. A crashed cache check for one bed blocks `fitneu`/`bestdfe` for all beds, while the DFE arrays still run. A bed skipped as already cached, whose old check file is missing, blocks permanently.

- **`DRY_RUN=1` only shows phase 1** (prep through bestdemo), not cache, DFE or summary, and the submit scripts ignore it if called directly.

- **Resource comments**: "~10 hours" for the cache contradicts the code's own 20-minute estimate and 2 h request. sfs and demo h_rt values are fixed and do not grow with the number of beds (up to 500 allowed), and demo has no deadline guard. The `DFE_MAX_TASKS` arithmetic can overshoot by up to (beds − 1).

- Minor: `worker_step.sh` logs `rc=` only on success; `submit_demo.sh` turns any directory entry into an array task; `RESUME` is not forwarded by three submit scripts (no current effect).

**Not checked.**

- qsub behaviour with `-hold_jid ""` and with spaces inside `-v` values.

- Whether `SGE_TASK_ID` can leak through `-V` from a launcher job, and whether `flock` exists on compute nodes (without it, appends to `fit_outcomes.tsv` race).

- The p0 grid contents on the cluster, and real check files on disk.
