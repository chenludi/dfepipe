# Round 2, scenario 1: expert, precise request

**Request.** An expert user asked for the DFE of phastCons 17-way 0-5 conserved sites from their own bed file, using the `three_epoch` demographic model, run on the cluster. The read-only agent had to report the exact commands, what to ask first, the config tiers, the science in brief, confusions, confidence, and how it found its way around.

**How it found its way in.** It started from `AGENTS.md` ("Start here"), which pointed it to `docs/PIPELINE.md` and gave it a map of the package.

## What the agent got right

- It split the work correctly: `sync/push.sh` with `CHECK=1` runs locally, then source the site profile, set the six run keys, check the reference data, preflight the bed, run `DRY_RUN=1` and submit on the cluster.
- It asked about assembly, chr prefix, autosomes and population (YRI, 108 individuals, projected to 160), and warned that results will not match the published table exactly (best fit, not last fit).
- It explained the three-phase self-submitting chain and resuming with `RESUME=1 FROM=dfe` plus a new `SUMMARY_RUN`, and warned against `FORCE=1 RESUME=1`.
- It gave a correct statistical summary: fit demography to local neutral controls, then fit a gamma or neutral-plus-gamma DFE over a precomputed cache.

## Problems it found

- The site profile's header said "Source AFTER pipeline.conf", but `config/pipeline.conf` asserts `OUT` and `REFS_DIR` and dies if sourced first. `AGENTS.md` said the opposite, and correctly.
- `AGENTS.md` §4 put the local push and the cluster `run_all.sh` call in one block, with no step that says "now move to the cluster".
- The `sync/push.sh` destination (host and remote path) was not documented. It could only be inferred from a file in `records/`.
- The manual preflight in `AGENTS.md` §2a used `$CHROMS`, which no profile sets, and it duplicated the preflight `run_all.sh` already runs.
- The profile's `OUT` pointed at a development output tree that `todo/final_clean_rerun.md` calls unfit, and that tree already held output for this same bed.
- `AGENTS.md` did not say the chain had never run clean end to end, though `todo/final_clean_rerun.md` did (for example, `fit_outcomes.tsv` had never been written by a completed run).
- Smaller issues: `CACHE_MARKER` is asserted by `run_all.sh` but not listed as required. `APPLY_CDF_INTEGRATION` sits in the profile but is read only by `setup/setup.sh`. `SGE_CACHE_EXTRA=highp` assumes access to the `highp` queue. `FORCE=1` has several meanings. There is no wall-clock or disk estimate for a whole run.
- Known caveats (the ~22% `no_bestfit`, 11 of 300 unusable cache spectra) were spread across several files, and none of those files was linked from `AGENTS.md`.

**Outcome.** A sound plan. The remaining doubt came from undocumented environment details, not from the commands.
