# Round 1, scenario 1: expert, precise request

**Request.** An experienced user asked to infer the DFE of a phastCons 17-way
conserved-site bin from their own bed file, using the `three_epoch` demographic
model, on the cluster. The report had to give exact commands, questions to ask
first, a config breakdown, a short statistical summary, and every confusion.

**What the agent got right.**
- The overall shape of a run: source a site profile, export the six run keys,
  `DRY_RUN=1` first, then one call to `scripts/chains/run_all.sh` that
  resubmits itself forward in three phases.
- Recognised that this exact bed was the project's validation bed, and asked
  whether the user wanted a fresh run or the existing answer.
- Noticed that a fresh run will not reproduce published numbers exactly,
  because the port takes the best converged fit where the original took the last.
- A correct statistical summary: demography fitted on local quiescent controls,
  held fixed, then a gamma or neutral-plus-gamma DFE fitted to the target SFS.

**What it could not find or had to guess.**
- No README or entry-point doc; `config/pipeline.conf` points at
  `docs/PIPELINE.md`, which did not exist. Invocation was reconstructed from the
  header of `run_all.sh` and a file in `records/`.
- `sync/` was empty: no way to push the package to the cluster. The agent
  invented its own rsync command.
- Two site profiles both claimed to be for the cluster and disagreed; the one
  labelled "the worked example" lacked `CONDA_ENV`, `VCF_TMPL`, `POPFILE`,
  `POP`, `NEUTRAL_CTRL_DIR` and could not run.
- `RESID_REFERENCE` defaulted to a `config/residual_reference.tsv` that did
  not exist and that nothing read.
- `MERGE_BEFORE_WINDOW`: no way to tell the right value for this bed family,
  and a wrong value silently changes the controls.
- Profile variables reached phases 2 and 3 only through `qsub -V`, so running
  without sourcing the profile in the same shell would fail late.
- `CACHE_MARKER` sits in a different position in cache and DFE directory names,
  so the agent had to cross-read four scripts to find check-file paths.
- `todo/final_clean_rerun.md` said the chain had never completed end to end;
  `DFE_TARGET_TASK_HOURS=19` was the agent's own arithmetic, not documented.

**Outcome.** Correct plan in outline, but the deploy step was invented, and
the agent rated its chance of a clean first run at about even.
