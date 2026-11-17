# Round 2, scenario 3: wrong assumptions

**Request.** The user wanted to run their bed file, which used hg19 coordinates and included chrX and chrY, said they were really interested in a mouse dataset, and asked for the `bottlegrowth` demographic model. The read-only agent had to list every incompatibility with the file that sets it, the commands to use once those were resolved, what to ask, confusions, and confidence.

**How it found its way in.** It started from `AGENTS.md`, whose §2a and §2b had already named all four of the user's problems.

## What the agent got right

- **hg19:** it called this blocking. All reference data is GRCh38, no liftover exists anywhere (it grepped for one), and the preflight's check of more than 1 Mb past the mask will rarely catch hg19.
- **chrX/chrY:** it called this blocking for those intervals. They are dropped silently at concat, and adding X/Y to `CHROMS` does not help, because the projection and the `Na` arithmetic assume diploid autosomes.
- **Mouse:** it called this a porting project, not a config change. Every reference input is human (VCF and panel, mask, GENCODE, ChromHMM, ENCODE-based controls, and the Roulette rate table with its human `ROULETTE_SCALE`).
- **bottlegrowth:** `config/pipeline.conf` accepts only `three_epoch` or `two_epoch`, and only those have p0 grids. It also flagged the dead `bottlegrowth` branches and a known parameter-count bug.

## Problems it found

- The real gap in the hg19 check: the preflight warns that it "could not check" only when `MASK_BED` cannot be read, not when the bed is simply on the wrong build.
- `bottlegrowth` branches in four scripts (`scripts/demo/run_demo.sh`, `scripts/cache/run_cache.sh`, `scripts/demo/run_bestfit.sh`, `scripts/dfe/run_fitneu.sh`) are unreachable and make the model look supported.
- No file said "this pipeline is for human data". The species had to be inferred from file paths.
- `ROULETTE_SCALE` has no source, as `docs/references/README.md` admits, so a non-human port has nothing to derive its own value from.
- The literal `chrs_1t22` is hard-coded in about 10 scripts. There was no guidance on whether a different `CHROMS` is just mislabelled or actually broken.
- Whether the p0 grids exist could only be taken from `todo/vendor_and_docs.md`, not checked.
- Smaller issues: `MODELSEL_PRIMARY=delta_50` decides the headline answer with no stated reason. `FORCE=1 RESUME=1` remains a trap.

**Outcome.** It found all four incompatibilities and confirmed each in the code.
