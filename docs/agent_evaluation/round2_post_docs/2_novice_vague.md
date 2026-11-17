# Round 2, scenario 2: novice, vague request

**Request.** A biologist with a bed file of enhancer regions, who is not a population geneticist, asked for help finding their DFE and did not know what the settings should be. The read-only agent had to report the commands, what to ask, the config tiers, a plain-language explanation of why control regions are needed, confusions, confidence, and how it found its way around.

**How it found its way in.** Nothing outside the package pointed anywhere. It found `AGENTS.md` by listing the package root, and from then on followed it as the briefing.

## What the agent got right

- It asked the questions that cause silent wrong answers: genome build (no liftover, and the preflight cannot confirm the build), chr prefix, intervals over 50 kb, YRI as the population, chrX/chrY being dropped, and whether to merge intervals.
- It spotted a bias for enhancer sets: prep subtracts quiescent chromatin from the target, and enhancers called from the user's own data could leave enhancer-like sites in the controls.
- It wrote a clear explanation for a non-specialist: population history shifts every spectrum, so local neutral controls fix the history before selection is inferred.
- It knew `MODELSEL_SMALL_L` might apply to a short enhancer bed and pointed the user to the `modelsel_level` column.

## Problems it found

- The preflight command in `AGENTS.md` was broken as written. `$CHROMS` is empty after sourcing only the profile, so `check_bed.py --chroms ""` fails every bed with a misleading "missing chr prefix" message.
- The profile header told users to source it after `config/pipeline.conf`, which fails. This was fixed while the agent was reading.
- Nothing said the chain must run on the cluster or how to get the bed there, since `sync/push.sh` syncs only the package.
- `DFEPIPE_REMOTE` was documented in `config/profiles/template.conf`, but `sync/push.sh` never read it.
- The header of `run_sfs_chrom.sh` said windows are computed on merged regions. The code merges only when `MERGE_BEFORE_WINDOW=1`, and `run_prep.sh` already merges, so "one control window per input interval" was in doubt.
- `DEMO_P0_SET` appeared in an `AGENTS.md` command, but no config tier listed it. `CACHE_MARKER` is required by `run_all.sh` but was listed in Tier 3, with no advice on when to use another marker.
- Neither `APPLY_CDF_INTEGRATION` (a build-time flag kept in the profile) nor a cost or wall-clock estimate for a whole run was explained anywhere.
- Smaller issues: caveats about published results were scattered across `docs/questions.md`. There is no manifest of reference-data formats for a new site. `check_cache.tsv` is stored apart from every other check file.

**Outcome.** A workable plan, but only for a user who can read the existing reference data. The agent listed the bed's location on the cluster, the ssh host and the remote path as guesses.
