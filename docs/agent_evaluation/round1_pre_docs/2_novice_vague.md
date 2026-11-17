# Round 1, scenario 2: novice, vague request

**Request.** A biologist who is not a population geneticist had a bed file of
enhancer regions, wanted their DFE, and did not know what the settings should
be. The report had to give commands, questions to ask, a config breakdown, a
plain-language explanation of why controls are needed, and every confusion.

**What the agent got right.**
- A clear plain-language explanation: selection shifts the frequency spectrum
  toward rare variants, demography does too, so demography is fitted on nearby
  quiescent controls and held fixed.
- Asked the right questions: genome build, `chr` prefix, interval length
  (intervals over 50 kb can never find a control), population (YRI only), and
  expected compute cost.
- Proposed `DRY_RUN=1` and a small two-chromosome shakedown before the full run.
- Identified `CTRL_FLANK`, `PROJECTION` and `MODELSEL_PRIMARY` as choices the
  user should know about, and left optimiser settings alone.

**What it could not find or had to guess.**
- No README or usage doc; `docs/PIPELINE.md` was referenced but missing, and
  scripts said "see the report" with no report in the package.
- The concat stage needed two legacy R scripts (`concatefs.R`,
  `concate_mut.R`) under `$REFS_DIR`, not shipped and not listed as
  prerequisites; the demographic p0 grid was also external.
- The profile labelled "the worked example" was broken (set `VCF` instead of
  `VCF_TMPL`, lacked `POPFILE`, `POP`, `NEUTRAL_CTRL_DIR`).
- No input validation of the bed: a bed with `1` instead of `chr1` ran for
  hours and died at concat with a confusing message.
- The six run keys were mandatory but undocumented; the naming convention was
  guessed from `records/`.
- Check-file paths had to be reconstructed by hand, and the one worked example
  in `records/` used a pre-`CACHE_MARKER` naming that was now wrong.
- The three-phase chain shows only the current phase in `qstat`, and nothing
  said this was normal; there was no single "did my run work" command.
- `APPLY_CDF_INTEGRATION` sat in a run profile but only affects environment setup.

**Outcome.** Sound advice for a novice, but the agent rated its chance of an
end-to-end first run at about even, mostly on unshipped external dependencies.
