# Round 1, scenario 5: parameter questions

**Request.** A user asked six questions: what `CTRL_FLANK` is and whether to
change it; how `QUI_BED` differs from `NEUTRAL_CTRL_DIR`; `two_epoch` or
`three_epoch`; what changes for a 200 kb target; what `PDFS` and `neugamma`
are; and what the DFE bins are. The report had to cite a source for each answer.

**What the agent got right.**
- `CTRL_FLANK` is the half-width of a window centred on each interval's
  midpoint; keep 25000 unless intervals exceed 50 kb, which can never find a control.
- `QUI_BED` is subtracted from the target; `NEUTRAL_CTRL_DIR` is a stricter
  control set (688 Mb against 968 Mb).
- `neugamma` is a neutral spike plus a gamma, from a patch to dadi applied at setup.
- The six DFE bins in |s|, why they sum to 1, and where the check enforces it.

**What it could not find or had to guess.**
- No criterion for choosing `two_epoch` vs `three_epoch`: no model comparison,
  no guidance, no record of a `two_epoch` run, p0 grid not visible.
- No bp-to-SNP guidance, so whether 200 kb clears `MIN_SNPS` was unanswerable.
- `MODELSEL_SMALL_L` was described as a mutational target but compared against
  raw bp length, and it switched short beds to a different model-selection cutoff.
- No README; `docs/PIPELINE.md` and `docs/tocheck.md` were cited but missing.
- Comments contradicted code: `run_demo.sh` said `bottlegrowth` was accepted;
  `docs/questions.md` described a config section that did not exist.
- The manuscript supplement described a 100 kb control window; every run used
  25 kb either side, i.e. 50 kb.
- `PDFS` was ignored by some scripts in favour of a hardcoded `gamma|neugamma`
  list, so adding a pdf needed edits in several places.
- Naming was misleading: `has_qui` is the test set; `CTRL_FLANK` is not a flank.

**Outcome.** Most questions answered well from the code; the demographic-model
choice and the minimum target size were left open.
