# To check — open questions, not blocking

## 1. Why does `has_qui` change when the test bed is merged before windowing?

**Observed (ELS, 2026-09-29).** Merging the test bed before centring the control
window changes `has_qui` on exactly two chromosomes out of 22, in opposite
directions:

| chromosome | merged | unmerged (pre-2024/11/04) | difference |
|---|---|---|---|
| chr12 | 10,123,537 | 10,123,385 | **+152 bp** |
| chr16 | 6,485,291 | 6,485,466 | **−175 bp** |
| all others | — | — | 0 |

Genome-wide this is −26 bp in 202.7 Mb (1.3 × 10⁻⁷), so it changes no
conclusion. The question is mechanistic, not numerical.

**Mechanism (traced 2026-09-29).** Merging *relocates* the window, it does not
enlarge it. Two bookended ELS intervals on chr12 (120230636-120230966 and
120230966-120231118) merge to midpoint 120230877 instead of 120231042, moving
the window 165 bp left onto a quiescent site it previously missed - so the
152 bp interval was gained, and nothing was lost on chr12. A separate mechanism
applies above 2 x CTRL_FLANK = 50 kb: the window falls entirely inside the
region, where controls cannot be (controls are candidate-neutral minus test),
so such regions can never find one. That is the one worth counting.

Magnitudes (ELS): control 602,593,295 merged vs 602,593,816 unmerged (merging
can only shrink the control set); test 205,209,325 vs 205,209,348.

**Default is now `MERGE_BEFORE_WINDOW=0`** (one window per input interval), so
this is off the default path - but the >50 kb question still applies to raw
intervals.

**Steps**

- [ ] Distribution of merged region lengths for ELS per chromosome: how many
      exceed 2 × `CTRL_FLANK` (50 kb), and what fraction of bp do they hold?
      Expect chr12 and chr16 to be unremarkable — if they are, the explanation
      is local, not a property of those chromosomes.
- [ ] Identify the specific regions that changed status. Diff the two `has_qui`
      beds for chr12 and chr16 (`MERGE_BEFORE_WINDOW=1` vs `0`) and list the
      intervals gained and lost, with their lengths.
- [ ] For each changed region, record: merged length, number of unmerged
      sub-intervals, distance to the nearest candidate-neutral site, and whether
      the merged window falls entirely inside the region. This distinguishes the
      two mechanisms above directly.
- [ ] Count genome-wide how many merged regions exceed 50 kb and therefore
      *cannot* find a control under the current rule. If that number is large,
      the rule is quietly discarding long regions, which matters well beyond
      this 26 bp — it would bias `has_qui` against long annotations.
- [ ] Check the same for a long-region annotation, not just ELS. Tx or a
      merged-promoter set would show the effect more strongly if it is real.
- [ ] With the default now unmerged, check the same question for raw intervals:
      how many single input intervals exceed 50 kb, so that their own window
      falls entirely inside them and they can never find a control? For ELS the
      intervals are a few hundred bp, but Tx or a long-segment annotation could
      hit this.
- [ ] Decide whether the window should scale with the region — e.g.
      `[start − FLANK, end + FLANK]` rather than `mid ± FLANK`. That is a
      different definition of "nearby", so it is a science decision, not a
      cleanup: it changes which controls pair with which test sites for every
      annotation. Do not change it without deciding deliberately.

**Decision, 2026-09-29: the default is now `MERGE_BEFORE_WINDOW=0`**, one window
per input interval. The intervals of an annotation are meaningful units - one
cCRE, one ChromHMM segment - and the pipeline should not silently redefine them;
anyone who wants merged regions can merge their input bed first. It also matches
the majority of published results, since merging only entered the original
scripts on 2024/11/04. The question below is still worth answering, but it is no
longer on the default path.

**Where the evidence is.** `config/pipeline.conf` (`CTRL_FLANK`,
`MERGE_BEFORE_WINDOW`), `scripts/sfs/run_sfs_chrom.sh` (the windowing), and the
original `../windctrl_splitchr.sh:41-42,68-76`. The published ELS SFS predates
the merge, which is why it matches the unmerged numbers.


## 2. How much did the worst-converged-row bug move the published results?

`plot_converge_dfe_0430.py:282-321` keyed its dict by filename and updated it
once per converged row, so each row overwrote the previous and the **last**
survived - the worst converged fit in that file, not the best. The convergence
tracking itself (`covfile`) was fine; the overwrite was incidental. Confirmed on
`178.InferDFE.bestfits` (phastCons17way 0-5): rows -610.1733440069229 /
-610.1733440099965 / -610.1733440139265, and the published value is the third.

**The pipeline now takes the best row (2026-09-30).** That is the only
deliberate numerical difference between the old and new pipelines.

- [ ] Quantify across every published fit: for each `*.InferDFE.bestfits`, the
      difference between the first and last converged row in ll, shape and
      scale, and how often a file has more than one converged row at all.
- [ ] Whether any bed's *winning file* changes when the best row is used - that
      is the only way the discretised DFE could move visibly.
- [ ] Re-derive the affected Data S2 rows if anything moves beyond the bootstrap
      spread (neutral bin 0.016, other bins ~0.008).

Expected to be tiny: for phastCons17way 0-5 the first-to-last gap is 7e-9 in ll
and the 6th decimal of shape.


## 3. Use measured fit times instead of DFE_SECONDS_PER_FIT

`submit_dfe.sh` sizes a task from `DFE_SECONDS_PER_FIT=330`, one shipped
constant. Measured 2026-09-30 on phastCons17way 0-5 gamma, timing each fit as
the gap between consecutive `.opts` mtimes inside one task:

| | n | median | p95 | max |
|---|---|---|---|---|
| produced a best fit | 165 | 1m15s | 5m12s | 7m53s |
| produced no best fit | 65 | 2m03s | 3m09s | 3m35s |

Two bestfit fits are recorded near 60m. At least one is an artifact: marker 54
has a `.bestfits`, so `timeout 60m` cannot have let it run 61m47s - the gap
includes marker 53's timeout, which immediately preceded it.

Note the shape of this, which is the opposite of what we first assumed: the
doomed small-shape seeds are the FASTEST group, not the slowest. They run all
100 optimisations but each bails at once with no likelihood, so they finish in
under 4 minutes. The 7-minute fits are genuine ones doing real work.

- [ ] Read a previous run's `fit_outcomes.tsv` and use its 95th percentile
      instead of the constant. The column is already there.
- [ ] 330 s is roughly the bestfit p95, so it is conservative by about 4x
      against the median. Fine for sizing, poor as a report.

Low priority: the deadline guard means a wrong estimate costs one extra
`RESUME=1` pass, never lost work.


## 4. Two caches got mixed inside one DFE run directory

**What happened (2026-09-30).** `DFE_RUN=dfe_260929` for phastCons17way 0-5
holds fits made against two different caches:

| run | reference cache `result_dadi/...-demo_three_epoch-1` | pipeline cache `demo_260929-1-cache_260929` |
|---|---|---|
| gamma | 171 | 59 |
| neugamma | 99 | 52 |

**How.** The validation run set `DFE_CACHE_FILE` to the existing reference cache
on purpose, to compare against the published result. Later resubmissions did not
set it, so `run_dfe.sh` fell back to its default path - which by then existed,
because the cache stage had run. `submit_dfe.sh` regenerates the packs on every
submission, so each resubmission rewrote the command lines with the new cache,
and `spec.txt` with them. Finished fits kept the old cache; new fits got the new
one. Nothing warned.

Worse, `spec.txt` now records only the later cache, so the provenance file
misstates what most of the fits used.

**Fixed:** `run_dfe.sh` refuses to resume when the resolved cache differs from
`spec.txt`'s `cachefile=`, naming both paths. `FORCE=1` or a new `DFE_RUN` is
then the way forward.

**Decision 2026-09-30: leave it.** `dfe_260929` stays as it is and development
continues on it. The two caches agree to the precision anything here is measured
at, so the mixture does not mislead a development check, and re-running it now
would cost hours for no new information. It is not a manuscript result.
Correctness is established instead by the clean end-to-end rerun, which is the
acceptance test for the finished package - see `todo/final_clean_rerun.md`.
- [ ] `spec.txt` is rewritten on every pack regeneration, so it describes the
      last submission rather than the run. Either stop rewriting it once fits
      exist, or version it.
- [ ] Packs are regenerated unconditionally by `submit_dfe.sh`. Stale packs from
      a larger PACKNUMBER survive if the new names do not collide with the old
      ones, which would run some starting points twice. Clear the old packs when
      regenerating.


## 5. RESUME=1 started a fresh run instead of resuming

**What happened (2026-09-30).** Adding `CACHE_MARKER` to the DFE run directory
name (`<CACHE_RUN>-<CACHE_MARKER>-<pdf>-<DFE_RUN>`, questions.md item 4) moved
where `run_dfe.sh` looks. The arrays of that morning ran a cluster copy from
before the rename, so their 151 neugamma fits landed in
`cache_260929-neugamma-dfe_260929`, while a later `RESUME=1` submission looked
for `cache_260929-1-neugamma-dfe_260929`, found nothing, and began from scratch.
`run_bestdfe.sh` then read those 9 orphan fits and reported ll=-6397 where the
real answer is -610.

Nothing warned at any point: an empty directory is indistinguishable from a new
run, and "resumed" output looks like success.

**Fixed:** `run_dfe.sh` now fails when `RESUME=1` is set and the run directory is
absent or empty, naming the keys to check. The real runs were renamed into the
marker convention and the 9-fit directory kept as
`...-dfe_260929.aborted-14984282`.

- [ ] Nothing reconciles a run directory against the config that made it. The
      cache guard (item 4) and this one are both after-the-fact patches for the
      same gap: `spec.txt` records the inputs but nothing checks it on re-entry.
      A single `verify_spec` step, run whenever a stage re-enters an existing run
      directory, would replace both.
