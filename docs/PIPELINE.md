# The pipeline, stage by stage

Per-stage reference: what each stage reads, what it writes, what its check means,
and why each threshold is the number it is. `../AGENTS.md` is the front door and
covers how to run the whole thing — read that first.

Six stages. Each writes into a directory named by its own run key, so a run can
never overwrite another, and each ends with a check file.

```
prep     bed            -> the analysable bases of that bed
sfs      prep bed, VCF  -> per-chromosome spectra for target and control
         concat         -> the genome-wide spectra and mutation tables
demo     control SFS    -> demographic fits
         bestfit        -> the winning demography
cache    demography     -> expected spectra over ~300 selection coefficients
dfe      cache, target  -> the DFE fits
         bestdfe        -> the winning fit, discretised
summary  everything     -> one row per bed
         bestpdf        -> which pdf the data supports
```

---

## prep — what of the bed can be analysed at all

`scripts/prep/run_prep.sh <bed>...`

Four set operations, in this order, though the order does not change the final
set of bases — only the intermediate counts:

1. subtract `HET_BED` (heterochromatin)
2. subtract `QUI_BED` (the quiescent annotation — so no target site is also a
   control)
3. subtract `CDS_BED` (coding sequence)
4. intersect `MASK_BED` (keep only accessible bases), then sort

Writes `s.nc_exHET_noq.<bed>.bed` and `lengths.tsv`, which records the bp
surviving each step. That file is how a user finds out where their bed went: a
bed of promoters can lose most of its bases to the quiescent subtraction alone.

**Check.** `PREP_MIN_RETAIN_FRAC=0.05` — a WARN if under 5% of the input bp
survives. Not a FAIL, because a legitimately small overlap is possible; it is
there so nobody discovers it three stages later.

---

## sfs — the frequency spectra

`scripts/sge/submit_sfs.sh <bed>...` — one array task per chromosome, then
`scripts/sfs/run_concat.sh <bed>...`.

Per chromosome, for the **target** set and the **control** set separately:
subset the VCF to the bed, build a folded spectrum with dadi, project to
`PROJECTION` chromosomes, and tabulate the mutational target from the Roulette
rate bed.

**The control window.** For each target interval, the midpoint is taken and the
window is `[mid − CTRL_FLANK, mid + CTRL_FLANK]`. Candidate-neutral sites inside
that window become controls. Two consequences worth knowing:

- The window width does **not** scale with the interval. An interval longer than
  `2 × CTRL_FLANK` has its window entirely inside itself, where controls cannot
  be, so it never finds one. `check_bed.py` warns about this before the run.
- `CTRL_FLANK=25000` means a 50 kb window. The manuscript supplement says
  100 kb; every run used 50. Recorded in `questions.md` item 19.

The target set that goes downstream is named `has_qui` — the intervals that
**found** a control. Intervals with no control in their window are not analysed.

**Why controls at all.** The demography is fitted to them. A spectrum shifted
toward rare variants is the signature of purifying selection, but population
history shifts every spectrum the same way; the controls carry the history
without the selection, so fitting them first separates the two. Keeping them
local also matches mutation rate, recombination and coverage.

**Checks.**

| threshold | default | why |
|---|---|---|
| `MIN_SNPS` | `500` | below this the spectrum is too thin to fit. **FAIL.** The smallest published bin had 1,560 SNPs, so 500 is well clear of anything real |
| `HASQUI_MIN_FRAC` | `0.90` | the fraction of prep bp that found a control. **FAIL** below it. ELS, the worst published case, was 0.9105 — so 0.90 passes everything that was published and catches a real collapse |
| `MIN_FS_BYTES` | `661` | a per-chromosome spectrum smaller than this holds no counts |
| `MIN_FSMUT_BYTES` | `500` | likewise for the mutation table |

A per-chromosome `hasqui_ratio` is reported against that chromosome's own prep
bp, so one bad chromosome is visible rather than averaged away.

`run_concat.sh` refuses to build a genome-wide spectrum if any chromosome with
intervals produced nothing — a partial spectrum that looks complete is worse
than no spectrum. A chromosome the bed never touches is not an error. The output
file name is the literal `chrs_1t22` whatever `CHROMS` was, because this pipeline
analyses autosomes; `concat_chroms.tsv` records what actually went in.

---

## demo — the demographic model

`scripts/sge/submit_demo.sh <bed>...` — one array task per starting-point file in
`DEMO_P0_DIR` — then `scripts/demo/run_bestfit.sh <bed>...`.

Fits `DEMO_MODEL` to the **control** spectrum with `dadi-cli InferDM`, from many
starting points, and picks the maximum-likelihood fit. Writes
`best_inference.txt`, whose `Na` and `theta` every later stage reads, and
`bestfit_model.fs`, the model spectrum.

`three_epoch` has four parameters (nuB nuF TB TF), `two_epoch` two (nu T). Only
these two are accepted, because only they have a starting-point grid.

**Checks.** Convergence, and the fit's residuals against the control spectrum.
`DEMO_LL_CUTOFF=500` bounds how far below the best a fit can be and still count.

**The residual thresholds**, which are the same for demo and dfe:

| threshold | default | what it means |
|---|---|---|
| `RESID_MAX_RARE_WARN` | `25` | the largest absolute Anscombe residual among the rare bins (1 to `RESID_RARE_MAX_COUNT=4`). A WARN above it |
| `RESID_MAX_RARE_FAIL` | `40` | beyond anything seen in past runs. Still only a WARN in effect — `classify()` FAILs only on structural problems |
| `RESID_BULK_BIAS_WARN` | `2.3` | the mean normalised residual over the non-rare bins. A systematic sign here means the model is wrong across the bulk of the spectrum, not just in one bin |
| `RESID_REST_OVER_RARE_MIN` | `5` | warns when the largest residual is outside the rare bins, which is unusual and worth seeing |
| `G2_PER_SNP_WARN` / `_FAIL` | `0.008` / `0.035` | Poisson deviance per SNP, `G² = 2Σ[d·ln(d/m) − (d−m)]`. Scaling by SNP count makes it comparable between a 1,500-SNP bin and a 500,000-SNP one |

These are percentiles of previously accepted fits, not statistical tests. They
say "unusual", not "wrong", which is why they warn and never block. **The figure
is the arbiter** — `sfs_rank1_poisson.pdf` and the demo plots.

**Sign convention:** dadi's residual is `(data_transformed − model_transformed)`,
so a **positive** residual in the rare bins means the model **over**estimates
them. The transform is Anscombe, hence the apparent inversion.

---

## cache — expected spectra over selection coefficients

`scripts/sge/submit_cache.sh <bed>...`

`dadi-cli GenerateCache1d` precomputes the expected spectrum for `GAMMA_PTS=300`
selection coefficients, spanning `2·GAMMA_S_LOW·Na` to `2·GAMMA_S_HIGH·Na`, on
the extrapolation grids `GRIDS="1500 1520 1540"`. Written as
`cache_<CACHE_MARKER>.spectra.bpkl`.

Cost is **data-independent** — grid size roughly quadratic, `GAMMA_PTS` linear,
sample size minor, sequence length not at all. Measured: 20 minutes on 10 cores.
`SGE_CACHE_RT=2:00:00` is generous for that. A cache killed part-way cannot be
resumed; raise the limit and rebuild.

**Checks.**

| threshold | default | why |
|---|---|---|
| `CACHE_MIN_SIZE_FRAC` | `0.5` | the file must be at least half of `gamma_pts × (projection+1) × 8 × 1.4`. A fixed byte floor was tried and would have failed every real cache: real ones range 335 KB to 1.8 MB |
| `CACHE_MIN_USABLE_FRAC` | `0.9` | the fraction of the 300 spectra that are usable. dadi warns that 11 of 300 come back with large negative values at these grids. That is **reproduced, not introduced** — the same 11 appear in the original cache — and larger grids were tested without benefit |

`CACHE_MARKER` names a row of `config/cache_markers.tsv`, which records the
bounds and grids a cache was built with, and appears in the directory name so two
caches cannot be confused.

---

## dfe — the fits

`scripts/dfe/run_fitneu.sh <bed>...` for the null, then
`scripts/sge/submit_dfe.sh <pdf> <bed>...` per pdf, then
`scripts/dfe/run_bestdfe.sh <pdf> <bed>...`.

**`neu`** is the null: the likelihood of the target spectrum under the demography
alone, with `theta_test = theta_demog × mutratio` and no selection. One number,
no optimisation.

**`gamma`** (shape, scale) and **`neugamma`** (shape, scale, pneu) are fitted by
`dadi-cli InferDFE` against the cache, with the demography held fixed.
`neugamma` adds a neutral spike: `pneu/1e-4` over `0 ≤ γ < 1e-4`, on top of
`(1−pneu) ×` the gamma. It is not stock dadi — `setup/patches/` adds it.

**One starting point = one pack line = one `dadi-cli` call** with one `--p0`,
which internally runs up to `OPTCOUNT=100` optimisations from perturbations of
that seed and stops early once `NUM_CONVERGENCE=5` agree. `PACKNUMBER` of those
go into one array task.

The starting-point grid log-spaces shape and scale over the pdf's bounds and
keeps combinations where `2e-6·Na ≤ shape·scale ≤ 2·0.5·Na`. **About 22% of the
grid is small-shape seeds that cannot produce a likelihood at all** — they run
all 100 optimisations, return `--`, and write no `.bestfits`. Expected, and the
fastest group: under 4 minutes each, against 7 m 53 s for the slowest fit that
worked.

`fit_outcomes.tsv` records every starting point as
`marker / outcome / converged / severity / seconds / note`. `seconds` is that
whole `dadi-cli` call's wall clock, covering all of its optimisations, and 0 when
the point was skipped rather than run.

| outcome | what came back | severity |
|---|---|---|
| `bestfit` | `.opts` and `.bestfits`. Whether that file carries a converged block is the stricter, separate question | `ok` if converged, else `fail_optimize` |
| `no_bestfit` | `.opts` only: it optimised and no row was usable — the likelihood came out `--`, or the parameters left the bounds | `fail_optimize` |
| `error` | nothing, not even `.opts`: dadi-cli never ran one optimisation. A technical fault — missing cache, unreadable `--demo-popt`, wrong environment | `broken` |
| `timeout` | a part-written `.opts`: killed by `DFE_TIMEOUT` | `fail_optimize` |
| `not_run` | nothing, and no log row: the array task never reached it | `rerun` |

The severity is the column to read first. `rerun` is the work a resume is for.

What `RESUME=1` skips, precisely, is a starting point that already has a
`.bestfits` and one logged `timeout`; everything else is attempted again. A
timeout is excluded deliberately — a point needing more than `DFE_TIMEOUT` is a
bad starting point, not a slow one, so retrying spends another `DFE_TIMEOUT`
reaching the same place. `no_bestfit` points are *not* excluded, so they repeat
on every pass; they are the fastest group, so this costs wall clock only.

`not_run` arrives two ways: no log row at all, because the array task died before
reaching that point, or a logged `not_run` from the deadline guard, which
declines to start a fit when under `DFE_TIMEOUT` of `h_rt` remains. Both are
`rerun`.

**Checks.**

| threshold | default | why |
|---|---|---|
| `DFE_MAX_NOTRUN_FRAC` | `0.20` | WARN when over 20% never ran — the parameter space was under-explored |
| `DFE_MIN_CONVERGED_FRAC` | `0.20` | and `DFE_MIN_CONVERGED_N=10`. Keyed on the **converged** count, not merely on fits that produced output: about a third of published `.bestfits` files carry no converged block. The count floor stays because a handful of converged fits usually means poor starting points or an unexplored grid, whatever the fraction says |

`run_bestdfe.sh` takes the **best** converged fit and discretises it over
`DFE_BINS` using the gamma CDF, scaled by `(1−pneu)` with `pneu` added to the
first bin. It draws `sfs_rank<k>_poisson.pdf` — the model at its fitted theta
against the data, with Anscombe Poisson residuals.

**Only the Poisson figure, not the multinomial one.** For a DFE the total SNP
count is a *prediction* — theta is fixed in advance from the neutral sites times
the mutation ratio — so rescaling the model to match the data's total would hide
the reduction in SNP count that selection causes, which is the thing being
measured. For the demography the opposite holds: theta is a free nuisance
parameter, so the model is optimally scaled and multinomial is correct. Measured
on phastCons17way 0-5: the demo model total equals the data total exactly, while
the DFE model over-predicts by 0.63%.

---

## summary — the table

`scripts/summary/run_summary.sh <bed>...` then
`scripts/summary/run_bestpdf.sh`.

`sumdfe.tsv` is one row per bed: `ll_neu`, each pdf's winning fit and parameters,
its six DFE bins, and `u`, `l`, `ul` from the genome-wide row of the mutation
table. Column names are prefixed by pdf (`gamma_shape`, `neugamma_shape`) so
every column can be addressed.

`bestpdf_sumdfe.tsv` adds model selection. Three likelihood differences, compared
against every confidence row of `config/cutvalue_modelsel.csv`:

```
gamma beats neutral  when ll_gamma    − ll_neu    > cut[gamma-neu]
then, from whichever won:
  from gamma     neugamma wins when ll_neugamma − ll_gamma > cut[neugamma-gamma]
  from neutral   neugamma wins when ll_neugamma − ll_neu   > cut[neugamma-neu]
```

Every level gets a `pdf_best_<level>` column, and `pdf_best_allconfidence` is
filled only when they all agree — a useful signal that the verdict is not an
artefact of where the cutoff sits.

`MODELSEL_PRIMARY=delta_50` decides the reported DFE. The `delta_N` rows are flat
log-likelihood margins and are what the published figures use; the `chi_*` rows
are the nominal likelihood-ratio cutoffs. **`delta_50` is an empirical choice,
not a derived statistic.**

Below `MODELSEL_SMALL_L=1500000` bp, `MODELSEL_SMALL_PRIMARY=chi_5%` decides
instead, because a flat 50-unit margin is a stricter test on a short sequence
than a long one. `modelsel_level` records which applied. Set
`MODELSEL_SMALL_L=0` for the original single-cutoff behaviour.

A `neutral` verdict reports the DFE as `[1, 0, 0, 0, 0, 0]` rather than the gamma
bins. `num_nonneutral = (1 − first bin) × ul`.

**Check.** The six bins must sum to 1 within `1e-6`; a pdf's likelihood no better
than `ll_neu` warns that selection buys that model nothing.

---

## How much a DFE actually moves

The yardstick for whether a difference matters, from 15 EnhA bootstrap
replicates (Data S2a, Fig. S14): the neutral bin spreads by **0.0161**, the other
bins by **0.0076–0.0089**. A difference smaller than that is noise.

For context, this pipeline reproduces the published phastCons17way 0-5 result to
~1e-5 (gamma) and ~4e-7 (neugamma) — four to six orders of magnitude inside that
spread.

---

## Reference tables in `config/`

| file | what it holds |
|---|---|
| `cache_markers.tsv` | one row per `CACHE_MARKER`: the gamma bounds, `GAMMA_PTS` and `GRIDS` a cache was built with |
| `cutvalue_modelsel.csv` | one row per confidence level: the three likelihood-difference cutoffs |
