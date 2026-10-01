# dfepipe — for an agent helping someone run this

**Start here.** You are helping a researcher infer a distribution of fitness
effects (DFE) for a set of genomic regions they supply as a bed file. This file
is the whole briefing: what the pipeline does, what to settle with the user
before spending their cluster time, and what to run. `docs/PIPELINE.md` has the
per-stage detail; read it when you need a specific stage.

Do not guess a configuration value. Everything the pipeline needs is either in
`config/pipeline.conf` with a default that reproduces the published analysis, or
is a path that only the user can tell you. If you find yourself inventing a
number, stop and ask.

---

**This pipeline is for human data.** Every reference input is human and GRCh38:
1000 Genomes VCFs and accessibility mask, GENCODE CDS, Ernst ChromHMM states,
ENCODE cCREs in the control definition, and Roulette mutation rates with a
human-specific scale constant. Another species is a porting project, not a
configuration change — see `todo/vendor_and_docs.md`.

## 1. What the pipeline does, in four sentences

A set of sequences under selection carries fewer common variants and more rare
ones, because deleterious mutations are removed before they spread — so the
**site frequency spectrum** of the target regions encodes selection. But
population history distorts every spectrum in the genome, selected or not, so
the pipeline first fits a **demographic model to nearby neutral control
regions** (quiescent chromatin within `CTRL_FLANK` of each target interval),
where the distortion is all there is to see. Holding that demography fixed, it
then asks what distribution of selection coefficients must be added to explain
the target spectrum, fitting a gamma (and optionally a neutral-plus-gamma
mixture) by maximum likelihood over a precomputed cache of expected spectra.

Keeping the controls **local** matters: it matches mutation rate, recombination
and sequencing coverage between target and control, so the residual difference
is attributable to selection rather than to where in the genome the targets sit.

**Two sentences for a non-specialist:** we count how many variants sit in your
regions and how common each one is, then compare that pattern against nearby
stretches of DNA we believe are not under selection. The comparison is what
separates "selection removed these mutations" from "human population history
makes every spectrum look like this."

---

## 2. The conversation to have first

Work through these with the user before running anything. Several cause silent
wrong answers rather than errors.

### 2a. The bed file — confirm, do not assume

| ask | why it matters |
|---|---|
| **Which assembly?** | All reference data is GRCh38. A bed from another build produces plausible, wrong numbers. Nothing lifts over, and **the preflight cannot detect it** - it only catches a coordinate running more than 1 Mb past the mask. This one is a human answer. |
| **`chr1` or `1`?** | Must be `chr`-prefixed. `scripts/checks/check_bed.py` now catches this in a second; before, it cost the whole sfs array. |
| **Autosomes only?** | The pipeline analyses chr1–chr22. chrX/chrY intervals survive prep and are **dropped silently at concat**, so the reported length is not the bed's length. The preflight warns with the bp cost. |
| **How long are the intervals?** | An interval longer than `2 × CTRL_FLANK` (50 kb) has its control window entirely inside itself, where controls cannot be, so it can **never** find one and is dropped. Fine for a few-hundred-bp cCRE set; serious for long domains. |
| **What do the intervals mean?** | They are kept as units by default (`MERGE_BEFORE_WINDOW=0`) — one control window per input interval. If the bed has adjacent or overlapping intervals the user may want them merged first. |
| **Which population?** | A DFE is population-specific. The configured data is 1000 Genomes YRI, 108 individuals, projected to 160 chromosomes. This is not a free choice without a different VCF and panel. |

Run the preflight and show the user its output — it answers most of the above
from their actual file:

```bash
python3 scripts/checks/check_bed.py --mask-bed "$MASK_BED" their.bed
```

It reads `CHROMS` and `CTRL_FLANK` from `config/pipeline.conf` itself, so do not
pass them: a profile does not export `CHROMS`, and passing an empty one used to
fail every bed while blaming the chromosome names.

### 2b. The demographic model — the one real science choice

`DEMO_MODEL` accepts `three_epoch` or `two_epoch`. Each needs a
starting-point grid at `$REFS_DIR/dadi_initialparams/<model>_<DEMO_P0_SET>/`;
those two are the only ones that have them.

Tell the user this: **`three_epoch` is what this project validated and what every
published result used. For an unfamiliar population, start with the simpler
`two_epoch`, then compare — fewer parameters is easier to fit and harder to
overfit, and if `three_epoch` does not improve the fit materially there is no
reason to prefer it.** Then let them decide. The pipeline will not choose for
them, and it has no automatic demographic model comparison — comparing means
running both and looking at the likelihoods and the residual figures.

### 2c. Reference data that must already exist

Not shipped, and the pipeline cannot start without them. Check they are there:

```bash
ls "$REFS_DIR/dadi_initialparams/${DEMO_MODEL}_${DEMO_P0_SET}/"   # the p0 grid
ls "$NEUTRAL_CTRL_DIR" | head                                    # qui2_chr1.bed ...
ls -l "$MASK_BED" "$CDS_BED" "$HET_BED" "$QUI_BED" "$MUTRATE_TAB"
ls "${VCF_TMPL/\{chrom\}/chr22}"
```

**`QUI_BED` and `NEUTRAL_CTRL_DIR` are different things** and the names invite
confusion. `QUI_BED` is the raw quiescent annotation, *subtracted from* the
target bed so no target site is quiescent. `NEUTRAL_CTRL_DIR` is the control
set, which is stricter — in this project it is quiescent **and** non-coding
**and** mask-passing **and** free of ENCODE enhancer/promoter annotation: 688 Mb
against `QUI_BED`'s 968 Mb.

---

## 3. Configuration, in three tiers

Everything lives in `config/pipeline.conf`. Nothing else carries defaults.

### Tier 1 — must be set, the run cannot start without them

The six **run keys**, which name the output directories so a run can never
overwrite another. Use one dated set across all six:

```bash
export PREP_RUN=prep_261001 FS_RUN=fs_261001 DEMO_RUN=demo_261001 \
       CACHE_RUN=cache_261001 DFE_RUN=dfe_261001 SUMMARY_RUN=sum_261001
```

Plus the site paths, which come from a profile — never set them by hand:

```bash
source config/profiles/hoffman_replicate.conf    # or your own, from template.conf
```

`config/profiles/template.conf` lists every key a profile must export and what
it is for. It is a template and will not run as-is.

### Tier 2 — the user should decide, defaults reproduce the manuscript

| key | default | what it decides |
|---|---|---|
| `DEMO_MODEL` | `three_epoch` | see 2b |
| `CTRL_FLANK` | `25000` | control window is midpoint ± this. **Changing it changes which controls pair with which targets for every annotation** — a science decision, not a tuning knob |
| `PROJECTION` | `160` | sites genotyped in fewer chromosomes are dropped; trades sites for sample size |
| `PDFS` | `neu gamma neugamma` | which DFE shapes to fit. `neu` is the null; `neugamma` adds a neutral spike |
| `MERGE_BEFORE_WINDOW` | `0` | **currently has no effect** - `run_prep.sh` already merges, so the bed reaching sfs is merged either way. Verified on chr22: re-merging changes neither the interval count nor the bp. `todo/tocheck.md` item 1 |
| `CHROMS` | `1`–`22` | only for a deliberate subset; concat refuses a partial set otherwise |

### Tier 3 — tell the user these exist, then leave them alone

Roughly sixty keys: the optimiser bounds and starting-point grid density per pdf,
`GAMMA_PTS=300` and `GRIDS="1500 1520 1540"` for the cache, `OPTCOUNT=100`,
`MAXEVAL=200`, `NUM_CONVERGENCE=5`, `DFE_BINS`, `MODELSEL_PRIMARY=delta_50`,
every warn/fail threshold in section 4 (they never change a number, only what
gets flagged), and the scheduler resources in section 5.

Two worth naming because a user may ask:

- **`MODELSEL_PRIMARY=delta_50`** decides which model-selection cutoff reports
  the headline DFE. It is an empirical choice, not a derived test statistic. The
  `chi_*` rows in `config/cutvalue_modelsel.csv` are the nominal
  likelihood-ratio cutoffs.
- **`GRIDS`** stays as it is. dadi warns that 11 of 300 cache spectra are
  unusable at these grids; larger grids were tested and did not help, so the
  warning is expected and reproduced rather than introduced.

Scheduling knobs that cost nothing to change: `PACKNUMBER` (fits per array
task), `DFE_TARGET_TASK_HOURS` and `DFE_MAX_TASKS` (which size `PACKNUMBER` for
you), `DFE_TIMEOUT` (30 m; the slowest real fit measured was 7 m 53 s).

---

## 4. Running it

**Two machines.** The push runs where you edit; everything else runs on the
cluster, because `run_all.sh` calls `qsub` and the bed and reference data are
there.

**On your machine**, from the package root — is the cluster copy current? A stale
script is the classic silent failure, and `CHECK=1` writes nothing:

```bash
CHECK=1 bash sync/push.sh <host>:<path>     # exits 1 and lists what differs
bash sync/push.sh <host>:<path>             # only if it reported differences
```

The destination for this project is recorded at the top of
`config/profiles/hoffman_replicate.conf`.

**Then on the cluster**, one shell:

```bash
cd <package>
source config/profiles/<site>.conf          # FIRST: pipeline.conf asserts OUT
export PREP_RUN=... FS_RUN=... DEMO_RUN=... CACHE_RUN=... DFE_RUN=... SUMMARY_RUN=...

DRY_RUN=1 bash scripts/chains/run_all.sh their.bed    # always, first
bash scripts/chains/run_all.sh their.bed
```

The chain preflights the bed itself and refuses to submit on a FAIL, so §2a's
manual `check_bed.py` call is only for showing the user the verdict before
committing their cluster time. If you do run it by hand, `source
config/pipeline.conf` too — `$CHROMS` and `$CTRL_FLANK` come from there, not
from a profile.

**The chain submits itself forward in three phases**, because two stages need the
previous stage's output at *submission* time (`submit_cache.sh` reads the
demography; `submit_dfe.sh` generates the command packs). So `qstat` shows only
the current phase, and the rest arrives later via launcher jobs:

```
phase 1  prep sfs concat demo bestfit   -> launcher
phase 2  cache                          -> launcher
phase 3  fitneu dfe bestdfe summary bestpdf
```

Tell the user that. Jobs appearing, vanishing and reappearing an hour later is
normal and not a failure.

`FROM=<stage>` starts part way — `prep sfs demo cache dfe summary` — when
earlier stages are already done under the same run keys.

**The answer** is `$OUT/summary/$SUMMARY_RUN/bestpdf_sumdfe.tsv`: one row per
bed, with each pdf's fit, the six discretised DFE bins, which pdf the data
supports, and `num_nonneutral`.

---

## 5. Reading the result

Every stage writes a check file. Four statuses, and the distinction matters:

- **PASS** — nothing to look at.
- **WARN** — unusual compared with past runs, not wrong. The residual and
  deviance thresholds are percentiles of previously accepted fits; they flag a
  fit worth looking at. **They never block.** Look at the figure.
- **FAIL** — structural: nothing survived prep, a chromosome with intervals but
  no spectrum, too few SNPs, an unusable cache, no fits at all. Inside a chain
  this **stops the next job**.
- **SKIP** — not applicable, e.g. a chromosome the bed does not touch.

```bash
cat $OUT/prep/check_prep.$PREP_RUN.tsv
cat $OUT/sfs/check_sfs.$FS_RUN.tsv
cat $OUT/demo/check_demo.$FS_RUN-$DEMO_MODEL-$DEMO_RUN.tsv
cat $OUT/cache/<bed>/$DEMO_RUN-$CACHE_MARKER-$CACHE_RUN/check_cache.tsv
cat $OUT/dfe/check_dfe.$CACHE_RUN-$CACHE_MARKER-<pdf>-$DFE_RUN.tsv
cat $OUT/summary/check_summary.$SUMMARY_RUN.tsv
```

`fit_outcomes.tsv` in each DFE run directory records what became of every
starting point, with a severity saying what to do:

| outcome | severity | what to do |
|---|---|---|
| `bestfit`, converged | `ok` | nothing |
| `bestfit`, not converged | `fail_optimize` | the optimiser did not settle — change the bounds or the starting points |
| `no_bestfit` | `fail_optimize` | it optimised and reached no usable likelihood. Expected at the edge of the grid |
| `timeout` | `fail_optimize` | a bad starting point, not a slow one. **Not worth re-running** |
| `error` | `broken` | dadi-cli never started: cache path, `--demo-popt`, environment. Fix it; re-running changes nothing |
| `not_run` | `rerun` | the array never reached it, or the deadline guard deferred it. **The work a resume is for** |

---

## 6. When it goes wrong

**Array tasks killed by the clock.**

```bash
SUMMARY_RUN=<new> RESUME=1 FROM=dfe bash scripts/chains/run_all.sh their.bed
```

What `RESUME=1` actually skips is **a starting point with a `.bestfits`, and one
logged `timeout`**. Everything else is attempted again — so the `no_bestfit`
points go round once more. They are the fastest group (under 4 minutes each), so
this costs some wall clock and nothing else, but a resume is not as cheap as the
`not_run` count alone suggests.

The new `SUMMARY_RUN` is **required**, not tidiness: `run_summary.sh` and
`run_bestpdf.sh` refuse an existing table and ignore `RESUME`.

**Expect `fitneu` to log a failure.** It refuses its existing run directory and
does not honour `RESUME`. Harmless — the neutral likelihood is already there and
the summary reads it from that directory.

**`FORCE=1` means three different things** — bypass an upstream gate, redo every
DFE fit, overwrite a finished table. So `FORCE=1 RESUME=1` *discards* the fits
`RESUME` was meant to keep. To get past a gate only, run that stage by hand: a
hand-run stage consults no gate.

**A gate stopped a job.** The message names the upstream stage and its check
file. Read the file before overriding anything.

**Resuming is guarded in two ways:** the pipeline refuses to resume with a
different cache than the run was built with, and refuses when `RESUME=1` has
nothing to resume. Both guards exist because both failures happened and were
nearly invisible.

---

## 7. Do not

- **Invent a config value.** Ask, or use the default and say which you used.
- **Use `FORCE=1` to get past a gate without telling the user what the gate
  said.** It is in the output; quote it.
- **Change `GRIDS`, `DFE_BINS`, `CTRL_FLANK`, `PROJECTION` or the optimiser
  bounds silently.** They change results. Say what you changed and why.
- **Submit without `DRY_RUN=1` first**, or without checking the cluster copy is
  current.
- **Promise a number that matches a published table exactly.** This pipeline
  takes the *best* converged fit; the original took the last, which was the
  worst. That is the one deliberate numerical difference, ~1e-7 in practice but
  unbounded for a poorly converged file. See `todo/tocheck.md` item 2.
- **Claim the assembly was checked when it was not.** The preflight says so
  explicitly when it could not check.

---

## 8. Where things are

```
AGENTS.md              this file
docs/PIPELINE.md       per-stage detail: inputs, outputs, what each check means
config/pipeline.conf   every tunable, grouped by purpose, with defaults
config/profiles/       site paths. template.conf to start a new one
scripts/chains/        run_all.sh, the entry point
scripts/<stage>/       prep sfs demo cache dfe summary
scripts/sge/           submit_*.sh and worker_*.sh, one pair per array stage
scripts/checks/        check_*.py, one per stage, plus check_bed.py and gate.py
sync/push.sh           copy the package to a cluster; CHECK=1 compares
setup/setup.sh         build the conda env with the patched dadi
todo/                  what is known to be unfinished, and why
records/               what was actually run, dated
```

Known-unfinished work is in `todo/`, not hidden: `tocheck.md` (open questions
about the pipeline's own behaviour), `final_clean_rerun.md` (the acceptance
test), `vendor_and_docs.md` (what is still outside the package, and the `FORCE=1`
overloading). **Read `todo/` before telling a user this package is finished** —
each file records a dated incident, and they are the reason most of the guards
exist.

`records/` holds what was actually run, dated. `records/260929_dfe_validation.sh`
is the most useful single file for seeing a real invocation end to end, though
note it predates the `CACHE_MARKER` directory naming.

### Three things that look like breakage and are not

- **~22% of DFE starting points return `no_bestfit`.** The grid log-spaces shape
  from 1e-10, and the smallest seeds cannot produce a likelihood. By design.
- **dadi warns that 11 of 300 cache spectra are unusable.** Reproduced from the
  original, not introduced; larger grids were tested and did not help.
- **A `RESUME` run logs a `fitneu` failure.** `run_fitneu.sh` refuses an existing
  run directory and ignores `RESUME`. The neutral likelihood is already there.
