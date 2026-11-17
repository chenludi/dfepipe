# Round 1, scenario 3: wrong assumptions

**Request.** A user wanted to run an hg19 bed that included chrX and chrY, said
they were really interested in mouse data, and asked for the `bottlegrowth`
demographic model. The report had to list every incompatibility with its source,
the commands if resolved, questions to ask, and every confusion.

**What the agent got right.**
- `bottlegrowth` is refused by the model gate in `config/pipeline.conf`
  because no starting-point grid exists for it.
- chrX/chrY are not analysed, and found they are dropped silently at concat,
  so the reported length is shorter than the input bed.
- hg19 vs hg38: every reference is GRCh38 and nothing checks the assembly, so
  a mismatch gives plausible, wrong output.
- Mouse would need a full replacement reference set (VCF, panel, mask, CDS,
  chromatin states, controls, mutation rates); it is a porting job.

**What it could not find or had to guess.**
- No README or `docs/PIPELINE.md`; everything came from script headers.
- `bottlegrowth` was half-implemented: several stages still mapped it, and a
  comment in `scripts/demo/run_demo.sh` claimed it was accepted.
- `ROULETTE_SCALE` was configurable in name only: hardcoded in three of five
  consumers, so changing it would give inconsistent scaling silently.
- `MUTRATE_PATTERN` could not be overridden from the environment despite the
  config header promising every value could be; `MUTRATE_TMPL` was undocumented.
- No validation of the input bed: assembly, contig naming, column count, or
  contigs outside `CHROMS`.
- The legacy R scripts for concat and the `dadi_initialparams/` grids lived
  outside the package with no manifest of what must exist.
- `PROJECTION` was never checked against the panel size; too high drops every
  site silently.
- The literal `chrs_1t22` in output names stays even when `CHROMS` is a subset.

**Outcome.** Found all four incompatibilities and two real bugs; the agent was
fairly confident it had found the hard constraints, less so for the rest.
