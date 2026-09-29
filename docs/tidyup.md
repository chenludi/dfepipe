# Disk policy

**Nothing is deleted automatically except each task's own temporaries.**

That is a deliberate choice, not an oversight. The optimisation files are the
record of how an estimate was reached, and are kept even after publication so a
result can be traced back to the fit that produced it.

## Deleted automatically

Each stage builds inside a `mktemp -d` under its destination and moves the
finished files into place, so an interrupted job cannot leave a truncated file
at its final name. The temp directory goes on exit, including:

- **the masked VCF copy**, which is the only large one. dadi's
  `make_data_dict_vcf()` takes a file path and reads the whole file, so
  restricting a spectrum to a bed means subsetting the VCF first
  (`bedtools intersect -a <vcf> -b <bed> -header | gzip`, as
  `get_fs_mut_bychr.sh:57` does). There is no way to hand dadi a filtered
  stream, so the copy is unavoidable — but it is per-task and short-lived.
  The original wrote it into the shared reference directory, where 22
  concurrent array tasks could collide over it.
- the per-chromosome intermediate beds of the windowing step (`sorted`,
  `merged`, `extest`, `window`, `q.window`).

## Kept

Everything else, indefinitely:

| what | why |
|---|---|
| `*.InferDFE.opts.*`, `*.bestfits` | how each fit was reached; kept for traceability after publication |
| cache `.bpkl` | 10 h on 10 cores to rebuild |
| genome-wide `.fs`, `fsmut.*.tab` | the inputs every later stage reads |
| `best_inference.txt` | demography feeding every downstream stage |
| `check_*.tsv`, `spec.txt`, `lengths.tsv` | what ran, with what config, and whether it looked right |
| per-chromosome `.fs` / `.tab` | 132 KB per annotation; deleting them saves nothing |
| per-chromosome beds | 67 MB per annotation; the only real candidate, and still small |

## If space ever runs short

Measured for ELS, 2026-09-29: the whole `sfs` run's per-chromosome beds are
67 MB per annotation, while `result_dadi` alone is 4.0 GB. So reclaiming `sfs`
intermediates is effort spent on a rounding error. Order of least to most
regret:

1. per-chromosome beds of a run whose next stage has already passed its check
   (one array job to rebuild, ~2 min per chromosome)
2. nothing else without a deliberate decision

Runs are keyed by `PREP_RUN`, `FS_RUN`, and so on, so old runs accumulate rather
than being overwritten. If a directory is genuinely obsolete, delete the whole
run directory rather than picking files out of it — a half-emptied run is worse
than either keeping or removing it.
