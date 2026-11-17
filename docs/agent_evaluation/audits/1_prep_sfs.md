# Audit 1: prep and sfs stages

**Scope.** `scripts/prep/run_prep.sh`, `scripts/sfs/run_sfs_chrom.sh` and `scripts/sfs/run_concat.sh`, checked claim by claim against the original scripts they port (`funsp_exhetqui.sh`, `filt_nc-unmask.sh`, `windctrl_splitchr.sh`, `get_fs_mut_bychr.sh`, `piperef_stddfe_2_concatefs.sh`, `concatefs.R`, `concate_mut.R`). Read-only; nothing was run.

Every line-number citation in the three headers resolved correctly. The problems are in the prose around them and in behaviour the headers do not mention.

**Main findings.**

- **Contradictory claim about `MERGE_BEFORE_WINDOW`** (`run_sfs_chrom.sh`). The header says pre-2024/11/04 results need `MERGE_BEFORE_WINDOW=0` to reproduce, but prep already merges its output, so the bed reaching sfs is merged either way. Runs from before that date cannot be reproduced at all. The same mistake appeared in the profile, `AGENTS.md` and `todo/tocheck.md`. **Could change a result:** merging moves the midpoints of adjacent intervals and so the ±25 kb control windows (measured +152 bp on chr12 and −175 bp on chr16 for ELS; not measured for long-interval annotations).

- **Wrong description of the control set** (`run_prep.sh`). The header says the quiescent sites `QUI_BED` "become the neutral control set". The control set is actually `NEUTRAL_CTRL_DIR`, a stricter subset, and `AGENTS.md` itself warns against mixing the two up. Documentation only.

- **concat now accepts fewer than 22 chromosomes** (`run_concat.sh`). The original wrote nothing unless all 22 were present. The port skips chromosomes with no sites and sums the rest under a filename that still says `chrs_1t22`. This is recorded in `concat_chroms.tsv` but missing from the Differences list. **Could change a result.**

- **False story about how a chromosome went missing** (`run_concat.sh`, repeated in `check_sfs.py`). The header says chr4 "went missing from one bin unnoticed". The original record says it was noticed and re-run, and the original concat cannot produce a genome-wide file with a chromosome missing.

- **Overstated check in concat**: `concatefs.R` only prints an error when headers differ and then silently leaves that chromosome out of the sum. The header says it checks the headers.

- **Vendored R scripts are not byte-identical**: each copy has an extra trailing blank line, so "Not modified" and the recorded hashes are not exactly right. A provenance line citation is also off (`:69` should be `:73`). Harmless.

- **Wrong scripts named as the source of prep**. The commands themselves are byte-identical to the published route, so behaviour is unaffected.

- **Undocumented deviations**: the prep merge reduces the bed to 3 columns; `MERGE_BEFORE_WINDOW=0` sorts (the original neither sorted nor merged); coverage is computed with awk instead of `bc -l` (differences below the 5th decimal); `.track`/`.log` bookkeeping files are dropped; stages refuse existing output and exit non-zero where the original overwrote and exited 0; a multi-match `.fs` glob is silently resolved with `head -1`.

- **Invented chain-blocking thresholds**: `MIN_SNPS=500` (its cited "1,560 SNPs" figure could not be found) and `HASQUI_MIN_FRAC=0.90` (cited to "the manuscript's runs", which could not be verified). Both stop the chain on FAIL. `check_sfs.py`'s docstring says the has_qui FAIL "never blocks" while the code returns FAIL. Smaller invented WARN-only values: `PREP_MIN_RETAIN_FRAC=0.05` and `OUT_OF_RANGE_SLACK=1e6`.

- **`MIN_FS_BYTES=661` only applies at `PROJECTION=160`**. Nothing ties the two together, so changing the projection silently invalidates the threshold.

**Not checked.**

- Anything on the cluster: the reference beds, the `HET_BED`/`QUI_BED` path (the profile uses a `temp/` subdirectory the originals do not; if that resolves to a different file, every number moves), and the bp totals.

- Nothing was executed, so the merge measurements and R-script equivalence were not reproduced.

- The manuscript tables behind `HASQUI_MIN_FRAC` and `MIN_SNPS`.
