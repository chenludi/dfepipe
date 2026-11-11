# Example: ENCODE enhancer-like signatures (GRCh38-ELS)

A complete chain, prep through summary, on one bed. It is the example to read
first because its result can be checked against a published number rather than
only against itself: this run reproduces Data S2 sheet `S2a.enhancers`, the
"DFEs of ELS" row, to 1e-5 or better on every DFE quantity.

`expected/` holds the small output files from the run of 2026-10-25. The inputs
are not here - the bed alone is ~25 MB and the VCFs and cache are far larger -
so this documents what to run and what you should get, not a dataset you can
re-run offline.

## The input

`GRCh38-ELS.bed`, the ENCODE cCRE enhancer-like set at the `raw` level from the
`old` ENCODE release, which is the release the manuscript used. It was built by
`record_260928_ELS-ChromHMMProm_length.sh`, which lives outside this package.
937,426 intervals, 255,869,536 bp over chr1-22.

The preflight WARNs on it, both reasons benign and worth seeing because they are
typical:

* chrX and chrY are present, 2.5% of bp. `CHROMS` is chr1-22, so those bases are
  dropped at concat and the reported `l` is not the bed's own length. The
  published ELS run did the same.
* 66 intervals are out of coordinate order. prep sorts anyway.

No interval exceeds `2 x CTRL_FLANK`; the longest is 917 bp. So the long-interval
control-window problem (`TODO.md`, open question 1) does not arise here.

## The command

```bash
source config/profiles/<yoursite>.conf
bash scripts/chains/run_all.sh /path/to/GRCh38-ELS.bed
```

Everything else is the default in `config/pipeline.conf`, which reproduces the
manuscript. The run keys used here were `prep_261025`, `fs_261025`,
`demo_261025`, `cache_261025`, `dfe_261025`, `sum_261025`; set your own with
those variables.

The real run was three beds in one chain (ELS and its two ChromHMM-promoter
subsets) - see `records/261025_els_prom_run.sh`, which is the record of it.

## What you should get

Demography, `expected/demo/best_inference.txt` - `three_epoch`, picked from 3
starting points by likelihood:

| | |
|---|---|
| ll | -517.0139615073276 |
| nuB, nuF | 1.7256242494267395, 1000.0 |
| TB, TF | 0.7495047456939522, 0.0056048427302586995 |
| theta | 470803.5407379631 |
| Na | 18364.008247474252 |

`nuF` sitting exactly at 1000.0 is the upper bound of its search range, not a
coincidence; it is what the published run also found.

The DFE, `expected/dfe/spec.bestdfe.txt` - gamma, from 152 fit files of which 66
converged, 4 kept above the likelihood cutoff, best `129.InferDFE.bestfits` at
ll -477.47922788115284.

The discrete DFE from those best-fit parameters,
`expected/dfe/dfe_exHET_noq.GRCh38-ELS_qui2.chrs_1t22_gamma.csv` and the bar
chart beside it, over `DFE_BINS = 0.0 1e-05 0.0001 0.001 0.01 0.5 inf`:

| \|s\| | proportion |
|---|---|
| 0 to 1e-05 | 0.8370603106807804 |
| 1e-05 to 1e-04 | 0.05515881035257475 |
| 1e-04 to 1e-03 | 0.05684737790114358 |
| 1e-03 to 1e-02 | 0.04423424443643853 |
| 1e-02 to 0.5 | 0.006699256629062789 |
| 0.5 and above | 0.0 |

These six numbers are also the `delta_50_*` columns of
`expected/summary/bestpdf_sumdfe.tsv`, identical to 12 significant figures,
because the delta-50 rule chose gamma for this bed
(`pdf_plotting_delta_50 = gamma`). Where the rule chooses a different pdf the two
will not match, and the summary column is the one the manuscript figures use.

## Reading the rest of `expected/`

One check file per stage, each the stage's own verdict on itself:
`prep/check_prep.*`, `sfs/check_sfs.*`, `demo/check_demo.*`,
`cache/check_cache.tsv`, `dfe/check_dfe.*`, `summary/check_summary.*`. A stage
that cannot write a PASS row stops the chain, so these are what to look at first
when a run stops.

Two spectra plots: `demo/sfs_bestfit.pdf`, the demographic fit against the
control spectrum, and `dfe/sfs_rank1_poisson.pdf`, the best DFE fit against the
test spectrum.

The absolute paths inside `best_inference.txt` and `spec.bestdfe.txt` are from
the machine this ran on. Yours will differ; nothing reads them back.
