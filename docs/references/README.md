# References

Verbatim quotes from primary sources, for any claim in this package that rests
on one. A document here carries the quote, the source URL and the section it
came from; the pipeline's own files cite it in short form and point here.

The rule is that a number or a method attributed to a paper, a tool's
documentation or a spec is verified against that source and quoted, rather than
paraphrased from memory or from another citation.

Claims in this package that still need a source quoted here:

- `ROULETTE_SCALE = 1.015e-7` (`config/pipeline.conf`), which converts a
  Roulette table value to a per-site diploid rate. Currently cited only to
  `dadi_cal_ll_neu_0723.py:160`, which is where the pipeline got it, not where
  it comes from. Needs the Roulette paper or its documentation.
- the `1e-4` threshold below which `neugamma` treats a selection coefficient as
  neutral (`setup/patches/dadi_DFE_PDFs.patch`). Cited to the patch itself.
- the `three_epoch` model's definition and parameter order, currently taken from
  dadi's own source.

Nothing here yet, which is why the list above exists.
