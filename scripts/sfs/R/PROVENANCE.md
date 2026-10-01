# Vendored R scripts

`concatefs.R` and `concate_mut.R`, copied verbatim on 2026-10-01 from
`/u/home/c/cdi/project-klohmuel/noncdoingdfe/test_0118/script/` on Hoffman2,
which is where the original pipeline kept them.

They are shipped because they were the last executable dependency outside the
package: dadi and dadi-cli are vendored in `setup/vendor/`, every path is in the
config, but `run_concat.sh` shelled out to these two, so the package could not
run anywhere but that one cluster.

| file | sha256 (first 12) | lines |
|---|---|---|
| concatefs.R | 042613de7578 | 80 |
| concate_mut.R | 06d00c097584 | 86 |

**Not modified.** Not even the obvious defects, because these produced every
published spectrum and the point of vendoring is that the behaviour does not
change:

- `concatefs.R:57` only *prints* "Error: head or tail text not equal" when a
  per-chromosome spectrum has a different projection or population, then carries
  on - so a mismatched chromosome is silently left out of the sum rather than
  stopping the run. The pipeline's own `check_sfs.py` catches the symptom.
- `concate_mut.R:75` writes `length_inputbed` twice in its output row, once at
  position 5 and again at position 11 where the length of `inputbed` belongs.
  That is why the genome-wide fsmut table has a `length_inputbed.1` column
  holding the same value. `run_summary.sh` reads column 5 by name.
- `concate_mut.R:69` computes `out_colnames` and never uses it.

`LEGACY_SCRIPT_DIR` overrides the vendored copies, which is how they were
checked: a concat run against the originals and one against these must produce
byte-identical output.
