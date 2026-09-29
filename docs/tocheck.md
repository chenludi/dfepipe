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

**Hypothesis (CD).** The 50 kb window is too small for a large interval. The
window is `[mid − 25000, mid + 25000]` centred on the midpoint, and its width
does not grow with the region. So for a merged region longer than 50 kb the
window lies entirely *inside* the region. Controls are candidate-neutral sites
*minus* the test sites, so there can be no control inside the region, and the
region is dropped from `has_qui` — whereas its unmerged sub-intervals each had
their own window reaching out into flanking sequence and could find one.

That predicts merging should mostly *lose* sequence. chr16 fits. chr12 gains,
so there is a second mechanism: a merged region's window is centred on the
merged midpoint and can reach sequence that no sub-interval's window reached,
picking up a control that none of them found. Both effects are real; which
dominates depends on the local interval structure.

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
- [ ] Decide whether the window should scale with the region — e.g.
      `[start − FLANK, end + FLANK]` rather than `mid ± FLANK`. That is a
      different definition of "nearby", so it is a science decision, not a
      cleanup: it changes which controls pair with which test sites for every
      annotation. Do not change it without deciding deliberately.

**Do not act on this yet.** The pipeline keeps the current behaviour
(`MERGE_BEFORE_WINDOW=1`, windows on merged regions) until this is understood.

**Where the evidence is.** `config/pipeline.conf` (`CTRL_FLANK`,
`MERGE_BEFORE_WINDOW`), `scripts/sfs/run_sfs_chrom.sh` (the windowing), and the
original `../windctrl_splitchr.sh:41-42,68-76`. The published ELS SFS predates
the merge, which is why it matches the unmerged numbers.
