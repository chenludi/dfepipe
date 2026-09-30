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

**Partly answered, 2026-09-29.** The chr12 gain is not the "region keeps its
sub-intervals as a unit" effect. Two bookended ELS intervals,

    chr12  120230636  120230966   330 bp
    chr12  120230966  120231118   152 bp   <- the interval that was gained

merge into [120230636, 120231118] with midpoint 120230877, where the 152 bp
interval's own midpoint is 120231042. The window does not grow, it **moves
165 bp left**: [120206042, 120256042] unmerged vs [120205877, 120255877]
merged. A quiescent site sits in the sliver that move gains, so the interval
finds a control and is kept. Nothing was lost on chr12 at all (1 interval
gained, 0 lost).

So the general rule is that merging **relocates** the window centre, and a
relocated window can gain or lose controls either way. The "long region swallows
its own window" effect is a separate mechanism that only bites above
2 x CTRL_FLANK = 50 kb, and is still the thing worth counting.

Measured magnitudes for ELS: control set 602,593,295 merged vs 602,593,816
unmerged (-521 bp, and merging can only shrink the control set, since the union
of per-interval windows can only get smaller); test set 205,209,325 vs
205,209,348 (-23 bp net, from +152 on chr12 and -175 on chr16).

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
