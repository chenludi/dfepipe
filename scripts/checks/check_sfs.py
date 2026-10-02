#!/usr/bin/env python3
"""Checkpoint for the sfs stage: did every chromosome that should have produced
an SFS produce one, is there enough data to fit, and how much of the test set
survived control pairing?

The stage runs as an SGE array job, one task per chromosome, so a task that is
killed or times out leaves a hole. The concatenation script that produced the
published results caught it - it demanded exactly 22 files and wrote nothing
otherwise - but did not say WHICH was absent, so finding it meant checking by
hand, as record_250624_243-40-45.sh records doing for chr4 of phyloP243 40-45.
This check therefore reports per chromosome and names the missing ones.

The check file carries one row per "<bed>:<chrom>" plus one "<bed>" summary row.
Items must be unique - checkfile.write_check enforces it - so a later row can
never hide an earlier failure.

Per chromosome:
    SKIP  the prep bed has no intervals on that chromosome. Normal, not an
          error: an annotation may simply be absent from a chromosome.
    FAIL  it HAS intervals but no .fs was produced. Either the array task was
          killed or timed out (look in the SGE log directory), or the VCF does
          not cover that chromosome.
    PASS  otherwise.

Summary row:
    FAIL  any chromosome FAILed, or fewer than MIN_SNPS SNPs in total. Below
          MIN_SNPS the spectrum is too thin to fit anything, so this blocks.
    FAIL  the has_qui bp is under HASQUI_MIN_FRAC of what prep produced: test
          sites are being lost for want of a candidate-neutral site inside
          [mid-CTRL_FLANK, mid+CTRL_FLANK], so whatever survives is a biased
          subset. This DOES block the chain, like every FAIL (gate.py); FORCE=1
          proceeds. An earlier version of this docstring said it never blocks,
          which was wrong.

FOR REFERENCE ONLY applies to the WARN, as everywhere in this pipeline: it means
"unusual compared with past runs", not "wrong".
"""
import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from conf import conf  # noqa: E402
# One definition of how prep names its output, shared rather than restated: if
# prep's naming ever changes, this check must not keep looking for the old name.
from check_prep import OUTPUT_TEMPLATE as PREP_BED_TEMPLATE  # noqa: E402

# The manuscript's runs kept at least 90% of test sites after control pairing.
MIN_HASQUI_FRAC = conf.get_float("HASQUI_MIN_FRAC", 0.90)
# Smallest published bin (PLS phastCons17 20-25) has 1,560 SNPs. Mirrors
MIN_SNPS = conf.get_int("MIN_SNPS", 500)

# Rows run_sfs_chrom.sh writes into bychr/lengths.<chrom>.tsv.
LENGTH_STAGES = ["prep_chr", "merged_chr", "has_qui", "quiecent_ctrl"]


def normalise_chrom(chrom):
    """'1' and 'chr1' both mean chr1.

    CHROMS in pipeline.conf is bare numbers ("1 2 3 ...") while every bed, VCF
    and output name uses the chr-prefixed form, so the two have to meet here.
    """
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else "chr" + chrom


# ------------------------------------------------------------------ verdicts --

def assess_chrom(n_intervals, has_fs, log_dir=None, chrom=""):
    """(status, reasons) for one <bed>:<chrom>. See the module docstring."""
    if not n_intervals:
        return "SKIP", ["no intervals in bed"]
    if not has_fs:
        where = log_dir or "the sfs SGE log directory"
        return "FAIL", [
            "the bed has %d interval(s) on %s but no .fs was produced: either the "
            "array task was killed or timed out (check the SGE logs in %s, and "
            "SGE_SFS_RT / SGE_SFS_MEM), or the VCF has no records for %s"
            % (n_intervals, chrom or "this chromosome", where, chrom or "it")]
    return "PASS", []


def assess_summary(n_snps, hasqui_ratio, failed_chroms,
                   min_snps=MIN_SNPS, min_hasqui_frac=MIN_HASQUI_FRAC):
    """(status, reasons) for the <bed> summary row.

    Every reason is collected, not just the first: a run that is both short of
    SNPs and short of controls should say so once, not across two re-runs.
    """
    reasons = []
    status = "PASS"

    if failed_chroms:
        status = "FAIL"
        reasons.append("no SFS for chromosome(s) that have intervals: %s"
                       % " ".join(failed_chroms))
    if n_snps is None:
        status = "FAIL"
        reasons.append("no SFS could be read, so the SNP count is unknown")
    elif n_snps < min_snps:
        status = "FAIL"
        reasons.append("only %.0f SNPs (minimum %d): too thin to fit"
                       % (n_snps, min_snps))
    if hasqui_ratio is not None and hasqui_ratio < min_hasqui_frac:
        # FAIL, not WARN: if a large share of the test sequence has no nearby
        # control it is dropped silently, and whatever survives is a biased
        # subset. Two causes worth stopping for - regions too long for a
        # midpoint-centred window to reach outside them, or too little
        # candidate-neutral sequence nearby. FORCE=1 proceeds anyway.
        status = "FAIL"
        reasons.append(
            "only %.1f%% of the prep'd test sequence has a nearby control "
            "(minimum %.1f%%): either the test regions are too long for a "
            "+/-CTRL_FLANK window centred on their midpoint, or there is too "
            "little candidate-neutral sequence near them"
            % (100 * hasqui_ratio, 100 * min_hasqui_frac))

    return status, reasons


# ------------------------------------------------------------------- readers --

def bed_intervals_by_chrom(path):
    """{chrom: (n_intervals, bp)} for a bed. A missing file reads as empty."""
    out = {}
    try:
        fh = open(path)
    except OSError:
        return out
    with fh:
        for line in fh:
            if not line.strip() or line[0] in "#t":   # blanks, comments, 'track'
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            try:
                span = abs(int(f[2]) - int(f[1]))
            except ValueError:
                continue
            n, bp = out.get(f[0], (0, 0))
            out[f[0]] = (n + 1, bp + span)
    return out


def fs_total(path):
    """sites_fsfile: the sum of a dadi .fs file's counts row, or None.

    Same quantity dadi_sfs_YRI.py recorded and the same one concatefs.R adds up
    across chromosomes, so the per-chromosome sum here equals what the
    genome-wide file holds.
    """
    try:
        with open(path) as fh:
            lines = fh.readlines()
    except OSError:
        return None
    if len(lines) < 2:
        return None
    try:
        return sum(float(x) for x in lines[1].split())
    except ValueError:
        return None


def read_lengths(path):
    """{stage: bp} from a lengths.<chrom>.tsv, or {} if unreadable."""
    try:
        with open(path) as fh:
            return {r["stage"]: int(r["bp"])
                    for r in csv.DictReader(fh, delimiter="\t") if r.get("stage")}
    except (OSError, KeyError, TypeError, ValueError):
        return {}


def _bed_base(bed):
    base = os.path.basename(bed)
    return base[:-4] if base.endswith(".bed") else base


def find_fs(bychr, prefix, chrom, bed_stem):
    """The per-chromosome .fs, whose name embeds the VCF's.

    get_fs_mut_bychr.sh named it <bedstem>.<vcf basename minus .vcf.gz>.fs, so
    the tail depends on VCF_TMPL and is globbed rather than reconstructed.
    """
    hits = sorted(glob.glob(os.path.join(
        bychr, "smer.%s.%s.%s.*.fs" % (prefix, chrom, bed_stem))))
    return hits[0] if hits else None


# ---------------------------------------------------------------------- rows --

def rows_for(out_dir, run, prep_out, prep_run, bed, chroms,
             log_dir=None, min_snps=MIN_SNPS, min_hasqui_frac=MIN_HASQUI_FRAC):
    """Every row for one bed: one per chromosome, then the summary."""
    base = _bed_base(bed)
    bed_stem = PREP_BED_TEMPLATE.format(base=base)[:-4]     # drop '.bed'
    prep_dir = os.path.join(prep_out, base, prep_run)
    prep_bed = os.path.join(prep_dir, PREP_BED_TEMPLATE.format(base=base))
    bychr = os.path.join(out_dir, base, run, "bychr")

    by_chrom = bed_intervals_by_chrom(prep_bed)
    prep_bp = read_lengths(os.path.join(prep_dir, "lengths.tsv")).get("after_mask")

    if not os.path.isfile(prep_bed):
        return [{"item": base, "status": "FAIL",
                 "reason": "no prep bed at %s: run the prep stage first, or "
                           "check PREP_RUN" % prep_bed}]

    rows = []
    failed, missing, with_data = [], [], 0
    n_snps, bp_has_qui, bp_ctrl = 0.0, 0, 0
    any_fs = False

    for chrom in chroms:
        n_intervals, bp_chrom = by_chrom.get(chrom, (0, 0))
        fs = find_fs(bychr, "has_qui", chrom, bed_stem)
        status, reasons = assess_chrom(n_intervals, fs is not None,
                                       log_dir=log_dir, chrom=chrom)
        lengths = read_lengths(os.path.join(bychr, "lengths.%s.tsv" % chrom))
        snps = fs_total(fs) if fs else None

        if status == "FAIL":
            failed.append(chrom)
        if fs is None:
            missing.append(chrom)
        else:
            with_data += 1
            any_fs = True
            if snps is not None:
                n_snps += snps
        bp_has_qui += lengths.get("has_qui", 0)
        bp_ctrl += lengths.get("quiecent_ctrl", 0)

        rows.append({
            "item": "%s:%s" % (base, chrom),
            "status": status,
            "reason": "; ".join(reasons) or "NA",
            "n_intervals_prep": n_intervals,
            "bp_prep": bp_chrom,
            "bp_has_qui": lengths.get("has_qui", "NA"),
            "bp_quiecent_ctrl": lengths.get("quiecent_ctrl", "NA"),
            # Per chromosome, against that chromosome's own prep bp - so one bad
            # chromosome is visible instead of being averaged away. The original
            # divided a per-chromosome numerator by the genome-wide bed length
            # (windctrl_splitchr.sh:122), giving ~1/22 of the truth.
            "hasqui_ratio": (round(lengths["has_qui"] / bp_chrom, 4)
                             if bp_chrom and "has_qui" in lengths else "NA"),
            "n_snps": "NA" if snps is None else snps,
        })

    hasqui_ratio = (bp_has_qui / prep_bp) if prep_bp else None
    status, reasons = assess_summary(n_snps if any_fs else None, hasqui_ratio, failed,
                                     min_snps=min_snps, min_hasqui_frac=min_hasqui_frac)
    rows.append({
        "item": base,
        "status": status,
        "reason": "; ".join(reasons) or "NA",
        "n_chroms_expected": len(chroms),
        "n_chroms_with_data": with_data,
        "missing_chroms": " ".join(missing) or "NA",
        "n_snps": n_snps if any_fs else "NA",
        "bp_prep": prep_bp if prep_bp is not None else "NA",
        "bp_has_qui": bp_has_qui,
        "bp_quiecent_ctrl": bp_ctrl,
        "hasqui_ratio": "NA" if hasqui_ratio is None else round(hasqui_ratio, 4),
    })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="the sfs output root, i.e. $OUT/sfs")
    ap.add_argument("--run", required=True, help="FS_RUN, the run key")
    ap.add_argument("--prep-out", required=True, help="the prep output root, $OUT/prep")
    ap.add_argument("--prep-run", required=True, help="PREP_RUN")
    ap.add_argument("--chroms", required=True,
                    help="space-separated chromosomes; '1' and 'chr1' both work")
    ap.add_argument("--log-dir", default=None,
                    help="where the SGE array logs are, named in a FAIL reason")
    ap.add_argument("--min-snps", type=int, default=MIN_SNPS)
    ap.add_argument("--min-hasqui-frac", type=float, default=MIN_HASQUI_FRAC)
    ap.add_argument("beds", nargs="+", help="the input beds prep was run on")
    a = ap.parse_args(argv)

    chroms = [normalise_chrom(c) for c in a.chroms.split()]
    rows = []
    for bed in a.beds:
        rows.extend(rows_for(a.out, a.run, a.prep_out, a.prep_run, bed, chroms,
                             log_dir=a.log_dir, min_snps=a.min_snps,
                             min_hasqui_frac=a.min_hasqui_frac))

    os.makedirs(a.out, exist_ok=True)
    dest = os.path.join(a.out, "check_sfs.%s.tsv" % a.run)
    write_check(dest, rows)
    print("wrote %s" % dest)
    for r in rows:
        if r["status"] not in ("PASS", "SKIP"):
            print("  %s %s: %s" % (r["status"], r["item"], r["reason"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
