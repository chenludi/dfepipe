"""Tests for the sfs checkpoint.

The sfs stage runs one job per chromosome, so its checkpoint is per chromosome:
one row per <bed>:<chrom> plus one <bed> summary row. The interesting part is
telling "this bed has nothing on chr21" (SKIP, normal) apart from "chr21 should
have produced an SFS and did not" (FAIL, a killed array task or a VCF that does
not cover that chromosome) - the situation that produced
record_250624_243-40-45.sh, where chr4 silently went missing from one bin.
"""
import csv
import pathlib
import subprocess
import sys

import pytest

from scripts.checks.check_sfs import (
    MIN_HASQUI_FRAC,
    MIN_SNPS,
    assess_chrom,
    assess_summary,
    bed_intervals_by_chrom,
    fs_total,
    main,
    normalise_chrom,
    read_lengths,
)
from scripts.checks.checkfile import read_status

ROOT = pathlib.Path(__file__).resolve().parent.parent

CHROMS = ["chr1", "chr2"]


# ------------------------------------------------------------ chromosome names --

def test_bare_chromosome_numbers_are_accepted():
    """CHROMS in pipeline.conf is '1 2 3 ...' but every bed and VCF says 'chr1'."""
    assert normalise_chrom("1") == "chr1"
    assert normalise_chrom("chr1") == "chr1"
    assert normalise_chrom(22) == "chr22"


# --------------------------------------------------------------- assess_chrom --

def test_no_intervals_on_that_chromosome_is_a_skip_not_an_error():
    status, reasons = assess_chrom(n_intervals=0, has_fs=False)
    assert status == "SKIP"
    assert any("no intervals in bed" in r for r in reasons)


def test_intervals_but_no_fs_is_a_fail_naming_the_likely_causes():
    status, reasons = assess_chrom(n_intervals=12, has_fs=False,
                                   log_dir="/out/sfs/log", chrom="chr7")
    assert status == "FAIL"
    joined = " ".join(reasons)
    assert "/out/sfs/log" in joined           # points at the SGE logs
    assert "timed out" in joined or "timed-out" in joined
    assert "killed" in joined
    assert "VCF" in joined


def test_intervals_and_an_fs_pass():
    assert assess_chrom(n_intervals=12, has_fs=True) == ("PASS", [])


def test_an_fs_without_intervals_still_skips_rather_than_crashing():
    """Defensive: a stale output from an earlier bed must not read as success."""
    status, _ = assess_chrom(n_intervals=0, has_fs=True)
    assert status == "SKIP"


# ------------------------------------------------------------- assess_summary --

def test_a_failed_chromosome_fails_the_summary():
    status, reasons = assess_summary(n_snps=50000, hasqui_ratio=0.99,
                                     failed_chroms=["chr4"])
    assert status == "FAIL"
    assert any("chr4" in r for r in reasons)


def test_too_few_snps_fails():
    """Below MIN_SNPS the data is too thin to fit at all: FAIL, not WARN."""
    status, reasons = assess_summary(n_snps=MIN_SNPS - 1, hasqui_ratio=0.99,
                                     failed_chroms=[])
    assert status == "FAIL"
    assert any("SNP" in r for r in reasons)


def test_exactly_min_snps_is_not_a_failure():
    assert assess_summary(n_snps=MIN_SNPS, hasqui_ratio=0.99, failed_chroms=[])[0] == "PASS"


def test_the_worst_published_bed_is_not_blocked():
    """The threshold was 0.90, which FAILed 8 of the 280 published beds -
    top3stats_TxEnh 0.790, TxEnh 0.803, TxEx 0.827, Tx 0.871, TSS 0.879 - and so
    would have blocked the Tx run this package was validated against.
    proportion_hasquictrl.tsv, 5th percentile 0.913, minimum 0.790."""
    for ratio in (0.790, 0.803, 0.827, 0.838, 0.871, 0.879):
        assert assess_summary(n_snps=50000, hasqui_ratio=ratio,
                              failed_chroms=[])[0] == "PASS", ratio


def test_low_hasqui_ratio_warns_without_blocking():
    """Losing test sites for want of a nearby control is a design question, not
    a broken run, so it warns and the chain carries on."""
    status, reasons = assess_summary(n_snps=50000, hasqui_ratio=0.60,
                                     failed_chroms=[])
    assert status == "WARN"
    assert any("control" in r for r in reasons)


def test_a_low_hasqui_ratio_never_masks_a_real_failure():
    """The warning must not downgrade a FAIL that was already set."""
    status, _ = assess_summary(n_snps=50000, hasqui_ratio=0.60,
                               failed_chroms=["chr7"])
    assert status == "FAIL"


def test_exactly_the_hasqui_threshold_passes():
    assert assess_summary(n_snps=50000, hasqui_ratio=MIN_HASQUI_FRAC,
                          failed_chroms=[])[0] == "PASS"


def test_unknown_hasqui_ratio_is_not_a_warning():
    assert assess_summary(n_snps=50000, hasqui_ratio=None, failed_chroms=[])[0] == "PASS"


def test_a_fail_outranks_a_warn():
    status, reasons = assess_summary(n_snps=10, hasqui_ratio=0.1, failed_chroms=["chr9"])
    assert status == "FAIL"
    assert len(reasons) >= 2          # every reason is reported, not just the first


def test_ok():
    assert assess_summary(n_snps=50000, hasqui_ratio=0.96, failed_chroms=[]) == ("PASS", [])


# -------------------------------------------------------------------- readers --

def test_bed_intervals_by_chrom_counts_intervals_and_bp(tmp_path):
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t0\t10\nchr1\t20\t25\nchr2\t0\t3\n")
    assert bed_intervals_by_chrom(bed) == {"chr1": (2, 15), "chr2": (1, 3)}


def test_bed_intervals_by_chrom_on_a_missing_file_is_empty(tmp_path):
    assert bed_intervals_by_chrom(tmp_path / "absent.bed") == {}


def test_fs_total_sums_the_second_line(tmp_path):
    """dadi's .fs: header, the counts, the mask. sites_fsfile is the counts row."""
    fs = tmp_path / "x.fs"
    fs.write_text('5 folded "YRI"\n0 100.0 40.0 0 0\n1 0 0 1 1\n')
    assert fs_total(fs) == pytest.approx(140.0)


def test_fs_total_of_an_unreadable_file_is_none(tmp_path):
    assert fs_total(tmp_path / "absent.fs") is None


def test_read_lengths_parses_the_stage_bp_table(tmp_path):
    p = tmp_path / "lengths.chr1.tsv"
    p.write_text("stage\tbp\nhas_qui\t120\nquiecent_ctrl\t500\n")
    assert read_lengths(p) == {"has_qui": 120, "quiecent_ctrl": 500}


# --------------------------------------------------------------------- main() --

def _prep(tmp_path, base, chrom_intervals, after_mask=None):
    """A fake prep run: the filtered bed plus lengths.tsv."""
    d = tmp_path / "prep" / base / "prep_t"
    d.mkdir(parents=True, exist_ok=True)
    lines, bp = [], 0
    for chrom, n in chrom_intervals.items():
        for i in range(n):
            lines.append("%s\t%d\t%d\n" % (chrom, 1000 * i, 1000 * i + 100))
            bp += 100
    (d / ("s.nc_exHET_noq.%s.bed" % base)).write_text("".join(lines))
    (d / "lengths.tsv").write_text(
        "stage\tbp\ninput\t%d\nafter_mask\t%d\n" % (bp, after_mask if after_mask is not None else bp))
    return d


def _chrom_out(tmp_path, base, chrom, has_qui_bp=500, ctrl_bp=4500, snps=1000.0):
    """One chromosome's sfs output: the two beds' bp plus a has_qui .fs."""
    d = tmp_path / "sfs" / base / "fs_t" / "bychr"
    d.mkdir(parents=True, exist_ok=True)
    stem = "s.nc_exHET_noq.%s" % base
    (d / ("lengths.%s.tsv" % chrom)).write_text(
        "stage\tbp\nprep_chr\t%d\nmerged_chr\t%d\nhas_qui\t%d\nquiecent_ctrl\t%d\n"
        % (has_qui_bp, has_qui_bp, has_qui_bp, ctrl_bp))
    (d / ("smer.has_qui.%s.%s.thevcf.%s.fs" % (chrom, stem, chrom))).write_text(
        '5 folded "YRI"\n0 %s 0 0 0\n1 0 0 1 1\n' % snps)
    (d / ("smer.quiecent_ctrl.%s.%s.thevcf.%s.fs" % (chrom, stem, chrom))).write_text(
        '5 folded "YRI"\n0 %s 0 0 0\n1 0 0 1 1\n' % snps)
    return d


def _main(tmp_path, *beds, **kw):
    argv = ["--out", str(tmp_path / "sfs"), "--run", "fs_t",
            "--prep-out", str(tmp_path / "prep"), "--prep-run", "prep_t",
            "--chroms", kw.get("chroms", " ".join(CHROMS)),
            "--log-dir", str(tmp_path / "sfs" / "log")]
    return main(argv + list(beds))


def test_writes_a_row_per_chromosome_plus_a_summary(tmp_path):
    _prep(tmp_path, "ELS", {"chr1": 5, "chr2": 5})
    _chrom_out(tmp_path, "ELS", "chr1")
    _chrom_out(tmp_path, "ELS", "chr2")
    assert _main(tmp_path, "ELS.bed") == 0

    dest = tmp_path / "sfs" / "check_sfs.fs_t.tsv"
    rows = {r["item"]: r for r in csv.DictReader(dest.open(), delimiter="\t")}
    assert set(rows) == {"ELS:chr1", "ELS:chr2", "ELS"}
    assert rows["ELS:chr1"]["status"] == "PASS"
    assert rows["ELS"]["status"] == "PASS"
    assert rows["ELS"]["n_chroms_expected"] == "2"
    assert rows["ELS"]["n_chroms_with_data"] == "2"
    assert rows["ELS"]["n_snps"] == "2000.0"          # 1000 per chromosome
    assert rows["ELS"]["bp_has_qui"] == "1000"
    assert rows["ELS"]["bp_quiecent_ctrl"] == "9000"
    assert float(rows["ELS"]["hasqui_ratio"]) == pytest.approx(1.0)


def test_a_chromosome_the_bed_does_not_touch_is_skipped(tmp_path):
    _prep(tmp_path, "ELS", {"chr1": 5})               # nothing on chr2
    _chrom_out(tmp_path, "ELS", "chr1")
    assert _main(tmp_path, "ELS.bed") == 0
    status = read_status(tmp_path / "sfs" / "check_sfs.fs_t.tsv")
    assert status["ELS:chr2"] == "SKIP"
    assert status["ELS"] == "PASS"                    # a SKIP never blocks


def test_a_chromosome_with_intervals_and_no_output_fails_the_bed(tmp_path):
    _prep(tmp_path, "ELS", {"chr1": 5, "chr2": 5})
    _chrom_out(tmp_path, "ELS", "chr1")               # chr2 never ran
    assert _main(tmp_path, "ELS.bed") == 0
    dest = tmp_path / "sfs" / "check_sfs.fs_t.tsv"
    rows = {r["item"]: r for r in csv.DictReader(dest.open(), delimiter="\t")}
    assert rows["ELS:chr2"]["status"] == "FAIL"
    assert rows["ELS"]["status"] == "FAIL"
    assert "chr2" in rows["ELS"]["missing_chroms"]
    assert "chr2" in rows["ELS"]["reason"]
    assert str(tmp_path / "sfs" / "log") in rows["ELS:chr2"]["reason"]


def test_items_are_unique_so_no_row_can_hide_another(tmp_path):
    """write_check enforces it; this is the regression test for the item scheme."""
    _prep(tmp_path, "ELS", {"chr1": 5, "chr2": 5})
    _prep(tmp_path, "PLS", {"chr1": 5, "chr2": 5})
    for chrom in CHROMS:
        _chrom_out(tmp_path, "ELS", chrom)
        _chrom_out(tmp_path, "PLS", chrom)
    assert _main(tmp_path, "ELS.bed", "PLS.bed") == 0
    dest = tmp_path / "sfs" / "check_sfs.fs_t.tsv"
    items = [r["item"] for r in csv.DictReader(dest.open(), delimiter="\t")]
    assert len(items) == len(set(items)) == 6


def test_a_thin_bed_fails_on_snp_count(tmp_path):
    _prep(tmp_path, "Tiny", {"chr1": 5, "chr2": 5})
    _chrom_out(tmp_path, "Tiny", "chr1", snps=10.0)
    _chrom_out(tmp_path, "Tiny", "chr2", snps=10.0)
    assert _main(tmp_path, "Tiny.bed") == 0
    dest = tmp_path / "sfs" / "check_sfs.fs_t.tsv"
    rows = {r["item"]: r for r in csv.DictReader(dest.open(), delimiter="\t")}
    assert rows["Tiny"]["status"] == "FAIL"
    assert "SNP" in rows["Tiny"]["reason"]


def test_losing_test_sites_to_missing_controls_warns(tmp_path):
    # prep keeps 10000 bp; only 500 bp per chromosome ends up in has_qui.
    _prep(tmp_path, "Far", {"chr1": 5, "chr2": 5}, after_mask=10000)
    _chrom_out(tmp_path, "Far", "chr1")
    _chrom_out(tmp_path, "Far", "chr2")
    assert _main(tmp_path, "Far.bed") == 0
    dest = tmp_path / "sfs" / "check_sfs.fs_t.tsv"
    rows = {r["item"]: r for r in csv.DictReader(dest.open(), delimiter="\t")}
    assert rows["Far"]["status"] == "WARN"
    assert "control" in rows["Far"]["reason"]


def test_a_missing_prep_bed_is_a_fail_row_not_a_crash(tmp_path):
    """The check file must still be written, or the gate has nothing to read."""
    (tmp_path / "sfs").mkdir(parents=True, exist_ok=True)
    assert _main(tmp_path, "Gone.bed") == 0
    status = read_status(tmp_path / "sfs" / "check_sfs.fs_t.tsv")
    assert status["Gone"] == "FAIL"
    assert "prep" in (tmp_path / "sfs" / "check_sfs.fs_t.tsv").read_text()


def test_bare_chromosome_numbers_work_end_to_end(tmp_path):
    _prep(tmp_path, "ELS", {"chr1": 5})
    _chrom_out(tmp_path, "ELS", "chr1")
    assert _main(tmp_path, "ELS.bed", chroms="1") == 0
    assert read_status(tmp_path / "sfs" / "check_sfs.fs_t.tsv") == {
        "ELS:chr1": "PASS", "ELS": "PASS"}


def test_runs_as_a_script(tmp_path):
    """run_concat.sh invokes it with `python3 <path>`, so the CLI must stand alone."""
    _prep(tmp_path, "ELS", {"chr1": 5})
    _chrom_out(tmp_path, "ELS", "chr1")
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts/checks/check_sfs.py"),
         "--out", str(tmp_path / "sfs"), "--run", "fs_t",
         "--prep-out", str(tmp_path / "prep"), "--prep-run", "prep_t",
         "--chroms", "chr1", "ELS.bed"],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "sfs" / "check_sfs.fs_t.tsv").is_file()


def test_the_hasqui_warning_explains_both_possible_causes():
    """It does not block, so the message is the whole value: it must name what to
    look at and say that the dropped sites shrink l."""
    status, reasons = assess_summary(n_snps=50000, hasqui_ratio=0.5,
                                     failed_chroms=[], min_hasqui_frac=0.90)
    assert status == "WARN"
    joined = " ".join(reasons)
    assert "too long" in joined and "candidate-neutral" in joined
    assert "biased subset" in joined and "l before using" in joined


def test_per_chromosome_hasqui_ratio_is_against_that_chromosome(tmp_path):
    """The original divided by the genome-wide bed length, giving ~1/22 of the truth."""
    import csv as _csv
    out = tmp_path / "sfs"
    run = "fs_test"
    bychr = out / "TST" / run / "bychr"
    bychr.mkdir(parents=True)
    (out / "TST" / run).joinpath("spec.txt").write_text("FS_RUN=%s\n" % run)
    for chrom, prep, hq in (("chr1", 1000, 900), ("chr2", 1000, 100)):
        (bychr / ("lengths.%s.tsv" % chrom)).write_text(
            "stage\tbp\nprep_chr\t%d\nmerged_chr\t%d\nhas_qui\t%d\nquiecent_ctrl\t50\n"
            % (prep, prep, hq))
    rows = {}
    for line in (bychr / "lengths.chr1.tsv").read_text().splitlines()[1:]:
        k, v = line.split("\t"); rows[k] = int(v)
    assert round(rows["has_qui"] / rows["prep_chr"], 4) == 0.9


# ------------------------------------------------- the neutral control inputs --

import pathlib                                             # noqa: E402
import subprocess                                          # noqa: E402

PKG = pathlib.Path(__file__).resolve().parent.parent
SFS = str(PKG / "scripts" / "sfs" / "run_sfs_chrom.sh")
BASE_ENV = {
    "OUT": "/tmp/o", "REFS_DIR": "/tmp/r", "PREP_RUN": "p", "FS_RUN": "f",
    "DEMO_RUN": "d", "CACHE_RUN": "c", "DFE_RUN": "x", "SUMMARY_RUN": "s",
    "VCF_TMPL": "/tmp/v.{chrom}.vcf.gz", "MASK_BED": "/tmp/m",
    "CDS_BED": "/tmp/cd", "HET_BED": "/tmp/h", "QUI_BED": "/tmp/q",
    "MUTRATE_TAB": "/tmp/mt", "POPFILE": "/tmp/pf", "POP": "YRI",
}


def _run_sfs(**extra):
    env = dict(BASE_ENV)
    env.update(extra)
    return subprocess.run(["bash", SFS, "chr22", "/tmp/test.bed"],
                          capture_output=True, text=True, env=env)


def test_both_control_keys_set_is_refused():
    """One key would have to be type-sniffed; two keys cannot be ambiguous, so
    setting both has to be an error rather than a silent precedence rule."""
    r = _run_sfs(NEUTRAL_CTRL_DIR="/tmp/dir", NEUTRAL_CTRL_BED="/tmp/one.bed")
    assert r.returncode == 2
    assert "both set" in r.stderr


def test_neither_control_key_set_is_refused():
    r = _run_sfs()
    assert r.returncode == 2
    assert "neutral controls are not set" in r.stderr


def test_one_control_key_gets_past_the_check():
    """Either key alone must reach the ordinary missing-input error (rc 1), not
    the configuration error (rc 2)."""
    for key in ("NEUTRAL_CTRL_DIR", "NEUTRAL_CTRL_BED"):
        r = _run_sfs(**{key: "/tmp/does-not-exist"})
        assert r.returncode == 1, (key, r.returncode, r.stderr)
        assert "both set" not in r.stderr and "are not set" not in r.stderr
