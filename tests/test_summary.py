"""run_summary.sh end to end against a synthetic output tree.

The point of this stage is the join, so the tests build a complete fake run and
then remove one input at a time: the failure that matters is a column shifting,
which only an end-to-end run can show.
"""
import csv
import os
import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = str(ROOT / "scripts" / "summary" / "run_summary.sh")

BED = "phastCons17way_0-5_nomiss"
OUTBASE = "exHET_noq.%s_qui2" % BED          # what short_base() produces
FSBASE = "%s.chrs_1t22" % OUTBASE
BINS = ["0.0 <= |s| < 1e-05", "1e-05 <= |s| < 0.0001", "0.0001 <= |s| < 0.001",
        "0.001 <= |s| < 0.01", "0.01 <= |s| < 0.5", "0.5 <= |s| < inf"]
PROPS = ["0.4491", "0.2088", "0.2518", "0.0899", "0.0001", "0.0003"]

LARGE_LL = {
    "gamma": ("filename\tLog(likelihood)\tshape\tscale\ttheta\tscale_div2Na\tsdadi\tNe\n"
              "153.InferDFE.bestfits\t-610.1731\t0.16899\t65.5901\t54890.8\t0.0017765\t0.00030022\t18460.05\n"),
    "neugamma": ("filename\tLog(likelihood)\tshape\tscale\tpneu\ttheta\tscale_div2Na\tsdadi\tNe\n"
                 "205.InferDFE.bestfits\t-610.1731\t0.16899\t65.5901\t0.0\t54890.8\t0.0017765\t0.00030022\t18460.05\n"),
}

FSMUT_HEADER = ("fsfile\tsites_fsfile\tsites_vcf\tratio_keptsites\tlength_inputbed\t"
                "length_adj_inputbed\tprojection_num\tmutbed\tlength_withmutinfo\t"
                "inputbed\tlength_inputbed.1\tratio_withmutinfo\t"
                "avepersite_diploi_Roulette_mutrate\n")
# The genome-wide row, carrying the .etal marker, then one per-chromosome row so
# the file looks like the real thing.
FSMUT_ROWS = ("genomewide.fs.etal\t372407\t375578\t0.9916\t77600048\t76941215\t160\t"
              "mut.etal\t77598417\tbed.etal\t77600048\t0.99998\t0.19037571\n"
              "chr1.fs.etal\t7592\t7600\t0.99\t7592972\t7516215\t160\tmut.etal\t"
              "7592000\tbed.etal\t7592972\t0.99\t0.18732675\n")

ROULETTE = 1.015e-7


def build_tree(tmp_path, pdfs=("neu", "gamma", "neugamma"), omit=()):
    """A complete fake OUT tree. `omit` drops named inputs."""
    out = tmp_path / "out"
    run = "cache_r-1-%s-dfe_r"
    for pdf in pdfs:
        d = out / "dfe" / BED / (run % pdf)
        d.mkdir(parents=True)
        if pdf == "neu":
            if "neu_log" not in omit:
                (d / ("has_qui.%s.dadi_cal_ll_neu_0723.log" % FSBASE)).write_text(
                    "mutratio: 0.109\nll: -28858.7066\n")
            continue
        if "large_ll_%s" % pdf not in omit:
            (d / "large_ll.txt").write_text(LARGE_LL[pdf])
        if "csv_%s" % pdf not in omit:
            with open(d / ("dfe_%s_%s.csv" % (FSBASE, pdf)), "w") as fh:
                w = csv.writer(fh)
                w.writerow([""] + BINS)
                w.writerow(["dfe_%s_%s" % (FSBASE, pdf)] + PROPS)

    sfs = out / "sfs" / BED / "fs_r"
    sfs.mkdir(parents=True)
    if "fsmut" not in omit:
        (sfs / ("fsmut.has_qui.%s.tab" % FSBASE)).write_text(FSMUT_HEADER + FSMUT_ROWS)
    return out


def run_summary(tmp_path, out, pdfs="neu gamma neugamma", summary_run="sum_r"):
    env = dict(os.environ)
    env.update({
        "OUT": str(out), "REFS_DIR": str(tmp_path), "PREP_RUN": "prep_r",
        "FS_RUN": "fs_r", "DEMO_RUN": "demo_r", "CACHE_RUN": "cache_r",
        "DFE_RUN": "dfe_r", "SUMMARY_RUN": summary_run, "PDFS": pdfs,
    })
    return subprocess.run(["bash", SCRIPT, "%s.bed" % BED],
                          capture_output=True, text=True, env=env)


def read_table(out, summary_run="sum_r"):
    path = pathlib.Path(out) / "summary" / summary_run / "sumdfe.tsv"
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        rows = [r for r in rd if r]
    return header, rows


# --------------------------------------------------------------- the happy path --

def test_complete_run_gives_one_row_with_every_column(tmp_path):
    out = build_tree(tmp_path)
    r = run_summary(tmp_path, out)
    assert r.returncode == 0, r.stderr
    header, rows = read_table(out)
    assert len(rows) == 1
    assert len(rows[0]) == len(header)
    row = dict(zip(header, rows[0]))
    assert row["outputbase"] == OUTBASE
    assert float(row["ll_neu"]) == pytest.approx(-28858.7066)
    assert float(row["ll_gamma"]) == pytest.approx(-610.1731)
    assert float(row["ll_neugamma"]) == pytest.approx(-610.1731)
    assert row["gamma_filename"] == "153.InferDFE.bestfits"
    assert row["neugamma_filename"] == "205.InferDFE.bestfits"


def test_u_l_and_ul_come_from_the_genome_wide_row(tmp_path):
    """Not from chr1, which is the row below it."""
    out = build_tree(tmp_path)
    run_summary(tmp_path, out)
    header, rows = read_table(out)
    row = dict(zip(header, rows[0]))
    assert float(row["u"]) == pytest.approx(ROULETTE * 0.19037571 / 2)
    assert float(row["l"]) == pytest.approx(76941215)
    assert float(row["l_inputbed"]) == pytest.approx(77600048)
    assert float(row["ul"]) == pytest.approx(float(row["u"]) * float(row["l"]))


def test_each_pdf_gets_its_own_prefixed_columns(tmp_path):
    """The original repeats bare `shape`/`scale` per pdf, which cannot be read."""
    out = build_tree(tmp_path)
    run_summary(tmp_path, out)
    header, _ = read_table(out)
    assert "gamma_shape" in header and "neugamma_shape" in header
    assert "neugamma_pneu" in header
    assert "shape" not in header            # no bare, ambiguous name
    assert len(header) == len(set(header))  # no duplicates at all


def test_bins_are_prefixed_and_all_six_present(tmp_path):
    out = build_tree(tmp_path)
    run_summary(tmp_path, out)
    header, _ = read_table(out)
    for pdf in ("gamma", "neugamma"):
        got = [h for h in header if h.startswith("%s_bin_" % pdf)]
        assert len(got) == 6, (pdf, got)


# ------------------------------------------------- a missing input must not shift --

@pytest.mark.parametrize("omit", ["neu_log", "large_ll_gamma", "csv_gamma",
                                  "fsmut"])
def test_a_missing_input_becomes_NA_and_does_not_shift_columns(tmp_path, omit):
    """The defect this stage exists to avoid: the original pastes blocks by line
    number, so one missing block slides every later column by one."""
    full = build_tree(tmp_path / "full")
    run_summary(tmp_path / "full", full)
    full_header, _ = read_table(full)

    out = build_tree(tmp_path / "part", omit=(omit,))
    r = run_summary(tmp_path / "part", out)
    assert r.returncode == 0, r.stderr
    header, rows = read_table(out)
    row = dict(zip(header, rows[0]))

    assert len(rows[0]) == len(header)
    # the columns that are still present must hold the same values as before
    if omit != "fsmut":
        assert float(row["u"]) == pytest.approx(ROULETTE * 0.19037571 / 2)
    if omit not in ("large_ll_gamma",):
        assert row["neugamma_filename"] == "205.InferDFE.bestfits"
    if omit == "neu_log":
        assert row["ll_neu"] == "NA"
        assert float(row["ll_gamma"]) == pytest.approx(-610.1731)
    if omit == "large_ll_gamma":
        assert row["ll_gamma"] == "NA" and row["gamma_filename"] == "NA"
        assert float(row["ll_neugamma"]) == pytest.approx(-610.1731)
    if omit == "csv_gamma":
        assert all(row[h] == "NA" for h in header if h.startswith("gamma_bin_"))
        assert float(row["ll_gamma"]) == pytest.approx(-610.1731)
    if omit == "fsmut":
        assert row["u"] == "NA" and row["ul"] == "NA"
        assert float(row["ll_gamma"]) == pytest.approx(-610.1731)


def test_the_header_is_the_same_whatever_is_missing(tmp_path):
    """A fixed column list is what makes a shift impossible."""
    full = build_tree(tmp_path / "a")
    run_summary(tmp_path / "a", full)
    header_full, _ = read_table(full)
    part = build_tree(tmp_path / "b", omit=("neu_log", "csv_gamma"))
    run_summary(tmp_path / "b", part)
    header_part, _ = read_table(part)
    assert header_full == header_part


# ------------------------------------------------------------------ guard rails --

def test_an_existing_table_is_not_overwritten(tmp_path):
    out = build_tree(tmp_path)
    assert run_summary(tmp_path, out).returncode == 0
    again = run_summary(tmp_path, out)
    assert again.returncode == 1
    assert "already exists" in again.stderr


def test_FORCE_replaces_it(tmp_path):
    out = build_tree(tmp_path)
    run_summary(tmp_path, out)
    env = dict(os.environ)
    env.update({"FORCE": "1"})
    r = subprocess.run(["bash", SCRIPT, "%s.bed" % BED], capture_output=True,
                       text=True, env={**env, **{
                           "OUT": str(out), "REFS_DIR": str(tmp_path),
                           "PREP_RUN": "p", "FS_RUN": "fs_r", "DEMO_RUN": "demo_r",
                           "CACHE_RUN": "cache_r", "DFE_RUN": "dfe_r",
                           "SUMMARY_RUN": "sum_r", "PDFS": "neu gamma neugamma",
                           "FORCE": "1"}})
    assert r.returncode == 0, r.stderr


def test_no_beds_is_a_usage_error(tmp_path):
    out = build_tree(tmp_path)
    env = dict(os.environ)
    env.update({"OUT": str(out), "REFS_DIR": str(tmp_path), "PREP_RUN": "p",
                "FS_RUN": "fs_r", "DEMO_RUN": "demo_r", "CACHE_RUN": "cache_r",
                "DFE_RUN": "dfe_r", "SUMMARY_RUN": "sum_r",
                "PDFS": "neu gamma"})
    r = subprocess.run(["bash", SCRIPT], capture_output=True, text=True, env=env)
    assert r.returncode == 2 and "usage" in r.stderr
