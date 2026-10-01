"""The chain: ordering, FROM, and that a gate actually stops a job.

The ordering is tested through DRY_RUN=1, which prints the qsub commands and
submits nothing. The gate is tested for real, by running worker_step.sh against
a check file holding a FAIL and asserting the step never ran.
"""
import os
import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHAIN = str(ROOT / "scripts" / "chains" / "run_all.sh")
WORKER_STEP = str(ROOT / "scripts" / "sge" / "worker_step.sh")

BASE = {
    "PREP_RUN": "p1", "FS_RUN": "f1", "DEMO_RUN": "d1", "CACHE_RUN": "c1",
    "DFE_RUN": "x1", "SUMMARY_RUN": "s1", "DRY_RUN": "1",
}


def dry(tmp_path, beds=("a.bed",), **extra):
    for b in beds:
        (tmp_path / b).write_text("chr1\t0\t100\n")
    env = dict(os.environ)
    env.update(BASE)
    env.update({"OUT": str(tmp_path / "out"), "REFS_DIR": str(tmp_path)})
    env.update(extra)
    r = subprocess.run(["bash", CHAIN] + [str(tmp_path / b) for b in beds],
                       capture_output=True, text=True, env=env)
    return r


def order(r):
    """The stage names in submission order, from the dry-run lines."""
    out = []
    for line in (r.stdout + r.stderr).splitlines():
        m = re.match(r"\s+(\S+)\s+(TERSE=1|qsub )", line)
        if m:
            out.append(m.group(1))
    return out


# ------------------------------------------------------------------- ordering --

def test_phase_1_stops_before_cache_and_launches_the_next(tmp_path):
    """submit_cache.sh reads best_inference.txt at submit time, so cache cannot
    go in before bestdemo has run. The phase ends with a launcher instead."""
    r = dry(tmp_path)
    assert r.returncode == 0, r.stderr
    assert order(r) == ["prep", "sfs", "concat", "demo", "bestdemo",
                        "launch:cache"]
    assert "phase 1 submitted" in r.stdout


def test_phase_2_is_cache_then_a_launcher(tmp_path):
    """submit_dfe.sh generates the packs, which needs the spectrum and Na, and
    the array size depends on how many packs result."""
    r = dry(tmp_path, FROM="cache")
    assert order(r) == ["cache", "launch:dfe"]
    assert "phase 2 submitted" in r.stdout


def test_phase_3_runs_to_the_end(tmp_path):
    r = dry(tmp_path, FROM="dfe")
    assert order(r) == ["fitneu", "dfe", "bestdfe_gamma", "dfe",
                        "bestdfe_neugamma", "summary", "bestpdf"]
    assert "launch" not in r.stdout + r.stderr


def test_the_launcher_passes_the_next_FROM_and_the_run_keys(tmp_path):
    """It re-runs this script, so it must carry enough to resume the same run."""
    r = dry(tmp_path)
    line = [l for l in (r.stdout + r.stderr).splitlines()
            if "launch:cache" in l and "qsub" in l][0]
    assert "FROM=cache" in line and "run_all.sh" in line


def test_nothing_is_submitted_in_dry_run(tmp_path):
    r = dry(tmp_path)
    assert "nothing will be submitted" in r.stdout
    assert "qsub -terse" not in r.stdout


def test_each_stage_holds_on_the_one_before(tmp_path):
    r = dry(tmp_path)
    text = r.stdout + r.stderr
    for stage, held in [("sfs", "DRY-prep"), ("concat", "DRY-sfs"),
                        ("demo", "DRY-concat"), ("bestdemo", "DRY-demo"),
                        ("launch:cache", "DRY-bestdemo")]:
        line = [l for l in text.splitlines()
                if re.match(r"\s+%s\s+(TERSE|qsub)" % re.escape(stage), l)]
        assert line, stage
        assert held in line[0], (stage, line[0])


def test_summary_waits_for_every_pdf_and_the_null(tmp_path):
    """It reads all of them, so it must not start when one is still running."""
    r = dry(tmp_path, FROM="dfe")
    line = [l for l in (r.stdout + r.stderr).splitlines()
            if re.match(r"\s+summary\s+qsub", l)][0]
    for held in ("DRY-fitneu", "DRY-bestdfe_gamma", "DRY-bestdfe_neugamma"):
        assert held in line, (held, line)


def test_every_stage_gates_on_its_upstream(tmp_path):
    """Including the arrays: a FAILed prep must not cost 22 sfs tasks."""
    text = dry(tmp_path).stdout + dry(tmp_path).stderr
    text += dry(tmp_path, FROM="cache").stderr + dry(tmp_path, FROM="dfe").stderr
    for stage, gate in [("sfs", "GATE=prep"), ("demo", "GATE=sfs"),
                        ("cache", "GATE=demo"), ("dfe", "GATE=cache")]:
        line = [l for l in text.splitlines()
                if re.match(r"\s+%s\s+TERSE" % re.escape(stage), l)]
        assert line, stage
        assert gate in line[0], (stage, line[0])


def test_prep_gates_on_nothing(tmp_path):
    r = dry(tmp_path)
    line = [l for l in (r.stdout + r.stderr).splitlines()
            if re.match(r"\s+prep\s+TERSE", l)][0]
    assert "GATE=none" in line and "HOLD=none" in line


# ----------------------------------------------------------------------- FROM --

@pytest.mark.parametrize("start,expected_first", [
    ("prep", "prep"), ("sfs", "sfs"), ("demo", "demo"),
    ("cache", "cache"), ("dfe", "fitneu"), ("summary", "summary"),
])  # noqa: E501
def test_FROM_starts_at_that_stage(tmp_path, start, expected_first):
    r = dry(tmp_path, FROM=start)
    assert r.returncode == 0, r.stderr
    got = order(r)
    assert got and got[0] == expected_first, (start, got)


def test_FROM_summary_skips_straight_to_the_table(tmp_path):
    r = dry(tmp_path, FROM="summary")
    assert order(r) == ["summary", "bestpdf"]


def test_an_unknown_FROM_is_refused(tmp_path):
    r = dry(tmp_path, FROM="nonsense")
    assert r.returncode == 2
    assert "is not one of" in r.stderr


def test_a_missing_bed_is_refused_before_anything_is_submitted(tmp_path):
    env = dict(os.environ)
    env.update(BASE)
    env.update({"OUT": str(tmp_path / "out"), "REFS_DIR": str(tmp_path)})
    r = subprocess.run(["bash", CHAIN, str(tmp_path / "absent.bed")],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 1 and "no such bed" in r.stderr


def test_no_beds_is_a_usage_error(tmp_path):
    env = dict(os.environ)
    env.update(BASE)
    env.update({"OUT": str(tmp_path / "out"), "REFS_DIR": str(tmp_path)})
    r = subprocess.run(["bash", CHAIN], capture_output=True, text=True, env=env)
    assert r.returncode == 2 and "usage" in r.stderr


# ------------------------------------------------- the gate, for real this time --

def run_worker_step(tmp_path, gate_files, script_args, force=None):
    env = dict(os.environ)
    env.update({"REPO_DIR": str(ROOT), "GATE_FILES": gate_files,
                "GATE_STAGE": "demo"})
    if force:
        env["FORCE"] = force
    return subprocess.run(["bash", WORKER_STEP] + script_args,
                          capture_output=True, text=True, env=env)


@pytest.fixture
def marker_script(tmp_path):
    """A stand-in step: it writes a file, so its absence proves it never ran."""
    path = tmp_path / "step.sh"
    path.write_text("#!/bin/bash\necho ran > %s/ran.txt\n" % tmp_path)
    return path


def test_a_FAIL_upstream_stops_the_step(tmp_path, marker_script):
    chk = tmp_path / "check_demo.tsv"
    chk.write_text("item\tstatus\treason\nEnhA\tFAIL\ttoo few SNPs\n")
    r = run_worker_step(tmp_path, str(chk), [str(marker_script)])
    assert r.returncode == 1
    assert not (tmp_path / "ran.txt").exists(), "the step ran despite a FAIL"
    assert "gate stopped this task" in r.stderr


def test_FORCE_lets_it_through(tmp_path, marker_script):
    chk = tmp_path / "check_demo.tsv"
    chk.write_text("item\tstatus\treason\nEnhA\tFAIL\ttoo few SNPs\n")
    r = run_worker_step(tmp_path, str(chk), [str(marker_script)], force="1")
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "ran.txt").exists()


def test_a_WARN_does_not_stop_the_step(tmp_path, marker_script):
    chk = tmp_path / "check_demo.tsv"
    chk.write_text("item\tstatus\treason\nEnhA\tWARN\tunusual residual\n")
    r = run_worker_step(tmp_path, str(chk), [str(marker_script)])
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "ran.txt").exists()


def test_a_missing_check_file_stops_the_step(tmp_path, marker_script):
    """It means the upstream check crashed, which is not the same as fine."""
    r = run_worker_step(tmp_path, str(tmp_path / "absent.tsv"),
                        [str(marker_script)])
    assert r.returncode == 1
    assert not (tmp_path / "ran.txt").exists()


def test_no_gate_asked_for_runs_the_step(tmp_path, marker_script):
    """A stage run by hand sets no GATE_FILES, and must not be blocked."""
    env = dict(os.environ)
    env.update({"REPO_DIR": str(ROOT)})
    r = subprocess.run(["bash", WORKER_STEP, str(marker_script)],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "ran.txt").exists()


def test_several_gate_files_are_all_consulted(tmp_path, marker_script):
    """cache writes one check file per bed, so the chain passes a list."""
    ok = tmp_path / "ok.tsv"
    ok.write_text("item\tstatus\treason\nA\tPASS\tNA\n")
    bad = tmp_path / "bad.tsv"
    bad.write_text("item\tstatus\treason\nB\tFAIL\tcache too small\n")
    r = run_worker_step(tmp_path, "%s %s" % (ok, bad), [str(marker_script)])
    assert r.returncode == 1
    assert not (tmp_path / "ran.txt").exists()
