import pathlib
import subprocess
from scripts.checks.checkfile import write_check, read_status

ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_roundtrip(tmp_path):
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "PASS", "max_abs_resid": 1.2},
                    {"item": "ELS", "status": "FAIL", "max_abs_resid": 19.0}])
    assert read_status(f) == {"EnhA": "PASS", "ELS": "FAIL"}
    header = f.read_text().splitlines()[0].split("\t")
    assert header[:2] == ["item", "status"]


def test_missing_columns_are_blank_not_dropped(tmp_path):
    f = tmp_path / "c.tsv"
    write_check(f, [{"item": "a", "status": "PASS", "x": 1},
                    {"item": "b", "status": "PASS", "y": 2}])
    lines = f.read_text().splitlines()
    assert lines[0].split("\t") == ["item", "status", "x", "y"]
    assert lines[1].split("\t")[3] == "NA"


GATE = str(ROOT / "scripts/checks/gate.py")


def run_gate(*args):
    return subprocess.run(["python3", GATE, *args], capture_output=True, text=True)


def test_gate_warns_but_exits_zero(tmp_path):
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "FAIL", "reason": "max|resid| 19"}])
    r = run_gate(str(f), "demo")
    assert r.returncode == 0
    assert "WARNING" in r.stderr and "EnhA" in r.stderr
    assert "demo" in r.stderr


def test_gate_reports_warn_rows(tmp_path):
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "PASS"},
                    {"item": "ELS", "status": "WARN", "reason": "few sites"}])
    r = run_gate(str(f), "demo")
    assert r.returncode == 0
    assert "ELS" in r.stderr
    assert "EnhA" not in r.stderr


def test_gate_reports_both_fail_and_warn(tmp_path):
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "FAIL", "reason": "max|resid| 19"},
                    {"item": "ELS", "status": "WARN", "reason": "few sites"},
                    {"item": "PLS", "status": "PASS"}])
    r = run_gate(str(f), "demo")
    assert r.returncode == 0
    assert "EnhA" in r.stderr
    assert "ELS" in r.stderr
    assert "PLS" not in r.stderr
    assert "reference only" in r.stderr


def test_gate_silent_when_all_pass(tmp_path):
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "PASS"},
                    {"item": "ELS", "status": "PASS"}])
    r = run_gate(str(f), "demo")
    assert r.returncode == 0
    assert r.stderr == ""
    assert r.stdout == ""


def test_gate_on_missing_file_exits_zero(tmp_path):
    r = run_gate(str(tmp_path / "absent.tsv"), "demo")
    assert r.returncode == 0
    assert "no check file" in r.stderr


def test_gate_on_empty_file_exits_zero(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text("")
    r = run_gate(str(f), "demo")
    assert r.returncode == 0
    assert "no check file" in r.stderr


def test_gate_survives_tab_and_newline_in_reason(tmp_path):
    # checkfile.write_check quotes a reason holding a tab or a newline, so the
    # row still has exactly two leading columns. A line-and-tab parser such as
    # `awk -F'\t' 'NR>1 && $2=="FAIL"'` does not honour that quoting: the
    # embedded newline splits the PLS row into a second physical line whose
    # first two tab-fields read as item=PLS_bin3, status=FAIL, so awk invents a
    # failing item that does not exist.
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "FAIL",
                     "reason": "max|resid|\t19.0 over\tthreshold"},
                    {"item": "PLS", "status": "PASS",
                     "reason": "per-bin resid:\nPLS_bin3\tFAIL\t0.4"}])
    # the library round-trips it correctly; only the naive parser does not
    assert read_status(f) == {"EnhA": "FAIL", "PLS": "PASS"}
    r = run_gate(str(f), "demo")
    assert r.returncode == 0
    assert "EnhA" in r.stderr
    assert "PLS_bin3" not in r.stderr


# --- write-time validation: a bad row must fail here, not silently later -------

def test_write_rejects_missing_status(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="status"):
        write_check(tmp_path / "c.tsv", [{"item": "EnhA"}])


def test_write_rejects_unknown_status(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="expected one of"):
        write_check(tmp_path / "c.tsv", [{"item": "EnhA", "status": "ok"}])


def test_write_rejects_duplicate_items(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="unique"):
        write_check(tmp_path / "c.tsv", [{"item": "EnhA", "status": "FAIL"},
                                         {"item": "EnhA", "status": "PASS"}])


def test_skip_is_a_valid_status(tmp_path):
    f = tmp_path / "c.tsv"
    write_check(f, [{"item": "EnhA:chr4", "status": "SKIP", "reason": "no intervals in bed"}])
    assert read_status(f) == {"EnhA:chr4": "SKIP"}


def test_duplicate_items_collapse_to_the_worst_status(tmp_path):
    """A hand-edited file must not hide a FAIL behind a later PASS."""
    f = tmp_path / "c.tsv"
    f.write_text("item\tstatus\nEnhA\tFAIL\nEnhA\tPASS\n")
    assert read_status(f) == {"EnhA": "FAIL"}


def test_gate_reports_unrecognised_status(tmp_path):
    import subprocess
    f = tmp_path / "c.tsv"
    f.write_text("item\tstatus\nEnhA\tNA\n")
    r = subprocess.run(["python3", str(ROOT / "scripts/checks/gate.py"), str(f), "demo"],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert "does not understand" in r.stderr and "EnhA" in r.stderr
