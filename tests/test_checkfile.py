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


def test_gate_warns_but_exits_zero(tmp_path):
    f = tmp_path / "check_demo.tsv"
    write_check(f, [{"item": "EnhA", "status": "FAIL", "reason": "max|resid| 19"}])
    r = subprocess.run(["bash", str(ROOT / "scripts/checks/gate.sh"), str(f), "demo"],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert "WARNING" in r.stderr and "EnhA" in r.stderr


def test_gate_on_missing_file_exits_zero(tmp_path):
    r = subprocess.run(["bash", str(ROOT / "scripts/checks/gate.sh"),
                        str(tmp_path / "absent.tsv"), "demo"],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert "no check file" in r.stderr
