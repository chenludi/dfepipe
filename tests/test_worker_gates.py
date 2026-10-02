"""Every worker must run its upstream gate unconditionally, before doing work.

`worker_dfe.sh` had the gate nested inside its `if [[ "$DFE_PDF" == "neu" ]]`
branch, and `run_all.sh:293` skips `neu`, so the branch was never taken from the
chain. The dfe array was therefore the only stage with no gate at all: on
2026-10-01 a missing `check_cache.tsv` stopped `fitneu` and both `bestdfe` jobs
while 347 fits ran happily against a cache nothing had verified.

`run_all.sh` and `submit_dfe.sh` both passed GATE_FILES correctly; the worker
ignored it. So the bug was invisible from the submitting side, which is why these
tests look at the worker scripts themselves.
"""
import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
SGE = ROOT / "scripts" / "sge"

WORKERS = sorted(p.name for p in SGE.glob("worker_*.sh"))


def test_the_workers_were_found():
    assert len(WORKERS) >= 6, WORKERS


@pytest.mark.parametrize("name", WORKERS)
def test_worker_calls_the_gate(name):
    text = (SGE / name).read_text()
    assert "dfepipe_gate_or_exit" in text, (
        "%s never calls the gate, so nothing stops it when its upstream failed" % name)


@pytest.mark.parametrize("name", WORKERS)
def test_the_gate_call_is_not_nested_in_a_conditional(name):
    """The defect was indentation, not a missing call: the gate was inside an
    `if`. A gate that runs only sometimes is worse than none, because the
    submitting side still reports it as gated."""
    for lineno, line in enumerate((SGE / name).read_text().splitlines(), 1):
        if "dfepipe_gate_or_exit" in line and not line.lstrip().startswith("#"):
            assert line == line.lstrip(), (
                "%s:%d the gate call is indented, so it is inside a conditional "
                "or loop and will not always run: %r" % (name, lineno, line))


@pytest.mark.parametrize("name", WORKERS)
def test_the_gate_runs_before_any_stage_script(name):
    """Gating after the work has started would not prevent anything."""
    lines = (SGE / name).read_text().splitlines()
    gate = next((i for i, l in enumerate(lines)
                 if "dfepipe_gate_or_exit" in l and not l.lstrip().startswith("#")),
                None)
    assert gate is not None, name
    # worker_step.sh runs `bash "$script"`, the others name the path inline
    runs = [i for i, l in enumerate(lines)
            if re.search(r'^\s*bash "\$(\{REPO_DIR\}/scripts/|script")', l)]
    assert runs, "%s never runs a stage script" % name
    assert gate < min(runs), (
        "%s calls the gate at line %d but starts work at line %d"
        % (name, gate + 1, min(runs) + 1))


def test_the_gate_blocks_when_a_check_file_is_missing(tmp_path):
    """The behaviour the whole mechanism exists for, exercised directly."""
    script = tmp_path / "probe.sh"
    script.write_text(
        '#!/bin/bash\nset -euo pipefail\n'
        'source "%s/gate_upstream.sh"\n'
        'dfepipe_gate_or_exit "probe"\n'
        'echo WORK_RAN\n' % SGE)
    env = {"PATH": "/usr/bin:/bin", "REPO_DIR": str(ROOT),
           "GATE_FILES": str(tmp_path / "absent.tsv"), "GATE_STAGE": "cache"}
    r = subprocess.run(["bash", str(script)], capture_output=True, text=True, env=env)
    assert "WORK_RAN" not in r.stdout
    assert r.returncode != 0


def test_the_gate_passes_a_clean_check_file(tmp_path):
    chk = tmp_path / "check_cache.tsv"
    chk.write_text("item\tstatus\treason\nbed1\tPASS\tNA\n")
    script = tmp_path / "probe.sh"
    script.write_text(
        '#!/bin/bash\nset -euo pipefail\n'
        'source "%s/gate_upstream.sh"\n'
        'dfepipe_gate_or_exit "probe"\n'
        'echo WORK_RAN\n' % SGE)
    env = {"PATH": "/usr/bin:/bin", "REPO_DIR": str(ROOT),
           "GATE_FILES": str(chk), "GATE_STAGE": "cache"}
    r = subprocess.run(["bash", str(script)], capture_output=True, text=True, env=env)
    assert "WORK_RAN" in r.stdout, r.stderr
    assert r.returncode == 0
