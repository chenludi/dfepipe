"""The sfs and demo arrays must cover every (outer, bed) pair exactly once.

Both stages used to be one task per chromosome / per p0 file, with the beds
looped inside. A task's runtime then scaled with the bed count while its h_rt
did not, so on 2026-10-25 a three-bed demo array needed about 5h12 against a
3h10 limit and had to be qaltered before it ran. The fits were independent all
along; only the packing serialised them.

Splitting introduces the one bug a fan-out can have and tests cannot be talked
out of: an index that maps two tasks to the same pair, or misses one. A missed
(chromosome, bed) is a spectrum silently absent from concat; a duplicated one is
two tasks writing the same output. So these walk every task id and assert the
mapping is a bijection onto the product.
"""
import os
import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
SGE = ROOT / "scripts" / "sge"


def fake_repo(tmp_path, stage_rel):
    """A REPO_DIR whose stage script just prints the arguments it was handed."""
    repo = tmp_path / "repo"
    (repo / "scripts" / "sge").mkdir(parents=True)
    (repo / stage_rel).parent.mkdir(parents=True, exist_ok=True)
    (repo / stage_rel).write_text('#!/bin/bash\necho "CALLED $*"\n')
    # the gate is exercised by test_worker_gates.py; here it must not interfere
    (repo / "scripts" / "sge" / "gate_upstream.sh").write_text(
        "dfepipe_gate_or_exit() { :; }\n")
    return repo


def run_task(worker, repo, task_id, beds, **env):
    e = dict(os.environ)
    e.pop("FORCE", None)
    e.update({"REPO_DIR": str(repo), "SGE_TASK_ID": str(task_id)})
    e.update(env)
    r = subprocess.run(["bash", str(SGE / worker)] + beds,
                       capture_output=True, text=True, env=e)
    return r


CHROMS = "chr1 chr2 chr3 chr4"
P0 = "p0_a.txt p0_b.txt p0_c.txt"
BEDS = ["/data/one.bed", "/data/two.bed", "/data/three.bed"]


@pytest.mark.parametrize("worker,stage_rel,outer_env,outer_list", [
    ("worker_sfs.sh", "scripts/sfs/run_sfs_chrom.sh", {"CHROMS": CHROMS}, CHROMS.split()),
    ("worker_demo.sh", "scripts/demo/run_demo.sh", {"DEMO_P0_FILES": P0}, P0.split()),
])
def test_every_pair_is_covered_exactly_once(tmp_path, worker, stage_rel,
                                            outer_env, outer_list):
    repo = fake_repo(tmp_path, stage_rel)
    total = len(outer_list) * len(BEDS)
    seen = []
    for task_id in range(1, total + 1):
        r = run_task(worker, repo, task_id, BEDS, **outer_env)
        assert r.returncode == 0, (task_id, r.stderr)
        called = [l for l in r.stdout.splitlines() if l.startswith("CALLED ")]
        assert len(called) == 1, (task_id, r.stdout)
        args = called[0].split()[1:]
        assert len(args) == 2, "a task must run exactly one bed: %r" % args
        seen.append(tuple(args))
    assert len(seen) == total
    assert len(set(seen)) == total, "a pair is handled twice"
    assert set(seen) == {(o, b) for o in outer_list for b in BEDS}


@pytest.mark.parametrize("worker,stage_rel,outer_env,n_outer", [
    ("worker_sfs.sh", "scripts/sfs/run_sfs_chrom.sh", {"CHROMS": CHROMS}, 4),
    ("worker_demo.sh", "scripts/demo/run_demo.sh", {"DEMO_P0_FILES": P0}, 3),
])
def test_a_task_id_past_the_end_fails_loudly(tmp_path, worker, stage_rel,
                                             outer_env, n_outer):
    """Off the end must not silently wrap onto a pair another task owns."""
    repo = fake_repo(tmp_path, stage_rel)
    r = run_task(worker, repo, n_outer * len(BEDS) + 1, BEDS, **outer_env)
    assert r.returncode != 0
    assert "CALLED" not in r.stdout


@pytest.mark.parametrize("worker,stage_rel,outer_env", [
    ("worker_sfs.sh", "scripts/sfs/run_sfs_chrom.sh", {"CHROMS": CHROMS}),
    ("worker_demo.sh", "scripts/demo/run_demo.sh", {"DEMO_P0_FILES": P0}),
])
def test_no_beds_is_an_error_not_an_empty_run(tmp_path, worker, stage_rel, outer_env):
    repo = fake_repo(tmp_path, stage_rel)
    r = run_task(worker, repo, 1, [], **outer_env)
    assert r.returncode != 0


@pytest.mark.parametrize("submit,outer_expr", [
    ("submit_sfs.sh", "n_chroms * ${#beds[@]}"),
    ("submit_demo.sh", "${#p0files[@]} * ${#beds[@]}"),
])
def test_the_submitter_sizes_the_array_over_the_product(submit, outer_expr):
    """If the submitter and the worker disagree about the array size, the extra
    tasks fail and the missing pairs are never run."""
    text = (SGE / submit).read_text()
    assert outer_expr in text, "%s does not size its array over beds" % submit


def test_prep_is_still_one_job_and_says_why():
    """prep writes its own check file, so an array would race on it. That is a
    deliberate exception and must stay documented, not look like an oversight."""
    text = (SGE / "submit_prep.sh").read_text()
    assert ' -t ' not in text, "prep became an array; check_prep.py races"
    assert "check_prep.py" in text
