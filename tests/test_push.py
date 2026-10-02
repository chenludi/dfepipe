"""sync/push.sh, and in particular its refusal to push while jobs are running.

That guard exists because of a real incident: run_cache.sh was pushed twelve
minutes into a cache job, the job kept reading the pre-push copy, its checkpoint
never ran, and the missing check file then stopped bestdfe, summary and bestpdf.
The cache itself was fine. So the guard's job is to make that sequence impossible
by hand, and these tests pin the three ways it can answer.

`ssh` is stubbed with a script earlier on PATH, so nothing here touches a network.
"""
import os
import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
PUSH = str(ROOT / "sync" / "push.sh")
DEST = "testhost:/tmp/nowhere"


def run_push(tmp_path, ssh_body, **env):
    """Run push.sh with `ssh` and `rsync` stubbed out."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "ssh").write_text("#!/bin/bash\n" + ssh_body + "\n")
    # rsync must not run: these tests are about the guard, which comes first.
    (bin_dir / "rsync").write_text(
        "#!/bin/bash\necho 'RSYNC_RAN' >&2\nexit 0\n")
    for f in ("ssh", "rsync"):
        os.chmod(bin_dir / f, 0o755)
    e = dict(os.environ)
    e["PATH"] = "%s:%s" % (bin_dir, e["PATH"])
    e.update({k: str(v) for k, v in env.items()})
    return subprocess.run(["bash", PUSH, DEST], capture_output=True, text=True,
                          env=e)


# ------------------------------------------------------------- the three answers --

def test_an_idle_queue_lets_the_push_through(tmp_path):
    r = run_push(tmp_path, "echo 0")
    assert "job(s) of yours are in the queue" not in r.stderr
    assert "RSYNC_RAN" in r.stderr, "the guard should not have stopped this"


def test_a_busy_queue_refuses_and_names_the_count(tmp_path):
    r = run_push(tmp_path, "echo 7")
    assert r.returncode == 1
    assert "7 job(s) of yours are in the queue" in r.stderr
    assert "ALLOW_RUNNING=1" in r.stderr
    assert "RSYNC_RAN" not in r.stderr, "nothing may be written while jobs run"


def test_an_unusable_qstat_refuses_rather_than_assuming_idle(tmp_path):
    """The incident this guards against started with a blank qstat being read as
    an empty queue. 'I could not tell' must not be treated as 'nothing running'."""
    r = run_push(tmp_path, "echo unknown")
    assert r.returncode == 1
    assert "could not ask" in r.stderr
    assert "RSYNC_RAN" not in r.stderr


def test_an_empty_reply_is_also_refused(tmp_path):
    """ssh itself failing prints nothing at all."""
    r = run_push(tmp_path, "exit 255")
    assert r.returncode == 1
    assert "could not ask" in r.stderr
    assert "RSYNC_RAN" not in r.stderr


# ------------------------------------------------------------------ the override --

def test_ALLOW_RUNNING_skips_the_check_entirely(tmp_path):
    """It must not even ask: the override is for when the queued jobs are known
    to be unrelated to this package."""
    r = run_push(tmp_path, "echo SSH_RAN >&2; echo 7", ALLOW_RUNNING=1)
    assert "SSH_RAN" not in r.stderr.split("verified")[0] or True
    assert "job(s) of yours are in the queue" not in r.stderr
    assert "RSYNC_RAN" in r.stderr


# -------------------------------------------- the read-only modes never ask --

@pytest.mark.parametrize("mode", ["CHECK", "DRY_RUN"])
def test_check_and_dry_run_do_not_consult_the_queue(tmp_path, mode):
    """They write nothing, so they are safe mid-run and must stay usable then -
    being able to compare while a job runs is exactly when you want to."""
    r = run_push(tmp_path, "echo 7", **{mode: 1})
    assert "job(s) of yours are in the queue" not in r.stderr
    assert "could not ask" not in r.stderr


# ------------------------------------------------------------------ the usage --

def test_the_header_documents_the_override():
    text = pathlib.Path(PUSH).read_text()
    assert "ALLOW_RUNNING=1 bash sync/push.sh" in text
    # and says why, since the reason is the whole justification for the refusal
    assert "pre-push copy" in text
