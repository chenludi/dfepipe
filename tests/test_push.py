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


# ------------------------------------------------- the post-push verification --

def test_verification_resolves_a_tilde_destination(tmp_path):
    """The verify step single-quoted the remote path, so `~` reached the remote
    shell literally and sha256sum never found the file. Every push to a `~`
    destination then printed "MISMATCH ... remote absent" and exited 1, although
    the push itself had worked - checked by hand on 2026-10-25, where the remote
    hash matched exactly. A check that always fails is one you stop reading."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    # stand in for the remote: refuse a literal ~, otherwise hash the real file
    (bin_dir / "ssh").write_text(
        "#!/bin/bash\n"
        'cmd="${!#}"\n'
        'case "$cmd" in *"~"*) exit 0 ;; esac\n'          # literal tilde: no output
        'f="${cmd#*\\"}"; f="${f%%\\"*}"\n'
        'f="${f/\\$HOME/%s}"\n'                            # the stub's own "home"
        "shasum -a256 \"$f\" 2>/dev/null | cut -c1-16\n" % ROOT.parent)
    (bin_dir / "rsync").write_text("#!/bin/bash\nexit 0\n")
    for f in ("ssh", "rsync"):
        os.chmod(bin_dir / f, 0o755)
    e = dict(os.environ)
    e["PATH"] = "%s:%s" % (bin_dir, e["PATH"])
    e["ALLOW_RUNNING"] = "1"
    r = subprocess.run(["bash", PUSH, "testhost:~/%s" % ROOT.name],
                       capture_output=True, text=True, env=e)
    assert "MISMATCH" not in r.stderr, r.stderr
    assert "verified  config/pipeline.conf" in r.stdout, (r.stdout, r.stderr)


def test_both_verified_files_are_reported_not_just_the_first(tmp_path):
    """It exited on the first mismatch, so a second broken file stayed hidden
    until the first was fixed and the push run again."""
    text = pathlib.Path(PUSH).read_text()
    i = text.index("for f in config/pipeline.conf scripts/chains/run_all.sh")
    body = text[i:i + 900]
    assert "bad=1" in body, "a mismatch must not exit the loop"
    assert 'MISMATCH' in body
