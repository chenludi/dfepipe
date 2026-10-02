"""The config file is the only place a cluster path may appear."""
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent


# The one allowed exception: a provenance record of where result files came
# from on the cluster is only meaningful with the real paths in it.
PATH_EXCEPTIONS = {
    "setup/patches/optional/reference_integration/PROVENANCE.log",
}

# Assembled at runtime: this file has to name the forbidden prefixes to look for
# them, and must not be exempted from its own check.
CLUSTER_PREFIXES = ("/u/" + "home/", "/u/" + "project/")


def test_no_hardcoded_cluster_paths_outside_config():
    offenders = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix not in {".sh", ".py", ".log"}:
            continue
        rel = str(path.relative_to(ROOT))
        if rel.startswith("config/") or rel in PATH_EXCEPTIONS:
            continue
        text = path.read_text(errors="ignore")
        if any(p in text for p in CLUSTER_PREFIXES):
            offenders.append(rel)
    assert offenders == [], f"cluster paths outside config: {offenders}"


def test_conf_requires_out():
    """Sourcing without OUT set must fail loudly, not default silently."""
    r = subprocess.run(
        ["bash", "-c", f"source {ROOT}/config/pipeline.conf"],
        capture_output=True, text=True,
    )
    assert r.returncode != 0
    assert "OUT not set" in r.stderr


def test_conf_rejects_unknown_demo_model():
    env = "OUT=/tmp REFS_DIR=/tmp VCF=/tmp/a MASK_BED=/tmp/b CDS_BED=/tmp/c HET_BED=/tmp/d "
    env += "MUTRATE_TAB=/tmp/e QUI_BED=/tmp/f "
    env += "PREP_RUN=r FS_RUN=r DEMO_RUN=r CACHE_RUN=r DFE_RUN=r SUMMARY_RUN=r "
    r = subprocess.run(
        ["bash", "-c", f"{env} DEMO_MODEL=nonsense source {ROOT}/config/pipeline.conf"],
        capture_output=True, text=True,
    )
    assert r.returncode != 0
    assert "not supported" in r.stderr


def _config_keys():
    import re
    text = (ROOT / "config" / "pipeline.conf").read_text()
    # some lines hold several assignments: `A="..."; B="..."`
    pat = re.compile(r'(?:^|;)\s*([A-Za-z_][A-Za-z0-9_]*)=', re.M)
    return set(pat.findall(text))


# Paths a stage derives from REPO_DIR or from another config key, not settings.
DERIVED = {"REPO_DIR", "LEGACY_SCRIPT_DIR", "MUTRATE_TMPL"}


def test_no_hidden_defaults_in_stage_scripts():
    """A stage may not invent a default for a tunable; pipeline.conf owns them.

    `NAME="${NAME:-value}"` at the top of a stage is a second copy of a number
    that nobody reading the config would see. Stages assert instead.
    """
    import re
    pat = re.compile(r'^([A-Z_][A-Z0-9_]*)="\$\{\1:-', re.M)
    offenders = []
    for path in sorted(ROOT.glob("scripts/*/*.sh")):
        for name in pat.findall(path.read_text()):
            if name not in DERIVED:
                offenders.append(f"{path.relative_to(ROOT)}:{name}")
    assert offenders == [], (
        "hidden defaults; move the value to config/pipeline.conf and assert "
        f'with : "${{NAME:?...}}": {offenders}')


def test_every_config_key_a_check_reads_exists():
    """conf.get_*("KEY", fallback) must name a real pipeline.conf key."""
    import re
    keys = _config_keys()
    pat = re.compile(r'conf\.get(?:_float|_int)?\("([A-Z_][A-Z0-9_]*)"')
    missing = []
    for path in sorted(ROOT.glob("scripts/checks/*.py")):
        for name in pat.findall(path.read_text()):
            if name not in keys:
                missing.append(f"{path.relative_to(ROOT)}:{name}")
    assert missing == [], f"not in pipeline.conf: {missing}"


def test_conf_loader_reads_the_real_config():
    import sys
    sys.path.insert(0, str(ROOT / "scripts" / "checks"))
    from conf import Conf
    c = Conf(str(ROOT / "config" / "pipeline.conf"))
    assert c.get_int("MIN_SNPS", -1) == 500
    assert c.get_float("HASQUI_MIN_FRAC", -1.0) == 0.75
    assert c.get("NEUGAMMA_UPPER", "") == "2 1e8 1"
    # a missing key falls back rather than raising
    assert c.get_int("NO_SUCH_KEY", 7) == 7


# Set by the caller or derived by the stage, not by pipeline.conf.
# Set by the scheduler, or by the submit script for its own worker.
RUNTIME_SET = {"REPO_DIR", "DFE_PACK", "FORCE", "RESUME", "SGE_TASK_ID",
               "CACHE_BEDS", "DEMO_P0_FILES", "DFE_PDF", "DFE_PACKS",
               "DFE_TASK_BEDS"}
# config/env.sh, which every stage sources after pipeline.conf.
ENV_SH_SET = {"BEDTOOLS_MOD", "R_MOD"}


def test_every_asserted_name_is_a_config_key():
    """`: "${X:?...}"` in a stage must name something pipeline.conf defines."""
    import re
    keys = _config_keys()
    pat = re.compile(r'^:\s*"\$\{([A-Z_][A-Z0-9_]*):\?', re.M)
    missing = []
    for path in sorted(ROOT.glob("scripts/*/*.sh")):
        for name in pat.findall(path.read_text()):
            if name not in keys | RUNTIME_SET | ENV_SH_SET:
                missing.append(f"{path.relative_to(ROOT)}:{name}")
    assert missing == [], f"asserted but never set by pipeline.conf: {missing}"


def test_scheduling_inputs_default_to_no_change():
    """Both preferences empty means PACKNUMBER stands, so a published run is
    reproduced byte for byte."""
    conf = (ROOT / "config" / "pipeline.conf").read_text()
    assert 'DFE_MAX_TASKS="${DFE_MAX_TASKS:-}"' in conf
    assert 'DFE_TARGET_TASK_HOURS="${DFE_TARGET_TASK_HOURS:-}"' in conf
    assert 'PACKNUMBER="${PACKNUMBER:-50}"' in conf


def test_submit_dfe_derives_packnumber_from_the_preferences():
    sub = (ROOT / "scripts" / "sge" / "submit_dfe.sh").read_text()
    # a target time sets PACKNUMBER before the packs are generated ...
    assert "DFE_TARGET_TASK_HOURS * 3600 / DFE_SECONDS_PER_FIT" in sub
    # ... and a task cap is checked once the array size is known, and refuses
    assert "DFE_MAX_TASKS" in sub and "Nothing was submitted" in sub


def test_conf_sources_cleanly_under_set_u():
    """Every derived path must resolve with `set -u`: a dropped helper variable
    (_conf_dir) only shows up here, not in a plain source."""
    env = ("OUT=/tmp REFS_DIR=/tmp PREP_RUN=r FS_RUN=r DEMO_RUN=r CACHE_RUN=r "
           "DFE_RUN=r SUMMARY_RUN=r ")
    script = (
        f"set -euo pipefail; {env} source {ROOT}/config/pipeline.conf; "
        'printf "%s\\n%s\\n" "$CACHE_MARKER_TABLE" "$MODELSEL_TABLE"')
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for line in r.stdout.split():
        assert line.startswith("/") and "/config/" in line, line


def test_no_stage_hardcodes_a_config_value():
    """ROULETTE_SCALE was a config key honoured by two stages and hardcoded in
    three, so overriding it would have given an Na from one constant and DFE
    scaling from another, silently. Any config value a stage needs must be
    passed in, never written as a literal beside the variable's own name."""
    import re
    offenders = []
    for path in sorted(ROOT.glob("scripts/*/*.sh")):
        text = path.read_text()
        for name in ("ROULETTE_SCALE",):
            for m in re.finditer(r"^%s\s*=\s*([0-9].*)$" % name, text, re.M):
                offenders.append("%s: %s = %s" % (path.relative_to(ROOT), name,
                                                  m.group(1)))
    assert offenders == [], offenders


def test_every_roulette_consumer_asserts_the_key():
    """A stage that APPLIES the scale must demand the value, not invent one.

    Keyed on `ROULETTE_SCALE *`, the multiplication, not on the column name:
    run_sfs_chrom.sh writes the raw averaged rate into the fsmut table and never
    scales it, so it is not a consumer.
    """
    missing = []
    for path in sorted(ROOT.glob("scripts/*/*.sh")):
        text = path.read_text()
        if "ROULETTE_SCALE *" in text and "ROULETTE_SCALE:?" not in text:
            missing.append(str(path.relative_to(ROOT)))
    assert missing == [], missing
