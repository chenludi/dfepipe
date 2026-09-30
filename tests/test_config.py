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
    env += "MUTRATE_TAB=/tmp/e CANDIDATE_NEUTRAL_BED=/tmp/f "
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
    assert c.get_float("HASQUI_MIN_FRAC", -1.0) == 0.90
    assert c.get("NEUGAMMAZERO_UPPER", "") == "3 1e8 1"
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
