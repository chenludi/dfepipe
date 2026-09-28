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
