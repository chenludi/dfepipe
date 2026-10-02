"""AGENTS.md and docs/PIPELINE.md must stay true as the code changes.

Five blind agents were handed this package and all five reconstructed the right
commands but had to guess at the same handful of things. These docs are the
answer to that, so a claim in them that has quietly gone stale is worse than no
doc at all: it is a confident wrong answer. Everything checkable is checked.
"""
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
DOCS = [ROOT / "AGENTS.md", ROOT / "docs" / "PIPELINE.md"]


def existing_docs():
    return [p for p in DOCS if p.exists()]


def config_keys():
    text = (ROOT / "config" / "pipeline.conf").read_text()
    text += (ROOT / "config" / "env.sh").read_text()
    return set(re.findall(r"(?:^|;)\s*([A-Za-z_]\w*)=", text, re.M))


def config_default(key):
    text = (ROOT / "config" / "pipeline.conf").read_text()
    m = re.search(r'^%s="\$\{%s:-([^}]*)\}"' % (key, key), text, re.M)
    return m.group(1) if m else None


# Names that look like config keys but are something else.
NOT_CONFIG = {
    "DRY_RUN", "FROM", "CHECK", "FORCE", "RESUME", "TERSE", "HOLD", "PNG",
    "DFE_GENERATE_ONLY", "DRY_RUN", "GATE_FILES", "GATE_STAGE", "JOB_ID",
    "SGE_TASK_ID", "DFE_PACK", "DFE_TASK_DEADLINE", "PASS", "WARN", "FAIL",
    "SKIP", "AGENTS", "README", "PIPELINE", "MEMORY", "NA",
}


@pytest.mark.parametrize("doc", existing_docs(), ids=lambda p: p.name)
def test_every_config_key_the_doc_names_exists(doc):
    """A doc naming a key that was renamed sends the reader to nothing."""
    keys = config_keys()
    text = doc.read_text()
    # ALL_CAPS tokens of 4+ chars with an underscore, as the docs write them
    named = set(re.findall(r"`(\$?\{?[A-Z][A-Z0-9_]{3,}\}?)`", text))
    missing = []
    for raw in named:
        key = raw.strip("${}")
        if key in NOT_CONFIG or "_" not in key:
            continue
        if key not in keys:
            missing.append(key)
    assert missing == [], "%s names keys that do not exist: %s" % (doc.name,
                                                                  missing)


@pytest.mark.parametrize("doc", existing_docs(), ids=lambda p: p.name)
def test_every_path_the_doc_names_exists(doc):
    """Including the one pipeline.conf used to point at and which did not."""
    text = doc.read_text()
    missing = []
    for rel in set(re.findall(r"`((?:scripts|config|setup|sync|docs|todo|records)/[\w./-]+)`",
                              text)):
        rel = rel.rstrip(".")
        if rel.endswith("/"):
            ok = (ROOT / rel).is_dir()
        elif "<" in rel or "*" in rel:
            continue
        else:
            ok = (ROOT / rel).exists()
        if not ok:
            missing.append(rel)
    assert missing == [], "%s names paths that do not exist: %s" % (doc.name,
                                                                   missing)


@pytest.mark.parametrize("doc", existing_docs(), ids=lambda p: p.name)
def test_quoted_defaults_match_the_config(doc):
    """A table saying CTRL_FLANK is 25000 must be right about it."""
    text = doc.read_text()
    wrong = []
    # | `KEY` | `value` |  or  **`KEY=value`**
    for key, shown in re.findall(r"\|\s*`([A-Z][A-Z0-9_]+)`\s*\|\s*`?([^|`]+?)`?\s*\|",
                                 text):
        actual = config_default(key)
        if actual is None:
            continue
        shown = shown.strip()
        if shown in ("", "-"):
            continue
        # a range like `1`-`22` is prose for the CHROMS list
        if shown.replace("`", "").replace("–", "-") not in actual and \
                actual not in shown:
            wrong.append("%s: doc says %r, config has %r" % (key, shown, actual))
    assert wrong == [], "%s: %s" % (doc.name, wrong)


@pytest.mark.parametrize("doc", existing_docs(), ids=lambda p: p.name)
def test_check_file_patterns_match_what_the_checks_write(doc):
    """The probes all had to reconstruct these from four variables and one got
    them wrong, so the doc spelling them out must be right."""
    text = doc.read_text()
    checks = {
        "check_prep": "check_prep.%s.tsv",
        "check_sfs": "check_sfs.%s.tsv",
        "check_demo": "check_demo.%s.tsv",
        "check_dfe": "check_dfe.%s.tsv",
        "check_summary": "check_summary.%s.tsv",
    }
    for module, pattern in checks.items():
        src = (ROOT / "scripts" / "checks" / ("%s.py" % module)).read_text()
        stem = pattern.split("%s")[0]
        assert stem in src, "%s no longer writes %s" % (module, stem)
        if stem in text:
            # the doc must show it keyed by a run, not bare
            assert re.search(re.escape(stem) + r"\$?\{?\w", text), \
                "%s shows %s without a run key" % (doc.name, stem)


def test_agents_md_names_the_real_entry_point():
    """PIPELINE.md is a per-stage reference and defers to AGENTS.md for how to
    run the thing, so only the front door has to carry the command."""
    text = (ROOT / "AGENTS.md").read_text()
    assert "scripts/chains/run_all.sh" in text
    assert "DRY_RUN=1" in text


@pytest.mark.parametrize("doc", existing_docs(), ids=lambda p: p.name)
def test_outcome_vocabulary_matches_the_code(doc):
    """fit_outcomes.tsv's words were renamed twice; the doc must track them."""
    from importlib import import_module
    import sys
    sys.path.insert(0, str(ROOT / "scripts" / "checks"))
    kinds = set(import_module("check_dfe").OUTCOME_KINDS)
    text = doc.read_text()
    if "fit_outcomes.tsv" not in text:
        pytest.skip("%s does not describe the outcome log" % doc.name)
    for kind in kinds:
        assert "`%s`" % kind in text, \
            "%s does not mention the outcome %r" % (doc.name, kind)
    # and must not mention one that no longer exists
    for gone in ("no_likelihood", "crashed", "fitted", "converged_block"):
        assert "`%s`" % gone not in text, \
            "%s still mentions the removed outcome %r" % (doc.name, gone)


def test_agents_md_exists_and_is_the_front_door():
    """The gap every probe reported first: no entry point."""
    p = ROOT / "AGENTS.md"
    assert p.exists(), "AGENTS.md is the documented entry point"
    head = p.read_text()[:400]
    assert "Start here" in head


def test_pipeline_conf_points_only_at_docs_that_exist():
    """config/pipeline.conf referenced docs/PIPELINE.md before it existed."""
    text = (ROOT / "config" / "pipeline.conf").read_text()
    missing = [m for m in re.findall(r"((?:docs|todo)/[\w.]+\.md)", text)
               if not (ROOT / m).exists()]
    assert missing == [], "pipeline.conf points at %s" % missing


def test_no_script_points_at_a_doc_that_does_not_exist():
    missing = []
    for p in sorted(ROOT.rglob("*.sh")) + sorted(ROOT.rglob("*.py")):
        if ".git" in p.parts or "vendor" in p.parts:
            continue
        if p.name == "test_docs.py":       # it names the docs in order to test them
            continue
        for m in re.findall(r"((?:docs|todo)/[\w.]+\.md)", p.read_text()):
            if not (ROOT / m).exists():
                missing.append("%s -> %s" % (p.relative_to(ROOT), m))
    assert missing == [], missing
