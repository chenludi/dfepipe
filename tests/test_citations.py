"""Every `file.ext:line` citation in this package must resolve.

The package explains itself by pointing at the original scripts it was ported
from, so a citation that points at the wrong line is not a cosmetic error: it is
the evidence for a claim about what the original did. An audit of this package
found nine citations that were off by anywhere from one line to a hundred, each
written from memory rather than checked, and each one made its surrounding claim
unverifiable.

This test cannot know whether a cited line *means* what the comment says. It can
know that the file exists and the line exists and is not blank, which is what
every one of those nine failed.

Originals live in the parent directory (the user's working scripts). Some
citations name dadi or dadi-cli files, which live inside setup/vendor/*.tar.gz;
those are resolved by extracting the tarballs once.
"""
import pathlib
import re
import sys
import tarfile
import tempfile

import pytest

PKG = pathlib.Path(__file__).resolve().parent.parent
ORIG = PKG.parent

# `name.ext:12` or `name.ext:12-34`, and the `,56` tails of a compound citation
# like `concate_mut.R:27,71` - which is one citation of two lines, not two.
CITATION = re.compile(
    r'([A-Za-z0-9_][A-Za-z0-9_.\-]*\.(?:sh|py|R)):(\d+)(?:-(\d+))?((?:,\d+(?:-\d+)?)*)')

SKIP_PARTS = {".git", "__pycache__", ".pytest_cache"}
READ_SUFFIXES = {".sh", ".py", ".R", ".md", ".conf"}


@pytest.fixture(scope="module")
def index():
    """basename -> path, over the originals, our own files, and the vendored deps."""
    out = {}
    tmp = tempfile.mkdtemp()
    for tar in sorted((PKG / "setup" / "vendor").glob("*.tar.gz")):
        try:
            with tarfile.open(tar) as t:
                t.extractall(tmp, filter="data")
        except (tarfile.TarError, OSError, TypeError):
            pass
    # lowest priority first, so an original shadows a vendored file of the same name
    for root in (pathlib.Path(tmp), PKG, ORIG):
        for p in root.rglob("*"):
            if p.is_file() and p.suffix in (".sh", ".py", ".R"):
                if any(s in p.parts for s in SKIP_PARTS):
                    continue
                out[p.name] = p
    return out


def package_files():
    for p in sorted(PKG.rglob("*")):
        if not p.is_file() or p.suffix not in READ_SUFFIXES:
            continue
        if any(s in p.parts for s in SKIP_PARTS):
            continue
        if "vendor" in p.parts or "tests" in p.parts:
            continue
        yield p


def citations():
    """(source file, source line, raw citation, target name, [line numbers])."""
    for p in package_files():
        for lineno, line in enumerate(p.read_text(errors="ignore").splitlines(), 1):
            for m in CITATION.finditer(line):
                name, a, b, tail = m.groups()
                nums = [int(a)] + ([int(b)] if b else [])
                for extra in re.findall(r'(\d+)', tail or ""):
                    nums.append(int(extra))
                yield p, lineno, m.group(0), name, nums


def test_some_citations_exist():
    """A regex that silently matches nothing would make every test below pass."""
    found = list(citations())
    assert len(found) > 80, "only found %d citations; the regex is probably wrong" % len(found)


def test_every_cited_file_can_be_found(index):
    missing = {}
    for p, lineno, raw, name, _ in citations():
        if name not in index:
            missing.setdefault(name, []).append("%s:%d" % (p.relative_to(PKG), lineno))
    assert missing == {}, "cited files not found anywhere: %s" % missing


def test_every_cited_line_is_within_its_file(index):
    bad = []
    for p, lineno, raw, name, nums in citations():
        tgt = index.get(name)
        if tgt is None:
            continue
        n = len(tgt.read_text(errors="ignore").splitlines())
        for want in nums:
            if want > n:
                bad.append("%s:%d cites %s but %s has only %d lines"
                           % (p.relative_to(PKG), lineno, raw, name, n))
    assert bad == [], "\n".join(bad)


def test_no_citation_points_at_a_blank_line(index):
    """A blank target is the signature of a line number that drifted."""
    bad = []
    for p, lineno, raw, name, nums in citations():
        tgt = index.get(name)
        if tgt is None:
            continue
        lines = tgt.read_text(errors="ignore").splitlines()
        for want in nums:
            if want <= len(lines) and not lines[want - 1].strip():
                bad.append("%s:%d cites %s -> %s:%d is blank"
                           % (p.relative_to(PKG), lineno, raw, name, want))
    assert bad == [], "\n".join(bad)


def test_citations_to_our_own_scripts_are_not_stale(index):
    """A citation into this package drifts every time the file is edited, and
    nothing else would notice. run_prep.sh:112 is cited from three files."""
    ours = {p.name for p in PKG.rglob("*")
            if p.is_file() and p.suffix in (".sh", ".py", ".R")
            and not any(s in p.parts for s in SKIP_PARTS)}
    checked = 0
    bad = []
    for p, lineno, raw, name, nums in citations():
        if name not in ours or (ORIG / name).is_file():
            continue          # an original of the same name wins; not ours
        tgt = PKG / "scripts"
        hit = [q for q in PKG.rglob(name) if q.is_file()]
        if not hit:
            continue
        lines = hit[0].read_text(errors="ignore").splitlines()
        checked += 1
        for want in nums:
            if want > len(lines) or not lines[want - 1].strip():
                bad.append("%s:%d cites %s, which is stale"
                           % (p.relative_to(PKG), lineno, raw))
    assert bad == [], "\n".join(bad)
