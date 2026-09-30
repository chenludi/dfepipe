#!/usr/bin/env python3
"""Checkpoint for the cache stage: is the spectra cache usable?

FOR REFERENCE ONLY, as everywhere in this pipeline: a WARN means "unusual
compared with previous runs", not "wrong". The cache itself remains the arbiter.

FAIL is reserved for a cache the DFE stage cannot work with, and does block a
chain (FORCE=1 overrides):
    the .bpkl is missing        GenerateCache did not finish
    the .bpkl will not unpickle the job was killed mid-write
    the file is below CACHE_MIN_SIZE_FRAC of the size gamma_pts and PROJECTION
                               imply. A cache is a gamma_pts x (projection+1)
                               float64 array, so its size follows from the
                               config and NOT from the data - see
                               expected_bytes(). Anything much smaller was
                               truncated, or was built with fewer gamma points
                               than this run asked for. A fixed byte count
                               cannot do this job: 600000, inherited from
                               dadi-cli_cache_simp.sh:129 where it only printed
                               a warning, is larger than a correct marker-1
                               cache (measured: 559,739 bytes) and would have
                               failed every standard run.
    fewer than CACHE_MIN_USABLE_FRAC of the gamma points produced a usable
                               spectrum. dadi extrapolates each gamma from three
                               grids; where that breaks down the row comes back
                               summing to <= 0, and the DFE can only be
                               evaluated over the gamma points that survived.
                               Widen or shift the gamma range, or raise the
                               grids. Fig. S16 of the manuscript (proportion of
                               mutations per |s| class across log-shape by
                               log-scale) is the reference for the range.

WARN is advisory: the covered gamma range topping out well short of
2*GAMMA_S_HIGH*Na, which leaves strong selection unexplorable.

Unlike check_prep/check_sfs/check_demo this writes ONE FILE PER BED, into the
bed's own run directory, because the cache stage is submitted as an array over
beds: a single stage-level file would have 10-hour tasks racing to rewrite it.
The chain points gate.py at the per-bed file.
"""
import argparse
import math
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from conf import conf  # noqa: E402

# The covered upper bound may fall no more than this fraction short of
# 2*GAMMA_S_HIGH*Na before it is worth remarking on. A factor of ten: the gammas
# are log-spaced, so anything tighter would fire on rounding.
RANGE_TOLERANCE = 0.1

# Pickling the float64 array costs about 40% on top of it. Measured on three real
# caches: 300x161 -> 559,739 B (payload 386,400, x1.45), 300x101 -> ~335 kB
# (242,400, x1.4), 1000x161 -> ~1.8 MB (1,288,000, x1.4).
PICKLE_OVERHEAD = 1.4

CHECK_NAME = "check_cache.tsv"
SPEC_NAME = "spec.txt"
# spec.txt keys that are numbers; everything else stays a string.
_FLOAT_KEYS = ("Na", "gamma_lo", "gamma_hi", "gamma_pts", "projection")


# ------------------------------------------------------------------- verdict --

def expected_bytes(gamma_pts, projection):
    """The size a correct cache should have, in bytes.

    A cache is a gamma_pts x (projection+1) float64 array, plus PICKLE_OVERHEAD.
    It depends only on the config, never on the data, which is what makes a
    fraction of it a usable floor: CACHE_MIN_SIZE_FRAC of this.
    """
    return int(int(gamma_pts) * (int(projection) + 1) * 8 * PICKLE_OVERHEAD)


def assess_cache(nbytes, usable_frac, gamma_lo, gamma_hi, na, min_bytes,
                 min_usable, s_high=0.5, killed=False):
    """(status, reasons) for one cache. See the module docstring.

    `nbytes` is None when the .bpkl is absent, `usable_frac` and the gamma
    bounds None when it could not be read. `na` is the Na the bounds were
    derived from, None when unknown; `s_high` is GAMMA_S_HIGH, so the range
    warning compares against the bound this run actually asked for.
    """
    if nbytes is None and killed:
        return "FAIL", ["no cache file, and the log stops mid-build: the job was "
                        "probably killed at its h_rt limit. Raise SGE_CACHE_RT "
                        "and re-run - a partial cache cannot be resumed."]
    if nbytes is None:
        return "FAIL", ["the cache file is missing: dadi-cli GenerateCache did "
                        "not finish. Check the SGE log for this bed"]
    if nbytes < min_bytes:
        return "FAIL", ["cache size %d below %d: the file is truncated, or was "
                        "built with fewer gamma points than this run asked for"
                        % (nbytes, min_bytes)]
    if usable_frac is None:
        return "FAIL", ["the cache file will not unpickle: the job was probably "
                        "killed while writing it. Re-run the stage"]
    if usable_frac < min_usable:
        return "FAIL", ["only %.0f%% of gamma points produced a spectrum: "
                        "extrapolation failed; adjust the gamma range"
                        % (100 * usable_frac)]

    status, reasons = "PASS", []
    if na is not None and gamma_hi is not None:
        expected_hi = 2 * s_high * na
        if gamma_hi < expected_hi * RANGE_TOLERANCE:
            status = "WARN"
            reasons.append(
                "gamma range tops out at %.3g, well below 2*%g*Na=%.3g: strong "
                "selection is unexplorable" % (gamma_hi, s_high, expected_hi))
    return status, reasons


# ---------------------------------------------------------------- statistics --

def cache_stats(gammas, spectra):
    """(n_total, n_usable, gamma_lo, gamma_hi) over a cache's gammas and spectra.

    dadi's Cache1D stores `gammas` as a negative, log-spaced float array and
    `spectra` as one row per gamma (a plain (gamma_pts x ns+1) float array, not
    Spectrum objects). A row that came back non-finite, or summing to <= 0, is
    where the extrapolation broke down and is not usable.

    The bounds are |gamma| over the USABLE rows only: that is the range a DFE
    can actually be integrated over. They are None when nothing is usable.

    Deliberately numpy-free, so it can be tested without dadi installed; a
    numpy array iterates row-wise and works unchanged.
    """
    usable_gammas = []
    n_total = 0
    for i, row in enumerate(spectra):
        n_total += 1
        if row is None:
            continue
        try:
            values = [float(v) for v in row]
        except (TypeError, ValueError):
            continue
        if not all(math.isfinite(v) for v in values):
            continue
        if sum(values) <= 0:
            continue
        if i < len(gammas):
            usable_gammas.append(abs(float(gammas[i])))

    if not usable_gammas:
        return n_total, 0, None, None
    return n_total, len(usable_gammas), min(usable_gammas), max(usable_gammas)


def load_cache_stats(path):
    """cache_stats() for a dadi cache pickle, or None if it cannot be read.

    Unpickling a real cache imports dadi (the pickle names
    dadi.DFE.Cache1D_mod.Cache1D), so this only ever runs inside the pipeline's
    conda environment.
    """
    try:
        with open(path, "rb") as fh:
            cache = pickle.load(fh)
    except Exception:            # a truncated or half-written pickle
        return None
    gammas = getattr(cache, "gammas", None)
    spectra = getattr(cache, "spectra", None)
    if gammas is None or spectra is None:
        return None
    return cache_stats(gammas, spectra)


# ------------------------------------------------------------------- readers --

def read_spec(path):
    """run_cache.sh's spec.txt as a dict, `key=value` per line; {} if absent.

    This is where Na comes from: run_cache.sh records the value it took from
    best_inference.txt and derived the gamma bounds with, so the check reports
    the number this cache was actually built on and needs no path to the demo
    stage. A value that should be a number but is not is dropped rather than
    kept as a string, so a caller can test `"Na" in spec`.
    """
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return {}
    out = {}
    for line in lines:
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        out[key.strip()] = value.strip()
    for key in _FLOAT_KEYS:
        if key in out:
            try:
                out[key] = float(out[key])
            except ValueError:
                del out[key]
    return out


def _bed_base(bed):
    base = os.path.basename(bed)
    return base[:-4] if base.endswith(".bed") else base


# ---------------------------------------------------------------------- rows --

def row_for(cache_dir, base, marker, min_size_frac, min_usable, s_high,
            gamma_pts, projection):
    """The single check row for one bed."""
    bpkl = os.path.join(cache_dir, "cache_%s.spectra.bpkl" % marker)
    spec = read_spec(os.path.join(cache_dir, SPEC_NAME))
    na = spec.get("Na")

    # A cache.log that exists while the .bpkl does not means the build started
    # and did not finish - almost always the h_rt kill, since GenerateCache
    # writes the pickle only at the end.
    killed = (not os.path.isfile(bpkl)
              and os.path.isfile(os.path.join(cache_dir, "cache.log")))

    nbytes = os.path.getsize(bpkl) if os.path.isfile(bpkl) else None
    stats = load_cache_stats(bpkl) if nbytes else None
    if stats is None:
        n_gamma, n_usable, lo, hi, frac = 0, 0, None, None, None
    else:
        n_gamma, n_usable, lo, hi = stats
        frac = (float(n_usable) / n_gamma) if n_gamma else 0.0

    expected = expected_bytes(gamma_pts, projection)
    min_bytes = int(min_size_frac * expected)
    status, reasons = assess_cache(nbytes, frac, lo, hi, na, min_bytes,
                                   min_usable, s_high=s_high, killed=killed)
    return {
        "item": base,
        "status": status,
        "reason": "; ".join(reasons) or "NA",
        "cache_file": bpkl,
        "file_bytes": nbytes if nbytes is not None else "NA",
        "expected_bytes": expected,
        "min_bytes": min_bytes,
        "n_gamma": n_gamma if stats is not None else "NA",
        "n_usable": n_usable if stats is not None else "NA",
        "usable_frac": round(frac, 4) if frac is not None else "NA",
        "gamma_lo": lo if lo is not None else "NA",
        "gamma_hi": hi if hi is not None else "NA",
        "Na": na if na is not None else "NA",
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="the cache output root, $OUT/cache")
    ap.add_argument("--run", required=True,
                    help="the run directory name, <DEMO_RUN>-<CACHE_MARKER>-<CACHE_RUN>")
    ap.add_argument("--marker", default=conf.get("CACHE_MARKER", "1"),
                    help="CACHE_MARKER, the <marker> in cache_<marker>.spectra.bpkl")
    ap.add_argument("--min-size-frac", type=float,
                    default=conf.get_float("CACHE_MIN_SIZE_FRAC", 0.5),
                    help="CACHE_MIN_SIZE_FRAC; fail below this fraction of the "
                         "expected gamma_pts x (projection+1) x 8 bytes")
    ap.add_argument("--gamma-pts", type=int, required=True, help="GAMMA_PTS")
    ap.add_argument("--projection", type=int, required=True, help="PROJECTION")
    ap.add_argument("--min-usable", type=float,
                    default=conf.get_float("CACHE_MIN_USABLE_FRAC", 0.9),
                    help="CACHE_MIN_USABLE_FRAC")
    ap.add_argument("--gamma-s-high", type=float,
                    default=conf.get_float("GAMMA_S_HIGH", 0.5),
                    help="GAMMA_S_HIGH; the cache should reach 2*GAMMA_S_HIGH*Na")
    ap.add_argument("beds", nargs="+", help="the input beds the stage was run on")
    a = ap.parse_args(argv)

    for bed in a.beds:
        base = _bed_base(bed)
        cache_dir = os.path.join(a.out, base, a.run)
        row = row_for(cache_dir, base, a.marker, a.min_size_frac, a.min_usable,
                      a.gamma_s_high, a.gamma_pts, a.projection)
        # One file per bed, in the bed's own run directory: array tasks must not
        # race. Written even when the directory is otherwise empty, so a missing
        # cache is recorded rather than silently absent.
        if not os.path.isdir(cache_dir):
            os.makedirs(cache_dir, exist_ok=True)
        dest = os.path.join(cache_dir, CHECK_NAME)
        write_check(dest, [row])
        print("wrote %s" % dest)
        print("  for reference only: these flag caches worth a second look, not "
              "an accept/reject test")
        print("  %s %s: %s gamma points, %s usable, |gamma| %s..%s, Na=%s"
              % (row["status"], base, row["n_gamma"], row["n_usable"],
                 row["gamma_lo"], row["gamma_hi"], row["Na"]))
        if row["status"] not in ("PASS", "SKIP"):
            print("  %s %s: %s" % (row["status"], base, row["reason"]))
    # Always 0: the verdict is the check file, and gate.py is what blocks a
    # chain. A non-zero exit here would look like the stage itself had failed.
    return 0


if __name__ == "__main__":
    sys.exit(main())
