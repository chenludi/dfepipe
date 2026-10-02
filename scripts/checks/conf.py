#!/usr/bin/env python3
"""Read default values out of config/pipeline.conf.

Every check takes its thresholds on the command line and the stage script
passes the configured value, so this is only the default for running a check by
hand. It reads pipeline.conf rather than repeating the numbers, so a threshold
exists in exactly one place.

    from conf import conf
    conf.get_float("HASQUI_MIN_FRAC", 0.90)

The fallback should match the configured value. If the two drift, a check run by
hand and the same check run by its stage use different numbers, which is the one
thing this module exists to prevent; tests/test_naming_honesty.py enforces it.

Only `KEY="${KEY:-value}"` lines are read; anything else in the file is
ignored. If the file cannot be read, the fallback passed at the call site is
used, so a check still runs outside a checkout.
"""
import os
import re

_LINE = re.compile(
    r'^([A-Za-z_][A-Za-z0-9_]*)="\$\{\1:-([^}]*)\}"\s*(?:#.*)?$')

DEFAULT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "config", "pipeline.conf")


class Conf(object):
    def __init__(self, path=None):
        self.path = path or os.environ.get("DFEPIPE_CONF") or DEFAULT_PATH
        self.values = {}
        try:
            with open(self.path) as fh:
                for line in fh:
                    # some lines hold several assignments: `A="..."; B="..."`
                    for part in line.rstrip("\n").split(";"):
                        m = _LINE.match(part.strip())
                        if m:
                            self.values[m.group(1)] = m.group(2)
        except (IOError, OSError):
            pass

    def get(self, key, fallback):
        return self.values.get(key, fallback)

    def get_float(self, key, fallback):
        try:
            return float(self.values[key])
        except (KeyError, ValueError):
            return fallback

    def get_int(self, key, fallback):
        try:
            return int(self.values[key])
        except (KeyError, ValueError):
            return fallback


conf = Conf()
