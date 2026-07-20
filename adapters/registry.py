"""Registry of parser adapters + the isolated runner that drives them.

Single source of truth shared by run.py (Phase 1 side-by-side view) and
compare/ (Phase 2 classification). Each adapter is an external command that
takes a VCF path and prints a ParseSummary as JSON. The runner isolates each in
its own subprocess with a timeout, so a parser that segfaults or hangs is caught
here instead of taking down the caller.
"""
from __future__ import annotations

import os
import subprocess
import sys

from adapters.common import (
    KIND_CRASH,
    KIND_PLUMBING,
    KIND_TIMEOUT,
    ParseSummary,
)

# registry.py lives in adapters/, so the repo root is two levels up.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NOODLES_BIN = os.path.join(
    REPO_ROOT, "adapters", "noodles_adapter", "target", "release", "noodles_adapter"
)

# (display name, argv prefix run as `<prefix...> <vcf_path>`).
# Two htslib-independent Python parsers plus the pure-Rust noodles binary.
ADAPTERS = [
    ("pysam (htslib)", [sys.executable, "-m", "adapters.pysam_adapter"]),
    ("vcfpy (pure-python)", [sys.executable, "-m", "adapters.vcfpy_adapter"]),
    ("noodles (rust)", [NOODLES_BIN]),
]
TIMEOUT_S = 10


def run_adapter(display: str, argv_prefix: list[str], vcf_path: str) -> ParseSummary:
    """Run one adapter as an isolated subprocess and classify the outcome."""
    try:
        proc = subprocess.run(
            argv_prefix + [vcf_path],
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
            cwd=REPO_ROOT,
        )
    except subprocess.TimeoutExpired:
        return ParseSummary(
            tool=display, kind=KIND_TIMEOUT,
            error=f"no result within {TIMEOUT_S}s (killed)",
        )
    except FileNotFoundError:
        return ParseSummary(
            tool=display, kind=KIND_PLUMBING,
            error=f"adapter not found: {argv_prefix[0]} "
                  f"(build it with `cargo build --release`)",
        )

    if proc.returncode != 0:
        where = f"signal {-proc.returncode}" if proc.returncode < 0 else f"exit {proc.returncode}"
        tail = " | ".join((proc.stderr or "").strip().splitlines()[-3:])
        return ParseSummary(
            tool=display, kind=KIND_CRASH,
            error=f"worker died ({where}): {tail}" if tail else f"worker died ({where})",
        )

    try:
        return ParseSummary.from_json(proc.stdout)
    except Exception as e:  # noqa: BLE001
        return ParseSummary(
            tool=display, kind=KIND_PLUMBING,
            error=f"unreadable worker output ({e}); stdout={proc.stdout!r}",
        )


def run_all(vcf_path: str) -> list[ParseSummary]:
    """Run every registered adapter on one file, in registry order."""
    abspath = os.path.abspath(vcf_path)  # subprocesses run with cwd=REPO_ROOT
    return [run_adapter(disp, prefix, abspath) for (disp, prefix) in ADAPTERS]
