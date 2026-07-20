#!/usr/bin/env python3
"""Phase 1 runner: feed one VCF to every adapter and print the summaries side by side.

The adapter machinery lives in adapters/registry.py (shared with the Phase 2
comparator). This script is just the human-facing view plus a coarse sanity
gate -- on a valid file every parser must agree exactly.
"""
from __future__ import annotations

import sys

from adapters.common import ParseSummary
from adapters.registry import run_all
from compare.compare import BUCKET_ALL_AGREE, BUCKET_ARTIFACT, classify


def _fmt_alts(alts) -> str:
    return ",".join(alts) if alts else "."


def print_summary(s: ParseSummary) -> None:
    print(f"== {s.tool} ==")
    if s.ok:
        print(f"   status : OK  ({s.num_records} records)")
        for (chrom, pos, ref, alts, info, samples) in s.records:
            info_str = f"  INFO={dict(info)}" if info else ""
            samp_str = f"  FMT={[dict(x) for x in samples]}" if samples else ""
            print(f"     {chrom}:{pos}  {ref} -> {_fmt_alts(alts)}{info_str}{samp_str}")
    else:
        print(f"   status : FAILED [{s.kind}]")
        print(f"   detail : {s.error}")
    print()


def sanity_gate(summaries: list[ParseSummary]) -> int:
    """Sanity gate on a valid file. Uses the same triage-aware classifier as the
    fuzzer, so a representation-only difference (e.g. float32 vs float64 on a
    numeric INFO value) passes -- it is quarantined, not a content disagreement."""
    print("---- sanity check ----")
    bucket, detail, _ = classify(summaries)
    if bucket == BUCKET_ALL_AGREE:
        print(f"PASS: {detail}")
        return 0
    if bucket == BUCKET_ARTIFACT:
        print(f"PASS (representation-only, quarantined): {detail}")
        return 0
    print(f"FAIL [{bucket}]: {detail}")
    return 1


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: python run.py <vcf_path>")
        return 2
    print(f"Input file: {argv[1]}\n")
    summaries = run_all(argv[1])
    for s in summaries:
        print_summary(s)
    return sanity_gate(summaries)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
