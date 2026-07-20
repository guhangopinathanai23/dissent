#!/usr/bin/env python3
"""Phase 2 demo: generate the named mutations, classify each, print a table.

For every mutation we run the file through all three adapters, bucket the
result, and (for the interesting buckets) save the offending file + all three
summaries under findings/. At the end we print a table and a detailed
walkthrough of every interesting catch.

Deliberately small -- a few dozen understandable files. Scaling is Phase 3.
"""
from __future__ import annotations

import sys
from collections import Counter

from compare.compare import (
    INTERESTING,
    compare_file,
    marker,
    render_text_report,
    save_finding,
    sides,
)
from generators.mutate import SEED_DEFAULT, generate_all


def main() -> int:
    pairs = generate_all(SEED_DEFAULT)
    results = []
    for name, path in pairs:
        comp = compare_file(path)
        saved = save_finding(comp, label=name) if comp.bucket in INTERESTING else None
        results.append((name, comp, saved))

    # ---- summary table ----
    print("=" * 82)
    print(f"    {'mutation':<22} {'bucket':<20} split")
    print("-" * 82)
    for name, comp, _ in results:
        print(f"{marker(comp)} {name:<22} {comp.bucket:<20} {sides(comp)}")
    print("=" * 82)
    print("buckets:", dict(Counter(comp.bucket for _, comp, _ in results)))
    print("legend: *** all parse but differ | ** split + acceptors also disagree "
          "| * accept/reject split")

    # ---- detailed walkthrough of interesting catches ----
    interesting = [(n, c, s) for (n, c, s) in results if c.bucket in INTERESTING]
    if not interesting:
        print("\nNo interesting catches (everything was all_agree or all_fail).")
        return 0

    print("\n" + "#" * 82)
    print(f"# {len(interesting)} INTERESTING CATCH(ES) -- detail below (also saved under findings/)")
    print("#" * 82)
    for name, comp, saved in interesting:
        print("\n" + render_text_report(comp, name))
        print(f"[saved: {saved}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
