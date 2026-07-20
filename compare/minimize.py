"""Test-case minimization by delta debugging.

When a file triggers a real disagreement, shrink it to the smallest file that
still triggers the SAME disagreement -- turning a big messy catch into a tiny,
shareable example.

"Same disagreement" is defined by compare.signature(): the bucket, how the
parsers that parsed partition by content, and which parsers failed. Requiring
the signature to stay identical means we never accidentally shrink our way into
a *different* bug. It also means structurally-required lines (e.g. #CHROM) are
kept automatically: removing them changes the signature, so the oracle rejects
the removal.
"""
from __future__ import annotations

import os

from compare.compare import compare_file, signature


def minimize(vcf_path: str, target_sig, workdir: str, tag: str = "min"):
    """Greedy line-level 1-minimization. Returns (minimized_text, num_lines)."""
    with open(vcf_path, newline="") as f:
        original = f.read()
    trailing_nl = original.endswith("\n")
    lines = original.split("\n")
    if trailing_nl and lines and lines[-1] == "":
        lines = lines[:-1]

    candidate = os.path.join(workdir, f"_min_{tag}.vcf")

    def reproduces(cand_lines) -> bool:
        text = "\n".join(cand_lines) + ("\n" if trailing_nl else "")
        with open(candidate, "w", newline="") as f:
            f.write(text)
        return signature(compare_file(candidate)) == target_sig

    # Guard: the original must reproduce under our oracle.
    if not reproduces(lines):
        return original, len([l for l in lines if l != ""])

    changed = True
    while changed and len(lines) > 1:
        changed = False
        i = 0
        while i < len(lines):
            trial = lines[:i] + lines[i + 1:]
            if trial and reproduces(trial):
                lines = trial
                changed = True  # re-try at same index (list shifted)
            else:
                i += 1

    minimized = "\n".join(lines) + ("\n" if trailing_nl else "")
    return minimized, len(lines)
