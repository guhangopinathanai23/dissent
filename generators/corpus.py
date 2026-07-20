"""Scaled corpus generator (Phase 3).

Stacks the Phase 2 idea into volume: for each file, apply 1-3 random mutation
operators at random record positions, starting from the valid seed. Every file
stays fully traceable -- a manifest records exactly which operator instances
produced it, so any catch is explainable ("f01234 = alt[r0='<DEL>'] + dup_record[r2]").

These are still understandable transformations, not random bytes; we just apply
them combinatorially and at scale.
"""
from __future__ import annotations

import json
import os
import random

from generators.mutate import _data_indices, _fields, _join, _lines, SEED_DEFAULT

# Deliberately nasty-but-namable field values.
BAD_POS = ["0", "-1", "-999", "12ABC", "14370.5", "999999999999", " 100", "100 ", "1e5", "0x1F", ""]
WEIRD_REF = ["", ".", "g", "N", "ACGT", "n", " A", "A ", "R", "U", "A.T"]
WEIRD_ALT = ["", ".", "A,", ",A", "a", "A,,T", "<DEL>", "N", "A;T", "*", "A/T", "A.T",
             "<INS:ME>", "G]17:198982]", "ACGTACGT", "A,<DEL>"]


class _NotApplicable(Exception):
    """Raised when an operator cannot apply to the current text (e.g. no data
    records remain after a prior CR-only line-ending mutation)."""


def _ord(di, i):
    return di.index(i)


def _pick(lines, rng):
    di = _data_indices(lines)
    if not di:
        raise _NotApplicable
    return rng.choice(di), di


def _require(fields, n):
    """A prior mutation may have shortened the record; skip ops that need a
    column that no longer exists."""
    if len(fields) < n:
        raise _NotApplicable


# Each operator: (text, rng) -> (new_text, label) or None if not applicable.

def op_pos_shift(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 2)
    try:
        p = int(f[1])
    except ValueError:
        return None
    d = rng.choice([-2, -1, 1, 2]); f[1] = str(p + d); lines[i] = _join(f)
    return "\n".join(lines), f"pos_shift[r{_ord(di, i)}{d:+d}]"


def op_pos_corrupt(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 2)
    v = rng.choice(BAD_POS); f[1] = v; lines[i] = _join(f)
    return "\n".join(lines), f"pos[r{_ord(di, i)}={v!r}]"


def op_ref_corrupt(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 4)
    v = rng.choice(WEIRD_REF); f[3] = v; lines[i] = _join(f)
    return "\n".join(lines), f"ref[r{_ord(di, i)}={v!r}]"


def op_alt_corrupt(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 5)
    v = rng.choice(WEIRD_ALT); f[4] = v; lines[i] = _join(f)
    return "\n".join(lines), f"alt[r{_ord(di, i)}={v!r}]"


def op_delete_field(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i])
    if len(f) <= 1:
        return None
    k = rng.randrange(len(f)); del f[k]; lines[i] = _join(f)
    return "\n".join(lines), f"del_col[r{_ord(di, i)},c{k}]"


def op_add_tab(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng)
    lines[i] = lines[i].replace("\t", "\t\t", 1)
    return "\n".join(lines), f"stray_tab[r{_ord(di, i)}]"


def op_dup_record(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng)
    lines.insert(i + 1, lines[i])
    return "\n".join(lines), f"dup_record[r{_ord(di, i)}]"


def op_truncate_record(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng)
    if len(lines[i]) < 2:
        return None
    lines[i] = lines[i][: rng.randint(1, len(lines[i]) - 1)]
    return "\n".join(lines), f"truncate[r{_ord(di, i)}]"


def op_blank_line(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng)
    lines.insert(i + 1, "")
    return "\n".join(lines), f"blank_line[after r{_ord(di, i)}]"


def op_drop_header(text, rng):
    lines = _lines(text)
    hdr = [j for j, l in enumerate(lines) if l.startswith("#")]
    if not hdr:
        return None
    j = rng.choice(hdr); dropped = lines[j][:24]; del lines[j]
    return "\n".join(lines), f"drop_header[{dropped!r}]"


def op_line_ending(text, rng):
    style = rng.choice(["crlf", "cr"])
    new = text.replace("\n", "\r\n") if style == "crlf" else text.replace("\n", "\r")
    return new, f"line_ending[{style}]"


def op_swap_ref_alt(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i])
    if len(f) < 5:
        return None
    f[3], f[4] = f[4], f[3]; lines[i] = _join(f)
    return "\n".join(lines), f"swap_ref_alt[r{_ord(di, i)}]"


def op_dup_pos(text, rng):
    lines = _lines(text); di = _data_indices(lines)
    if len(di) < 2:
        return None
    a = rng.randrange(len(di) - 1); i, j = di[a], di[a + 1]
    fi, fj = _fields(lines[i]), _fields(lines[j])
    if len(fi) < 2 or len(fj) < 2:
        return None
    fj[1] = fi[1]; lines[j] = _join(fj)
    return "\n".join(lines), f"dup_pos[r{a}->r{a + 1}]"


OPERATORS = [
    op_pos_shift, op_pos_corrupt, op_ref_corrupt, op_alt_corrupt, op_delete_field,
    op_add_tab, op_dup_record, op_truncate_record, op_blank_line, op_drop_header,
    op_line_ending, op_swap_ref_alt, op_dup_pos,
]


# ---- INFO value-mutation corpus (headers PRESERVED) -----------------------
# A fair test of whether the parsers agree on DECLARED INFO content: keep every
# header line intact and mutate only the INFO field VALUES, so all three parsers
# have the Type declaration and any disagreement is real content interpretation.

INFO_VALUES = ["0", "-1", "999999999", "1.5", "abc", ".", "0.5,0.6", "1e400",
               "NaN", "0x10", "", "3.141592653589793", "+5", "00", "1,2,3"]


def _parse_info(s):
    if s == ".":
        return []
    parts = []
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            parts.append((k, v))
        else:
            parts.append((kv, None))  # flag
    return parts


def _fmt_info(parts):
    if not parts:
        return "."
    return ";".join(k if v is None else f"{k}={v}" for k, v in parts)


def op_info_mutate(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 8)
    parts = _parse_info(f[7])
    if not parts:
        return None
    j = rng.randrange(len(parts)); k = parts[j][0]
    v = rng.choice(INFO_VALUES)
    parts[j] = (k, v); f[7] = _fmt_info(parts); lines[i] = _join(f)
    return "\n".join(lines), f"info[r{_ord(di, i)},{k}={v!r}]"


def op_info_extra_value(text, rng):
    # append an extra comma-value to a field (stress Number cardinality)
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 8)
    parts = _parse_info(f[7]); valued = [p for p in parts if p[1] is not None]
    if not valued:
        return None
    j = parts.index(rng.choice(valued)); k, v = parts[j]
    parts[j] = (k, f"{v},99"); f[7] = _fmt_info(parts); lines[i] = _join(f)
    return "\n".join(lines), f"info_extra[r{_ord(di, i)},{k}]"


INFO_OPERATORS = [op_info_mutate, op_info_mutate, op_info_extra_value]


def generate_info_corpus(n: int, out_dir: str = "data/corpus_info",
                         seed_path: str = SEED_DEFAULT, rng_seed: int = 4242) -> dict:
    """Header-preserving corpus: mutate only INFO values, keep all header lines,
    so the comparison tests agreement on DECLARED INFO content."""
    rng = random.Random(rng_seed)
    with open(seed_path) as f:
        seed = f.read()
    os.makedirs(out_dir, exist_ok=True)
    manifest = {}
    for idx in range(n):
        k = rng.choice([1, 1, 2])
        text = seed
        labels = []
        attempts = 0
        while len(labels) < k and attempts < 10:
            attempts += 1
            try:
                res = rng.choice(INFO_OPERATORS)(text, rng)
            except _NotApplicable:
                continue
            if res is None:
                continue
            text, label = res
            labels.append(label)
        if not labels:
            res = op_info_mutate(seed, rng)
            text, label = res if res else (seed, "noop")
            labels = [label]
        fname = f"i{idx:05d}.vcf"
        with open(os.path.join(out_dir, fname), "w", newline="") as f:
            f.write(text)
        manifest[fname] = labels
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f)
    return manifest


# ---- genotype value-mutation corpus (headers PRESERVED) -------------------
# Fair test of whether the parsers agree on per-sample FORMAT/genotype content:
# keep the header intact and mutate GT / FORMAT subfield VALUES. Exercises the
# interpret-not-passthrough territory (phasing, ploidy, allele order, missing).

GT_VALUES = ["0/1", "1/0", "0|1", "1|0", "1", "0", "1/1", "2/2", "1|1|1", "./.",
             ".", "0|.", "1/2", "3/0", "-1/0", "00/1", "0/0/0"]
FMT_VALUES = ["999", "-1", "1.5", "abc", ".", "1,2", "", "0/0", "NaN"]


def op_gt_mutate(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 10)
    si = rng.randrange(9, len(f))
    sub = f[si].split(":")
    v = rng.choice(GT_VALUES); sub[0] = v  # GT is always the first FORMAT subfield
    f[si] = ":".join(sub); lines[i] = _join(f)
    return "\n".join(lines), f"gt[r{_ord(di, i)},s{si - 9}={v!r}]"


def op_fmt_mutate(text, rng):
    lines = _lines(text); i, di = _pick(lines, rng); f = _fields(lines[i]); _require(f, 10)
    si = rng.randrange(9, len(f)); sub = f[si].split(":")
    if len(sub) < 2:
        return None
    j = rng.randrange(1, len(sub))  # a non-GT subfield
    v = rng.choice(FMT_VALUES); sub[j] = v
    f[si] = ":".join(sub); lines[i] = _join(f)
    return "\n".join(lines), f"fmt[r{_ord(di, i)},s{si - 9},sub{j}={v!r}]"


GT_OPERATORS = [op_gt_mutate, op_gt_mutate, op_fmt_mutate]


def generate_gt_corpus(n: int, out_dir: str = "data/corpus_gt",
                       seed_path: str = SEED_DEFAULT, rng_seed: int = 9999) -> dict:
    """Header-preserving corpus that mutates GT / FORMAT subfield values only."""
    rng = random.Random(rng_seed)
    with open(seed_path) as f:
        seed = f.read()
    os.makedirs(out_dir, exist_ok=True)
    manifest = {}
    for idx in range(n):
        k = rng.choice([1, 1, 2])
        text = seed
        labels = []
        attempts = 0
        while len(labels) < k and attempts < 10:
            attempts += 1
            try:
                res = rng.choice(GT_OPERATORS)(text, rng)
            except _NotApplicable:
                continue
            if res is None:
                continue
            text, label = res
            labels.append(label)
        if not labels:
            res = op_gt_mutate(seed, rng)
            text, label = res if res else (seed, "noop")
            labels = [label]
        fname = f"g{idx:05d}.vcf"
        with open(os.path.join(out_dir, fname), "w", newline="") as f:
            f.write(text)
        manifest[fname] = labels
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f)
    return manifest


def generate_corpus(n: int, out_dir: str = "data/corpus", seed_path: str = SEED_DEFAULT,
                    rng_seed: int = 1337) -> dict:
    """Write n stacked-mutation files + a provenance manifest. Deterministic
    given rng_seed, so a run is fully reproducible."""
    rng = random.Random(rng_seed)
    with open(seed_path) as f:
        seed = f.read()
    os.makedirs(out_dir, exist_ok=True)
    manifest = {}
    for idx in range(n):
        k = rng.choice([1, 1, 1, 2, 2, 3])  # mostly 1-2, some 3-stacks
        text = seed
        labels = []
        attempts = 0
        while len(labels) < k and attempts < 12:
            attempts += 1
            try:
                res = rng.choice(OPERATORS)(text, rng)
            except _NotApplicable:
                continue
            if res is None:
                continue
            text, label = res
            labels.append(label)
        if not labels:  # guarantee at least one mutation
            text, label = op_pos_corrupt(seed, rng)
            labels = [label]
        fname = f"f{idx:05d}.vcf"
        with open(os.path.join(out_dir, fname), "w", newline="") as f:
            f.write(text)
        manifest[fname] = labels
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f)
    return manifest


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    m = generate_corpus(n)
    print(f"generated {len(m)} files in data/corpus/ (+ manifest.json)")
    for name in list(m)[:8]:
        print(f"  {name}: {m[name]}")
