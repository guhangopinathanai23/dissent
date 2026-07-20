"""Comparator: run a VCF through all adapters, classify, and triage.

Correctness lives here. Two ideas layered on top of each other:

1. CONTENT-DISAGREEMENT-FIRST classification. A crash is loud -- someone
   notices. A *silent content disagreement* (two parsers both succeed but
   report different variants) corrupts an analysis with nobody the wiser. So
   that is the top-priority bucket, above accept/reject splits.

2. A NORMALIZATION-ARTIFACT TRIAGE GATE. Before we call a content disagreement
   a real finding, we re-compare the parsers' output under a deeper
   "canonicalized-output" normalization (missing-token + case + whitespace). If
   the disagreement collapses under that, it was only a representation
   difference -- bucketed as `normalization_artifact` and NOT reported as a
   finding. This is the automated version of the .-vs-'' investigation, so we
   never have to re-litigate representation noise by hand at scale.

Buckets, highest priority first:
   content_disagreement   >=2 parse, differ, SURVIVES triage      (GOLD, a finding)
   normalization_artifact >=2 parse, differ, collapses under triage (quarantined)
   split_decision         some parse some reject; parsers agree     (lower value)
   all_fail               none parse
   all_agree              all parse and agree
   harness_error          our own plumbing failed
"""
from __future__ import annotations

import dataclasses
import json
import os
import shutil
import struct
from typing import Dict, List, Tuple

from adapters.common import (
    KIND_CRASH,
    KIND_PLUMBING,
    KIND_TIMEOUT,
    ParseSummary,
    canonicalize_alts,
)
from adapters.registry import REPO_ROOT, run_all

BUCKET_CONTENT = "content_disagreement"
BUCKET_ARTIFACT = "normalization_artifact"
BUCKET_SPLIT = "split_decision"
BUCKET_ALL_FAIL = "all_fail"
BUCKET_ALL_AGREE = "all_agree"
BUCKET_HARNESS_ERROR = "harness_error"

# Real findings the fuzzer reports + minimizes.
FINDING_BUCKETS = {BUCKET_CONTENT}
# Saved for inspection but NOT headline findings.
QUARANTINE_BUCKETS = {BUCKET_ARTIFACT}
# Everything worth saving/eyeballing (used by the Phase 2 demo).
INTERESTING = {BUCKET_CONTENT, BUCKET_ARTIFACT, BUCKET_SPLIT}

PRIORITY = {
    BUCKET_CONTENT: 0,
    BUCKET_ARTIFACT: 1,
    BUCKET_SPLIT: 2,
    BUCKET_ALL_FAIL: 3,
    BUCKET_ALL_AGREE: 4,
    BUCKET_HARNESS_ERROR: 5,
}


@dataclasses.dataclass
class Comparison:
    vcf_path: str
    bucket: str
    detail: str
    summaries: List[ParseSummary]
    all_parsed: bool = False  # True => a *silent* disagreement (no parser failed)


def short(tool: str) -> str:
    return tool.split()[0]


# ---- grouping keys --------------------------------------------------------

def _content_groups(ok: List[ParseSummary]) -> Dict[Tuple, List[str]]:
    """Group successful parsers by their (already base-canonicalized) records."""
    groups: Dict[Tuple, List[str]] = {}
    for s in ok:
        groups.setdefault(tuple(s.records), []).append(s.tool)
    return groups


def _triage_alleles(alts: Tuple[str, ...]) -> Tuple[str, ...]:
    # Deeper representation collapse: base missing-token canon + case-fold +
    # whitespace strip. If a disagreement vanishes here, it was representation.
    return canonicalize_alts(tuple(a.strip().upper() for a in alts))


def _f32(x: float) -> float:
    """Narrow a float to float32 precision and back."""
    return struct.unpack("<f", struct.pack("<f", x))[0]


def _triage_num(v):
    # NARROW float64 to float32 for the triage key. This collapses ONLY the
    # htslib/noodles(f32)-vs-vcfpy(f64) storage difference: two values that
    # agree to float32 precision map to the same key. A genuine numeric
    # difference (e.g. 0.017 vs 0.018) differs at float32 and survives as a
    # finding. Integers/strings/flags/None are untouched.
    if isinstance(v, bool):
        return v
    if isinstance(v, float):
        return _f32(v)
    if isinstance(v, tuple):
        return tuple(_triage_num(x) for x in v)
    return v


def _triage_info(info) -> Tuple:
    return tuple((k, _triage_num(v)) for (k, v) in info)


def _triage_samples(samples) -> Tuple:
    # same numeric (float32) narrowing on sample subfield values; GT strings pass
    # through unchanged (their equality rule is literal string equality).
    return tuple(
        tuple((k, _triage_num(v)) for (k, v) in sample) for sample in samples
    )


def _triage_key(records) -> Tuple:
    # Deeper representation collapse: REF/ALT (case + whitespace) and INFO
    # numerics narrowed to float32 (the approved f32-vs-f64 rule). A content
    # difference that vanishes here was representation-only -> quarantined.
    return tuple(
        (c, p, (r or "").strip().upper(), _triage_alleles(a),
         _triage_info(info), _triage_samples(samples))
        for (c, p, r, a, info, samples) in records
    )


def _triage_groups(ok: List[ParseSummary]) -> Dict[Tuple, List[str]]:
    groups: Dict[Tuple, List[str]] = {}
    for s in ok:
        groups.setdefault(_triage_key(s.records), []).append(s.tool)
    return groups


# ---- classification -------------------------------------------------------

def classify(summaries: List[ParseSummary]) -> Tuple[str, str, bool]:
    """Return (bucket, detail, all_parsed)."""
    if any(s.kind == KIND_PLUMBING for s in summaries):
        bad = [s for s in summaries if s.kind == KIND_PLUMBING]
        return (
            BUCKET_HARNESS_ERROR,
            "harness failure -- " + "; ".join(f"{short(s.tool)}: {s.error}" for s in bad),
            False,
        )

    ok = [s for s in summaries if s.ok]
    fail = [s for s in summaries if not s.ok]
    all_parsed = not fail

    # Content axis first: needs >= 2 parsers that succeeded.
    if len(ok) >= 2:
        cgroups = _content_groups(ok)
        if len(cgroups) > 1:
            triaged = _triage_groups(ok)
            artifact = len(triaged) == 1  # collapses under deeper canon
            bucket = BUCKET_ARTIFACT if artifact else BUCKET_CONTENT
            return bucket, _describe_content(cgroups, fail, artifact), all_parsed

    # No content disagreement.
    if not ok:
        return BUCKET_ALL_FAIL, _describe_failures(fail), False
    if not fail:
        return BUCKET_ALL_AGREE, f"all {len(ok)} parsed and agree ({ok[0].num_records} records)", True
    return BUCKET_SPLIT, _describe_split(ok, fail), False


def _describe_content(cgroups, fail, artifact: bool) -> str:
    parts = sorted(cgroups.values(), key=len, reverse=True)
    split = " vs ".join("+".join(short(t) for t in p) for p in parts)
    tag = "representation-only (collapses under triage)" if artifact else "REAL content disagreement"
    scope = "all parsers parsed" if not fail else f"partial; failed: {', '.join(short(s.tool) for s in fail)}"
    return f"{tag} [{scope}]: {split}"


def _describe_failures(fail: List[ParseSummary]) -> str:
    parts = [f"{short(s.tool)}[{s.kind}]: {s.error}" for s in fail]
    hard = [short(s.tool) for s in fail if s.kind in (KIND_CRASH, KIND_TIMEOUT)]
    note = f"  NOTE: {', '.join(hard)} crashed/hung rather than cleanly rejecting." if hard else ""
    return "all rejected. " + " | ".join(parts) + note


def _describe_split(ok: List[ParseSummary], fail: List[ParseSummary]) -> str:
    ok_names = ", ".join(short(s.tool) for s in ok)
    fail_desc = ", ".join(f"{short(s.tool)}[{s.kind}]" for s in fail)
    base = f"parsed: {ok_names}; rejected/crashed: {fail_desc}."
    if len(ok) == 2 and len(fail) == 1:
        base += f" 2-vs-1: {short(fail[0].tool)} is the outlier."
    elif len(ok) == 1 and len(fail) == 2:
        base += f" 1-vs-2: {short(ok[0].tool)} accepted alone."
    return base


# ---- signature (for minimization oracle) ----------------------------------

_ABSENT = object()


def _typename(v) -> str:
    if v is _ABSENT:
        return "absent"
    if v is None:
        return "none"
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, float):
        return "float"
    if isinstance(v, int):
        return "int"
    if isinstance(v, str):
        return "str"
    if isinstance(v, tuple):
        return "tuple"
    return type(v).__name__


def _diff_descriptor(ok: List[ParseSummary]) -> frozenset:
    """WHICH fields differ across the parsers, and the value-TYPES involved.
    Computed on TRIAGE-normalized records, so representation-only differences
    (float32-vs-float64 etc.) do NOT pollute it -- only the real, triage-
    surviving differences appear. This separates distinct kinds of disagreement
    that share a tool partition -- 'INFO:DP as float-vs-int' (a coercion) is a
    different finding from 'INFO:AF as none-vs-str' (a drop)."""
    tgroups = {}
    for s in ok:
        tgroups.setdefault(_triage_key(s.records), []).append(s)
    if len(tgroups) <= 1:
        return frozenset()
    reps = list(tgroups.keys())  # each is a tuple of triage-normalized records
    descr = set()
    if len({len(r) for r in reps}) > 1:
        descr.add(("NUM_RECORDS", tuple(sorted(len(r) for r in reps))))
    for ri in range(min(len(r) for r in reps)):
        recs = [r[ri] for r in reps]  # each = (chrom, pos, ref, alts, info, samples)
        for fi, lab in enumerate(("CHROM", "POS", "REF", "ALT")):
            vals = [rec[fi] for rec in recs]
            if len(set(vals)) > 1:
                descr.add((lab, tuple(sorted(_typename(v) for v in vals))))
        infos = [dict(rec[4]) for rec in recs]
        for k in set().union(*[set(i) for i in infos]):
            vals = [i.get(k, _ABSENT) for i in infos]
            if len(set(vals)) > 1:
                descr.add((f"INFO:{k}", tuple(sorted(_typename(v) for v in vals))))
        # per-sample FORMAT fields
        samp = [rec[5] for rec in recs]
        if len({len(sl) for sl in samp}) > 1:
            descr.add(("NUM_SAMPLES", tuple(sorted(len(sl) for sl in samp))))
        for si in range(max((len(sl) for sl in samp), default=0)):
            sdicts = [dict(sl[si]) if si < len(sl) else {} for sl in samp]
            for k in set().union(*[set(sd) for sd in sdicts]):
                vals = [sd.get(k, _ABSENT) for sd in sdicts]
                if len(set(vals)) > 1:
                    descr.add((f"FMT[{si}]:{k}", tuple(sorted(_typename(v) for v in vals))))
    return frozenset(descr)


def signature(comp: "Comparison") -> Tuple:
    """A structural fingerprint of the outcome: bucket + how the parsers
    partition by TRIAGE content (representation-only diffs collapsed) + which
    parsers failed + WHAT differs (fields + value types). Minimization keeps this
    identical, and the reporter groups genuinely-distinct findings apart."""
    ok = [s for s in comp.summaries if s.ok]
    tgroups = {}
    for s in ok:
        tgroups.setdefault(_triage_key(s.records), []).append(short(s.tool))
    partition = frozenset(frozenset(tools) for tools in tgroups.values())
    failed = frozenset(short(s.tool) for s in comp.summaries if not s.ok)
    return (comp.bucket, partition, failed, _diff_descriptor(ok))


# ---- compact reporting ----------------------------------------------------

def sides(comp: "Comparison") -> str:
    ss = comp.summaries
    b = comp.bucket
    if b == BUCKET_ALL_AGREE:
        return "all agree"
    if b == BUCKET_HARNESS_ERROR:
        return "HARNESS ERROR"
    if b == BUCKET_ALL_FAIL:
        notes = [f"{short(s.tool)}:{s.kind}!" for s in ss if s.kind in (KIND_CRASH, KIND_TIMEOUT)]
        return "all reject" + (" (" + ",".join(notes) + ")" if notes else "")
    if b in (BUCKET_CONTENT, BUCKET_ARTIFACT):
        parts = sorted(_content_groups([s for s in ss if s.ok]).values(), key=len, reverse=True)
        g = " != ".join("+".join(short(t) for t in p) for p in parts)
        failed = [short(s.tool) for s in ss if not s.ok]
        return g + (f" (+{','.join(failed)} failed)" if failed else "")
    if b == BUCKET_SPLIT:
        ok = "+".join(short(s.tool) for s in ss if s.ok)
        bad = "+".join(f"{short(s.tool)}:{s.kind}" for s in ss if not s.ok)
        return f"[{ok}]parsed vs [{bad}]"
    return "?"


def marker(comp: "Comparison") -> str:
    return {
        BUCKET_CONTENT: "###",
        BUCKET_ARTIFACT: "~~~",
        BUCKET_SPLIT: "  *",
        BUCKET_HARNESS_ERROR: "!!!",
    }.get(comp.bucket, "   ")


def render_text_report(comp: "Comparison", label: str) -> str:
    lines = [
        f"CATCH: {label}",
        f"bucket: {comp.bucket}" + ("  (SILENT: every parser succeeded)" if comp.all_parsed and comp.bucket == BUCKET_CONTENT else ""),
        f"detail: {comp.detail}",
        f"input : {os.path.basename(comp.vcf_path)}",
        "",
    ]
    for s in comp.summaries:
        lines.append(f"[{short(s.tool)}] {s.tool}  kind={s.kind}")
        if s.ok:
            lines.append(f"    {s.num_records} records:")
            for (c, p, r, a, info, samples) in s.records:
                info_str = f"  INFO={dict(info)}" if info else ""
                samp_str = f"  FMT={[dict(x) for x in samples]}" if samples else ""
                lines.append(f"      {c}:{p}  REF={r!r}  ALT={list(a)!r}{info_str}{samp_str}")
        else:
            lines.append(f"    error: {s.error}")
        for note in getattr(s, "notes", []) or []:
            lines.append(f"    NOTE: {note}")
        lines.append("")
    return "\n".join(lines)


# ---- entry points ---------------------------------------------------------

def compare_file(vcf_path: str) -> Comparison:
    summaries = run_all(vcf_path)
    bucket, detail, all_parsed = classify(summaries)
    return Comparison(os.path.abspath(vcf_path), bucket, detail, summaries, all_parsed)


def save_finding(comp: Comparison, label: str, findings_dir: str | None = None,
                 extra: dict | None = None) -> str:
    findings_dir = findings_dir or os.path.join(REPO_ROOT, "findings")
    dest = os.path.join(findings_dir, f"{comp.bucket}__{label}")
    os.makedirs(dest, exist_ok=True)
    shutil.copyfile(comp.vcf_path, os.path.join(dest, "input.vcf"))
    report = {
        "label": label,
        "bucket": comp.bucket,
        "all_parsed": comp.all_parsed,
        "detail": comp.detail,
        "input_vcf": os.path.basename(comp.vcf_path),
        "parsers": [
            {"tool": s.tool, "kind": s.kind, "num_records": s.num_records,
             "records": s.records, "error": s.error}
            for s in comp.summaries
        ],
    }
    if extra:
        report.update(extra)
    with open(os.path.join(dest, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(dest, "report.txt"), "w") as f:
        f.write(render_text_report(comp, label))
    return dest


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("usage: python -m compare.compare <vcf_path>")
        raise SystemExit(2)
    c = compare_file(sys.argv[1])
    print(f"bucket : {c.bucket}")
    print(f"detail : {c.detail}")
    print(f"split  : {sides(c)}")
    print()
    print(render_text_report(c, os.path.basename(sys.argv[1])))
