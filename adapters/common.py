"""Shared contract every VCF parser adapter must speak.

The whole project depends on two adapters looking at the SAME file producing
directly comparable output. That comparability lives here.

Canonical per-record key:  (chrom, pos, ref, alts)
    chrom : str            -- contig name, verbatim
    pos   : int            -- 1-based VCF POS (NOT any tool's 0-based coord)
    ref   : str            -- reference allele, verbatim
    alts  : tuple[str,...] -- alternate alleles as VCF strings;
                              the EMPTY tuple means "no ALT" (the '.' record)

Design principle: we keep only the *raw content* of each site. We deliberately
exclude every tool's *interpretation* -- allele-type classification
(SNV/INDEL/...), 0-based coordinates, parsed genotypes -- because those are
opinions, not facts about the file, and two tools with different opinions but
identical content would look like a fake disagreement.

Failure is a first-class outcome. A parser fed a broken file may:
  - raise a catchable exception  -> KIND_EXCEPTION  (it rejected the input)
  - crash the process (segfault) -> KIND_CRASH      (detected by the runner)
  - hang forever                 -> KIND_TIMEOUT    (detected by the runner)
Only the first can be caught inside the worker; the other two are caught by the
parent because they kill or freeze the worker process.
"""
from __future__ import annotations

import dataclasses
import json
import math
import sys
from typing import List, Optional, Tuple

# A single normalized variant record.
Record = Tuple[str, int, str, Tuple[str, ...]]

# The VCF spec's "missing value" marker.
MISSING_ALLELE = "."


def canonicalize_alts(alts: Tuple[str, ...]) -> Tuple[str, ...]:
    """Collapse representation-only differences in how parsers render a
    missing/absent alternate allele, so the comparator does not flag them as
    content disagreements.

    Two equivalences, both validated by round-tripping through each tool's own
    writer (pysam and noodles emit identical bytes for these once collapsed):
      1. an empty-string allele ''  ==  the VCF missing marker '.'
         (e.g. noodles yields 'A,'->['A',''] where pysam yields ['A','.']).
      2. an ALT that is ENTIRELY missing (all '.' / '' / empty list) means
         "no alternate alleles" -> ()  (e.g. an empty ALT field: pysam holds
         ['.'], noodles holds [], but both WRITE '.').

    This only ever collapses missing-vs-missing / missing-vs-absent. A missing
    allele is never merged with a real one, so genuine disagreements survive
    (e.g. ('.',) vs ('G',) still differ: () vs ('G',)).
    """
    norm = tuple(MISSING_ALLELE if a == "" else a for a in alts)
    if all(a == MISSING_ALLELE for a in norm):  # all() on () is True -> ()
        return ()
    return norm


def canonicalize_ref(ref: str) -> str:
    """Same missing-token equivalence for the REF field. An empty REF ('') and
    the '.' marker are two renderings of an absent/missing reference (pysam
    emits '.', noodles emits '' -- validated by round-trip). Collapse them so
    that is not mistaken for a content disagreement. Never touches a real base,
    so ('.' vs 'G') style genuine differences survive.
    """
    return MISSING_ALLELE if ref == "" else ref


_NAN = "__nan__"  # canonical NaN marker (nan != nan breaks equality/hashing)


def _canon_info_value(v):
    # Container normalization (list -> tuple) + two approved equalities:
    #  - MISSING value: '', None, an empty array, or an all-missing array are
    #    all renderings of "no value" -> canonical None. Only collapses
    #    missing-vs-missing; a real value vs missing still differs.
    #  - NaN: mapped to a marker so nan == nan (Python's nan != nan otherwise).
    # Numbers are otherwise EXACT (no rounding) -- real value differences survive.
    if isinstance(v, list):
        t = tuple(_canon_info_value(x) for x in v)
        if len(t) == 0 or all(x is None for x in t):
            return None
        return t
    if v == "":
        return None
    if isinstance(v, float) and math.isnan(v):
        return _NAN
    return v


def canonicalize_info(info) -> tuple:
    """Canonical form of a record's INFO: an UNORDERED mapping key -> value,
    represented as a sorted tuple of (key, value) pairs (hashable, order-
    independent -- INFO is semantically unordered). Values are compared exactly
    after container normalization; nothing numeric is rounded here.
    """
    if not info:
        return ()
    return tuple(sorted(
        ((k, _canon_info_value(v)) for k, v in info.items()),
        key=lambda kv: kv[0],
    ))


def canonicalize_samples(samples) -> tuple:
    """Canonical form of the per-sample FORMAT data: an ORDERED tuple (samples
    are positional) of per-sample mappings (unordered key -> value).

    Per the approved genotype rule:
      - GT is a STRING ('a|b' phased / 'a/b' unphased). We keep it verbatim, so
        string equality implements "do NOT sort unphased" ('1/0' != '0/1') and
        "ploidy literal" ('1' != '1|1') -- no interpretation here.
      - a present-but-ALL-MISSING subfield (e.g. HQ=.,.) is dropped, so it equals
        an absent subfield. (Only collapses missing-vs-missing; a real value vs
        missing still differs.)
      - other subfield values reuse the INFO value canon (container norm, missing
        collapse, NaN marker; numerics exact).

    The drop is ROUND-TRIP VALIDATED (watch-point 1): parsers agree on an explicit
    present-but-missing subfield (all give HQ=[None,None]); they differ internally
    only on an *omitted* trailing subfield (pysam pads HQ=[None], vcfpy/noodles
    leave it absent) -- but pysam and vcfpy both re-emit the identical bytes
    (`0/0:61:2:.`), so that difference is representation, not content. Collapsing
    it is therefore correct, not hiding a finding.
    """
    if not samples:
        return ()
    out = []
    for sample in samples:
        items = []
        for k, v in sample.items():
            cv = _canon_info_value(v)
            if cv is None:  # present-but-all-missing subfield == absent -> drop
                continue
            items.append((k, cv))
        out.append(tuple(sorted(items, key=lambda kv: kv[0])))
    return tuple(out)


# Outcome classes.
KIND_OK = "ok"                # parser read the file and produced records
KIND_EXCEPTION = "exception"  # parser raised (a clean rejection of the input)
KIND_CRASH = "crash"          # worker process died (signal / abort / segfault)
KIND_TIMEOUT = "timeout"      # parser hung and was killed by the runner
KIND_PLUMBING = "plumbing"    # OUR harness failed -- not a verdict about a parser


@dataclasses.dataclass
class ParseSummary:
    """What one parser thinks one file contains (or how it failed)."""

    tool: str
    kind: str
    num_records: int = 0
    records: List[Record] = dataclasses.field(default_factory=list)
    error: Optional[str] = None
    # display-only annotations (NOT part of the compared content); used e.g. to
    # surface pysam's raw-record GT alongside its reinterpreted .samples GT.
    notes: List[str] = dataclasses.field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.kind == KIND_OK

    def to_json(self) -> str:
        return json.dumps(dataclasses.asdict(self))

    @classmethod
    def from_json(cls, text: str) -> "ParseSummary":
        d = json.loads(text)
        # Single choke point for ALL adapter output (pysam/vcfpy/noodles all
        # arrive here as JSON): rebuild tuples AND apply the shared missing-token
        # canonicalization, so every tool is normalized identically.
        d["records"] = [
            (r[0], r[1], canonicalize_ref(r[2]), canonicalize_alts(tuple(r[3])),
             canonicalize_info(r[4] if len(r) > 4 else {}),
             canonicalize_samples(r[5] if len(r) > 5 else []))
            for r in d.get("records", [])
        ]
        return cls(**d)


def run_worker(tool: str, parse_fn) -> None:
    """Subprocess entry point for an adapter module.

    Reads a VCF path from argv, runs parse_fn(path) -> list[Record], and writes
    a ParseSummary as JSON to stdout. Any exception from the parser is caught
    and reported as KIND_EXCEPTION. Crashes and hangs cannot be caught here --
    they are the runner's job -- and that separation is intentional.
    """
    if len(sys.argv) != 2:
        sys.stderr.write("usage: python -m adapters.<name> <vcf_path>\n")
        raise SystemExit(2)
    path = sys.argv[1]
    try:
        result = parse_fn(path)
        records, notes = result if isinstance(result, tuple) else (result, [])
        records = list(records)
        summary = ParseSummary(
            tool=tool, kind=KIND_OK, num_records=len(records), records=records,
            notes=list(notes),
        )
    except Exception as e:  # noqa: BLE001 -- any parser error is data we want
        summary = ParseSummary(
            tool=tool, kind=KIND_EXCEPTION, error=f"{type(e).__name__}: {e}"
        )
    sys.stdout.write(summary.to_json())
    sys.stdout.flush()
