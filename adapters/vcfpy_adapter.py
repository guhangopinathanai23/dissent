"""Adapter for vcfpy -- an independent, pure-Python VCF parser.

vcfpy shares no code with htslib, which is exactly why it is a meaningful
cross-check against pysam.
"""
from __future__ import annotations

from adapters.common import Record, run_worker

TOOL = "vcfpy (pure-python)"


def _alt_to_str(alt) -> str:
    # vcfpy wraps each ALT allele in an AltRecord object (Substitution,
    # SymbolicAllele, BreakEnd, SingleBreakEnd). .serialize() returns the exact
    # VCF string form ('A', '<DEL>', 'G]17:198982]'), which is what pysam
    # reports as a plain string -- so serializing keeps the two comparable and
    # discards vcfpy's SNV/INDEL type-classification, which we don't want.
    return alt.serialize()


def _info(rec, declared):
    # Only DECLARED INFO fields (see pysam_adapter for the rationale). vcfpy
    # types declared INFO per the header: Flag -> True, Number>1 -> list, else
    # scalar. Floats are Python float64 -- the f32-vs-f64 split with the other
    # two parsers is exactly what the triage gate's float32 rule quarantines.
    return {k: v for k, v in rec.INFO.items() if k in declared}


def _samples(rec):
    # vcfpy gives each call's FORMAT data as an OrderedDict; GT is already the
    # VCF string ('0|1'). Emit verbatim -- no interpretation.
    return [dict(call.data) for call in rec.calls]


def parse(path: str):
    import vcfpy

    records: list[Record] = []
    reader = vcfpy.Reader.from_path(path)
    try:
        declared = set(reader.header.info_ids())
        for rec in reader:
            # rec.ALT is [] when there is no ALT ('.') -> normalizes to ().
            alts = tuple(_alt_to_str(a) for a in rec.ALT)
            records.append((rec.CHROM, rec.POS, rec.REF, alts,
                            _info(rec, declared), _samples(rec)))
    finally:
        reader.close()
    return records, []


if __name__ == "__main__":
    run_worker(TOOL, parse)
