"""Adapter for pysam -- the representative of the htslib parser family.

(bcftools, cyvcf2 and samtools share this same htslib backend, so we keep only
ONE of them in the comparison.)
"""
from __future__ import annotations

from adapters.common import Record, run_worker

TOOL = "pysam (htslib)"


def _info(rec, declared):
    # Only compare INFO fields DECLARED in this file's header. Undeclared fields
    # are excluded because the parsers legitimately guess their type/structure
    # differently -- comparing those guesses is not a fair content test.
    # htslib types declared INFO: Flag -> True, Number>1 -> tuple, else scalar.
    # Floats come back as float32 widened to Python float; kept exact.
    out = {}
    for k, v in rec.info.items():
        if k not in declared:
            continue
        if v is True:
            out[k] = True
        elif isinstance(v, tuple):
            out[k] = list(v)
        else:
            out[k] = v
    return out


def _gt_string(sample):
    # Reconstruct the VCF GT string LITERALLY from pysam's index tuple + phased
    # flag: 'a|b' if phased else 'a/b', '.' for a missing allele. No sorting, no
    # interpretation -- equality decisions live in the central canonicalizer.
    gt = sample.get("GT")
    if gt is None:
        return None
    sep = "|" if sample.phased else "/"
    return sep.join("." if a is None else str(a) for a in gt)


def _samples(rec, ri, notes):
    # The GT we COMPARE is what pysam's high-level .samples['GT'] API returns --
    # that is how real pysam users read genotypes. But that API reinterprets
    # invalid allele indices (e.g. GT=1 at a no-ALT site) as missing '.'. To make
    # any such divergence self-demonstrating, we ALSO read the RAW genotype from
    # pysam's own serialization (str(rec), which preserves the file value) and, if
    # it differs from the .samples value, emit a display-only note. The note is
    # never compared -- it just proves the divergence enters at the reader API.
    raw_cols = str(rec).rstrip("\n").split("\t")
    out = []
    for si, name in enumerate(rec.samples):
        s = rec.samples[name]
        d = {}
        for k, v in s.items():
            if k == "GT":
                continue  # reconstructed as a string below
            d[k] = list(v) if isinstance(v, tuple) else v
        gt = _gt_string(s)
        if gt is not None:
            d["GT"] = gt
            col = 9 + si
            if col < len(raw_cols):
                raw_gt = raw_cols[col].split(":")[0]
                if raw_gt != gt:
                    notes.append(
                        f"{rec.chrom}:{rec.pos} sample[{si}] GT: file says {raw_gt!r} "
                        f"but pysam .samples API returns {gt!r} "
                        f"(vcfpy/noodles report the file value)"
                    )
        out.append(d)
    return out


def parse(path: str):
    import pysam

    records: list[Record] = []
    notes: list[str] = []
    with pysam.VariantFile(path) as vf:
        declared = set(vf.header.info)
        for ri, rec in enumerate(vf):
            # rec.alts is a tuple of allele strings, or None when ALT is '.'.
            alts = tuple(rec.alts) if rec.alts is not None else ()
            # rec.pos is the 1-based VCF POS (rec.start would be 0-based).
            records.append((rec.chrom, rec.pos, rec.ref, alts,
                            _info(rec, declared), _samples(rec, ri, notes)))
    return records, notes


if __name__ == "__main__":
    run_worker(TOOL, parse)
