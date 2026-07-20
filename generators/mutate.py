"""Small, deliberate, NAMED mutations of the seed VCF (Phase 2).

Each mutation is ONE understandable transformation applied to the valid seed,
so any catch is traceable straight back to its cause. These are not random
bytes -- every entry has a name that says exactly what it does. Kept small on
purpose; scaling to a large corpus is Phase 3.

Field indices in a VCF data line (tab-separated):
    0 CHROM  1 POS  2 ID  3 REF  4 ALT  5 QUAL  6 FILTER  7 INFO  8 FORMAT  9+ samples
Most POS/REF/ALT mutations target the FIRST data record (20:14370 G->A).
"""
from __future__ import annotations

import os

SEED_DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "seed.vcf")


# ---- line helpers ---------------------------------------------------------

def _lines(text: str) -> list[str]:
    return text.split("\n")


def _data_indices(lines: list[str]) -> list[int]:
    return [i for i, ln in enumerate(lines) if ln and not ln.startswith("#")]


def _fields(line: str) -> list[str]:
    return line.split("\t")


def _join(fields: list[str]) -> str:
    return "\t".join(fields)


def _edit_field(text: str, rec: int, field: int, value: str) -> str:
    """Set field `field` of the `rec`-th data record to `value`."""
    lines = _lines(text)
    i = _data_indices(lines)[rec]
    f = _fields(lines[i])
    f[field] = value
    lines[i] = _join(f)
    return "\n".join(lines)


# ---- mutations: POS (field 1, first record = 14370) -----------------------

def m_control_valid(t):          return t                                  # sanity control
def m_pos_shift_plus1(t):        return _edit_field(t, 0, 1, "14371")
def m_pos_shift_minus1(t):       return _edit_field(t, 0, 1, "14369")
def m_pos_zero(t):               return _edit_field(t, 0, 1, "0")          # 1-based -> 0 is edge
def m_pos_negative(t):           return _edit_field(t, 0, 1, "-1")
def m_pos_non_integer(t):        return _edit_field(t, 0, 1, "12ABC")
def m_pos_float(t):              return _edit_field(t, 0, 1, "14370.5")
def m_pos_huge(t):               return _edit_field(t, 0, 1, "999999999999")  # > contig length
def m_pos_leading_space(t):      return _edit_field(t, 0, 1, " 14370")     # stray whitespace


# ---- mutations: REF / ALT (fields 3, 4; first record G->A) ----------------

def m_ref_lowercase(t):          return _edit_field(t, 0, 3, "g")
def m_alt_lowercase(t):          return _edit_field(t, 0, 4, "a")
def m_alt_empty(t):              return _edit_field(t, 0, 4, "")           # empty field, not '.'
def m_alt_to_dot(t):             return _edit_field(t, 0, 4, ".")          # monomorphic (valid)
def m_alt_trailing_comma(t):     return _edit_field(t, 0, 4, "A,")         # empty trailing allele
def m_alt_leading_comma(t):      return _edit_field(t, 0, 4, ",A")         # empty leading allele


def m_split_multiallelic(t):
    """Split the multi-allelic record (A -> G,T) into two single-ALT lines."""
    lines = _lines(t)
    for i in _data_indices(lines):
        f = _fields(lines[i])
        if "," in f[4]:
            expanded = []
            for allele in f[4].split(","):
                g = f[:]
                g[4] = allele
                expanded.append(_join(g))
            lines[i:i + 1] = expanded
            break
    return "\n".join(lines)


# ---- mutations: structural ------------------------------------------------

def m_delete_col_qual(t):
    """Remove the QUAL column (field 5) from the first record -> column shift."""
    lines = _lines(t)
    i = _data_indices(lines)[0]
    f = _fields(lines[i])
    del f[5]
    lines[i] = _join(f)
    return "\n".join(lines)


def m_add_stray_tab(t):
    """Insert an extra tab after CHROM in the first record -> empty column."""
    lines = _lines(t)
    i = _data_indices(lines)[0]
    lines[i] = lines[i].replace("\t", "\t\t", 1)
    return "\n".join(lines)


def m_duplicate_record(t):
    """Duplicate the first data record line."""
    lines = _lines(t)
    i = _data_indices(lines)[0]
    lines.insert(i + 1, lines[i])
    return "\n".join(lines)


def m_extra_blank_line(t):
    """Insert a blank line between the first and second records."""
    lines = _lines(t)
    i = _data_indices(lines)[0]
    lines.insert(i + 1, "")
    return "\n".join(lines)


def m_truncate_last_line(t):
    """Cut the final data record in half and drop the trailing newline."""
    lines = _lines(t)
    i = _data_indices(lines)[-1]
    lines[i] = lines[i][: max(1, len(lines[i]) // 2)]
    return "\n".join(lines[: i + 1])


def m_remove_final_newline(t):
    return t[:-1] if t.endswith("\n") else t


# ---- mutations: header ----------------------------------------------------

def m_drop_fileformat(t):
    return "\n".join(l for l in _lines(t) if not l.startswith("##fileformat"))


def m_drop_chrom_header(t):
    return "\n".join(l for l in _lines(t) if not l.startswith("#CHROM"))


# ---- mutations: line endings ----------------------------------------------

def m_crlf(t):      return t.replace("\n", "\r\n")   # Windows
def m_cr_only(t):   return t.replace("\n", "\r")     # classic Mac


# ---- registry -------------------------------------------------------------

MUTATIONS = {
    "control_valid": m_control_valid,
    "pos_shift_plus1": m_pos_shift_plus1,
    "pos_shift_minus1": m_pos_shift_minus1,
    "pos_zero": m_pos_zero,
    "pos_negative": m_pos_negative,
    "pos_non_integer": m_pos_non_integer,
    "pos_float": m_pos_float,
    "pos_huge": m_pos_huge,
    "pos_leading_space": m_pos_leading_space,
    "ref_lowercase": m_ref_lowercase,
    "alt_lowercase": m_alt_lowercase,
    "alt_empty": m_alt_empty,
    "alt_to_dot": m_alt_to_dot,
    "alt_trailing_comma": m_alt_trailing_comma,
    "alt_leading_comma": m_alt_leading_comma,
    "split_multiallelic": m_split_multiallelic,
    "delete_col_qual": m_delete_col_qual,
    "add_stray_tab": m_add_stray_tab,
    "duplicate_record": m_duplicate_record,
    "extra_blank_line": m_extra_blank_line,
    "truncate_last_line": m_truncate_last_line,
    "remove_final_newline": m_remove_final_newline,
    "drop_fileformat": m_drop_fileformat,
    "drop_chrom_header": m_drop_chrom_header,
    "crlf": m_crlf,
    "cr_only": m_cr_only,
}


def generate_all(seed_path: str = SEED_DEFAULT, out_dir: str = "data/mutations") -> list[tuple[str, str]]:
    """Write one file per mutation. Returns [(name, path), ...]."""
    with open(seed_path) as f:
        seed_text = f.read()
    os.makedirs(out_dir, exist_ok=True)
    pairs = []
    for name, fn in MUTATIONS.items():
        mutated = fn(seed_text)
        path = os.path.join(out_dir, f"{name}.vcf")
        # newline="" so \r\n / \r survive exactly as written (no translation).
        with open(path, "w", newline="") as f:
            f.write(mutated)
        pairs.append((name, path))
    return pairs


if __name__ == "__main__":
    made = generate_all()
    print(f"wrote {len(made)} mutation files to data/mutations/")
    for name, path in made:
        print(f"  {name:<22} -> {path}")
