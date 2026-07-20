# dissent — a differential fuzzer for VCF parsers

*Complete write-up. Covers the shallow fields, typed INFO, and per-sample
genotypes; indels are a deliberate deferral (see Future work), not a gap.*

**Three production VCF parsers silently disagree on typed field content — INFO
values, and how one of them hands genotypes to its callers — on ordinary-looking
inputs. This tool finds and minimizes the exact files where they do.**

## What the tool does

Feed the same (often malformed) VCF file to three **independently-built** VCF
parsers, then catch where they silently disagree or crash. A disagreement proves
at least one parser has a bug **without needing to know the "correct" answer in
advance** — that is the point of differential testing.

The three voices are genuinely separate codebases (this was verified, not
assumed — many popular VCF tools are all htslib underneath and would just echo
each other):

| Parser | Language | Backend | Independent? |
|--------|----------|---------|--------------|
| **pysam** | Python/Cython | htslib (C) | representative of the htslib family |
| **vcfpy** | pure Python | own | yes |
| **noodles** | Rust | own — pure Rust, no C bindings (verified via `cargo tree`) | yes |

Pipeline: mutate a valid seed → run each parser in an **isolated subprocess**
(so a crash or hang is caught, not fatal) → reduce each to a normalized summary
→ classify (with a triage gate that quarantines representation-only differences)
→ minimize any real catch to a tiny file.

## The result — the deep-field hypothesis, confirmed and closed

An early version of this project could only report a **shallow negative**: on the
four simplest fields the parsers agreed, and the real bugs — if any — had to live
in the typed fields we hadn't yet compared. We then compared them. The complete
arc:

- **Shallow fields `(chrom, pos, ref, alt)` → the parsers agree.** Zero silent
  content disagreements, zero crashes across a 3,000-file corpus. The only
  divergence is *accept/reject* — what each parser is willing to accept.
- **Typed `INFO` content → they diverge silently.** A header-preserving corpus
  surfaced genuine silent content disagreements, each round-trip-confirmed real
  (the parsers re-emit different bytes).
- **Per-sample genotypes → one finding the round-trip test cannot even see.** A
  divergence on pysam's *reader* API that is invisible to any comparison of
  output files.

These parsers agree on the rigid, simple fields and **diverge exactly where the
format asks them to interpret** — typed INFO values and genotype semantics. Every
deep field the format is known to be fragile in has now been looked at.

## Headline — a silent `INFO` content disagreement

**Given `AA=T,99` in a `Number=1` (single-value) String field, pysam splits it
into two values; vcfpy and noodles keep it as one.** All three accept the file —
nothing errors — yet they report different content. This is the only case where
*all three parse and diverge*: a fully-silent disagreement with no safety net.

Minimized (`findings/content_disagreement__i00384.vcf_min`):

```
##fileformat=VCFv4.2
##INFO=<ID=AA,Number=1,Type=String,...>
#CHROM  POS      ID  REF  ALT  ...  INFO
20      1230237  .   T    .    ...  NS=3;DP=13;AA=T,99
```

| parser | parsed `AA` | re-emits |
|--------|-------------|----------|
| **pysam** | `('T', '99')` — **two values** | `AA=T,99` |
| vcfpy | `'T,99'` — one value | `AA=T%2C99` (comma escaped) |
| noodles | `'T,99'` — one value | `AA=T,99` |

The round-trip is the proof this is **content, not representation**: pysam and
vcfpy re-emit **different bytes** (`AA=T,99` vs `AA=T%2C99`). They genuinely
disagree — pysam reads two values; vcfpy reads one value containing a literal
comma, which the spec requires escaping as `%2C`. A downstream tool asking "how
many values does `AA` have?" gets **2** via pysam and **1** via the others,
silently, with no error anywhere.

### Supporting — pysam silently coerces malformed `INFO` values

A family of pysam silent coercions on declared fields. Here noodles **rejects**
the malformed value and vcfpy keeps the raw token, so pysam is the outlier: it
alone accepts the value and reshapes it into plausible-but-wrong content. Each row
is confirmed against the saved findings and round-trips to a value different from
what was read (pysam's re-emit shown):

| input (declared field) | pysam (re-emits) | vcfpy | noodles |
|---|---|---|---|
| `NS=1.5` (Integer) | **`1`** (`NS=1`) — silently truncates | `'1.5'` | rejects |
| `DP=abc` (Integer) | **missing** (`DP=.`) — silently drops | `'abc'` | rejects |
| `NS=3,99` (Number=1) | **`(3,99)`** (`NS=3,99`) — keeps both | `'3,99'` | rejects |
| `DB=0` (Flag) | **`('0',)`** (`DB=0`) — keeps value | flag `True` | rejects |

The pattern: **pysam (htslib) is the most permissive — it silently reshapes
malformed INFO into plausible wrong content; noodles is strictest (rejects); vcfpy
keeps the raw token.** `NS=1.5 → 1` is the most dangerous single value, but because
noodles rejects it, it is not a fully-silent three-way disagreement the way the
comma-split is.

> **The discipline that makes these trustworthy.** Two *earlier* apparent content
> disagreements (empty `ALT`/`REF`: pysam `.`, noodles `''`) turned out to be
> **representation artifacts** — round-tripping through each parser's writer
> produced identical bytes. A **triage gate** quarantines that whole class (and
> absorbs the float32-vs-float64 numeric difference), so it never masquerades as a
> finding. Every disagreement above survived that gate and was round-trip-checked.

## Genotypes — a serialization-invisible finding

Extending the comparison into per-sample genotypes surfaced exactly **one** new
finding — and it is one the round-trip test **cannot catch**, which is precisely
what makes it notable.

**pysam's `.samples['GT']` reader API reinterprets an invalid allele index as a
missing call, where vcfpy and noodles report the file's value.** For `GT=1` at a
site with no ALT allele (allele index 1 does not exist), pysam's `.samples['GT']`
— the standard way pysam users read genotypes — returns a missing call `.`; vcfpy
and noodles return `1`.

**Why the round-trip test is blind to it — and why that matters.** Every other
finding in this document is proven real by "the parsers re-emit different bytes."
This one re-emits **identical** bytes: pysam's *writer* faithfully preserves the
raw genotype `1`, so **no comparison of output files could ever detect it**. The
divergence lives entirely on the *reader* surface — what pysam hands a program
that calls `.samples['GT']` — which the writer never touches. The round-trip test
worked for comma-split and the missing-subfield case because, there, "what the
parser serializes" was a faithful proxy for "what the parser hands its consumer."
Here those two surfaces have come apart, and the round-trip only sees the writer.
The question that actually matters — *do downstream consumers diverge?* — is
unambiguously yes. This is a signature example of what this tool finds that
nothing else can: a divergence invisible to any file diff.

**The proof, instead, is a three-way display.** The fuzzer surfaces pysam's own
raw-record value next to its `.samples` value:

```
raw file value        : 1
pysam .samples['GT']  : .      <- reinterpreted (invalid index -> missing)
vcfpy / noodles       : 1
```

pysam's own raw value agreeing with vcfpy/noodles proves the divergence enters at
pysam's `.samples` API reinterpretation, **not at parse time**. Confirmed specific
to invalid indices: a valid genotype (`GT=1/1` at a biallelic site) parses
identically in all three, with no note.

**Scope, precisely.** This is a finding about pysam's `.samples['GT']` API
reinterpreting invalid allele indices as missing — *not* a claim that pysam parses
the genotype wrong (its writer is faithful). The honest statement:

> A program reading genotypes via pysam's standard `.samples` API gets a different
> genotype than vcfpy/noodles, even though the underlying file data is identical.

**Genotypes contribute this one finding, not a family.** FORMAT subfield value
mutations (`GQ=1.5`, etc.) reproduce the already-documented INFO silent-coercion
class — the same pysam mechanism, no new one.

## The shallow-field result — parsers agree on `(chrom, pos, ref, alt)`

Bucket breakdown of the 3,000-file shallow corpus (~49s):

| count | bucket | meaning |
|------:|--------|---------|
| 1213 | `all_agree` | all parsed, identical content |
| 1038 | `split_decision` | some accept, some reject — the accept/reject signal |
| 749 | `all_fail` | all rejected (agreement) |
| 0 | `content_disagreement` | none — parsers agree on these fields |

*"Crash or hang" means a parser that **segfaulted or hung past the timeout** — a
failure the isolated-subprocess harness cannot catch from inside. A parser that
raises a **catchable** error and cleanly rejects a file is recorded as a
**reject**, not a crash.* Zero crashes occurred in any run.

The divergence on the shallow fields is entirely **accept/reject**. Four minimized
examples, ordered by value:

- **① pysam silently accepts an empty `POS` as position 0** — emits a variant at
  `20:0`; vcfpy and noodles reject. The shallow-field analogue of the INFO
  coercions: a plausible-looking wrong value that flows downstream undetected.
- **② noodles accepts a record with a missing column** (11 fields, one sample
  dropped); pysam and vcfpy reject `invalid number of fields`. Two robust parsers
  on opposite sides.
- **③ vcfpy accepts a trailing blank line** (0 records); pysam and noodles reject.
- **④ vcfpy errors on undeclared FORMAT fields**; pysam and noodles tolerate them.
  (A clean rejection, not a crash.)

## Future work — optional extensions of a complete project

These are deliberate choices about experiment design, not unfinished business:

- **Metamorphic indel mode.** In a *differential-read* test, indel normalization is
  a non-issue: all three parsers report indels **verbatim** (verified — none
  normalize on read), so there is nothing to reconcile. The classic "same variant,
  two valid spellings" class is real but belongs to a **different experiment**:
  generate equivalent spellings of one variant and check each parser is
  *self-consistent* across them. That is metamorphic testing, not differential —
  it would live as its own clearly-labeled mode, not folded into this pipeline.
- **A fourth independent parser** (e.g. htsjdk). Three voices already localize most
  disagreements (2-vs-1 tells you the outlier); a fourth strengthens
  majority-voting on accept/reject splits and would further isolate which parser
  is the odd one out.

## Reproducing

```
python run.py data/seed.vcf         # 3-way sanity on the valid seed
python run_mutations.py             # the 26 named single mutations
python run_fuzzer.py 3000           # shallow scaled corpus        (rng seed 1337)
python run_fuzzer.py 3000 info      # header-preserving INFO corpus (rng seed 4242)
python run_fuzzer.py 3000 gt        # header-preserving genotype corpus (rng seed 9999)
```
Catches are saved under `findings/<bucket>__<name>/` as `input.vcf`,
`minimized.vcf`, `report.txt`, and `report.json`.
