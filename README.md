# dissent: a differential fuzzer for VCF parsers

Scientists store information about DNA in files called VCF
files, and a great deal of software exists to read those files. This project
checks whether three of the most widely used VCF-reading programs actually agree
on what a file says. It does that by feeding them thousands of tricky files and
watching for the moments where they quietly give different answers without any of
them reporting an error.

You can think of a person's DNA as an extremely long
string of letters. When researchers study someone's genome, they usually do not
write down the entire string. Instead they record only the spots where that
person's DNA differs from a standard reference, and those differences are called
*variants*. The standard way to store a list of variants is a text-based file
format called **VCF** (Variant Call Format). VCF files are everywhere in genetics
and medical research, and a huge amount of scientific software exists just to
read them, filter them, and analyze them.

A VCF file is, in the end, just text arranged into
columns and fields. Reading it sounds simple, but the format has many small rules
and awkward edge cases, and real-world files are often slightly messy or even
malformed. Every program that reads VCF, each written by different people in
different programming languages, has to decide how to interpret every field. Most
of the time they interpret them identically. But once in a while, on an unusual or
slightly broken input, two programs read the *same bytes* and come away with
*different answers*. Here is the dangerous part: neither program raises an error.
Everything looks fine. The programs simply, silently, disagree.

In science, a silent disagreement like this is worse than a
crash. A crash is loud, so a researcher notices it and fixes it. A silent
difference slips right through. If your analysis pipeline reads a variant one way,
and a colleague's pipeline reads the very same file a different way, your results
can quietly diverge, and no error message ever warns either of you. A bug like
that can sit undetected inside real, published research.

dissent hunts for exactly those silent disagreements,
automatically, in a handful of steps:

1. Start from one known-good VCF file.
2. Generate thousands of variations of it, each deliberately altered in a small,
   understandable way (shift a number by one, delete a field, add a stray
   character, put a value where a different type is expected, and so on).
3. Hand each generated file to three different, independently written VCF-reading
   programs.
4. Boil each program's answer down to a simple, shared summary and compare them.
5. Whenever the programs disagree, or one of them crashes, save the offending file.
6. Automatically shrink that file down to the smallest version that still triggers
   the disagreement, so the result is a tiny, shareable example instead of a large
   messy one.

he comparison only means something
if the three programs were built separately. Many popular VCF tools are actually
thin wrappers around the *same* underlying code, which means they would happily
agree with each other even when all of them are wrong. dissent deliberately uses
three genuinely independent programs, one written in C, one in pure Python, and
one in pure Rust. That way, when they disagree, it points to a real bug in at
least one of them, and nobody needs to know ahead of time what the "correct"
answer was. That last point is the whole trick. You never need an answer key. You
only need independent programs that ought to agree.

> **An analogy.** Imagine handing the same contract to three professional
> translators and comparing their translations. If all three produce the same
> result, wonderful. If two agree and the third quietly produces a *different*
> translation, not by throwing up their hands and saying "I cannot read this," but
> by confidently handing you a different answer, then you have found a mistake in
> one of them. dissent does this automatically, thousands of times over, for the
> software that biologists actually rely on.

**What it found.** Run at scale, dissent turned up real cases where these three
trusted programs silently disagree on the contents of a file. The most striking
one is invisible even if you compare the files the programs write back out. One
program's *reader* hands your code a different genotype than the others do, yet
that same program's *writer* still produces the original bytes, so no
file-to-file comparison could ever catch it. Only by comparing what each program's
code actually returns to its caller does the difference show up. The full catalog
is in the report linked below.

---

## For the technical reader

In one line: three production VCF parsers silently disagree on typed field content
(INFO values, and how one of them hands genotypes to its callers) on
ordinary-looking inputs, and this tool finds and minimizes the exact files where
they do.

Feed one (often malformed) VCF to three independently built parsers (pysam/htslib,
vcfpy, and noodles in pure Rust) and catch where they silently disagree or crash.
A disagreement proves at least one of them has a bug without needing to know the
"correct" answer in advance, which is the whole point of differential testing.

## The headline: a disagreement all three parsers accept

Give every parser this one line, a valid-looking `AA=T,99` in a field the header
declares single-valued (`Number=1`):

```
##INFO=<ID=AA,Number=1,Type=String,...>
20   1230237   .   T   .   ...   NS=3;DP=13;AA=T,99
```

Nobody errors. Yet:

- **pysam** reads `AA` as **two** values, `('T', '99')`
- **vcfpy** and **noodles** read it as **one** value, `'T,99'`

A downstream tool asking "how many values does `AA` have?" gets **2** from pysam
and **1** from the others, silently, with no error anywhere. Proof that this is
real content and not merely a difference in representation: pysam re-emits
`AA=T,99` while vcfpy re-emits `AA=T%2C99`, which are different bytes.

## And one you could never catch with a file diff

Subtler still: pysam's `.samples['GT']` API, which is how most pysam users read
genotypes, reinterprets an *invalid* allele index as **missing**, where vcfpy and
noodles report the file's value. pysam's *writer* faithfully preserves the raw
genotype, so no comparison of output files could ever detect it. Only comparing
what each parser hands its caller does.

**See every finding, with the round-trip discipline that separates real bugs from
representation noise: [findings/REPORT.md](findings/REPORT.md)**

## Why it's useful

- **No answer key required.** Differential testing sidesteps the oracle problem.
  You do not need to know the "right" parse, only that two independent parsers
  should agree. When they do not, at least one of them is wrong.
- **It targets the dangerous class.** A crash is loud, so someone notices. A
  *silent* content disagreement corrupts an analysis with nobody the wiser,
  including one finding here that is invisible even to a diff of the parsers'
  output files.
- **Every catch is tiny and reproducible.** Findings are automatically minimized
  to a few lines, so a bug report is a 3-line file rather than a 500-line mess.

## The three voices (why independence matters)

Many popular VCF tools (pysam, bcftools, cyvcf2, samtools) are thin wrappers over
the **same** C library, [htslib](https://github.com/samtools/htslib), so they
would agree even when all of them are wrong. This project only compares genuinely
independent codebases:

| Parser | Language | Backend |
|--------|----------|---------|
| pysam | Python/Cython | htslib (C), the htslib representative |
| vcfpy | pure Python | own |
| noodles | Rust | own, pure Rust, no C bindings |

noodles' independence is verified from its dependency graph (`cargo tree`): no
`hts` / `*-sys` / `bindgen` crate appears anywhere, so it shares no code with
htslib and cannot merely echo pysam.

## Getting started

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cargo build --release --manifest-path adapters/noodles_adapter/Cargo.toml

python run.py data/seed.vcf        # sanity: all 3 agree on a valid file
python run_fuzzer.py 3000 info     # the INFO hunt (the comma-split lives here)
python run_fuzzer.py 3000 gt       # the genotype hunt
```

Each parser runs in an isolated subprocess, so a crash or hang is caught rather
than fatal, and catches are minimized to a tiny file under `findings/`.

### Layout

```
generators/  mutate a valid seed into pathological / value-mutated files
adapters/    one adapter per parser: file in -> normalized summary out (noodles is a Rust binary)
compare/     classify agree / disagree / crash; triage representation vs content; minimize catches
findings/    saved catches + REPORT.md (the full write-up)
data/        seed VCF
```

## Getting help

- **A question, or a disagreement of your own?** Open an issue on the GitHub repo,
  or email the maintainer (below).
- **The full methodology and every finding:** [findings/REPORT.md](findings/REPORT.md).
- **Reproduce any catch:** each `findings/<bucket>__<name>/` directory holds the
  minimized `input.vcf`, all three parsers' outputs (`report.txt` and
  `report.json`), and the original file it was shrunk from.

## Maintainers & contributing

Maintained by **[@guhangopinathanai23](https://github.com/guhangopinathanai23)**
(guhangopinathanai@gmail.com). Contributions and bug reports are welcome; open an
issue or a pull request. Two good places to start, both described in the
[report](findings/REPORT.md)'s Future Work:

- add a **fourth independent parser** (for example htsjdk) as a new `adapters/`
  adapter, which strengthens majority-voting on which parser is the outlier;
- build the **metamorphic indel mode**, the honest way to chase the classic
  "same variant, two valid spellings" bug class.

*Versions under test:* **pysam 0.24.0** (htslib 1.23.1), **vcfpy 0.13.8**,
**noodles-vcf 0.88.0**. Python 3.9, rustc 1.96.

## License

MIT license, see [LICENSE](LICENSE). Contributions are accepted under the same
license.
