#!/usr/bin/env python3
"""Phase 3: scaled differential fuzzer.

Pipeline:
  1. generate N stacked-mutation files (traceable via manifest)
  2. compare each through all 3 parsers, in parallel
  3. classify (content-disagreement-first + normalization-artifact triage gate)
  4. minimize a tiny representative of every distinct disagreement signature
  5. report: file count, bucket breakdown, crashes/hangs, and for each
     disagreement the minimized file + all three parsers' outputs

Content disagreements (all parse, differ) are the top-priority class. Where the
run finds none, the accept/reject SPLITS are the findings -- and the ones where
the two robust parsers (pysam vs noodles) land on opposite sides are the most
interesting, since that is two serious independent implementations disagreeing
on what is even valid.

Usage: python run_fuzzer.py [N]     (default N=3000)
"""
from __future__ import annotations

import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from adapters.common import KIND_CRASH, KIND_TIMEOUT
from adapters.registry import REPO_ROOT
from compare.compare import (
    BUCKET_ARTIFACT,
    BUCKET_CONTENT,
    BUCKET_SPLIT,
    compare_file,
    render_text_report,
    save_finding,
    short,
    sides,
    signature,
)
from compare.minimize import minimize
from generators.corpus import generate_corpus, generate_gt_corpus, generate_info_corpus

MAX_WORKERS = 12

# mode -> (generator, corpus dir). "info"/"gt" keep headers and mutate values.
MODES = {
    "standard": (generate_corpus, os.path.join(REPO_ROOT, "data", "corpus")),
    "info": (generate_info_corpus, os.path.join(REPO_ROOT, "data", "corpus_info")),
    "gt": (generate_gt_corpus, os.path.join(REPO_ROOT, "data", "corpus_gt")),
}
OUT_DIR = MODES["standard"][1]  # set by main() per mode
MAX_MINIMIZE_PER_BUCKET = 10  # cap minimized representatives per bucket


def _minimize_one(rep, comp, paths, manifest, nfiles):
    sig = signature(comp)
    mtext, mlines = minimize(paths[rep], sig, OUT_DIR, tag=rep)
    minpath = os.path.join(OUT_DIR, f"_final_{rep}")
    with open(minpath, "w", newline="") as f:
        f.write(mtext)
    mcomp = compare_file(minpath)
    dest = save_finding(mcomp, f"{rep}_min", extra={
        "provenance": manifest[rep],
        "files_with_this_signature": nfiles,
        "minimized_lines": mlines,
        "minimized_vcf": mtext,
    })
    with open(os.path.join(dest, "minimized.vcf"), "w", newline="") as f:
        f.write(mtext)
    return (rep, nfiles, manifest[rep], mtext, mcomp, dest)


def minimize_bucket(results, paths, manifest, bucket):
    """Group all results in `bucket` by signature; minimize a small
    representative of each distinct signature, in parallel."""
    members = {f: c for f, c in results.items() if c.bucket == bucket}
    by_sig = defaultdict(list)
    for f, c in members.items():
        by_sig[signature(c)].append(f)
    # order signatures by frequency (most common first), cap the count
    ordered = sorted(by_sig.items(), key=lambda kv: -len(kv[1]))[:MAX_MINIMIZE_PER_BUCKET]
    reps = [(min(fs, key=lambda f: os.path.getsize(paths[f])), len(fs)) for _, fs in ordered]

    out = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(_minimize_one, rep, results[rep], paths, manifest, n): rep
                for rep, n in reps}
        for fut in as_completed(futs):
            out.append(fut.result())
    out.sort(key=lambda r: -r[1])
    return out, len(by_sig)


def _print_findings(title, reports, total_sigs):
    print("\n" + "#" * 84)
    extra = "" if total_sigs <= len(reports) else f" (showing top {len(reports)} of {total_sigs})"
    print(f"# {title}: {total_sigs} distinct signature(s){extra}")
    print("#" * 84)
    for rep, nfiles, prov, mtext, mcomp, dest in reports:
        nlines = len(mtext.rstrip("\n").split("\n"))
        print(f"\n--- {rep}  (x{nfiles} files, same signature)  minimized to {nlines} line(s) ---")
        print(f"    provenance: {prov}")
        print(f"    split: {sides(mcomp)}   silent(all parsed)={mcomp.all_parsed}")
        print("    minimized file:")
        for ln in mtext.rstrip("\n").split("\n"):
            print(f"        {ln!r}")
        print("    parser outputs:")
        for s in mcomp.summaries:
            if s.ok:
                recs = "; ".join(f"{c}:{p} {r}->{','.join(a) if a else '.'}{('  INFO='+str(dict(info))) if info else ''}{('  FMT='+str([dict(x) for x in samples])) if samples else ''}" for (c, p, r, a, info, samples) in s.records)
                print(f"        [{short(s.tool):8}] OK  {recs}")
            else:
                print(f"        [{short(s.tool):8}] {s.kind}: {s.error}")
            for note in getattr(s, "notes", []) or []:
                print(f"                   NOTE: {note}")
        print(f"    [saved: {os.path.relpath(dest, REPO_ROOT)}]")


def main(argv) -> int:
    global OUT_DIR
    n = int(argv[1]) if len(argv) > 1 else 3000
    mode = argv[2] if len(argv) > 2 else "standard"
    generator, OUT_DIR = MODES[mode]
    t0 = time.time()

    print(f"[1/4] generating {n} files (mode={mode}) ...")
    manifest = generator(n, OUT_DIR)
    files = sorted(manifest.keys())
    paths = {f: os.path.join(OUT_DIR, f) for f in files}

    print(f"[2/4] comparing {len(files)} files x 3 parsers (parallel) ...")
    results = {}
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(compare_file, paths[f]): f for f in files}
        for fut in as_completed(futs):
            results[futs[fut]] = fut.result()
    t_cmp = time.time()

    counts = Counter(c.bucket for c in results.values())

    # crash / hang scan across every parser on every file
    crashes = []
    for f, c in results.items():
        for s in c.summaries:
            if s.kind in (KIND_CRASH, KIND_TIMEOUT):
                crashes.append((f, short(s.tool), s.kind, s.error))

    print("[3/4] minimizing representatives of content disagreements + splits (parallel) ...")
    content_reports, n_content_sigs = minimize_bucket(results, paths, manifest, BUCKET_CONTENT)
    split_reports, n_split_sigs = minimize_bucket(results, paths, manifest, BUCKET_SPLIT)
    t1 = time.time()

    # ---- report ----
    print("\n" + "=" * 84)
    print(f"SCALED RUN: {n} files  |  compare {t_cmp - t0:.0f}s, minimize {t1 - t_cmp:.0f}s, total {t1 - t0:.0f}s")
    print("-" * 84)
    print("bucket breakdown (priority order):")
    for b in [BUCKET_CONTENT, BUCKET_ARTIFACT, BUCKET_SPLIT, "all_fail", "all_agree", "harness_error"]:
        if counts.get(b):
            print(f"   {counts[b]:6}  {b}")
    print("=" * 84)

    print(f"\n[4/4] ROBUSTNESS (crashes/hangs -- loudest bugs): {len(crashes)} occurrence(s)")
    if crashes:
        by_tool = Counter((t, k) for _, t, k, _ in crashes)
        for (tool, kind), ct in by_tool.most_common():
            example = next(e for e in crashes if e[1] == tool and e[2] == kind)
            print(f"   {tool} {kind} x{ct}  e.g. {example[0]}: {example[3][:70]}")
    else:
        print("   none -- no parser segfaulted or hung on any file.")

    print(f"\nCONTENT DISAGREEMENTS (all/some parse but differ, survived triage): "
          f"{counts.get(BUCKET_CONTENT, 0)} file(s)")
    print(f"NORMALIZATION ARTIFACTS quarantined by triage gate: {counts.get(BUCKET_ARTIFACT, 0)} file(s)")
    if content_reports:
        _print_findings("REAL CONTENT DISAGREEMENTS", content_reports, n_content_sigs)
    else:
        print("   -> none. Wherever >=2 parsers both accept a file, they agree on content.")

    if split_reports:
        _print_findings("ACCEPT/REJECT SPLITS (minimized examples)", split_reports, n_split_sigs)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
