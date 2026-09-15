#!/usr/bin/env python3
"""Reproduce the tie typology from a real BWA-MEM alignment and compare it with
this repository's own ungapped aligner on the SAME reads.

THE QUESTION THIS ANSWERS
  The project's sharpest claim is that mis-assignment splits into exact score
  ties and unopposed wins, and that the tie half is free to detect because
  `AS == XS` is an equality test. That was measured with an ungapped aligner.
  Gapped local alignment may not preserve the exactness: if a gap buys a base
  somewhere, a tie becomes a near-tie and the equality test becomes a threshold.

  So the number to look at is not only the tied FRACTION. It is the distribution
  of `AS - XS` among false calls. Under the ungapped aligner it is a spike at
  exactly 0. If BWA spreads it to 1 or 2, the claim needs rewording.

ONE DIFFERENCE TO EXPECT, BEFORE THE NUMBERS ARRIVE
  The two aligners break ties differently, and that is not a bug in either. The
  first-party aligner sorts by score alone, so a tie resolves to whichever
  reference entered the panel first -- in practice the exogenous one. BWA
  resolves a tie by picking a primary essentially arbitrarily, so roughly half
  of the tied reads will be reported against the host and will never become
  exogenous calls at all.

  So the `award_ties_to_exo` call count MUST differ between the two, and the
  first-party number is the pessimistic bound. The `three_bin_ambiguous` numbers
  should agree, because that policy excludes ties under either tie-break. If
  they do agree, that is one more argument that the three-bin policy is the
  robust one: it is the only policy whose output does not depend on an arbitrary
  implementation choice inside the aligner.

INPUT
  The directory written by export_for_bwa.py, plus a SAM produced from the same
  reads by `bwa mem -a` (all alignments, not just the primary -- the category
  tie structure cannot be recovered from primary records alone).

  Run inside WSL, where bwa and samtools live:
    bash run_bwa.sh bwa_compare
  then
    python compare_with_bwa.py bwa_compare
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_reference_map(path: Path) -> dict[str, str]:
    cats = {}
    with path.open(encoding="ascii") as fh:
        next(fh)
        for line in fh:
            ref_id, category, _len = line.rstrip("\n").split("\t")
            cats[ref_id] = category
    return cats


def read_truth(path: Path) -> dict[str, dict]:
    truth = {}
    with path.open(encoding="ascii") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("read_name"):
                continue
            name, source, locus, start, ldiv = line.rstrip("\n").split("\t")
            truth[name] = {"source": source, "locus": locus,
                           "frag_start": int(start),
                           "local_div": float(ldiv) if ldiv else float("nan")}
    return truth


def softclip_len(cigar: str) -> int:
    """Total soft-clipped bases at both ends."""
    total, num = 0, ""
    ops = []
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num or 0), ch))
            num = ""
    if ops and ops[0][1] == "S":
        total += ops[0][0]
    if len(ops) > 1 and ops[-1][1] == "S":
        total += ops[-1][0]
    return total


def parse_sam(path: Path, categories: dict[str, str]) -> dict[str, dict]:
    """Collect, per read END, the best alignment score on each reference.

    `bwa mem -a` emits one record per alignment, primary and secondary. Grouping
    by reference and keeping the maximum AS is what makes this comparable to the
    first-party aligner, which returns one best hit per reference by
    construction.
    """
    per_end: dict[tuple[str, int], dict[str, int]] = defaultdict(dict)
    clip: dict[tuple[str, int], int] = {}
    opened = sys.stdin if str(path) == "-" else path.open(encoding="ascii")
    try:
        for line in opened:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            name, flag, ref, cigar = f[0], int(f[1]), f[2], f[5]
            if ref == "*" or flag & 0x4:
                continue
            end = 2 if flag & 0x80 else 1
            score = None
            for tag in f[11:]:
                if tag.startswith("AS:i:"):
                    score = int(tag[5:])
                    break
            if score is None:
                continue
            key = (name, end)
            prev = per_end[key].get(ref)
            if prev is None or score > prev:
                per_end[key][ref] = score
                if not (flag & 0x100):
                    clip[key] = softclip_len(cigar)
    finally:
        if opened is not sys.stdin:
            opened.close()

    out: dict[str, dict] = {}
    for (name, end), by_ref in per_end.items():
        if end != 1:
            continue                     # features are defined on the first end
        best_score = max(by_ref.values())
        tied_refs = [r for r, s in by_ref.items() if s == best_score]
        tied_cats = {categories.get(r, "?") for r in tied_refs}
        others = sorted((s for r, s in by_ref.items() if r not in tied_refs),
                        reverse=True)
        xs = others[0] if others else (best_score if len(tied_refs) > 1 else 0)
        best_ref = tied_refs[0]
        out[name] = {
            "best_ref": best_ref,
            "best_cat": categories.get(best_ref, "?"),
            "AS": best_score,
            "XS": xs,
            "as_minus_xs": best_score - xs,
            "n_tied_top": len(tied_refs),
            "n_tied_top_categories": len(tied_cats),
            "tied_categories": sorted(tied_cats),
            "softclip_len": clip.get((name, 1), 0),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirname")
    ap.add_argument("--sam", default=None,
                    help="SAM path; defaults to <dirname>/aln.sam, '-' reads stdin")
    args = ap.parse_args()
    d = Path(args.dirname)
    sam = Path(args.sam) if args.sam else d / "aln.sam"

    categories = read_reference_map(d / "reference_map.tsv")
    truth = read_truth(d / "truth.tsv")
    feats = parse_sam(sam, categories)

    # the same selection rule as run_slice: reads a pipeline would COUNT
    calls = [(n, f) for n, f in feats.items() if f["best_cat"] == "EXO"]
    if not calls:
        print("no read was assigned to the exogenous reference")
        return

    y = np.array([truth[n]["source"] == "EXO" for n, _ in calls])
    gap = np.array([f["as_minus_xs"] for _, f in calls])
    amb = np.array([f["n_tied_top_categories"] > 1 for _, f in calls])

    n_exo_reads = sum(1 for t in truth.values() if t["source"] == "EXO")
    print(f"reads in truth        : {len(truth):,}  (exogenous {n_exo_reads:,})")
    print(f"aligned somewhere     : {len(feats):,}")
    print(f"assigned to EXO       : {len(calls):,}  "
          f"(true {int(y.sum()):,}, false {int((~y).sum()):,})")
    print()

    for label, mask in (("FALSE calls", ~y), ("TRUE calls", y)):
        if not mask.any():
            continue
        g = gap[mask]
        print(f"{label}: n={mask.sum():,}")
        print(f"  cross-category tie at the top : {100*amb[mask].mean():.1f}%")
        print(f"  AS - XS == 0                  : {100*(g == 0).mean():.1f}%")
        print(f"  AS - XS <= 2                  : {100*(g <= 2).mean():.1f}%")
        print(f"  AS - XS distribution          : "
              f"min {g.min()} p50 {int(np.median(g))} p95 {int(np.percentile(g, 95))} max {g.max()}")
        # THE number: if the ungapped spike at 0 has spread, it shows up here
        nonzero = g[g > 0]
        if nonzero.size:
            print(f"  among AS-XS > 0               : n={nonzero.size:,}, "
                  f"p50 {int(np.median(nonzero))}, max {nonzero.max()}")
        print()

    tp_all, fp_all = int(y.sum()), int((~y).sum())
    keep = ~amb
    tp, fp = int((y & keep).sum()), int(((~y) & keep).sum())
    print("policy comparison, as in run_slice:")
    print(f"  award_ties_to_exo   : calls {tp_all+fp_all:,}  false "
          f"{100*fp_all/(tp_all+fp_all):.1f}%  sensitivity "
          f"{tp_all/n_exo_reads:.3f}")
    print(f"  three_bin_ambiguous : calls {tp+fp:,}  false "
          f"{(100*fp/(tp+fp)) if (tp+fp) else float('nan'):.1f}%  sensitivity "
          f"{tp/n_exo_reads:.3f}")


if __name__ == "__main__":
    main()
