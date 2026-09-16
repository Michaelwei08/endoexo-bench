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

WHAT BWA ACTUALLY EMITS, learned the hard way 2026-09-16
  `bwa mem -a` does NOT emit secondary alignments for properly paired reads. On
  the first real run there was exactly one record per read end -- 17,246 records
  for 8,623 pairs -- and zero records carrying flag 0x100 or 0x800. A parser that
  reconstructs the second-best score by comparing records therefore saw a single
  reference per read, concluded XS = 0, and reported an AS - XS gap of 150 for
  every false call with 0.0 percent ties. Those numbers were a parser artefact
  and were nearly reported as a finding.

  The information is in the SAM, just not there. Every record carries BWA's own
  `XS:i:` second-best score, and 1,286 records carried MAPQ 0, which is BWA's
  signal for an ambiguous placement. Both are used now.

  ONE THING REMAINS UNOBSERVABLE, and it matters. BWA's XS tag gives the
  second-best SCORE but not the reference it was on, so the CATEGORY of the
  competitor -- the thing this project's whole tie typology is defined on -- is
  not recoverable from default output. It needs `-h`, which emits XA:Z listing
  alternative hits by reference name. That is a practical limitation on whether
  the three-bin rule is implementable in a standard pipeline at all, and it is
  reported rather than papered over: `competitor_category_known` says whether the
  question could be answered for each read.

INPUT
  The directory written by export_for_bwa.py, plus a SAM produced from the same
  reads by `bwa mem -a -h 200`.

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
    primary: dict[tuple[str, int], dict] = {}
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
            tags = {t[:2]: t for t in f[11:]}
            score = int(tags["AS"][5:]) if "AS" in tags else None
            if score is None:
                continue
            key = (name, end)
            prev = per_end[key].get(ref)
            if prev is None or score > prev:
                per_end[key][ref] = score
            if not (flag & 0x900):
                primary[key] = {
                    "ref": ref,
                    "AS": score,
                    # BWA's OWN second-best score, which is the number to use.
                    "XS_tag": int(tags["XS"][5:]) if "XS" in tags else None,
                    "mapq": int(f[4]),
                    "softclip_len": softclip_len(cigar),
                    # XA:Z lists alternative hits as ref,pos,CIGAR,NM; it is the
                    # only place BWA names the COMPETING reference, and it is
                    # emitted only when -h is given.
                    "xa": tags.get("XA"),
                }
    finally:
        if opened is not sys.stdin:
            opened.close()

    out: dict[str, dict] = {}
    for (name, end), prim in primary.items():
        if end != 1:
            continue                     # features are defined on the first end
        by_ref = per_end[(name, end)]
        best_score = prim["AS"]
        best_ref = prim["ref"]

        # Prefer BWA's XS tag. Fall back to cross-record comparison only when
        # secondary records actually exist, which for properly paired reads they
        # do not -- see the module docstring.
        others = sorted((s for r, s in by_ref.items() if r != best_ref),
                        reverse=True)
        xs = prim["XS_tag"] if prim["XS_tag"] is not None else (
            others[0] if others else 0)

        # Which CATEGORY is the competitor in? Recoverable only from XA:Z or from
        # genuine secondary records. Without either, it is unobservable and says
        # so rather than defaulting to "no tie".
        # Only alternatives TIED AT THE TOP count. Taking every alternative
        # regardless of score was a bug: a reference scoring far below the
        # primary is not a competitor, and counting its category inflated the
        # ambiguity. Caught by a unit test, 2026-09-16.
        alt_refs = [r for r, s in by_ref.items()
                    if r != best_ref and s == best_score]
        if prim["xa"]:
            # XA entries carry ref,pos,CIGAR,NM and NO score. BWA emits them only
            # for hits within 80 percent of the best, so an XA reference is
            # near-best but not verifiably tied. Included, with that caveat.
            for hit in prim["xa"][5:].rstrip(";").split(";"):
                if hit:
                    alt_refs.append(hit.split(",")[0])
        alt_cats = {categories.get(r, "?") for r in alt_refs}
        competitor_known = bool(alt_refs)

        out[name] = {
            "best_ref": best_ref,
            "best_cat": categories.get(best_ref, "?"),
            "AS": best_score,
            "XS": xs,
            "as_minus_xs": best_score - xs,
            "mapq": prim["mapq"],
            "softclip_len": prim["softclip_len"],
            "competitor_category_known": competitor_known,
            "competitor_categories": sorted(alt_cats),
            # Two observable tie proxies, since the exact category-level tie is
            # not recoverable from default BWA output.
            "tie_by_score": int(best_score == xs),
            "tie_by_mapq0": int(prim["mapq"] == 0),
            "n_tied_top_categories": (
                len(alt_cats | {categories.get(best_ref, "?")})
                if competitor_known and best_score == xs else 1),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirname")
    ap.add_argument("--sam", default=None,
                    help="SAM path; defaults to <dirname>/aln.sam, '-' reads stdin")
    ap.add_argument("--out", default=None,
                    help="write the measured quantities as JSON, so the paper's "
                         "numbers trace to a file like every other result")
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

    mapq0 = np.array([f["tie_by_mapq0"] for _, f in calls], dtype=bool)
    known = np.array([f["competitor_category_known"] for _, f in calls], dtype=bool)
    print(f"competitor category recoverable from the SAM for "
          f"{100*known.mean():.1f}% of calls (needs XA:Z, i.e. -h)")
    print()

    for label, mask in (("FALSE calls", ~y), ("TRUE calls", y)):
        if not mask.any():
            continue
        g = gap[mask]
        print(f"{label}: n={mask.sum():,}")
        print(f"  AS - XS == 0  (BWA's own XS tag) : {100*(g == 0).mean():.1f}%")
        print(f"  AS - XS <= 2                     : {100*(g <= 2).mean():.1f}%")
        print(f"  MAPQ == 0 (BWA's ambiguity flag)  : {100*mapq0[mask].mean():.1f}%")
        print(f"  cross-category tie, where knowable: "
              f"{100*amb[mask & known].mean():.1f}%" if (mask & known).any()
              else "  cross-category tie, where knowable: n/a")
        print(f"  AS - XS distribution             : "
              f"min {g.min()} p50 {int(np.median(g))} p95 {int(np.percentile(g, 95))} max {g.max()}")
        # THE number: if the ungapped spike at 0 has spread, it shows up here
        near = g[(g > 0) & (g <= 10)]
        print(f"  in the near-tie band 1..10       : {near.size:,} "
              f"({100*near.size/max(1, g.size):.1f}%)")
        print()

    record = {
        "n_reads": len(truth), "n_exogenous_reads": n_exo_reads,
        "n_calls": len(calls), "n_true": int(y.sum()), "n_false": int((~y).sum()),
        "competitor_category_knowable_frac": round(float(known.mean()), 4),
        "false_calls": {
            "tie_by_score_frac": round(float((gap[~y] == 0).mean()), 4),
            "tie_by_mapq0_frac": round(float(mapq0[~y].mean()), 4),
            "near_tie_band_1_10": int(((gap[~y] > 0) & (gap[~y] <= 10)).sum()),
        },
        "true_calls": {
            "tie_by_score_frac": round(float((gap[y] == 0).mean()), 4),
            "tie_by_mapq0_frac": round(float(mapq0[y].mean()), 4),
            "as_minus_xs_min": int(gap[y].min()),
        },
        "policies": {},
    }
    print("policy comparison. Sensitivity is against ALL simulated exogenous reads.")
    policies = [
        ("award_ties_to_exo   ", np.ones(len(y), dtype=bool)),
        # The implementable rule. Category is unobservable from standard output,
        # but BWA's own MAPQ 0 marks the same reads -- and it is one field that
        # every pipeline already has.
        ("discard_mapq0       ", ~mapq0),
        # The rule as this project defines it, which needs the COMPETITOR's
        # category. Reported only where that is knowable, so its absence is
        # visible rather than silently collapsing into the first row.
        ("three_bin_by_category", ~amb if known.any() else None),
    ]
    for name, keep in policies:
        if keep is None:
            print(f"  {name}: NOT COMPUTABLE -- no XA:Z in the SAM, so the "
                  f"competitor's category is unknown for every call")
            continue
        tp = int((y & keep).sum())
        fp = int(((~y) & keep).sum())
        frac = (100 * fp / (tp + fp)) if (tp + fp) else float("nan")
        print(f"  {name}: calls {tp+fp:,}  false {frac:.1f}%  "
              f"sensitivity {tp/n_exo_reads:.3f}")
        record["policies"][name.strip()] = {
            "calls": tp + fp, "false_calls": fp,
            "false_call_frac": round(frac / 100, 4) if frac == frac else None,
            "sensitivity_vs_all_exo_reads": round(tp / n_exo_reads, 4)}

    if args.out:
        import json
        versions = (d / "aligner_versions.txt")
        record["aligner"] = (versions.read_text(encoding="ascii").strip()
                             if versions.exists() else None)
        with open(args.out, "w", encoding="ascii") as fh:
            json.dump(record, fh, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
