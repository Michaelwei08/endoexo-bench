#!/usr/bin/env python3
"""D: cross-individual junction recurrence -- a discriminator that lives in the
COHORT and not in any sample.

THE IDEA
  An endogenous insertion sits at the same host coordinate in every individual
  who inherited it. An exogenous integration sits at a coordinate private to the
  individual, because it happened in their lifetime. So sort the host-side
  junction coordinates by how many samples carry them: recurrent means
  endogenous, private means exogenous.

  This does not solve the tie / unopposed-win dichotomy. It sidesteps it. No
  sequence has to be told apart -- only positions counted -- which is why it can
  work in the regime where F024 showed competition tests are a no-op and F050
  showed no per-sample statistic beats counting reads.

WHAT IS OBSERVABLE AND WHAT IS NOT
  Observable: a read pair whose one end best-hits the exogenous reference while
  its mate best-hits the host reference, and the host mate's aligned position.
  Nothing else is used. The true integration coordinates are read back only to
  check the detector, never to build it.

RECURRENCE IS COUNTED LEAVE-ONE-OUT
  For sample s and bin b, "private" means no OTHER sample has a junction in b.
  Counting s itself would make every bin it carries non-private by construction
  and the statistic would be vacuous. That is the same class of error as
  read-level cross-validation and it is avoided the same way.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict

import numpy as np
from sklearn.metrics import roc_auc_score

from endoexo.align import PanelIndex, align_read
from endoexo.cohort import build_cohort, build_sample, locate
from endoexo.evaluate import fpr_at_sensitivity
from endoexo.simulate import add_errors


def _sample_pairs(template, segments, n_frags, rng, read_len, frag_mean, frag_sd,
                  error_rate):
    """Fragments from one template, with the SOURCE of each end recorded.

    A fragment straddling an insertion boundary is what produces a junction pair,
    and nothing special is done to create them -- they fall out of sampling a
    genome that has insertions in it.
    """
    out = []
    for _ in range(n_frags):
        flen = int(rng.normal(frag_mean, frag_sd))
        flen = max(read_len + 20, min(flen, len(template)))
        start = int(rng.integers(0, max(1, len(template) - flen)))
        end = start + flen
        r1 = add_errors(template[start:start + read_len].copy(), error_rate, rng)
        r2 = add_errors(template[end - read_len:end].copy(), error_rate, rng)
        src1, coord1 = locate(segments, start + read_len // 2)
        src2, coord2 = locate(segments, end - read_len // 2)
        out.append((r1, r2, src1, coord1, src2, coord2))
    return out


def one_sample(spec, infected, sample_index, cohort_seed, *, host_depth,
               clonal_fraction, read_len, frag_mean, frag_sd, error_rate,
               seed_k, index):
    """Detected junction bins for one individual, plus the truth to check them."""
    (with_prov, without, seg_with, seg_without, carried, exo_coord,
     d_strain) = build_sample(spec, infected=infected,
                              sample_seed=cohort_seed * 1000 + sample_index)
    rng = np.random.default_rng(cohort_seed * 1000 + sample_index + 7)

    n_total = int(host_depth * len(without) / (2 * read_len))
    pairs = []
    if infected and clonal_fraction > 0:
        n_clonal = int(round(n_total * clonal_fraction))
        pairs += _sample_pairs(with_prov, seg_with, n_clonal, rng, read_len,
                               frag_mean, frag_sd, error_rate)
        pairs += _sample_pairs(without, seg_without, n_total - n_clonal, rng,
                               read_len, frag_mean, frag_sd, error_rate)
    else:
        pairs = _sample_pairs(without, seg_without, n_total, rng, read_len,
                              frag_mean, frag_sd, error_rate)

    # OBSERVABLE detection only: one end on the exogenous reference, its mate on
    # the host reference, and the host mate's aligned position.
    junctions: list[int] = []
    n_exo_calls = 0
    for r1, r2, _s1, _c1, _s2, _c2 in pairs:
        h1 = align_read(r1, index)
        h2 = align_read(r2, index)
        if not h1 or not h2:
            continue
        a, b = h1[0], h2[0]
        if a.category == "EXO":
            n_exo_calls += 1
        if b.category == "EXO":
            n_exo_calls += 1
        if a.category == "EXO" and b.category == "HOST":
            junctions.append(b.ref_start)
        elif b.category == "EXO" and a.category == "HOST":
            junctions.append(a.ref_start)

    return {
        "infected": bool(infected),
        "carried_loci": carried,
        "true_exo_coord": exo_coord,
        "realized_strain_divergence": None if d_strain != d_strain else round(d_strain, 5),
        "n_pairs": len(pairs),
        "n_exo_calls": n_exo_calls,
        "junction_positions": junctions,
    }


def evaluate(samples, spec, bin_size: int, min_support: int) -> dict:
    """Bin the junctions, count recurrence leave-one-out, and score two rules."""
    bins_per_sample = []
    for s in samples:
        counts = defaultdict(int)
        for p in s["junction_positions"]:
            counts[p // bin_size] += 1
        bins_per_sample.append({b for b, n in counts.items() if n >= min_support})

    bin_sample_count = defaultdict(int)
    for bs in bins_per_sample:
        for b in bs:
            bin_sample_count[b] += 1

    y = np.array([s["infected"] for s in samples], dtype=int)
    # Rule A: does this sample hold a junction bin that no OTHER sample holds?
    private = np.array([
        sum(1 for b in bs if bin_sample_count[b] - 1 == 0) for bs in bins_per_sample
    ], dtype=float)
    # Rule B: the naive rule -- any junction at all
    any_junction = np.array([len(bs) for bs in bins_per_sample], dtype=float)

    out = {
        "bin_size": bin_size, "min_support": min_support,
        "n_samples": len(samples), "n_infected": int(y.sum()),
        "n_junction_bins_total": len(bin_sample_count),
    }
    if len(np.unique(y)) == 2:
        for name, score in (("n_private_junction_bins", private),
                            ("n_junction_bins", any_junction)):
            entry = {
                "roc_auc": round(float(roc_auc_score(y, score)), 4),
                "mean_infected": round(float(score[y == 1].mean()), 3),
                "mean_uninfected": round(float(score[y == 0].mean()), 3),
            }
            # FPR at 95 percent sensitivity is the WRONG metric for this rule and
            # is reported only so the mistake is visible. The score is a small
            # integer count, mostly 0 or 1, so demanding 95 percent sensitivity
            # forces the threshold down to "everything" and the FPR is 1.0 by
            # construction. This discriminator is high specificity and limited
            # sensitivity, so the operating point that means anything is a
            # threshold of one.
            entry["fpr_at_95_sens_DEGENERATE"] = round(
                float(fpr_at_sensitivity(y, score, 0.95)), 4)
            pred = score >= 1
            entry["at_threshold_1"] = {
                "sensitivity": round(float(pred[y == 1].mean()), 4),
                "fpr": round(float(pred[y == 0].mean()), 4),
                "specificity": round(float(1 - pred[y == 0].mean()), 4),
            }
            out[name] = entry

    # Does the detector find the truth? Checked, never used to build the rule.
    hits, misses = 0, 0
    for s, bs in zip(samples, bins_per_sample):
        if not s["infected"]:
            continue
        true_bin = s["true_exo_coord"] // bin_size
        if any(abs(b - true_bin) <= 1 for b in bs):
            hits += 1
        else:
            misses += 1
    out["true_integration_recovered"] = {"hit": hits, "miss": misses}

    # F064a: DECOMPOSE the non-private integrations by what they collide with.
    # Without this the sensitivity figure is uninterpretable: collision with
    # another integration is a simulation-scale artefact that a real genome does
    # not have, collision with an endogenous locus is a real failure mode, and a
    # spurious recurrent bin is a third thing entirely.
    true_bins = {i: s["true_exo_coord"] // bin_size
                 for i, s in enumerate(samples) if s["infected"]}
    endo_bins = {spec.locus_coord[n] // bin_size for n in spec.locus_coord}
    reasons = {"private_ok": 0, "collided_with_integration": 0,
               "collided_with_endogenous": 0, "collided_with_spurious": 0,
               "not_detected": 0}
    for i, s in enumerate(samples):
        if not s["infected"]:
            continue
        tb = true_bins[i]
        if not any(abs(b - tb) <= 1 for b in bins_per_sample[i]):
            reasons["not_detected"] += 1
            continue
        if bin_sample_count.get(tb, 0) <= 1:
            reasons["private_ok"] += 1
            continue
        others_with_integration = any(
            j != i and abs(true_bins[j] - tb) <= 1 for j in true_bins)
        near_endo = any(abs(b - tb) <= 1 for b in endo_bins)
        if others_with_integration:
            reasons["collided_with_integration"] += 1
        elif near_endo:
            reasons["collided_with_endogenous"] += 1
        else:
            reasons["collided_with_spurious"] += 1
    out["non_private_decomposition"] = reasons

    # A true integration can be DETECTED and still not counted as private, if its
    # bin collides with another sample's integration or with an endogenous locus.
    # That loss is a function of genome size: with a backbone of B bases and bin
    # size b there are B/b bins, and k integrations collide at the birthday rate.
    # It is a SIMULATION-SCALE artefact -- a real 3 Gb genome at this bin size has
    # about six million bins, where 25 integrations essentially never collide --
    # so the number is recorded to be discounted, not to be quoted.
    n_bins_available = len(spec.backbone) // bin_size
    k = int(y.sum())
    # The +/-1 bin tolerance in the recovery check makes each integration occupy
    # THREE bins, not one. Omitting that underestimated the expected collisions by
    # a factor of three and made the observed gap look partly unexplained.
    effective = 3.0 / n_bins_available
    out["collision_scale"] = {
        "bins_available": n_bins_available,
        "n_integrations": k,
        "bins_per_integration_with_tolerance": 3,
        "expected_colliding_pairs": round(k * (k - 1) / 2 * effective, 3),
        "expected_samples_involved": round(k * (k - 1) * effective, 2),
        "note": "artefact of the simulated backbone size, not of the method; a "
                "3 Gb genome at this bin size gives about 6e6 bins",
    }

    # And do the endogenous loci recur at roughly their population frequency?
    endo = []
    for name, coord in spec.locus_coord.items():
        b = coord // bin_size
        seen = max((bin_sample_count.get(b + d, 0) for d in (-1, 0, 1)), default=0)
        endo.append({"locus": name,
                     "population_frequency": round(spec.locus_freq[name], 3),
                     "samples_with_junction": seen,
                     "observed_fraction": round(seen / len(samples), 3)})
    out["endogenous_loci"] = endo

    rec = sorted(bin_sample_count.values())
    if rec:
        out["recurrence_distribution"] = {
            "n_bins": len(rec),
            "private_bins": int(sum(1 for r in rec if r == 1)),
            "private_fraction": round(sum(1 for r in rec if r == 1) / len(rec), 4),
            "max_samples_in_one_bin": rec[-1],
            "median": int(np.median(rec)),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", type=float, default=0.10)
    ap.add_argument("--n-samples", type=int, default=40)
    ap.add_argument("--infected-frac", type=float, default=0.5)
    ap.add_argument("--n-loci", type=int, default=10)
    ap.add_argument("--backbone-bp", type=int, default=100_000)
    ap.add_argument("--host-depth", type=float, default=6.0)
    ap.add_argument("--clonal-fraction", type=float, nargs="+", default=[0.5])
    ap.add_argument("--bin-size", type=int, default=500)
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--cohort-seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = []
    for cf in args.clonal_fraction:
        t0 = time.time()
        spec = build_cohort(args.divergence, n_loci=args.n_loci,
                            backbone_bp=args.backbone_bp, seed=args.cohort_seed)
        index = PanelIndex(spec.panel(), k=19)
        n_inf = int(round(args.n_samples * args.infected_frac))
        print(f"[divergence={args.divergence:.2f} clonal_fraction={cf:.2f} "
              f"n={args.n_samples} infected={n_inf}]")
        print(f"  locus divergence "
              f"{min(spec.realized_locus_divergence.values()):.3f}-"
              f"{max(spec.realized_locus_divergence.values()):.3f}, "
              f"frequencies {min(spec.locus_freq.values()):.2f}-"
              f"{max(spec.locus_freq.values()):.2f}")

        samples = []
        for i in range(args.n_samples):
            samples.append(one_sample(
                spec, i < n_inf, i, args.cohort_seed, host_depth=args.host_depth,
                clonal_fraction=cf, read_len=args.read_len, frag_mean=350,
                frag_sd=50, error_rate=0.002, seed_k=19, index=index))
            if (i + 1) % 10 == 0:
                print(f"    {i+1}/{args.n_samples} ({time.time()-t0:.0f}s)")

        ev = evaluate(samples, spec, args.bin_size, args.min_support)
        ev["divergence"] = args.divergence
        ev["clonal_fraction"] = cf
        ev["cohort_seed"] = args.cohort_seed
        ev["seconds"] = round(time.time() - t0, 1)
        print(json.dumps({k: v for k, v in ev.items() if k != "endogenous_loci"},
                         indent=2))
        results.append(ev)

    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            json.dump(results, fh, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
