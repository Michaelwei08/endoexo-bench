#!/usr/bin/env python3
"""The SAMPLE-level detection task: is this sample infected, or is it carrying
endogenous cross-mapping only?

WHY THIS EXISTS
  F026 established a hard negative result. When the endogenous locus is absent
  from the panel and the elements are young, no read-pair-level feature beats
  chance -- grouped PR AUC 0.138 against a 0.173 prevalence floor. The three-bin
  rule is a no-op there (F024) and the mate carries nothing (F025).

  But a detection assay does not report reads. It reports SAMPLES. And the one
  thing that distinguishes an infection from cross-mapping is not visible in any
  single read: an infection puts reads across the WHOLE provirus, while
  cross-mapping puts them only where the endogenous copy happens to be conserved
  relative to the exogenous reference. That is blocky, and it differs from
  sample to sample because individuals carry different insertions.

  So this asks whether an aggregate statistic succeeds where the read-level one
  is below chance -- and compares the statistic in production use (a read count
  over a threshold) against breadth and uniformity.

WHAT IS COMPARED
  Each candidate is treated as a continuous SCORE over samples, so no threshold
  has to be chosen and the comparison is threshold-free:

    n_calls              the count a pipeline reports today
    n_calls_unopposed    the same count after the three-bin rule
    breadth_1x           fraction of the exogenous reference covered at all
    neg_max_bin_frac     uniformity: how much depth sits in the single worst bin
    neg_depth_cv         uniformity: coefficient of variation across bins
    learned model        logistic regression and a GBM over all of the above

  Viral load VARIES across infected samples, drawn log-uniformly, because the
  discriminating power of a read count collapses near the detection limit and
  that is the regime the question matters in.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from endoexo.align import PanelIndex, features_for_pair
from endoexo.evaluate import _make_model, fpr_at_sensitivity
from endoexo.simulate import simulate_reads, simulate_world
from run_slice import PANEL_MODES, restrict_panel

SAMPLE_FEATURES = [
    "n_calls", "n_calls_unopposed", "tie_frac",
    "breadth_1x", "breadth_5x", "depth_cv", "max_bin_frac",
    "mean_AS", "mean_gap", "mean_softclip", "mean_aligned_len",
    "frac_in_ltr", "frac_mate_same_ref", "frac_mate_exo",
]

# Single statistics evaluated on their own, and the sign that makes "higher is
# more likely infected" true for each.
SINGLE_SCORES = {
    "n_calls": +1,
    "n_calls_unopposed": +1,
    "breadth_1x": +1,
    "max_bin_frac": -1,        # concentrated depth means cross-mapping
    "depth_cv": -1,
}


def coverage_profile(starts: np.ndarray, lengths: np.ndarray, ref_len: int,
                     bin_size: int = 250) -> dict:
    """Per-base depth on the exogenous reference, reduced to breadth and shape."""
    depth = np.zeros(ref_len, dtype=np.int32)
    for s, ln in zip(starts, lengths):
        lo, hi = max(0, int(s)), min(ref_len, int(s) + int(ln))
        if hi > lo:
            depth[lo:hi] += 1
    n_bins = max(1, ref_len // bin_size)
    binned = np.array([depth[i * bin_size:(i + 1) * bin_size].sum()
                       for i in range(n_bins)], dtype=float)
    total = binned.sum()
    return {
        "breadth_1x": float((depth >= 1).mean()),
        "breadth_5x": float((depth >= 5).mean()),
        "depth_cv": float(binned.std() / binned.mean()) if binned.mean() > 0 else 0.0,
        "max_bin_frac": float(binned.max() / total) if total > 0 else 0.0,
    }


def one_sample(divergence: float, panel_mode: str, cohort_seed: int,
               sample_index: int, exo_depth: float, host_depth: float,
               n_herv_loci: int, host_filler_bp: int, read_len: int,
               seed_k: int) -> dict | None:
    panel, host, herv_spans, exo_strain, meta = simulate_world(
        divergence, n_herv_loci=n_herv_loci, host_filler_bp=host_filler_bp,
        seed=cohort_seed, sample_seed=cohort_seed * 1000 + sample_index)

    ref_len = len(panel.refs["EXO_REF"])
    restrict_panel(panel, panel_mode)

    reads = simulate_reads(exo_strain, host, herv_spans, meta,
                           exo_depth=exo_depth, host_depth=host_depth,
                           read_len=read_len, seed=cohort_seed * 1000 + sample_index + 7)
    index = PanelIndex(panel, k=seed_k)

    cols = {k: [] for k in ("AS", "gap", "softclip", "aligned_len", "in_ltr",
                            "mate_same_ref", "mate_exo", "tied", "start")}
    for t in reads:
        f = features_for_pair(t.read, t.mate, index, read_len)
        if f is None or f.best_cat != "EXO":
            continue
        cols["AS"].append(f.AS)
        cols["gap"].append(f.as_minus_xs)
        cols["softclip"].append(f.softclip_len)
        cols["aligned_len"].append(f.aligned_len)
        cols["in_ltr"].append(f.in_ltr)
        cols["mate_same_ref"].append(f.mate_same_ref)
        cols["mate_exo"].append(f.mate_best_cat_is_exo)
        cols["tied"].append(int(f.n_tied_top_categories > 1))
        cols["start"].append(f.best_ref_start)

    n = len(cols["AS"])
    if n == 0:
        # A sample with no calls at all is a legitimate observation, not a
        # failure: it is what a clean uninfected sample looks like.
        row = {k: 0.0 for k in SAMPLE_FEATURES}
        row["breadth_1x"] = row["breadth_5x"] = 0.0
        return row

    tied = np.array(cols["tied"])
    prof = coverage_profile(np.array(cols["start"]), np.array(cols["aligned_len"]),
                            ref_len)
    return {
        "n_calls": float(n),
        "n_calls_unopposed": float((tied == 0).sum()),
        "tie_frac": float(tied.mean()),
        "mean_AS": float(np.mean(cols["AS"])),
        "mean_gap": float(np.mean(cols["gap"])),
        "mean_softclip": float(np.mean(cols["softclip"])),
        "mean_aligned_len": float(np.mean(cols["aligned_len"])),
        "frac_in_ltr": float(np.mean(cols["in_ltr"])),
        "frac_mate_same_ref": float(np.mean(cols["mate_same_ref"])),
        "frac_mate_exo": float(np.mean(cols["mate_exo"])),
        **prof,
    }


def run_cohort(divergence: float, panel_mode: str, *, n_samples: int = 40,
               infected_frac: float = 0.5, load_lo: float = 0.5, load_hi: float = 8.0,
               host_depth: float = 8.0, n_herv_loci: int = 8,
               host_filler_bp: int = 120_000, read_len: int = 150,
               seed_k: int = 19, cohort_seed: int = 42, verbose: bool = True) -> dict:
    t0 = time.time()
    rng = np.random.default_rng(cohort_seed)
    n_inf = int(round(n_samples * infected_frac))
    labels = np.array([1] * n_inf + [0] * (n_samples - n_inf))
    # log-uniform viral load, so some infected samples sit near the detection limit
    loads = np.where(labels == 1,
                     np.exp(rng.uniform(np.log(load_lo), np.log(load_hi), n_samples)),
                     0.0)

    rows = []
    for i in range(n_samples):
        rows.append(one_sample(divergence, panel_mode, cohort_seed, i,
                               float(loads[i]), host_depth, n_herv_loci,
                               host_filler_bp, read_len, seed_k))
        if verbose and (i + 1) % 10 == 0:
            print(f"    {i+1}/{n_samples} samples ({time.time()-t0:.0f}s)")

    X = np.array([[r[c] for c in SAMPLE_FEATURES] for r in rows], dtype=float)
    y = labels
    out = {"divergence": divergence, "panel_mode": panel_mode,
           "n_samples": n_samples, "n_infected": int(y.sum()),
           "cohort_seed": cohort_seed,
           "viral_load_range": [round(float(loads[y == 1].min()), 3),
                                round(float(loads[y == 1].max()), 3)],
           "single_scores": {}, "model": {}}

    col = {c: i for i, c in enumerate(SAMPLE_FEATURES)}
    for name, sign in SINGLE_SCORES.items():
        s = sign * X[:, col[name]]
        out["single_scores"][name] = {
            "roc_auc": round(float(roc_auc_score(y, s)), 4),
            "fpr_at_95_sens": round(float(fpr_at_sensitivity(y, s, 0.95)), 4),
        }

    # Cross-validated over SAMPLES. Each sample carries its own redrawn loci, so
    # there is no within-sample unit to leak across.
    n_splits = int(min(5, np.bincount(y).min()))
    out["n_splits"] = n_splits
    if n_splits < 2:
        out["model"] = {"note": "too few samples per class to cross-validate"}
        out["seconds"] = round(time.time() - t0, 1)
        return out
    for model in ("logreg", "gbm"):
        oof = np.full(n_samples, np.nan)
        for tr, te in StratifiedKFold(n_splits, shuffle=True, random_state=0).split(X, y):
            clf = _make_model(model, 0)
            clf.fit(X[tr], y[tr])
            oof[te] = clf.predict_proba(X[te])[:, 1]
        out["model"][model] = {
            "roc_auc": round(float(roc_auc_score(y, oof)), 4),
            "fpr_at_95_sens": round(float(fpr_at_sensitivity(y, oof, 0.95)), 4),
        }
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", type=float, nargs="+", default=[0.02, 0.10])
    ap.add_argument("--panel-mode", nargs="+",
                    default=["poly_no_decoy", "no_decoy"], choices=list(PANEL_MODES))
    ap.add_argument("--n-samples", type=int, default=40)
    ap.add_argument("--load-lo", type=float, default=0.5)
    ap.add_argument("--load-hi", type=float, default=8.0)
    ap.add_argument("--host-depth", type=float, default=8.0)
    ap.add_argument("--n-herv-loci", type=int, default=8)
    ap.add_argument("--host-filler-bp", type=int, default=120_000)
    ap.add_argument("--cohort-seed", type=int, nargs="+", default=[42])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = []
    for mode in args.panel_mode:
        for d in args.divergence:
            for cs in args.cohort_seed:
                print(f"[sample-level  panel={mode}  divergence={d:.2f}  seed={cs}]")
                r = run_cohort(d, mode, n_samples=args.n_samples,
                               load_lo=args.load_lo, load_hi=args.load_hi,
                               host_depth=args.host_depth,
                               n_herv_loci=args.n_herv_loci,
                               host_filler_bp=args.host_filler_bp,
                               cohort_seed=cs)
                print(json.dumps({k: r[k] for k in
                                  ("single_scores", "model", "viral_load_range",
                                   "seconds")}, indent=2))
                results.append(r)
    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            json.dump(results, fh, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
