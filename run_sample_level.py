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
    max_bin_frac         uniformity: how much depth sits in the single worst bin
    depth_cv             uniformity: coefficient of variation across bins
    prop_variable_sites  intra-host variation over the pileup; see below
    mean_maf             the same, as a mean minor allele fraction
    nucleotide_diversity the same, as per-site pi
    learned model        logistic regression and a GBM over all of the above

  Viral load VARIES across infected samples, drawn log-uniformly, because the
  discriminating power of a read count collapses near the detection limit and
  that is the regime the question matters in.

THE PILEUP STATISTICS ARE HERE FOR A SPECIFIC REASON
  Every rule measured in this repository so far is a COMPETITION test, and F024
  showed that a competition test is worthless when the read's true source is
  absent from the panel. A pileup statistic is computed from the reads that
  mapped, full stop, so it needs no competitor and could work in exactly the
  regime where the alignment-level rules fail. That is why the prior-art check
  flagged it (F049) as the outstanding item that could change the project's
  conclusion rather than its wording.

  Its direction had to be reversed against the literature, and the reversal is
  the finding. See the comment on SINGLE_SCORES.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from endoexo.align import PanelIndex, features_for_pair, revcomp
from endoexo.evaluate import _make_model, fpr_at_sensitivity
from endoexo.simulate import simulate_reads, simulate_world
from run_slice import PANEL_MODES, restrict_panel

SAMPLE_FEATURES = [
    "n_calls", "n_calls_unopposed", "tie_frac",
    "breadth_1x", "breadth_5x", "depth_cv", "max_bin_frac",
    "mean_AS", "mean_gap", "mean_softclip", "mean_aligned_len",
    "frac_in_ltr", "frac_mate_same_ref", "frac_mate_exo",
    # Hayward's discriminator. See diversity_stats.
    "prop_variable_sites", "mean_maf", "nucleotide_diversity",
]

DIVERSITY_FEATURES = ("prop_variable_sites", "mean_maf", "nucleotide_diversity")

# Single statistics evaluated on their own, and the sign that makes "higher is
# more likely infected" true for each.
SINGLE_SCORES = {
    "n_calls": +1,
    "n_calls_unopposed": +1,
    "breadth_1x": +1,
    "max_bin_frac": -1,        # concentrated depth means cross-mapping
    "depth_cv": -1,
    # SIGN REVERSED AGAINST THE LITERATURE, 2026-09-15, and the reversal is the
    # finding rather than a correction. Hayward's direction -- endogenous reads
    # high variation, exogenous low -- was measured with the literature's sign
    # first and gave ROC AUC 0.111, i.e. 0.889 reversed. The direction does not
    # transfer because the QUESTION is different. Hayward ask whether a set of
    # reads that all come from one candidate element is endogenous or exogenous,
    # so coalescence age is what the pileup reflects. Here the pileup is a
    # MIXTURE: an infected sample contributes both a viral strain, which differs
    # systematically from the reference at its own sites, and endogenous
    # cross-mappings, which differ at theirs. Two populations carry more variation
    # than one. So in this task the statistic detects a mixture, not an age, and
    # higher variation means infected.
    "prop_variable_sites": +1,
    "mean_maf": +1,
    "nucleotide_diversity": +1,
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


def diversity_stats(counts: np.ndarray, min_depth: int = 5,
                    maf_threshold: float = 0.05) -> dict:
    """Intra-host genetic variation over the exogenous-reference pileup.

    This is Hayward et al.'s discriminator, and the reason it is worth testing
    here is that it needs NO COMPETITOR IN THE PANEL. Every rule this repository
    has measured so far is a competition test, and F024 showed those are
    worthless when the read's true source is absent from the panel. A pileup
    statistic is computed from the reads that mapped, full stop, so it is the one
    published discriminator that could work in exactly the regime where the
    alignment-level rules fail.

    The direction: an endogenous family is a mixture of loci that diverged from
    one another over millions of years, so at a site where one locus differs the
    pileup shows a minor allele near 1/k for k contributing loci. An exogenous
    infection is one recently coalesced quasispecies, so its systematic
    differences from the reference are FIXED rather than variable, and the only
    residual variation is replication error and sequencing error. Endogenous
    therefore reads as high variation and exogenous as low.

    `counts` is (reference length, 4) of observed base counts. Sites below
    `min_depth` are excluded rather than imputed: a site seen twice carries no
    usable allele frequency.
    """
    depth = counts.sum(axis=1)
    keep = depth >= min_depth
    n_sites = int(keep.sum())
    if n_sites == 0:
        return {"prop_variable_sites": 0.0, "mean_maf": 0.0,
                "nucleotide_diversity": 0.0, "n_pileup_sites": 0}

    c = counts[keep].astype(float)
    n = depth[keep].astype(float)
    p = c / n[:, None]
    minor = 1.0 - p.max(axis=1)

    # Unbiased nucleotide diversity per site: n/(n-1) * (1 - sum p_i^2). The
    # correction matters here because depths are small by construction.
    homozygosity = (p ** 2).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        pi_site = np.where(n > 1, (n / (n - 1.0)) * (1.0 - homozygosity), 0.0)

    return {
        "prop_variable_sites": float((minor >= maf_threshold).mean()),
        "mean_maf": float(minor.mean()),
        "nucleotide_diversity": float(pi_site.mean()),
        "n_pileup_sites": n_sites,
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
    pileup = np.zeros((ref_len, 4), dtype=np.int32)
    for t in reads:
        f = features_for_pair(t.read, t.mate, index, read_len)
        if f is None or f.best_cat != "EXO":
            continue
        # Accumulate the pileup on the exogenous reference. The bases are
        # recovered from the stored strand and read offset, so nothing is
        # re-aligned; ungapped alignment makes the reference-to-base mapping one
        # to one, so no CIGAR walk is needed.
        seq = revcomp(t.read) if f.best_reverse else t.read
        bases = seq[f.best_read_start:f.best_read_start + f.aligned_len]
        if len(bases):
            lo = f.best_ref_start
            idx = np.arange(lo, lo + len(bases))
            inside = (idx >= 0) & (idx < ref_len)
            np.add.at(pileup, (idx[inside], bases[inside].astype(np.int64)), 1)
        cols["AS"].append(f.AS)
        cols["gap"].append(f.as_minus_xs)
        cols["softclip"].append(f.softclip_len)
        cols["aligned_len"].append(f.aligned_len)
        cols["in_ltr"].append(f.in_ltr)
        cols["mate_same_ref"].append(f.mate_same_ref)
        cols["mate_exo"].append(f.mate_best_cat_is_exo)
        cols["tied"].append(int(f.n_tied_top_categories > 1))
        cols["start"].append(f.best_ref_start)

    div = diversity_stats(pileup)
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
        **{k: v for k, v in div.items() if k in SAMPLE_FEATURES},
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

    # Is the diversity statistic independent signal, or just a read-count proxy?
    # If it correlates with n_calls as strongly as with the label, it is adding
    # nothing that counting does not already give.
    from scipy.stats import spearmanr
    out["diversity_diagnostics"] = {}
    for name in ("prop_variable_sites", "mean_maf", "nucleotide_diversity"):
        v = X[:, col[name]]
        rho_calls = spearmanr(v, X[:, col["n_calls"]]).statistic
        rho_load = spearmanr(v[y == 1], loads[y == 1]).statistic if (y == 1).sum() > 2 else float("nan")
        out["diversity_diagnostics"][name] = {
            "spearman_vs_n_calls": round(float(rho_calls), 4),
            "spearman_vs_load_within_infected": round(float(rho_load), 4),
        }

    # Cross-validated over SAMPLES. Each sample carries its own redrawn loci, so
    # there is no within-sample unit to leak across.
    n_splits = int(min(5, np.bincount(y).min()))
    out["n_splits"] = n_splits
    if n_splits < 2:
        out["model"] = {"note": "too few samples per class to cross-validate"}
        out["seconds"] = round(time.time() - t0, 1)
        return out
    # Ablate the pileup statistics, so any model improvement they bring is
    # ATTRIBUTABLE. The earlier comparison against a pre-F028 run could not
    # separate the new features from the rate-normalisation fix.
    featsets = {
        "all": list(range(len(SAMPLE_FEATURES))),
        "no_diversity": [i for i, c in enumerate(SAMPLE_FEATURES)
                         if c not in DIVERSITY_FEATURES],
    }
    for featset, cols in featsets.items():
        for model in ("logreg", "gbm"):
            oof = np.full(n_samples, np.nan)
            for tr, te in StratifiedKFold(n_splits, shuffle=True, random_state=0).split(X, y):
                clf = _make_model(model, 0)
                clf.fit(X[np.ix_(tr, cols)], y[tr])
                oof[te] = clf.predict_proba(X[np.ix_(te, cols)])[:, 1]
            out["model"][f"{model}/{featset}"] = {
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
