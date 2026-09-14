#!/usr/bin/env python3
"""One divergence point of the benchmark, end to end: simulate -> competitive
align -> select the reads a detection pipeline would count -> compare the
production rules against a learned discriminator under grouped and read-level CV.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np

from endoexo.align import (FEATURE_COLUMNS, MATE_FEATURES, PanelIndex,
                           features_for_pair)
from endoexo.evaluate import (confusion_rates, cv_scores, evaluate_oof,
                              rule_baselines)
from endoexo.simulate import simulate_reads, simulate_world

PANEL_MODES = {
    # host assembly CONTAINS the endogenous loci (they are reference insertions)
    "full": ("EXO_REF", "HOST", "HERV_CONSENSUS"),
    "no_decoy": ("EXO_REF", "HOST"),
    "viral_only": ("EXO_REF",),
    # host assembly LACKS them (insertionally polymorphic, non-reference loci)
    "poly_no_decoy": ("EXO_REF", "HOST_NOLOCI"),
    "poly_full": ("EXO_REF", "HOST_NOLOCI", "HERV_CONSENSUS"),
    # ... plus a COMPLETE catalogue of the polymorphic loci themselves, which is
    # the upper bound on what panel design can buy. HERV_LOCI is a wildcard for
    # every HERV_LOCUS_nn reference.
    "poly_catalogue": ("EXO_REF", "HOST_NOLOCI", "HERV_LOCI"),
}


def restrict_panel(panel, panel_mode: str) -> None:
    """Drop every reference the mode does not select, in place."""
    keep = PANEL_MODES[panel_mode]
    wildcard = "HERV_LOCI" in keep
    for ref_id in list(panel.refs):
        if ref_id in keep or (wildcard and ref_id.startswith("HERV_LOCUS_")):
            continue
        del panel.refs[ref_id]
        del panel.categories[ref_id]
        panel.ltr_spans.pop(ref_id, None)


def run_point(divergence: float, *, panel_mode: str = "full", seed: int = 42,
              exo_depth: float = 20.0, host_depth: float = 10.0,
              n_herv_loci: int = 20, host_filler_bp: int = 330_000,
              read_len: int = 150, seed_k: int = 19, rate_shape: float = 0.5,
              rate_block_len: int = 400, rate_jitter_shape: float = 4.0,
              sigma_locus: float = 0.7, sigma_strain: float = 0.4,
              verbose: bool = True) -> dict:
    t0 = time.time()
    panel, host, herv_spans, exo_strain, meta = simulate_world(
        divergence, n_herv_loci=n_herv_loci, host_filler_bp=host_filler_bp,
        rate_shape=rate_shape, rate_block_len=rate_block_len,
        rate_jitter_shape=rate_jitter_shape, sigma_locus=sigma_locus,
        sigma_strain=sigma_strain, seed=seed)

    restrict_panel(panel, panel_mode)

    reads = simulate_reads(exo_strain, host, herv_spans, meta, exo_depth=exo_depth,
                           host_depth=host_depth, read_len=read_len, seed=seed + 1)
    index = PanelIndex(panel, k=seed_k)
    if verbose:
        print(f"  panel={panel_mode} refs={list(panel.refs)} reads={len(reads)} "
              f"index_kmers={len(index.index)} "
              f"locus_div={meta['locus_divergence_min']:.3f}/"
              f"{meta['locus_divergence_median']:.3f}/{meta['locus_divergence_max']:.3f} "
              f"strain_div={meta['realized_strain_divergence']:.3f} "
              f"({time.time()-t0:.1f}s)")

    rows, truth_src, groups, local_div, ambiguous = [], [], [], [], []
    n_exo_reads = sum(1 for t in reads if t.source == "EXO")
    n_exo_lost_to_host = 0
    for t in reads:
        f = features_for_pair(t.read, t.mate, index, read_len)
        if f is None:
            continue
        if f.best_cat != "EXO":
            # not a call: the pipeline never counts it. For a genuine exogenous
            # read that is a sensitivity loss nobody sees, so it is counted.
            if t.source == "EXO":
                n_exo_lost_to_host += 1
            continue
        rows.append([getattr(f, c) for c in FEATURE_COLUMNS])
        truth_src.append(t.source)
        local_div.append(t.local_div)
        ambiguous.append(f.n_tied_top_categories > 1)
        # grouping unit: HERV/HOST locus, or a positional bin along the provirus
        groups.append(t.locus if t.source != "EXO"
                      else f"EXO_BIN{t.frag_start // 1000}")

    div_meta = {
        "realized_strain_divergence": round(meta["realized_strain_divergence"], 4),
        "locus_divergence_min": round(meta["locus_divergence_min"], 4),
        "locus_divergence_median": round(meta["locus_divergence_median"], 4),
        "locus_divergence_max": round(meta["locus_divergence_max"], 4),
        "rate_shape": meta["rate_shape"],
        "sigma_locus": meta["sigma_locus"],
        "sigma_strain": meta["sigma_strain"],
    }

    if not rows:
        return {"divergence": divergence, "panel_mode": panel_mode, "n_calls": 0,
                "divergence_model": div_meta,
                "note": "no read was assigned to the exogenous reference"}

    X = np.array(rows, dtype=float)
    y = np.array([s == "EXO" for s in truth_src], dtype=int)
    groups = np.array(groups)
    local_div = np.array(local_div, dtype=float)
    feats = {c: X[:, i] for i, c in enumerate(FEATURE_COLUMNS)}

    out = {
        "divergence": divergence,
        "panel_mode": panel_mode,
        "divergence_model": div_meta,
        "n_calls": int(len(y)),
        "n_true_exo": int(y.sum()),
        "n_false_exo": int((1 - y).sum()),
        "false_sources": {s: int(sum(1 for u in truth_src if u == s))
                          for s in sorted(set(truth_src)) if s != "EXO"},
        "rules": {}, "model": {},
    }

    # THE TWO FAILURE MECHANISMS, separated. A tie means the panel holds an equally
    # good explanation, which is free to detect; an unopposed win means the read
    # beat everything present and no competition-based rule can touch it.
    amb = np.array(ambiguous, dtype=bool)
    out["n_exo_reads_simulated"] = int(n_exo_reads)
    out["n_exo_lost_to_host_assignment"] = int(n_exo_lost_to_host)
    out["tie_typology"] = {
        "false_tied_frac": round(float(amb[y == 0].mean()), 4) if (y == 0).any() else None,
        "true_tied_frac": round(float(amb[y == 1].mean()), 4) if (y == 1).any() else None,
        "n_false_unopposed": int((~amb & (y == 0)).sum()),
        "n_true_unopposed": int((~amb & (y == 1)).sum()),
    }
    # What does the three-bin rule buy? Policy A awards a tie to the virus; policy
    # B routes it to an ambiguous bin and counts only unopposed wins. Sensitivity
    # is against ALL simulated exogenous reads, which is the denominator an assay
    # actually has.
    out["policies"] = {}
    for name, keep in (("award_ties_to_exo", np.ones(len(y), dtype=bool)),
                       ("three_bin_ambiguous", ~amb)):
        tp = int(((y == 1) & keep).sum())
        fp = int(((y == 0) & keep).sum())
        out["policies"][name] = {
            "calls": tp + fp,
            "false_calls": fp,
            "false_call_frac": round(fp / (tp + fp), 4) if (tp + fp) else None,
            "sensitivity_vs_all_exo_reads": round(tp / n_exo_reads, 4) if n_exo_reads else None,
        }

    # How hard were the reads that got through? local_div is unobservable at
    # inference time and is reported as a diagnostic only.
    fp_div = local_div[(y == 0) & ~np.isnan(local_div)]
    if fp_div.size:
        out["false_positive_local_divergence"] = {
            "min": round(float(fp_div.min()), 4),
            "p25": round(float(np.percentile(fp_div, 25)), 4),
            "median": round(float(np.median(fp_div)), 4),
            "max": round(float(fp_div.max()), 4),
            "frac_below_0.05": round(float((fp_div < 0.05).mean()), 4),
        }

    for name, pred in rule_baselines(feats).items():
        sens, fpr = confusion_rates(y, pred)
        out["rules"][name] = {"sensitivity": round(sens, 4), "fpr": round(fpr, 4)}

    # Feature sets. `no_mate` ablates everything the MATE contributes, which is the
    # only read-pair-level route past the conserved-window ceiling of F016. If the
    # mate is doing the work, dropping it should hurt most in the hardest cells.
    col_index = {c: i for i, c in enumerate(FEATURE_COLUMNS)}
    featsets = {
        "all": list(range(len(FEATURE_COLUMNS))),
        "no_mate": [col_index[c] for c in FEATURE_COLUMNS if c not in MATE_FEATURES],
    }

    if out["n_false_exo"] > 0 and out["n_true_exo"] > 0:
        n_groups = len(np.unique(groups))
        for featset, cols in featsets.items():
            for model in ("logreg", "gbm"):
                # read-level CV is run only on the full feature set; the leakage
                # contrast is established and repeating it per ablation is waste
                schemes = ("grouped", "read_level") if featset == "all" else ("grouped",)
                for scheme in schemes:
                    key = f"{model}/{scheme}/{featset}"
                    if scheme == "grouped" and n_groups < 5:
                        out["model"][key] = {"note": f"only {n_groups} groups"}
                        continue
                    oof = cv_scores(X[:, cols], y, groups, scheme=scheme, model=model)
                    out["model"][key] = {
                        k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in evaluate_oof(y, oof).items()}
                    out["model"][key]["prevalence"] = round(float(y.mean()), 4)
    out["n_groups"] = int(len(np.unique(groups)))
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", type=float, nargs="+", default=[0.10])
    ap.add_argument("--panel-mode", nargs="+", default=["full"], choices=list(PANEL_MODES))
    ap.add_argument("--exo-depth", type=float, default=20.0)
    ap.add_argument("--host-depth", type=float, default=10.0)
    ap.add_argument("--host-filler-bp", type=int, default=330_000)
    ap.add_argument("--n-herv-loci", type=int, default=20)
    ap.add_argument("--seed-k", type=int, default=19)
    ap.add_argument("--rate-shape", type=float, default=0.5,
                    help="gamma shape for REGIONAL across-site rate heterogeneity; "
                         "lower is more heterogeneous")
    ap.add_argument("--rate-block-len", type=int, default=400,
                    help="length of a constant-rate region; the autocorrelation scale")
    ap.add_argument("--rate-jitter", type=float, default=4.0,
                    help="gamma shape for per-site jitter on top of the regional "
                         "rate; 0 disables it")
    ap.add_argument("--sigma-locus", type=float, default=0.7,
                    help="lognormal sigma for per-locus divergence; 0 fixes all loci")
    ap.add_argument("--sigma-strain", type=float, default=0.4)
    ap.add_argument("--seed", type=int, nargs="+", default=[42])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = []
    for mode in args.panel_mode:
        for d in args.divergence:
          for sd in args.seed:
            print(f"[divergence={d:.2f} panel={mode} seed={sd}]")
            r = run_point(d, panel_mode=mode, seed=sd, exo_depth=args.exo_depth,
                          host_depth=args.host_depth, n_herv_loci=args.n_herv_loci,
                          host_filler_bp=args.host_filler_bp, seed_k=args.seed_k,
                          rate_shape=args.rate_shape,
                          rate_block_len=args.rate_block_len,
                          rate_jitter_shape=args.rate_jitter,
                          sigma_locus=args.sigma_locus,
                          sigma_strain=args.sigma_strain)
            r["seed"] = sd
            print(json.dumps(r, indent=2))
            results.append(r)
    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            json.dump(results, fh, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
