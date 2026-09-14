#!/usr/bin/env python3
"""Isolate the read-level-CV leakage seen at low divergence, and test its mechanism.

The hypothesis: a read-level split lets the model memorise WHERE on the
exogenous reference a given endogenous locus lands, because reads from one
locus pile onto overlapping reference positions. If that is the mechanism, then
removing the position features should shrink the grouped-vs-read-level gap.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from endoexo.align import FEATURE_COLUMNS, PanelIndex, features_for_pair
from endoexo.evaluate import cv_scores, evaluate_oof
from endoexo.simulate import simulate_reads, simulate_world

POSITION_FEATURES = ("dist_to_terminus", "in_ltr")


def build_matrix(divergence, n_herv_loci, host_filler_bp, exo_depth, host_depth, seed, read_len=150):
    panel, host, herv_spans, exo_strain = simulate_world(
        divergence, n_herv_loci=n_herv_loci, host_filler_bp=host_filler_bp, seed=seed)
    for ref_id in ("HOST", "HERV_CONSENSUS"):          # viral_only panel
        del panel.refs[ref_id]; del panel.categories[ref_id]; panel.ltr_spans.pop(ref_id, None)
    reads = simulate_reads(exo_strain, host, herv_spans, exo_depth=exo_depth,
                           host_depth=host_depth, read_len=read_len, seed=seed + 1)
    index = PanelIndex(panel, k=19)
    rows, y, groups = [], [], []
    for t in reads:
        f = features_for_pair(t.read, t.mate, index, read_len)
        if f is None or f.best_cat != "EXO":
            continue
        rows.append([getattr(f, c) for c in FEATURE_COLUMNS])
        y.append(int(t.source == "EXO"))
        groups.append(t.locus if t.source != "EXO" else f"EXO_BIN{t.frag_start // 1000}")
    return np.array(rows, float), np.array(y), np.array(groups)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", type=float, default=0.02)
    ap.add_argument("--n-herv-loci", type=int, default=30)
    ap.add_argument("--host-filler-bp", type=int, default=200_000)
    ap.add_argument("--exo-depth", type=float, default=20.0)
    ap.add_argument("--host-depth", type=float, default=8.0)
    ap.add_argument("--seeds", type=int, nargs="+", default=[42, 142, 242])
    ap.add_argument("--out", default="results_leakage_probe.json")
    args = ap.parse_args()

    keep_all = list(range(len(FEATURE_COLUMNS)))
    keep_nopos = [i for i, c in enumerate(FEATURE_COLUMNS) if c not in POSITION_FEATURES]

    records = []
    for seed in args.seeds:
        X, y, groups = build_matrix(args.divergence, args.n_herv_loci,
                                    args.host_filler_bp, args.exo_depth,
                                    args.host_depth, seed)
        print(f"seed {seed}: calls={len(y)} true={int(y.sum())} false={int((1-y).sum())} "
              f"groups={len(np.unique(groups))}")
        for featset, cols in (("all", keep_all), ("no_position", keep_nopos)):
            for model in ("logreg", "gbm"):
                row = {"seed": seed, "featset": featset, "model": model}
                for scheme in ("grouped", "read_level"):
                    oof = cv_scores(X[:, cols], y, groups, scheme=scheme, model=model)
                    m = evaluate_oof(y, oof)
                    row[f"{scheme}_pr_auc"] = m["pr_auc"]
                    row[f"{scheme}_fpr"] = m["fpr_at_sens"]
                row["leak_pr_auc"] = row["read_level_pr_auc"] - row["grouped_pr_auc"]
                records.append(row)
                print(f"  {featset:11s} {model:6s} "
                      f"grouped PR {row['grouped_pr_auc']:.3f} / "
                      f"read-level PR {row['read_level_pr_auc']:.3f} / "
                      f"leak +{row['leak_pr_auc']:.3f}")

    with open(args.out, "w", encoding="ascii") as fh:
        json.dump(records, fh, indent=2)

    print("\n=== mean over seeds ===")
    for featset in ("all", "no_position"):
        for model in ("logreg", "gbm"):
            sel = [r for r in records if r["featset"] == featset and r["model"] == model]
            g = np.mean([r["grouped_pr_auc"] for r in sel])
            rl = np.mean([r["read_level_pr_auc"] for r in sel])
            print(f"{featset:11s} {model:6s} grouped {g:.3f}  read_level {rl:.3f}  "
                  f"leak +{rl-g:.3f}")


if __name__ == "__main__":
    main()
