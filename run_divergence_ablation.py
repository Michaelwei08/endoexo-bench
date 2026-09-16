#!/usr/bin/env python3
"""Measure what the distributional divergence model buys, rather than asserting it.

The first version of this benchmark used a POINT divergence -- one value for the
endogenous family, one for the strain -- and reported perfect separation at every
divergence above 0.10. That was an artefact: with divergence fixed and site rates
independent, every 150 bp read window carries nearly the family average, so
alignment score alone separates the classes.

This script quantifies the difference by the only measure that matters for
difficulty: how often an ENDOGENOUS read window falls inside the range of
EXOGENOUS read windows. If that overlap is near zero the benchmark is trivial
regardless of what any model reports.

Existed as an ad-hoc check first. Promoted to a script because the paper quotes
its numbers, and every number the paper quotes has to trace to a file.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from endoexo.simulate import simulate_reads, simulate_world

POINT_MODEL = dict(rate_shape=1e6, rate_jitter_shape=0.0,
                   sigma_locus=0.0, sigma_strain=0.0)
DISTRIBUTIONAL = dict()          # the module defaults


def measure(divergence: float, *, n_herv_loci: int, host_filler_bp: int,
            exo_depth: float, host_depth: float, seed: int, **kw) -> dict:
    panel, host, spans, strain, meta = simulate_world(
        divergence, n_herv_loci=n_herv_loci, host_filler_bp=host_filler_bp,
        seed=seed, **kw)
    reads = simulate_reads(strain, host, spans, meta, exo_depth=exo_depth,
                           host_depth=host_depth, seed=seed + 1)
    exo = np.array([r.local_div for r in reads if r.source == "EXO"])
    herv = np.array([r.local_div for r in reads if r.source == "HERV"])
    herv = herv[~np.isnan(herv)]
    divs = np.array(list(meta["realized_locus_divergence"].values()))
    if exo.size == 0 or herv.size == 0:
        return {"divergence": divergence, "note": "a class produced no reads"}
    ceiling = float(np.percentile(exo, 95))
    return {
        "divergence": divergence,
        "locus_divergence_min": round(float(divs.min()), 4),
        "locus_divergence_median": round(float(np.median(divs)), 4),
        "locus_divergence_max": round(float(divs.max()), 4),
        "exo_window_divergence_p50": round(float(np.median(exo)), 4),
        "exo_window_divergence_p95": round(ceiling, 4),
        "herv_window_divergence_p5": round(float(np.percentile(herv, 5)), 4),
        "herv_window_divergence_p50": round(float(np.median(herv)), 4),
        # THE number: endogenous windows inside the exogenous range
        "overlap_frac": round(float((herv <= ceiling).mean()), 4),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", type=float, nargs="+",
                    default=[0.05, 0.10, 0.20, 0.30])
    ap.add_argument("--n-herv-loci", type=int, default=12)
    ap.add_argument("--host-filler-bp", type=int, default=200_000)
    ap.add_argument("--exo-depth", type=float, default=20.0)
    ap.add_argument("--host-depth", type=float, default=8.0)
    ap.add_argument("--seed", type=int, nargs="+", default=[42])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = []
    print(f"{'model':16s} {'seed':>5s} {'div':>5s} {'locus range':>16s} "
          f"{'HERV p5':>9s} {'overlap':>8s}")
    for name, kw in (("distributional", DISTRIBUTIONAL), ("point", POINT_MODEL)):
        for seed in args.seed:
            for d in args.divergence:
                r = measure(d, n_herv_loci=args.n_herv_loci,
                            host_filler_bp=args.host_filler_bp,
                            exo_depth=args.exo_depth, host_depth=args.host_depth,
                            seed=seed, **kw)
                r["model"] = name
                r["seed"] = seed
                results.append(r)
                if "overlap_frac" in r:
                    print(f"{name:16s} {seed:5d} {d:5.2f} "
                          f"{r['locus_divergence_min']:.3f}-{r['locus_divergence_max']:.3f}"
                          f"{'':>5s} {r['herv_window_divergence_p5']:9.4f} "
                          f"{100*r['overlap_frac']:7.1f}%")

    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            json.dump(results, fh, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
