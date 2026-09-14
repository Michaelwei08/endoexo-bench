#!/usr/bin/env python3
"""Build a reference panel and a read pool in which exogenous/endogenous homology
is controlled BY CONSTRUCTION, so every read carries its true source.

Geometry follows HTLV-1: a provirus of two identical 755 bp LTRs flanking an
internal region. The endogenous family is the same ancestral provirus carried
forward at substitution divergence and inserted at several host loci, which is
what makes a read from an endogenous locus able to align to the exogenous
reference at all.

DIVERGENCE IS A DISTRIBUTION, NOT A POINT (revised 2026-09-14)
  Three sources of spread, because the first version of this file had none and
  that made the high-divergence regime separable by alignment score alone:

  1 ACROSS LOCI. Each endogenous insertion gets its own divergence drawn
    lognormally, because insertions in a real genome are of different ages.
  2 ACROSS SITES. Relative substitution rates are gamma distributed and SHARED
    between lineages, because rate is a property of a site functional
    constraint, not of the lineage. This is the load-bearing one: under gamma
    rates an old endogenous locus still contains slowly-evolving windows, and a
    150 bp read from such a window is nearly indistinguishable from a read off
    the exogenous strain no matter how diverged the locus is overall.
  3 WITHIN THE STRAIN. The viral strain in the sample is itself drawn from a
    distribution rather than fixed, because the sample virus is never the
    reference.

  Every realized divergence is MEASURED from the generated sequences and
  reported, never assumed from the parameter that produced it.

Nucleotides are encoded A=0 C=1 G=2 T=3 as uint8 throughout.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

LTR_LEN = 755
INTERNAL_LEN = 7000
SATURATION = 0.75          # four bases: p above this is indistinguishable from random


def gamma_site_rates(length: int, shape: float, rng: np.random.Generator,
                     block_len: int = 400, jitter_shape: float = 4.0) -> np.ndarray:
    """Relative substitution rates across sites, gamma distributed with mean 1 and
    AUTOCORRELATED along the sequence.

    Low `shape` means strong heterogeneity: many nearly-invariant regions plus a
    few fast ones, which is the standard model for functional constraint.

    The autocorrelation is load-bearing and was missing from the first version.
    Functional constraint applies to REGIONS -- a gene, a domain, an LTR
    subregion -- not to independent sites. With iid site rates the mean of 150
    of them concentrates on 1 (coefficient of variation 1/sqrt(150*shape), about
    0.115 at shape 0.5), so every 150 bp read window carried nearly the family
    average divergence and alignment score alone separated the classes. Drawing
    the rate at block resolution instead lets a read land wholly inside a
    conserved region, which is what makes a read off an old endogenous locus
    genuinely indistinguishable from a read off the exogenous strain.

    Set `jitter_shape` to 0 and `shape` very high to recover flat rates.
    """
    n_blocks = int(np.ceil(length / block_len))
    regional = np.repeat(rng.gamma(shape, 1.0 / shape, size=n_blocks), block_len)[:length]
    if jitter_shape <= 0:
        return regional
    fine = rng.gamma(jitter_shape, 1.0 / jitter_shape, size=length)
    return regional * fine


def mutate(seq: np.ndarray, divergence: float, rng: np.random.Generator,
           rates: np.ndarray | None = None) -> np.ndarray:
    """Substitute sites to a DIFFERENT base, at per-site probability
    `divergence * rates[i]`, so the expected overall divergence is `divergence`."""
    if divergence <= 0:
        return seq.copy()
    if rates is None:
        rates = np.ones(len(seq))
    p = np.clip(divergence * rates, 0.0, SATURATION)
    mask = rng.random(len(seq)) < p
    out = seq.copy()
    n = int(mask.sum())
    if n:
        # +1..+3 mod 4 guarantees the new base differs from the old one
        out[mask] = (out[mask] + rng.integers(1, 4, size=n, dtype=np.uint8)) % 4
    return out


def lognormal_around(median: float, sigma: float, rng: np.random.Generator,
                     size: int | None = None):
    """Lognormal with the given median on the natural scale."""
    return rng.lognormal(np.log(median), sigma, size=size)


def _make_locus(exo_ref: np.ndarray, consensus: np.ndarray, family_mask: np.ndarray,
                rates: np.ndarray, target_div: float,
                rng: np.random.Generator) -> np.ndarray:
    """One endogenous insertion at a TARGET divergence from the exogenous reference.

    Divergence is drawn per locus and applied directly, not added on top of the
    family consensus. In the first version every locus inherited the whole family
    divergence as a hard floor, so no locus could be younger than the family and
    the across-loci spread was an artefact of the top-up term alone. Here a locus
    below the family divergence REVERTS a random subset of the family-specific
    substitutions -- which is what a young element looks like -- and a locus above
    it acquires private ones.
    """
    element = consensus.copy()
    family_div = len(family_mask) / len(exo_ref)
    target_div = float(np.clip(target_div, 0.0, SATURATION))

    if target_div < family_div and len(family_mask):
        keep_frac = target_div / family_div if family_div > 0 else 0.0
        n_revert = int(round((1.0 - keep_frac) * len(family_mask)))
        if n_revert > 0:
            revert = rng.choice(family_mask, size=n_revert, replace=False)
            element[revert] = exo_ref[revert]
    elif target_div > family_div:
        element = mutate(element, target_div - family_div, rng, rates)
    return element


def build_provirus(rng: np.random.Generator) -> np.ndarray:
    ltr = rng.integers(0, 4, LTR_LEN, dtype=np.uint8)
    internal = rng.integers(0, 4, INTERNAL_LEN, dtype=np.uint8)
    return np.concatenate([ltr, internal, ltr])


@dataclass
class Panel:
    """References available to the aligner, plus where the LTRs sit."""

    refs: dict[str, np.ndarray] = field(default_factory=dict)
    categories: dict[str, str] = field(default_factory=dict)
    ltr_spans: dict[str, list[tuple[int, int]]] = field(default_factory=dict)

    def add(self, ref_id: str, seq: np.ndarray, category: str,
            ltr_spans: list[tuple[int, int]] | None = None) -> None:
        self.refs[ref_id] = seq
        self.categories[ref_id] = category
        self.ltr_spans[ref_id] = ltr_spans or []


@dataclass
class Truth:
    """One simulated fragment and the source it actually came from.

    `local_div` is the realized divergence from the exogenous reference over the
    window this read came from. It is the ground-truth difficulty of the read and
    is a DIAGNOSTIC, never a feature: nothing observable at inference time
    reveals it.
    """

    read: np.ndarray
    mate: np.ndarray
    source: str           # "EXO" or "HERV" or "HOST"
    locus: str            # grouping unit for grouped CV
    frag_start: int
    local_div: float = float("nan")


def simulate_world(
    herv_divergence: float,
    *,
    exo_strain_divergence: float = 0.02,
    n_herv_loci: int = 20,
    host_filler_bp: int = 330_000,
    rate_shape: float = 0.5,
    rate_block_len: int = 400,
    rate_jitter_shape: float = 4.0,
    sigma_locus: float = 0.7,
    sigma_strain: float = 0.4,
    seed: int = 42,
):
    """Return (panel, host, herv_locus_spans, exo_strain, meta).

    `herv_divergence` is the MEDIAN of the per-locus distribution, not a fixed
    value. `meta` carries realized divergences and the per-locus difference masks
    the read simulator needs in order to report per-read difficulty.
    """
    rng = np.random.default_rng(seed)
    exo_ref = build_provirus(rng)            # the reference defines the coordinates
    site_rates = gamma_site_rates(len(exo_ref), rate_shape, rng,
                                  block_len=rate_block_len,
                                  jitter_shape=rate_jitter_shape)

    d_strain = float(lognormal_around(exo_strain_divergence, sigma_strain, rng))
    exo_strain = mutate(exo_ref, d_strain, rng, site_rates)

    herv_consensus = mutate(exo_ref, herv_divergence, rng, site_rates)

    pieces: list[np.ndarray] = []
    herv_spans: dict[str, tuple[int, int]] = {}
    locus_diff: dict[str, np.ndarray] = {}
    locus_div: dict[str, float] = {}
    filler_chunk = host_filler_bp // (n_herv_loci + 1)
    cursor = 0
    family_mask = np.flatnonzero(herv_consensus != exo_ref)
    for i in range(n_herv_loci):
        pieces.append(rng.integers(0, 4, filler_chunk, dtype=np.uint8))
        cursor += filler_chunk
        element = _make_locus(exo_ref, herv_consensus, family_mask, site_rates,
                              float(lognormal_around(herv_divergence, sigma_locus, rng)),
                              rng)
        name = f"HERV_L{i:02d}"
        herv_spans[name] = (cursor, cursor + len(element))
        diff = (element != exo_ref)
        locus_diff[name] = diff
        locus_div[name] = float(diff.mean())
        pieces.append(element)
        cursor += len(element)
    pieces.append(rng.integers(0, 4, filler_chunk, dtype=np.uint8))
    host = np.concatenate(pieces)
    # The same background with every endogenous element deleted. A locus that is
    # insertionally polymorphic is carried by the sample but absent from the
    # reference assembly, which is the only situation in which an explicit
    # endogenous decoy can earn its place in the panel.
    host_noloci = np.concatenate([p for i, p in enumerate(pieces) if i % 2 == 0])

    provirus_ltrs = [(0, LTR_LEN), (LTR_LEN + INTERNAL_LEN, 2 * LTR_LEN + INTERNAL_LEN)]

    panel = Panel()
    panel.add("EXO_REF", exo_ref, "EXO", provirus_ltrs)
    panel.add("HOST", host, "HOST")
    panel.add("HOST_NOLOCI", host_noloci, "HOST")
    panel.add("HERV_CONSENSUS", herv_consensus, "HERV", provirus_ltrs)

    strain_diff = (exo_strain != exo_ref)
    divs = list(locus_div.values())
    meta = {
        "realized_strain_divergence": float(strain_diff.mean()),
        "realized_locus_divergence": locus_div,
        "locus_divergence_min": float(min(divs)) if divs else float("nan"),
        "locus_divergence_median": float(np.median(divs)) if divs else float("nan"),
        "locus_divergence_max": float(max(divs)) if divs else float("nan"),
        "rate_shape": rate_shape,
        "rate_block_len": rate_block_len,
        "rate_jitter_shape": rate_jitter_shape,
        "sigma_locus": sigma_locus,
        "sigma_strain": sigma_strain,
        "_locus_diff": locus_diff,
        "_strain_diff": strain_diff,
    }
    return panel, host, herv_spans, exo_strain, meta


def _sample_fragments(template, n_frags, rng, read_len, frag_mean, frag_sd):
    out = []
    for _ in range(n_frags):
        flen = int(rng.normal(frag_mean, frag_sd))
        flen = max(read_len + 10, min(flen, len(template)))
        start = 0 if len(template) - flen <= 0 else int(rng.integers(0, len(template) - flen))
        frag = template[start:start + flen]
        out.append((start, frag[:read_len].copy(), frag[-read_len:].copy()))
    return out


def add_errors(read: np.ndarray, err: float, rng: np.random.Generator) -> np.ndarray:
    if err <= 0:
        return read
    mask = rng.random(len(read)) < err
    n = int(mask.sum())
    if n == 0:
        return read
    out = read.copy()
    out[mask] = (out[mask] + rng.integers(1, 4, size=n, dtype=np.uint8)) % 4
    return out


def _window_div(diff: np.ndarray, lo: int, hi: int) -> float:
    lo, hi = max(0, lo), min(len(diff), hi)
    return float(diff[lo:hi].mean()) if hi > lo else float("nan")


def simulate_reads(
    exo_strain: np.ndarray,
    host: np.ndarray,
    herv_spans: dict[str, tuple[int, int]],
    meta: dict,
    *,
    exo_depth: float = 20.0,
    host_depth: float = 10.0,
    read_len: int = 150,
    frag_mean: int = 350,
    frag_sd: int = 50,
    error_rate: float = 0.002,
    seed: int = 43,
) -> list[Truth]:
    """Reads from the exogenous strain (true positives) and from the host
    background including its endogenous loci (the false-positive reservoir)."""
    rng = np.random.default_rng(seed)
    reads: list[Truth] = []
    strain_diff = meta["_strain_diff"]
    locus_diff = meta["_locus_diff"]

    n_exo = int(exo_depth * len(exo_strain) / (2 * read_len))
    for start, r1, r2 in _sample_fragments(exo_strain, n_exo, rng, read_len, frag_mean, frag_sd):
        reads.append(Truth(add_errors(r1, error_rate, rng),
                           add_errors(r2, error_rate, rng),
                           "EXO", "EXO_STRAIN", start,
                           _window_div(strain_diff, start, start + read_len)))

    spans = sorted((s, e, name) for name, (s, e) in herv_spans.items())
    n_host = int(host_depth * len(host) / (2 * read_len))
    for start, r1, r2 in _sample_fragments(host, n_host, rng, read_len, frag_mean, frag_sd):
        end = start + read_len
        source, locus, ldiv = "HOST", "HOST_FILLER", float("nan")
        for s, e, name in spans:
            if start < e and end > s:
                source, locus = "HERV", name
                ldiv = _window_div(locus_diff[name], start - s, end - s)
                break
        reads.append(Truth(add_errors(r1, error_rate, rng),
                           add_errors(r2, error_rate, rng),
                           source, locus, start, ldiv))
    return reads
