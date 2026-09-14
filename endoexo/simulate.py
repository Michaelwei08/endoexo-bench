#!/usr/bin/env python3
"""Build a reference panel and a read pool in which exogenous/endogenous homology
is controlled BY CONSTRUCTION, so every read carries its true source.

Geometry follows HTLV-1: a provirus of two identical 755 bp LTRs flanking an
internal region. The endogenous family is the same ancestral provirus carried
forward at a higher substitution divergence and inserted at several loci in a
host background, which is what makes a read from a HERV locus able to align to
the exogenous reference at all.

Nucleotides are encoded A=0 C=1 G=2 T=3 as uint8 throughout.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
LTR_LEN = 755
INTERNAL_LEN = 7000


def mutate(seq: np.ndarray, divergence: float, rng: np.random.Generator) -> np.ndarray:
    """Substitute a fraction `divergence` of sites uniformly to a DIFFERENT base."""
    out = seq.copy()
    n = int(round(divergence * len(seq)))
    if n == 0:
        return out
    sites = rng.choice(len(seq), size=n, replace=False)
    # +1..+3 mod 4 guarantees the new base differs from the old one
    out[sites] = (out[sites] + rng.integers(1, 4, size=n, dtype=np.uint8)) % 4
    return out


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
    """One simulated fragment and the source it actually came from."""

    read: np.ndarray
    mate: np.ndarray
    source: str           # "EXO" or "HERV" or "HOST"
    locus: str            # grouping unit for grouped CV
    frag_start: int


def simulate_world(
    herv_divergence: float,
    *,
    exo_strain_divergence: float = 0.02,
    n_herv_loci: int = 20,
    host_filler_bp: int = 330_000,
    seed: int = 42,
) -> tuple[Panel, np.ndarray, dict[str, tuple[int, int]], np.ndarray]:
    """Return (panel, host_background, herv_locus_spans, exo_strain).

    The exogenous REFERENCE and the exogenous STRAIN that reads come from are
    separated by `exo_strain_divergence`, because in a real cohort the sample
    virus is never the reference. `herv_divergence` separates the endogenous
    family from the same ancestral provirus and is the axis of the benchmark.
    """
    rng = np.random.default_rng(seed)
    ancestral = build_provirus(rng)

    exo_ref = mutate(ancestral, exo_strain_divergence / 2.0, rng)
    exo_strain = mutate(exo_ref, exo_strain_divergence, rng)

    herv_consensus = mutate(ancestral, herv_divergence, rng)

    # Host background: random filler with HERV elements inserted at known spans.
    pieces: list[np.ndarray] = []
    herv_spans: dict[str, tuple[int, int]] = {}
    filler_chunk = host_filler_bp // (n_herv_loci + 1)
    cursor = 0
    for i in range(n_herv_loci):
        filler = rng.integers(0, 4, filler_chunk, dtype=np.uint8)
        pieces.append(filler)
        cursor += filler_chunk
        # each locus drifts a little further from the family consensus
        element = mutate(herv_consensus, herv_divergence * 0.25, rng)
        pieces.append(element)
        herv_spans[f"HERV_L{i:02d}"] = (cursor, cursor + len(element))
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

    return panel, host, herv_spans, exo_strain


def _sample_fragments(
    template: np.ndarray,
    n_frags: int,
    rng: np.random.Generator,
    read_len: int,
    frag_mean: int,
    frag_sd: int,
) -> list[tuple[int, np.ndarray, np.ndarray]]:
    out = []
    for _ in range(n_frags):
        flen = int(rng.normal(frag_mean, frag_sd))
        flen = max(read_len + 10, min(flen, len(template)))
        if len(template) - flen <= 0:
            start = 0
        else:
            start = int(rng.integers(0, len(template) - flen))
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


def simulate_reads(
    exo_strain: np.ndarray,
    host: np.ndarray,
    herv_spans: dict[str, tuple[int, int]],
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
    background including its HERV loci (the false-positive reservoir)."""
    rng = np.random.default_rng(seed)
    reads: list[Truth] = []

    n_exo = int(exo_depth * len(exo_strain) / (2 * read_len))
    for start, r1, r2 in _sample_fragments(exo_strain, n_exo, rng, read_len, frag_mean, frag_sd):
        reads.append(Truth(add_errors(r1, error_rate, rng),
                           add_errors(r2, error_rate, rng),
                           "EXO", "EXO_STRAIN", start))

    # Interval lookup so a host fragment is attributed to a HERV locus if it
    # overlaps one; otherwise it is plain host filler.
    spans = sorted((s, e, name) for name, (s, e) in herv_spans.items())
    n_host = int(host_depth * len(host) / (2 * read_len))
    for start, r1, r2 in _sample_fragments(host, n_host, rng, read_len, frag_mean, frag_sd):
        end = start + read_len
        source, locus = "HOST", "HOST_FILLER"
        for s, e, name in spans:
            if start < e and end > s:
                source, locus = "HERV", name
                break
        reads.append(Truth(add_errors(r1, error_rate, rng),
                           add_errors(r2, error_rate, rng),
                           source, locus, start))
    return reads
