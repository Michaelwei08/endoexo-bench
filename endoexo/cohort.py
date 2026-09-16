#!/usr/bin/env python3
"""A cohort of sample genomes sharing one host backbone, so that a host-side
junction coordinate MEANS THE SAME THING in every sample.

WHY THIS MODULE EXISTS
  Every discriminator measured so far is computed within one sample, and F050
  closed the last of them: when the read's true source is absent from the panel,
  nothing per-sample beats counting reads. What remains is a signal that does not
  live inside any sample at all.

  An endogenous insertion sits at the SAME host coordinate in every individual
  who carries it, because they inherited it. An exogenous integration sits at a
  coordinate PRIVATE to the individual, because it happened in their lifetime.
  So the discriminator is the cross-individual RECURRENCE of the host-side
  junction position, and it needs no sequence discrimination whatsoever. It
  sidesteps the tie / unopposed-win dichotomy entirely rather than solving it.

WHAT HAD TO CHANGE FROM simulate.py
  1 THE BACKBONE IS COHORT-LEVEL. The host filler is drawn once for the cohort,
    not per sample, because it is the reference genome. Without that a host
    coordinate is not comparable across samples and recurrence is meaningless.
    This was the one modelling error that would have silently invalidated the
    whole experiment.
  2 LOCI ARE PRESENCE-POLYMORPHIC AT FIXED COORDINATES. Each candidate
    endogenous locus has one cohort-level sequence, one cohort-level coordinate,
    and a population frequency. A sample either carries it or does not. In
    simulate.py every sample carried every locus, which would have made every
    endogenous junction recur in 100 percent of samples -- true but trivially so,
    and it would not have tested the discriminator.
  3 THE EXOGENOUS PROVIRUS IS INTEGRATED, at a private coordinate, in a clonal
    FRACTION of cells. Modelled by sampling that fraction of reads from a
    template that carries the provirus and the rest from one that does not, which
    is what a clonal integration actually looks like at the read level.

  The alignment panel's host reference is the BARE BACKBONE, with no insertions.
  That is what a reference assembly is with respect to an insertionally
  polymorphic locus, and it means every junction coordinate a read reports is
  already in the shared coordinate system.
"""
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, field

import numpy as np

from .simulate import (SATURATION, Panel, _make_locus, build_provirus,
                       gamma_site_rates, lognormal_around, mutate)


@dataclass
class Segment:
    """One stretch of a sample genome and where it came from."""

    start: int              # offset in the sample genome
    length: int
    source: str             # "BACKBONE", "HERV_Lxx" or "EXO"
    source_start: int       # offset within that source


@dataclass
class CohortSpec:
    """Everything shared across the cohort. Drawn once."""

    exo_ref: np.ndarray
    site_rates: np.ndarray
    backbone: np.ndarray
    locus_seq: dict[str, np.ndarray] = field(default_factory=dict)
    locus_coord: dict[str, int] = field(default_factory=dict)
    locus_freq: dict[str, float] = field(default_factory=dict)
    realized_locus_divergence: dict[str, float] = field(default_factory=dict)

    def panel(self) -> Panel:
        """The alignment panel: the exogenous reference and the BARE backbone.

        No insertions, which is the point -- an insertionally polymorphic locus is
        by definition absent from the assembly.
        """
        p = Panel()
        ltr = 755
        internal = len(self.exo_ref) - 2 * ltr
        p.add("EXO_REF", self.exo_ref, "EXO",
              [(0, ltr), (ltr + internal, len(self.exo_ref))])
        p.add("HOST", self.backbone, "HOST")
        return p


def build_cohort(
    herv_divergence: float,
    *,
    n_loci: int = 10,
    backbone_bp: int = 150_000,
    exo_strain_divergence: float = 0.02,
    rate_shape: float = 0.5,
    rate_block_len: int = 400,
    rate_jitter_shape: float = 4.0,
    sigma_locus: float = 0.7,
    freq_lo: float = 0.15,
    freq_hi: float = 0.75,
    seed: int = 42,
) -> CohortSpec:
    """Draw the cohort-level world once: reference, backbone, loci, frequencies."""
    rng = np.random.default_rng(seed)
    exo_ref = build_provirus(rng)
    site_rates = gamma_site_rates(len(exo_ref), rate_shape, rng,
                                  block_len=rate_block_len,
                                  jitter_shape=rate_jitter_shape)
    herv_consensus = mutate(exo_ref, herv_divergence, rng, site_rates)
    family_mask = np.flatnonzero(herv_consensus != exo_ref)

    backbone = rng.integers(0, 4, backbone_bp, dtype=np.uint8)

    spec = CohortSpec(exo_ref=exo_ref, site_rates=site_rates, backbone=backbone)
    # Coordinates spaced far enough apart that no two insertions share a bin.
    margin = backbone_bp // (n_loci + 2)
    for i in range(n_loci):
        name = f"HERV_L{i:02d}"
        target = float(lognormal_around(herv_divergence, sigma_locus, rng))
        element = _make_locus(exo_ref, herv_consensus, family_mask, site_rates,
                              min(target, SATURATION), rng)
        spec.locus_seq[name] = element
        spec.locus_coord[name] = margin * (i + 1)
        spec.locus_freq[name] = float(rng.uniform(freq_lo, freq_hi))
        spec.realized_locus_divergence[name] = float((element != exo_ref).mean())
    return spec


def _assemble(spec: CohortSpec, insertions: list[tuple[int, str, np.ndarray]]
              ) -> tuple[np.ndarray, list[Segment]]:
    """Splice insertions into the backbone, returning the genome and a map back."""
    insertions = sorted(insertions, key=lambda t: t[0])
    pieces: list[np.ndarray] = []
    segments: list[Segment] = []
    cursor_out = 0
    cursor_bb = 0
    for coord, name, seq in insertions:
        if coord > cursor_bb:
            n = coord - cursor_bb
            pieces.append(spec.backbone[cursor_bb:coord])
            segments.append(Segment(cursor_out, n, "BACKBONE", cursor_bb))
            cursor_out += n
            cursor_bb = coord
        pieces.append(seq)
        segments.append(Segment(cursor_out, len(seq), name, 0))
        cursor_out += len(seq)
    if cursor_bb < len(spec.backbone):
        n = len(spec.backbone) - cursor_bb
        pieces.append(spec.backbone[cursor_bb:])
        segments.append(Segment(cursor_out, n, "BACKBONE", cursor_bb))
    return np.concatenate(pieces), segments


def build_sample(spec: CohortSpec, *, infected: bool, sample_seed: int,
                 strain_divergence: float = 0.02, sigma_strain: float = 0.4):
    """One individual's two templates and the truth about their insertions.

    Returns (with_provirus, without_provirus, segments_with, segments_without,
    carried, exo_coord, realized_strain_divergence). The two templates exist so a
    CLONAL fraction can be modelled by sampling reads from each in proportion.
    """
    rng = np.random.default_rng(sample_seed)
    carried = [name for name, f in spec.locus_freq.items() if rng.random() < f]

    base = [(spec.locus_coord[n], n, spec.locus_seq[n]) for n in carried]
    without, seg_without = _assemble(spec, base)

    exo_coord = None
    d_strain = float("nan")
    with_prov, seg_with = without, seg_without
    if infected:
        d_strain = float(lognormal_around(strain_divergence, sigma_strain, rng))
        strain = mutate(spec.exo_ref, d_strain, rng, spec.site_rates)
        # a private coordinate: uniform on the backbone, away from every locus
        taken = set()
        for n in spec.locus_coord.values():
            taken.update(range(max(0, n - 3000), n + 3000))
        while True:
            c = int(rng.integers(2000, len(spec.backbone) - 2000))
            if c not in taken:
                exo_coord = c
                break
        with_prov, seg_with = _assemble(spec, base + [(exo_coord, "EXO", strain)])

    return (with_prov, without, seg_with, seg_without, carried, exo_coord,
            d_strain)


def locate(segments: list[Segment], pos: int) -> tuple[str, int]:
    """(source, backbone coordinate) for a sample-genome position.

    For a BACKBONE position the second value is its reference coordinate, which
    is the shared coordinate system recurrence is counted in. For a position
    inside an insertion it is the insertion's own offset and is not comparable
    across samples -- which is exactly why the HOST side of a junction is the
    informative one.
    """
    starts = [s.start for s in segments]
    i = bisect_right(starts, pos) - 1
    if i < 0:
        i = 0
    s = segments[i]
    off = pos - s.start
    return s.source, s.source_start + off
