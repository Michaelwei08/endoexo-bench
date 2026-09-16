#!/usr/bin/env python3
"""Competitive ungapped local alignment against a reference panel, in pure
Python + numpy, and the per-read features the discrimination task is built on.

WHY NOT BWA
  This is a design-validation harness: no aligner binary is available on the
  development machine, and the question under study is which alignment-level
  features separate a true exogenous read from an endogenous cross-mapping, not
  how fast the alignment runs. Scoring follows BWA-MEM conventions (match +1,
  mismatch -4, local alignment so the non-matching tail is soft-clipped) and the
  default seed length matches BWA-MEM's -k 19, so the seed-sensitivity limit is
  the same one BWA has.

  LIMITATION, stated: alignment is UNGAPPED. Indels between the endogenous copy
  and the exogenous reference are not modelled, so AS is a lower bound for any
  read whose true alignment needs a gap. A gapped rewrite, or a real BWA run on
  the same simulated FASTQs, is required before any published number.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

MATCH = 1
MISMATCH = -4


def revcomp(seq: np.ndarray) -> np.ndarray:
    """A=0 C=1 G=2 T=3, so the complement is simply 3 - base."""
    return (3 - seq)[::-1]


def _pack_kmers(seq: np.ndarray, k: int) -> np.ndarray:
    """Rolling 2-bit packing of every k-mer start in `seq`."""
    if len(seq) < k:
        return np.empty(0, dtype=np.int64)
    s = seq.astype(np.int64)
    windows = np.lib.stride_tricks.sliding_window_view(s, k)
    weights = (4 ** np.arange(k - 1, -1, -1, dtype=np.int64))
    return windows @ weights


class PanelIndex:
    def __init__(self, panel, k: int = 19, max_postings: int = 64):
        self.panel = panel
        self.k = k
        self.max_postings = max_postings
        self.index: dict[int, list[tuple[str, int]]] = {}
        for ref_id, seq in panel.refs.items():
            codes = _pack_kmers(seq, k)
            for pos, code in enumerate(codes.tolist()):
                bucket = self.index.setdefault(code, [])
                if len(bucket) < max_postings:
                    bucket.append((ref_id, pos))

    def candidates(self, read: np.ndarray, top: int = 6) -> list[tuple[str, int, int]]:
        """Seed hits collapsed onto (ref_id, diagonal), best-supported first."""
        codes = _pack_kmers(read, self.k)
        support: dict[tuple[str, int], int] = {}
        for read_pos, code in enumerate(codes.tolist()):
            for ref_id, ref_pos in self.index.get(code, ()):
                key = (ref_id, ref_pos - read_pos)
                support[key] = support.get(key, 0) + 1
        ranked = sorted(support.items(), key=lambda kv: -kv[1])[:top]
        return [(ref_id, diag, n) for (ref_id, diag), n in ranked]


def _best_local(read: np.ndarray, ref: np.ndarray,
                diag: int) -> tuple[int, int, int, int]:
    """Max-scoring ungapped local segment on one diagonal.

    Returns (score, ref_start, aligned_len, read_start). Vectorised max-subarray:
    max_j (prefix[j] - min_{i<=j} prefix[i]).

    `read_start` exists so a caller can recover the aligned BASES, not just the
    score -- which is what a pileup needs, and a pileup is what the intra-host
    diversity statistic is computed from.
    """
    lo = max(0, -diag)
    hi = min(len(read), len(ref) - diag)
    if hi - lo < 1:
        return 0, 0, 0, 0
    sub_read = read[lo:hi]
    sub_ref = ref[lo + diag:hi + diag]
    scores = np.where(sub_read == sub_ref, MATCH, MISMATCH).astype(np.int32)
    prefix = np.concatenate([[0], np.cumsum(scores)])
    running_min = np.minimum.accumulate(prefix[:-1])
    gains = prefix[1:] - running_min
    end = int(np.argmax(gains))
    best = int(gains[end])
    if best <= 0:
        return 0, 0, 0, 0
    start = int(np.argmin(prefix[:end + 1]))
    return best, lo + diag + start, end - start + 1, lo + start


@dataclass
class Hit:
    ref_id: str
    category: str
    score: int
    ref_start: int
    aligned_len: int
    # Enough to recover the aligned BASES: which strand won, and where on that
    # strand the aligned segment starts. A pileup needs bases, not scores.
    read_start: int = 0
    reverse: bool = False


def align_read(read: np.ndarray, index: PanelIndex, top: int = 6) -> list[Hit]:
    """Best hit per reference, both strands, sorted by score descending."""
    best_per_ref: dict[str, Hit] = {}
    for reverse, strand_seq in ((False, read), (True, revcomp(read))):
        for ref_id, diag, _n in index.candidates(strand_seq, top=top):
            score, ref_start, alen, read_start = _best_local(
                strand_seq, index.panel.refs[ref_id], diag)
            prev = best_per_ref.get(ref_id)
            if prev is None or score > prev.score:
                best_per_ref[ref_id] = Hit(ref_id, index.panel.categories[ref_id],
                                           score, ref_start, alen,
                                           read_start, reverse)
    return sorted(best_per_ref.values(), key=lambda h: -h.score)


def aligned_bases(read: np.ndarray, hit: Hit) -> tuple[int, np.ndarray]:
    """(ref_start, bases) for the segment of `read` that `hit` aligned.

    Ungapped, so the mapping from reference position to base is one to one and
    needs no CIGAR walk. The strand is taken from the hit, which is why Hit
    carries it.
    """
    seq = revcomp(read) if hit.reverse else read
    return hit.ref_start, seq[hit.read_start:hit.read_start + hit.aligned_len]


def _in_ltr(index: PanelIndex, ref_id: str, start: int, alen: int) -> int:
    for s, e in index.panel.ltr_spans.get(ref_id, ()):
        if start < e and start + alen > s:
            return 1
    return 0


@dataclass
class ReadFeatures:
    best_ref: str
    best_cat: str
    AS: int
    XS_any: int
    XS_other_cat: int
    as_minus_xs: int
    as_minus_xs_other_cat: int
    aligned_len: int
    softclip_len: int
    dist_to_terminus: int
    in_ltr: int
    n_alt_within5: int
    mate_same_ref: int
    mate_AS: int
    # Tie structure. The two ways a read can be mis-assigned are completely
    # different problems and these separate them: an EXACT TIE means the panel
    # contains something that explains the read equally well, which is free to
    # detect and removable at a cost in sensitivity; an UNOPPOSED WIN means the
    # read beat everything in the panel and is indistinguishable from a true
    # positive by any competition-based test.
    n_tied_top: int = 1
    n_tied_top_categories: int = 1
    # Pair-level discrimination. A read inside a conserved window carries no
    # information, but its mate sits a fragment away and may land outside that
    # window, where the competition is decidable.
    mate_best_cat_is_exo: int = 0
    mate_as_minus_xs: int = 0
    mate_n_tied_top_categories: int = 1
    pair_max_gap: int = 0
    pair_min_gap: int = 0
    # Position on the best reference. NOT a model feature -- it was the leakage
    # source in F004 -- but it is what a SAMPLE-level coverage profile is built
    # from, which is the one route past the read-pair ceiling of F026.
    best_ref_start: int = 0
    # Enough to recover the aligned BASES without re-aligning, which is what an
    # intra-host diversity statistic needs. Also not model features.
    best_read_start: int = 0
    best_reverse: bool = False


def _tie_structure(hits: list[Hit]) -> tuple[int, int]:
    """(number of references tied at the top score, number of distinct categories
    among them). A read whose top score is held by more than one CATEGORY cannot
    be assigned by competition at all.

    NOTE ON TIE-BREAKING. `align_read` sorts by score only, so a tie resolves to
    whichever reference was added to the panel first -- in practice the exogenous
    one. That is the OPTIMISTIC policy and it is kept as the raw behaviour on
    purpose, with the tie exposed here so the caller decides. A real aligner
    breaks the tie arbitrarily and reports mapping quality 0; treating the tie as
    an assignment is what silently inflates an exogenous call count.
    """
    if not hits:
        return 0, 0
    top = hits[0].score
    tied = [h for h in hits if h.score == top]
    return len(tied), len({h.category for h in tied})


def features_for_pair(read: np.ndarray, mate: np.ndarray, index: PanelIndex,
                      read_len: int) -> ReadFeatures | None:
    hits = align_read(read, index)
    if not hits or hits[0].score <= 0:
        return None
    best = hits[0]
    others = [h for h in hits[1:] if h.score > 0]
    xs_any = others[0].score if others else 0
    diff_cat = [h.score for h in others if h.category != best.category]
    xs_other = max(diff_cat) if diff_cat else 0
    n_tied, n_tied_cats = _tie_structure(hits)

    mate_hits = align_read(mate, index)
    mate_best = mate_hits[0] if mate_hits else None
    mate_others = [h for h in mate_hits[1:] if h.score > 0] if mate_hits else []
    mate_gap = (mate_best.score - (mate_others[0].score if mate_others else 0)) \
        if mate_best else 0
    _, mate_tied_cats = _tie_structure(mate_hits)
    read_gap = best.score - xs_any

    ref_len = len(index.panel.refs[best.ref_id])
    dist_term = min(best.ref_start, ref_len - (best.ref_start + best.aligned_len))

    return ReadFeatures(
        best_ref_start=best.ref_start,
        best_read_start=best.read_start,
        best_reverse=best.reverse,
        n_tied_top=n_tied,
        n_tied_top_categories=n_tied_cats,
        mate_best_cat_is_exo=int(mate_best is not None and mate_best.category == "EXO"),
        mate_as_minus_xs=mate_gap,
        mate_n_tied_top_categories=mate_tied_cats,
        pair_max_gap=max(read_gap, mate_gap),
        pair_min_gap=min(read_gap, mate_gap),
        best_ref=best.ref_id,
        best_cat=best.category,
        AS=best.score,
        XS_any=xs_any,
        XS_other_cat=xs_other,
        as_minus_xs=best.score - xs_any,
        as_minus_xs_other_cat=best.score - xs_other,
        aligned_len=best.aligned_len,
        softclip_len=read_len - best.aligned_len,
        dist_to_terminus=dist_term,
        in_ltr=_in_ltr(index, best.ref_id, best.ref_start, best.aligned_len),
        n_alt_within5=sum(1 for h in others if h.score >= best.score - 5),
        mate_same_ref=int(mate_best is not None and mate_best.ref_id == best.ref_id),
        mate_AS=mate_best.score if mate_best else 0,
    )


FEATURE_COLUMNS = [
    "AS", "XS_any", "XS_other_cat", "as_minus_xs", "as_minus_xs_other_cat",
    "aligned_len", "softclip_len", "dist_to_terminus", "in_ltr",
    "n_alt_within5", "mate_same_ref", "mate_AS",
    "n_tied_top", "n_tied_top_categories",
    "mate_best_cat_is_exo", "mate_as_minus_xs", "mate_n_tied_top_categories",
    "pair_max_gap", "pair_min_gap",
]

# Features whose value comes from the MATE rather than the read itself. F016
# showed that a read inside a conserved window is information-free, so these are
# the only read-pair-level route past that ceiling and are ablated separately.
MATE_FEATURES = ("mate_same_ref", "mate_AS", "mate_best_cat_is_exo",
                 "mate_as_minus_xs", "mate_n_tied_top_categories",
                 "pair_max_gap", "pair_min_gap")
