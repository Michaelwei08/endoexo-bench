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


def _best_local(read: np.ndarray, ref: np.ndarray, diag: int) -> tuple[int, int, int]:
    """Max-scoring ungapped local segment on one diagonal.

    Returns (score, ref_start, aligned_len). Vectorised max-subarray:
    max_j (prefix[j] - min_{i<=j} prefix[i]).
    """
    lo = max(0, -diag)
    hi = min(len(read), len(ref) - diag)
    if hi - lo < 1:
        return 0, 0, 0
    sub_read = read[lo:hi]
    sub_ref = ref[lo + diag:hi + diag]
    scores = np.where(sub_read == sub_ref, MATCH, MISMATCH).astype(np.int32)
    prefix = np.concatenate([[0], np.cumsum(scores)])
    running_min = np.minimum.accumulate(prefix[:-1])
    gains = prefix[1:] - running_min
    end = int(np.argmax(gains))
    best = int(gains[end])
    if best <= 0:
        return 0, 0, 0
    start = int(np.argmin(prefix[:end + 1]))
    return best, lo + diag + start, end - start + 1


@dataclass
class Hit:
    ref_id: str
    category: str
    score: int
    ref_start: int
    aligned_len: int


def align_read(read: np.ndarray, index: PanelIndex, top: int = 6) -> list[Hit]:
    """Best hit per reference, both strands, sorted by score descending."""
    best_per_ref: dict[str, Hit] = {}
    for strand_seq in (read, revcomp(read)):
        for ref_id, diag, _n in index.candidates(strand_seq, top=top):
            score, ref_start, alen = _best_local(strand_seq, index.panel.refs[ref_id], diag)
            prev = best_per_ref.get(ref_id)
            if prev is None or score > prev.score:
                best_per_ref[ref_id] = Hit(ref_id, index.panel.categories[ref_id],
                                           score, ref_start, alen)
    return sorted(best_per_ref.values(), key=lambda h: -h.score)


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

    mate_hits = align_read(mate, index)
    mate_best = mate_hits[0] if mate_hits else None

    ref_len = len(index.panel.refs[best.ref_id])
    dist_term = min(best.ref_start, ref_len - (best.ref_start + best.aligned_len))

    return ReadFeatures(
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
]
