#!/usr/bin/env python3
"""Tests for the invariants the findings actually rest on.

Run from the repository root:

    python -m unittest discover -s tests -t .

Two of these exist because the corresponding defect happened during development
and would have been invisible otherwise:

  test_point_model_flattens_rates -- the README claimed a flag combination
    recovered the original point-divergence model. It did not, because a
    hardcoded per-site jitter survived it. That claim was wrong in a shipped
    document for about twenty minutes.
  test_cohort_seed_fixes_the_virus -- the cohort design requires one viral
    reference across samples and different host loci per sample. The first
    version derived both from a single seed, which would have given every
    sample its own virus and quietly invalidated the sample-level task.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from endoexo import simulate
from endoexo.align import (MATCH, MISMATCH, Hit, PanelIndex, _best_local,
                           _tie_structure, align_read, aligned_bases, revcomp)
from endoexo.cohort import build_cohort, build_sample, locate
from endoexo.evaluate import fpr_at_sensitivity
from compare_with_bwa import parse_sam, softclip_len
from run_sample_level import coverage_profile


class TestSequenceOps(unittest.TestCase):
    def test_revcomp_is_an_involution(self):
        rng = np.random.default_rng(0)
        seq = rng.integers(0, 4, 500, dtype=np.uint8)
        np.testing.assert_array_equal(revcomp(revcomp(seq)), seq)

    def test_revcomp_complements_correctly(self):
        # A=0 C=1 G=2 T=3, so A<->T and C<->G
        seq = np.array([0, 1, 2, 3], dtype=np.uint8)
        np.testing.assert_array_equal(revcomp(seq), np.array([0, 1, 2, 3], np.uint8))
        seq = np.array([0, 0, 1], dtype=np.uint8)            # A A C
        np.testing.assert_array_equal(revcomp(seq), np.array([2, 3, 3], np.uint8))

    def test_mutate_hits_the_target_divergence(self):
        rng = np.random.default_rng(1)
        seq = rng.integers(0, 4, 20_000, dtype=np.uint8)
        for target in (0.01, 0.05, 0.20):
            out = simulate.mutate(seq, target, rng)
            realized = float((out != seq).mean())
            self.assertAlmostEqual(realized, target, delta=0.01 + 0.1 * target)

    def test_mutate_never_keeps_the_same_base(self):
        rng = np.random.default_rng(2)
        seq = np.zeros(5000, dtype=np.uint8)
        out = simulate.mutate(seq, 1.0, rng)         # every site targeted
        changed = out != seq
        self.assertGreater(changed.mean(), 0.7)      # 1.0 clips to SATURATION
        self.assertTrue((out[changed] != 0).all())

    def test_mutate_zero_divergence_is_a_copy(self):
        rng = np.random.default_rng(3)
        seq = rng.integers(0, 4, 100, dtype=np.uint8)
        np.testing.assert_array_equal(simulate.mutate(seq, 0.0, rng), seq)


class TestSiteRates(unittest.TestCase):
    def test_mean_rate_is_one(self):
        rng = np.random.default_rng(4)
        rates = simulate.gamma_site_rates(40_000, 0.5, rng)
        self.assertAlmostEqual(float(rates.mean()), 1.0, places=9)  # normalised exactly

    def test_rates_are_autocorrelated(self):
        """The whole point of F012: adjacent sites must share a regional rate."""
        rng = np.random.default_rng(5)
        rates = simulate.gamma_site_rates(40_000, 0.5, rng, block_len=400,
                                          jitter_shape=0.0)
        r = float(np.corrcoef(rates[:-1], rates[1:])[0, 1])
        self.assertGreater(r, 0.9, "block rates should be nearly constant locally")

    def test_point_model_flattens_rates(self):
        """Guards the documented ablation. See the module docstring."""
        rng = np.random.default_rng(6)
        rates = simulate.gamma_site_rates(20_000, 1e6, rng, jitter_shape=0.0)
        self.assertLess(float(rates.std()), 0.01)

    def test_jitter_is_what_the_flag_says(self):
        rng = np.random.default_rng(7)
        flat = simulate.gamma_site_rates(20_000, 1e6, rng, jitter_shape=0.0)
        jittered = simulate.gamma_site_rates(20_000, 1e6, rng, jitter_shape=4.0)
        self.assertGreater(float(jittered.std()), 10 * float(flat.std()) + 0.1)


class TestLocusConstruction(unittest.TestCase):
    def _setup(self, family_div):
        rng = np.random.default_rng(8)
        ref = rng.integers(0, 4, 8_000, dtype=np.uint8)
        rates = np.ones(len(ref))
        consensus = simulate.mutate(ref, family_div, rng, rates)
        mask = np.flatnonzero(consensus != ref)
        return ref, consensus, mask, rates, rng

    def test_locus_below_family_divergence_reverts(self):
        """A locus younger than the family must be able to exist at all -- the
        bug F011 recorded was that the family divergence acted as a hard floor."""
        ref, consensus, mask, rates, rng = self._setup(0.20)
        locus = simulate._make_locus(ref, consensus, mask, rates, 0.05, rng)
        realized = float((locus != ref).mean())
        self.assertAlmostEqual(realized, 0.05, delta=0.01)

    def test_locus_above_family_divergence_adds(self):
        ref, consensus, mask, rates, rng = self._setup(0.10)
        locus = simulate._make_locus(ref, consensus, mask, rates, 0.25, rng)
        realized = float((locus != ref).mean())
        self.assertAlmostEqual(realized, 0.25, delta=0.03)

    def test_reverted_locus_keeps_family_sites_only(self):
        """Reversion may only undo family substitutions, never invent new ones."""
        ref, consensus, mask, rates, rng = self._setup(0.20)
        locus = simulate._make_locus(ref, consensus, mask, rates, 0.05, rng)
        differing = set(np.flatnonzero(locus != ref).tolist())
        self.assertTrue(differing.issubset(set(mask.tolist())))


class TestCohortDecomposition(unittest.TestCase):
    def test_cohort_seed_fixes_the_virus(self):
        """Same cohort seed, different sample seed: one reference, different
        host loci and different strains. See the module docstring."""
        a = simulate.simulate_world(0.10, n_herv_loci=4, host_filler_bp=20_000,
                                    seed=42, sample_seed=1)
        b = simulate.simulate_world(0.10, n_herv_loci=4, host_filler_bp=20_000,
                                    seed=42, sample_seed=2)
        np.testing.assert_array_equal(a[0].refs["EXO_REF"], b[0].refs["EXO_REF"])
        np.testing.assert_array_equal(a[0].refs["HERV_CONSENSUS"],
                                      b[0].refs["HERV_CONSENSUS"])
        self.assertFalse(np.array_equal(a[1], b[1]), "host must differ per sample")
        self.assertFalse(np.array_equal(a[3], b[3]), "strain must differ per sample")

    def test_default_sample_seed_reproduces_single_sample(self):
        a = simulate.simulate_world(0.10, n_herv_loci=4, host_filler_bp=20_000,
                                    seed=42)
        b = simulate.simulate_world(0.10, n_herv_loci=4, host_filler_bp=20_000,
                                    seed=42, sample_seed=42)
        np.testing.assert_array_equal(a[1], b[1])
        np.testing.assert_array_equal(a[3], b[3])

    def test_host_without_loci_is_shorter_by_the_elements(self):
        panel, host, spans, _, _ = simulate.simulate_world(
            0.10, n_herv_loci=4, host_filler_bp=20_000, seed=42)
        removed = sum(e - s for s, e in spans.values())
        self.assertEqual(len(panel.refs["HOST_NOLOCI"]), len(host) - removed)


class TestAligner(unittest.TestCase):
    def _brute_force(self, read, ref, diag):
        lo, hi = max(0, -diag), min(len(read), len(ref) - diag)
        best = 0
        for i in range(lo, hi + 1):
            run = 0
            for j in range(i, hi):
                run += MATCH if read[j] == ref[j + diag] else MISMATCH
                best = max(best, run)
        return best

    def test_best_local_matches_brute_force(self):
        rng = np.random.default_rng(9)
        ref = rng.integers(0, 4, 300, dtype=np.uint8)
        for trial in range(25):
            read = simulate.mutate(ref[50:200].copy(), 0.1,
                                   np.random.default_rng(trial))
            score, _, alen, _ = _best_local(read, ref, 50)
            self.assertEqual(score, self._brute_force(read, ref, 50))
            self.assertLessEqual(alen, len(read))

    def test_perfect_match_scores_read_length(self):
        rng = np.random.default_rng(10)
        ref = rng.integers(0, 4, 200, dtype=np.uint8)
        score, start, alen, read_start = _best_local(ref[20:120].copy(), ref, 20)
        self.assertEqual(score, 100)
        self.assertEqual(start, 20)
        self.assertEqual(alen, 100)
        self.assertEqual(read_start, 0)

    def test_read_start_indexes_the_aligned_segment(self):
        """A read whose first 30 bases are junk must report read_start 30."""
        rng = np.random.default_rng(11)
        ref = rng.integers(0, 4, 400, dtype=np.uint8)
        read = np.concatenate([rng.integers(0, 4, 30, dtype=np.uint8),
                               ref[100:220].copy()])
        # diagonal such that read position 30 maps to reference position 100
        score, ref_start, alen, read_start = _best_local(read, ref, 70)
        self.assertEqual(alen, 120)
        self.assertEqual(ref_start, 100)
        self.assertEqual(read_start, 30)

    def test_tie_structure_counts_categories(self):
        hits = [Hit("EXO_REF", "EXO", 100, 0, 100),
                Hit("HOST", "HOST", 100, 5, 100),
                Hit("HERV_CONSENSUS", "HERV", 90, 0, 95)]
        self.assertEqual(_tie_structure(hits), (2, 2))

    def test_tie_structure_single_winner(self):
        hits = [Hit("EXO_REF", "EXO", 100, 0, 100),
                Hit("HOST", "HOST", 94, 5, 100)]
        self.assertEqual(_tie_structure(hits), (1, 1))

    def test_tie_structure_same_category_is_not_ambiguous(self):
        """Two references of the SAME category tied at the top is not an
        assignment problem -- the category is still decided."""
        hits = [Hit("HOST", "HOST", 100, 0, 100),
                Hit("HOST_NOLOCI", "HOST", 100, 0, 100)]
        self.assertEqual(_tie_structure(hits), (2, 1))


class TestAlignedBases(unittest.TestCase):
    """aligned_bases is correctness-critical: the diversity statistic is computed
    from the pileup it produces, so an off-by-one here silently corrupts every
    allele frequency."""

    def _index(self, ref):
        panel = simulate.Panel()
        panel.add("EXO_REF", ref, "EXO")
        return PanelIndex(panel, k=19)

    def test_forward_read_recovers_its_own_bases(self):
        rng = np.random.default_rng(20)
        ref = rng.integers(0, 4, 3000, dtype=np.uint8)
        read = ref[500:650].copy()
        hit = align_read(read, self._index(ref))[0]
        start, bases = aligned_bases(read, hit)
        self.assertEqual(start, 500)
        np.testing.assert_array_equal(bases, ref[500:650])

    def test_reverse_read_recovers_reference_orientation(self):
        rng = np.random.default_rng(21)
        ref = rng.integers(0, 4, 3000, dtype=np.uint8)
        read = revcomp(ref[800:950].copy())
        hit = align_read(read, self._index(ref))[0]
        self.assertTrue(hit.reverse)
        start, bases = aligned_bases(read, hit)
        self.assertEqual(start, 800)
        np.testing.assert_array_equal(bases, ref[800:950])

    def test_mismatches_are_preserved_not_corrected(self):
        """The pileup must see the READ's base, including where it disagrees."""
        rng = np.random.default_rng(22)
        ref = rng.integers(0, 4, 3000, dtype=np.uint8)
        read = ref[1000:1150].copy()
        read[40] = (read[40] + 1) % 4
        read[90] = (read[90] + 2) % 4
        hit = align_read(read, self._index(ref))[0]
        start, bases = aligned_bases(read, hit)
        self.assertEqual(start, 1000)
        self.assertEqual(int(bases[40]), int(read[40]))
        self.assertNotEqual(int(bases[40]), int(ref[1040]))
        self.assertEqual(int(bases[90]), int(read[90]))

    def test_bases_length_matches_aligned_len(self):
        rng = np.random.default_rng(23)
        ref = rng.integers(0, 4, 3000, dtype=np.uint8)
        read = np.concatenate([rng.integers(0, 4, 25, dtype=np.uint8),
                               ref[1500:1625].copy()])
        hit = align_read(read, self._index(ref))[0]
        _start, bases = aligned_bases(read, hit)
        self.assertEqual(len(bases), hit.aligned_len)


class TestCohortGenomes(unittest.TestCase):
    """The junction-recurrence experiment rests on one modelling invariant: the
    host backbone must be COHORT-level, so a host coordinate means the same thing
    in every sample. Getting that wrong would not crash anything -- it would
    silently make recurrence meaningless, which is why it is pinned here."""

    def _spec(self, **kw):
        return build_cohort(0.10, n_loci=4, backbone_bp=40_000, seed=42, **kw)

    def test_backbone_is_shared_across_the_cohort(self):
        a, b = self._spec(), self._spec()
        np.testing.assert_array_equal(a.backbone, b.backbone)
        for name in a.locus_seq:
            np.testing.assert_array_equal(a.locus_seq[name], b.locus_seq[name])
            self.assertEqual(a.locus_coord[name], b.locus_coord[name])

    def test_panel_host_reference_carries_no_insertions(self):
        """An insertionally polymorphic locus is absent from the assembly. If the
        panel's host reference contained the elements, the whole poly case would
        silently become the reference-insertion case."""
        spec = self._spec()
        panel = spec.panel()
        self.assertEqual(len(panel.refs["HOST"]), len(spec.backbone))
        np.testing.assert_array_equal(panel.refs["HOST"], spec.backbone)

    def test_presence_is_polymorphic_across_samples(self):
        spec = self._spec()
        carried = [set(build_sample(spec, infected=False, sample_seed=s)[4])
                   for s in range(1, 40)]
        self.assertGreater(len({frozenset(c) for c in carried}), 1,
                           "samples must differ in which loci they carry")

    def test_infected_sample_gets_a_private_coordinate(self):
        spec = self._spec()
        coords = set()
        for s in range(1, 25):
            out = build_sample(spec, infected=True, sample_seed=s)
            coords.add(out[5])
            self.assertIsNotNone(out[5])
        self.assertGreater(len(coords), 20, "integration coordinates must vary")

    def test_uninfected_sample_has_no_provirus(self):
        spec = self._spec()
        out = build_sample(spec, infected=False, sample_seed=3)
        self.assertIsNone(out[5])
        np.testing.assert_array_equal(out[0], out[1])   # the two templates agree

    def test_locate_maps_backbone_positions_to_reference_coordinates(self):
        spec = self._spec()
        with_prov, without, seg_with, _sw, carried, exo_coord, _d = build_sample(
            spec, infected=True, sample_seed=5)
        # the base immediately before the first insertion must report its own
        # backbone coordinate
        first = min([spec.locus_coord[n] for n in carried] + [exo_coord])
        src, coord = locate(seg_with, first - 1)
        self.assertEqual(src, "BACKBONE")
        self.assertEqual(coord, first - 1)

    def test_locate_identifies_the_provirus_interior(self):
        spec = self._spec()
        _wp, _wo, seg_with, _sw, carried, exo_coord, _d = build_sample(
            spec, infected=True, sample_seed=6)
        offset = sum(len(spec.locus_seq[n]) for n in carried
                     if spec.locus_coord[n] < exo_coord)
        src, _coord = locate(seg_with, exo_coord + offset + 100)
        self.assertEqual(src, "EXO")

    def test_genome_length_is_backbone_plus_insertions(self):
        spec = self._spec()
        wp, wo, _sw, _so, carried, _c, _d = build_sample(
            spec, infected=True, sample_seed=7)
        inserted = sum(len(spec.locus_seq[n]) for n in carried)
        self.assertEqual(len(wo), len(spec.backbone) + inserted)
        self.assertEqual(len(wp), len(wo) + len(spec.exo_ref))


class TestMetrics(unittest.TestCase):
    def test_fpr_at_sensitivity_perfect_separation(self):
        y = np.array([0, 0, 0, 1, 1, 1])
        self.assertEqual(fpr_at_sensitivity(y, np.array([0., 0., 0., 1., 1., 1.]), 0.95), 0.0)

    def test_fpr_at_sensitivity_inverted_scores(self):
        y = np.array([0, 0, 1, 1])
        self.assertEqual(fpr_at_sensitivity(y, np.array([1., 1., 0., 0.]), 0.95), 1.0)

    def test_fpr_at_sensitivity_needs_both_classes(self):
        y = np.ones(4, dtype=int)
        self.assertTrue(np.isnan(fpr_at_sensitivity(y, np.arange(4.0), 0.95)))


class TestCoverageProfile(unittest.TestCase):
    def test_full_uniform_coverage(self):
        starts = np.arange(0, 1000, 10)
        lengths = np.full(len(starts), 10)
        p = coverage_profile(starts, lengths, 1000, bin_size=250)
        self.assertAlmostEqual(p["breadth_1x"], 1.0)
        self.assertAlmostEqual(p["max_bin_frac"], 0.25, delta=0.01)
        self.assertLess(p["depth_cv"], 0.01)

    def test_single_pileup_is_detected(self):
        starts = np.zeros(50, dtype=int)
        lengths = np.full(50, 100)
        p = coverage_profile(starts, lengths, 1000, bin_size=250)
        self.assertAlmostEqual(p["breadth_1x"], 0.1, delta=0.01)
        self.assertAlmostEqual(p["max_bin_frac"], 1.0, delta=0.01)

    def test_empty_coverage_does_not_divide_by_zero(self):
        p = coverage_profile(np.array([]), np.array([]), 1000)
        self.assertEqual(p["breadth_1x"], 0.0)
        self.assertEqual(p["max_bin_frac"], 0.0)
        self.assertEqual(p["depth_cv"], 0.0)


class TestBwaSamParsing(unittest.TestCase):
    """The comparison script must be tested BEFORE it is handed to someone whose
    sudo is needed to produce its input. An untested parser waiting on an
    install is a way to waste another person's time."""

    CATS = {"EXO_REF": "EXO", "HOST": "HOST", "HOST_NOLOCI": "HOST"}

    def _sam(self, records):
        import tempfile
        nl, tab = chr(10), chr(9)
        fh = tempfile.NamedTemporaryFile("w", suffix=".sam", delete=False,
                                         encoding="ascii", newline=nl)
        fh.write(f"@HD{tab}VN:1.6{nl}@SQ{tab}SN:EXO_REF{tab}LN:8510{nl}")
        for r in records:
            fh.write(tab.join(str(x) for x in r) + nl)
        fh.close()
        return Path(fh.name)

    @staticmethod
    def _rec(name, flag, ref, cigar="150M", tags=("AS:i:150",)):
        return [name, flag, ref, 1, 60, cigar, "=", 200, 350,
                "A" * 150, "I" * 150, *tags]

    def test_softclip_both_ends(self):
        self.assertEqual(softclip_len("10S130M10S"), 20)
        self.assertEqual(softclip_len("150M"), 0)
        self.assertEqual(softclip_len("20S130M"), 20)
        self.assertEqual(softclip_len("130M20S"), 20)
        self.assertEqual(softclip_len("10S60M5I65M10S"), 20)

    def test_cross_category_tie_is_flagged(self):
        sam = self._sam([
            self._rec("r1", 65, "EXO_REF", tags=("AS:i:150",)),
            self._rec("r1", 65 | 0x100, "HOST", tags=("AS:i:150",)),
        ])
        f = parse_sam(sam, self.CATS)["r1"]
        self.assertEqual(f["n_tied_top_categories"], 2)
        self.assertEqual(f["as_minus_xs"], 0)
        self.assertEqual(f["AS"], 150)

    def test_clear_winner_gives_a_positive_gap(self):
        sam = self._sam([
            self._rec("r2", 65, "EXO_REF", tags=("AS:i:150",)),
            self._rec("r2", 65 | 0x100, "HOST", tags=("AS:i:120",)),
        ])
        f = parse_sam(sam, self.CATS)["r2"]
        self.assertEqual(f["n_tied_top_categories"], 1)
        self.assertEqual(f["as_minus_xs"], 30)
        self.assertEqual(f["best_cat"], "EXO")

    def test_same_category_tie_is_not_ambiguous(self):
        """Two host references tied at the top still decide the CATEGORY."""
        sam = self._sam([
            self._rec("r3", 65, "HOST", tags=("AS:i:150", "XS:i:150")),
            self._rec("r3", 65 | 0x100, "HOST_NOLOCI", tags=("AS:i:150",)),
            self._rec("r3", 65 | 0x100, "EXO_REF", tags=("AS:i:100",)),
        ])
        f = parse_sam(sam, self.CATS)["r3"]
        self.assertEqual(f["n_tied_top_categories"], 1)
        self.assertTrue(f["competitor_category_known"])
        self.assertEqual(f["tie_by_score"], 1)

    def test_as_comes_from_the_primary_not_the_best_record(self):
        """SPECIFICATION CHANGE, 2026-09-16, and the reason matters.

        This used to assert the maximum AS across records for a reference. The
        primary record's AS is the right one: it is the alignment BWA chose and
        therefore the one a downstream tool consumes. BWA optimises the PAIR, so
        the primary can legitimately NOT be the individually best-scoring
        alignment -- which is exactly why AS - XS reaches -5 on real data
        (F073). Taking a maximum across records would paper over that.
        """
        sam = self._sam([
            self._rec("r4", 65, "EXO_REF", tags=("AS:i:120", "XS:i:100")),
            self._rec("r4", 65 | 0x100, "EXO_REF", tags=("AS:i:145",)),
            self._rec("r4", 65 | 0x100, "HOST", tags=("AS:i:100",)),
        ])
        f = parse_sam(sam, self.CATS)["r4"]
        self.assertEqual(f["AS"], 120)
        self.assertEqual(f["as_minus_xs"], 20)

    def test_xs_tag_is_preferred_over_cross_record_reconstruction(self):
        """The tag is BWA's own answer. Reconstructing it from records is what
        produced the artefact recorded in F071."""
        sam = self._sam([
            self._rec("r6", 65, "EXO_REF", tags=("AS:i:150", "XS:i:150")),
        ])
        f = parse_sam(sam, self.CATS)["r6"]
        self.assertEqual(f["XS"], 150)
        self.assertEqual(f["as_minus_xs"], 0)
        self.assertEqual(f["tie_by_score"], 1)

    def test_missing_xs_tag_falls_back_to_records(self):
        sam = self._sam([
            self._rec("r7", 65, "EXO_REF", tags=("AS:i:150",)),
            self._rec("r7", 65 | 0x100, "HOST", tags=("AS:i:110",)),
        ])
        f = parse_sam(sam, self.CATS)["r7"]
        self.assertEqual(f["XS"], 110)

    def test_mapq_zero_is_recorded_as_its_own_tie_proxy(self):
        """MAPQ 0 is the only tie signal a standard pipeline actually has, since
        the competitor's category is not in default BWA output (F069)."""
        rec = self._rec("r8", 65, "EXO_REF", tags=("AS:i:150", "XS:i:150"))
        rec[4] = 0
        f = parse_sam(self._sam([rec]), self.CATS)["r8"]
        self.assertEqual(f["tie_by_mapq0"], 1)
        self.assertFalse(f["competitor_category_known"])

    def test_xa_tag_supplies_the_competitor_category(self):
        sam = self._sam([
            self._rec("r9", 65, "EXO_REF",
                      tags=("AS:i:150", "XS:i:150", "XA:Z:HOST,+5000,150M,0;")),
        ])
        f = parse_sam(sam, self.CATS)["r9"]
        self.assertTrue(f["competitor_category_known"])
        self.assertEqual(f["competitor_categories"], ["HOST"])
        self.assertEqual(f["n_tied_top_categories"], 2)

    def test_unmapped_missing_tag_and_second_end_are_excluded(self):
        sam = self._sam([
            self._rec("unmapped", 77, "*"),
            self._rec("no_as", 65, "EXO_REF", tags=("NM:i:0",)),
            self._rec("second_end_only", 129, "EXO_REF", tags=("AS:i:150",)),
            self._rec("kept", 65, "EXO_REF", tags=("AS:i:150",)),
        ])
        got = parse_sam(sam, self.CATS)
        self.assertEqual(set(got), {"kept"})

    def test_softclip_comes_from_the_primary_record(self):
        sam = self._sam([
            self._rec("r5", 65, "EXO_REF", cigar="20S130M", tags=("AS:i:130",)),
            self._rec("r5", 65 | 0x100, "HOST", cigar="150M", tags=("AS:i:100",)),
        ])
        self.assertEqual(parse_sam(sam, self.CATS)["r5"]["softclip_len"], 20)


if __name__ == "__main__":
    unittest.main()
