# CONTINUITY.md

## Snapshot

- Goal: 2026-09-14 [USER]: An independent benchmark, publishable and public, for
  separating exogenous retroviral reads from endogenous retroelement cross-mappings
  at the alignment level. Deliverables are BOTH a benchmark/dataset paper and a
  public GitHub artifact.
- Now: 2026-09-14 [CODE]: Vertical slice runs end to end on simulated data:
  panel construction, read simulation, competitive ungapped alignment, feature
  extraction, production-rule baselines, learned discriminator, grouped vs
  read-level CV.
- Next: 2026-09-14 [CODE]: Re-run the prior-art check that motivates the project
  before any novelty wording. Replace the two-point divergence model with a
  distribution (F007). Replace the pure-Python aligner with real BWA-MEM on the same
  simulated FASTQs. Then re-test F009, which is the candidate headline result.

## Invariants / Constraints

- 2026-09-14 [USER]: CLEAN ROOM. This project inherits METHOD from the Alizadeh Lab
  viral-sequencing work and inherits NO DATA and NO UNPUBLISHED NUMBER from it.
  Never copy a cohort statistic out of that repo's modules or docstrings. Inputs are
  public reference sequence and simulation only.
- 2026-09-14 [CODE]: No patient data, no controlled-access data, ever. Nothing here
  requires a data use agreement, which is what makes the project independently owned.
- 2026-09-14 [CODE]: Ground truth comes from construction, not from review. If a
  label needs a human to adjudicate it, it does not belong in the core benchmark.
- 2026-09-14 [CODE]: Grouped cross-validation by locus is mandatory. A read-level
  split lets the model memorise a locus and is reported only as the leakage contrast.
- 2026-09-14 [CODE]: Do not claim novel / first without a fresh adversarial prior-art
  check. The motivating gap is five months old as of today.

## Decisions

- D001 ACTIVE 2026-09-14 [CODE]: Frame the task as a decision about reads a pipeline
  WOULD COUNT -- those whose best competitive hit is the exogenous reference -- not as
  read classification in general. That is where the false positive is actually born.
- D002 ACTIVE 2026-09-14 [CODE]: Baselines are the rules in production use, not
  strawmen: unique-best (AS > XS), an AS-XS gap floor, and an aligned-length floor.
- D003 ACTIVE 2026-09-14 [CODE]: Ship a pure-Python ungapped aligner for the design
  slice because no aligner binary is installable on this machine (no admin rights),
  and treat it as a harness validator, NOT as a source of publishable numbers. Scoring
  follows BWA-MEM (match +1, mismatch -4, local) and the seed length defaults to
  BWA-MEM's -k 19 so the seed-sensitivity limit matches.
- D004 ACTIVE 2026-09-14 [CODE]: Report FPR at a fixed sensitivity as the headline
  metric, not ROC AUC. A detection assay runs at one operating point.
- D005 ACTIVE 2026-09-14 [CODE]: Panel composition is an experimental axis, not a
  setting: viral_only / no_decoy / full. The first slice reproduced the known effect
  that adding host and decoy references drives cross-mapping to zero.
- D006 ACTIVE 2026-09-14 [USER]: Divergence is a DISTRIBUTION, in three respects:
  across loci (lognormal, drawn directly, with reversion of family substitutions
  for loci younger than the family), within the strain (lognormal), and across
  sites with AUTOCORRELATION (gamma at block resolution, shared between lineages
  because constraint is a property of the site and not of the lineage). The
  point-divergence model stays reachable via --rate-shape 1e6 --sigma-locus 0
  --sigma-strain 0, so the effect of the distributional model can be measured
  rather than asserted.
- D007 ACTIVE 2026-09-14 [CODE]: Report the realized per-read local divergence as a
  DIAGNOSTIC and never as a feature. It is the ground-truth difficulty of a read and
  nothing observable at inference time reveals it; using it would be a label leak.

## Findings from the first slice (2026-09-14)

All from simulation with the pure-Python ungapped aligner, so these are harness
validation, NOT publishable numbers. Config: 200 kb host filler, exogenous depth
20x, host depth 8x, read length 150, error 0.002, seed length 19.

- F001 2026-09-14 [TOOL]: Panel composition dominates divergence. With 12
  endogenous loci at divergence 0.10, a viral-only panel made 2175 of 2742
  exogenous calls false (79 percent); adding the host reference took it to ZERO.
  At divergence 0.02 the host reference took 2772/3339 (83 percent) down to
  70/637 (11 percent). This is the first ground-truth measurement of an effect
  the source lab had observed but could not quantify.
- F002 2026-09-14 [TOOL]: The unique-best rule (AS > XS) and the no-close-alternative
  rule are STRUCTURALLY VACUOUS on a single-reference panel: with no competitor,
  XS is 0, so every read passes and the measured FPR is 1.0. They are not weak
  specificity filters, they are readouts of the competition, and they measure
  nothing unless the competitor is in the panel. This reframes what those filters
  are for.
- F003 2026-09-14 [TOOL]: Read-level cross-validation inflates the gradient-boosted
  model by +0.309 PR AUC (sd 0.009 over seeds 42/142/242; grouped 0.337 vs
  read-level 0.646). At the operating point it is worse: FPR at 95 percent
  sensitivity reads 0.380 under a read-level split when the honest grouped answer
  is 0.706 -- an understatement of nearly a factor of two.
- F004 2026-09-14 [TOOL]: The leakage mechanism is positional memorisation.
  Dropping dist_to_terminus and in_ltr cuts the GBM inflation from +0.309 to
  +0.026, a 92 percent reduction. The model was learning where each endogenous
  locus lands on the viral reference.
- F005 2026-09-14 [TOOL]: Read-level CV does not only inflate a number, it SELECTS
  THE WRONG MODEL. Under grouped CV all four configurations land at 0.298 to 0.344
  PR AUC, indistinguishable. Under a read-level split the GBM with position
  features jumps to 0.646 and looks like the clear winner. Its entire apparent
  advantage was leakage.
- F006 2026-09-14 [TOOL]: In the hard regime nothing works well, and that is the
  honest headline. At divergence 0.02 on a viral-only panel the production rules
  cannot fire at all, and the best grouped model reaches 0.344 PR AUC against a
  0.077 prevalence floor -- about four times chance.
- F007 2026-09-14 [CODE]: DESIGN DEFECT found by the trial, now FIXED (see D006).
  Divergence was a two-point model (all loci at d, the strain at 0.02), so at
  d >= 0.10 alignment score alone separated the classes and every model reported
  PR AUC 1.0 / FPR 0.0. That was a simulation artefact, not a result.
- F011 2026-09-14 [TOOL]: The first attempted fix -- iid gamma site rates plus a
  lognormal per-locus top-up -- did NOT fix it, and measuring instead of assuming
  is what caught that. Endogenous read windows falling inside the exogenous range
  were 37.7 percent at median divergence 0.05 but only 0.8 percent at 0.10 and
  0.1 percent at 0.20. Two reasons, both structural:
    Averaging. The mean of 150 iid rates concentrates on 1 with a coefficient of
    variation of 1/sqrt(150*shape), about 0.115 at shape 0.5, so every 150 bp
    window carried nearly the family average no matter how heterogeneous the
    sites were.
    A hard floor. Each locus was built by adding a top-up to the family
    consensus, so no locus could be younger than the family and the across-loci
    spread came entirely from the top-up term.
- F012 2026-09-14 [TOOL]: Both fixed -- rates drawn at block resolution
  (autocorrelated, 400 bp default) and per-locus divergence drawn directly with
  reversion of family substitutions below the family value. Realized per-locus
  divergence at median 0.10 now spans 0.024 to 0.154, a six-fold range, and the
  5th percentile of endogenous read-window divergence is 0.000 at EVERY median
  divergence tested up to 0.30. A read off a 27 percent diverged locus can be
  identical to the exogenous reference over its whole length, which is the real
  phenomenon the first model could not produce. Class overlap is now 75.8 / 57.8
  / 44.2 / 34.5 percent at median divergence 0.05 / 0.10 / 0.20 / 0.30.
- F008 2026-09-14 [CODE]: DESIGN GAP found by the trial. no_decoy and full scored
  identically at every divergence, which is correct rather than a bug: the
  endogenous elements sit inside the host reference and the host copy is closer to
  the read than the decoy consensus, so an explicit decoy can never win. The
  simulation therefore had no way to show a decoy helping. The case where it can
  is an insertionally polymorphic, non-reference locus -- carried by the sample,
  absent from the assembly. Added panel modes poly_no_decoy and poly_full.

- F009 2026-09-14 [TOOL]: PARTLY SUPERSEDED by F013 the same day -- every zero in
  this entry was an artefact of the point-divergence model (F007). The conditional
  itself survived and was strengthened; the claim that any panel drives
  cross-mapping to ZERO did not. Kept unedited below for the audit trail.
  FIRST CONCLUSIVE RESULT, and it is a conditional. The
  value of an explicit endogenous decoy depends entirely on whether the endogenous
  locus is present in the reference assembly.
    Reference insertion (host assembly contains it): the host reference alone takes
    false calls to zero at divergence >= 0.05, and the explicit decoy adds nothing --
    0 versus 0 at every divergence.
    Non-reference / insertionally polymorphic locus (carried by the sample, absent
    from the assembly): the host reference is worth almost exactly NOTHING. False
    calls stay at 82.9 / 82.7 / 79.2 / 39.8 percent for divergence
    0.02 / 0.05 / 0.10 / 0.20, within a fraction of a percent of having no host
    reference at all (83.0 / 82.9 / 79.3 / 39.9). The explicit decoy is what rescues
    it: 18.1 percent at divergence 0.02 and ZERO at 0.05 and above.
  So panel specificity against endogenous cross-mapping is not a property of the
  panel alone. It is a joint property of the panel and of the cohort's insertion
  polymorphism, and that is a testable, practically actionable claim.
- F010 2026-09-14 [TOOL]: The gap-plus-length rule does real work only at high
  divergence. On a polymorphic-locus panel without a decoy its FPR falls from 0.995
  at divergence 0.02 to 0.112 at 0.20, while unique-best stays vacuous at 1.0
  throughout.

### F009 re-run under the distributional divergence model (5 panels x 4 divergences x 3 seeds)

- F013 2026-09-14 [TOOL]: F009 CORRECTED. What holds, strengthened:
    When the endogenous locus is ABSENT from the assembly, adding the entire host
    genome to the panel is worth about 0.1 percentage points. False call fraction
    82.6 / 81.6 / 80.0 / 76.2 percent against a viral-only panel's 82.7 / 81.8 /
    80.1 / 76.3 at divergence 0.02 / 0.05 / 0.10 / 0.20. Confirmed across three
    seeds under the harder model.
    The explicit decoy is exactly redundant when the loci ARE in the assembly:
    no_decoy and full agree to three digits in every one of the twelve cells.
    The decoy does help in the polymorphic case: 82.6 to 69.3 percent at
    divergence 0.02, 76.2 to 50.9 percent at 0.20.
  What is WITHDRAWN as a point-model artefact: the claim that any panel drives
  cross-mapping to zero. Under distributional divergence NO panel composition
  eliminates it at any divergence tested. The best cell is a host-containing panel
  at divergence 0.20, and 31.1 percent of pre-filter calls are still false there.
  The point model also made high divergence look uniformly easy: it put viral_only
  at 39.9 percent false at divergence 0.20 where the distributional model puts it
  at 76.3.
- F014 2026-09-14 [TOOL]: F007 is resolved. No cell reports PR AUC 1.0 or FPR 0.0
  any more, and grouped performance now rises smoothly with divergence rather than
  saturating (viral_only grouped PR AUC 0.138 / 0.215 / 0.345 / 0.704).
- F015 2026-09-14 [TOOL]: The production rules fail on SENSITIVITY, not specificity,
  and that is the opposite of what was assumed. With the host in the panel,
  unique-best reaches FPR 0.000 -- but keeps only 34.7 percent of true exogenous
  reads at divergence 0.02 and 78.1 percent at 0.20, and gap-plus-length keeps only
  16.6 percent. An assay on those rules loses most of its genuine signal on a
  divergent strain.
- F016 2026-09-14 [TOOL]: MECHANISM, measured rather than argued. In the
  reference-insertion case the false positives that survive are ENTIRELY
  conserved-window reads: median realized window divergence 0.0000, max 0.0067,
  100 percent below 0.05, at both divergence extremes. In a conserved window the
  exogenous reference and the endogenous copy are the SAME SEQUENCE, so the two
  hypotheses are indistinguishable from the read alone; the unique-best rule
  resolves the tie by rejecting both, which is exactly where its sensitivity goes.
  No read-level feature can beat this, which is why the grouped model sits at or
  below chance in the hardest cells. The way out has to be information the read
  does not carry: the mate, coverage breadth across NON-conserved windows, and
  host-virus junction evidence.
  This independently derives the rubric that the motivating prior-art audit could
  only label "Supported + inference" -- require unique viral sequence or junction
  evidence. The benchmark now supplies the reason.
  The polymorphic case fails differently: there the surviving false positives are
  only 71.5 percent conserved-window reads at divergence 0.20, with window
  divergence up to 0.6118, because a genuinely divergent read can still win a
  competition that does not contain its true source.
- F017 2026-09-14 [TOOL]: The leakage result is confirmed, larger, and now has its
  sharpest form. Read-level CV inflates the GBM by up to +0.410 PR AUC (viral_only
  at divergence 0.10: grouped 0.345, read-level 0.750), and the inflation tracks
  task difficulty -- it collapses to +0.012 in the easiest cell. At the operating
  point, read-level reports FPR at 95 percent sensitivity of 0.273 where the honest
  grouped answer is 0.767, a 2.8-fold understatement.
  The sharpest form: in the hardest cell the grouped model is BELOW CHANCE
  (PR AUC 0.138 against a 0.173 prevalence floor, lift 0.80x) while the read-level
  split reports 2.43x chance. A leaky split does not merely inflate a real signal
  here -- it manufactures one where the honest answer is "worse than random".
- F018 2026-09-14 [CODE]: PR AUC is not comparable across these cells. Prevalence
  ranges from 0.173 to 0.689 because panel composition changes how many reads are
  called at all, and the PR AUC floor IS the prevalence. Always quote lift over
  prevalence alongside it, or quote FPR at fixed sensitivity instead.

## State

### Done

- 2026-09-14 [CODE]: Wrote endoexo/simulate.py, endoexo/align.py, endoexo/evaluate.py,
  run_slice.py and run_leakage_probe.py; 842 lines, one point runs in 4 to 27 s.
- 2026-09-14 [TOOL]: Ran the divergence x panel sweep (results_sweep1.json) and the
  three-seed leakage probe (results_leakage_probe.json). See Findings above.
- 2026-09-14 [CODE]: Added poly_no_decoy / poly_full panel modes for non-reference
  insertions, per F008, and ran them (results_poly.json). See F009 -- the gap that
  the trial exposed turned out to hold the project's first real result.

### Next

- 2026-09-14 [CODE]: Fresh adversarial prior-art check; record it as a claim matrix.
- 2026-09-14 [CODE]: Real BWA-MEM + samtools path, run on a machine that has them.
- 2026-09-14 [CODE]: Add indels to the divergence model; the current aligner is ungapped.
- 2026-09-14 [CODE]: Replace the two-point divergence model with a distribution, per F007.
  Until that is done, treat every divergence >= 0.10 result as uninformative.
- 2026-09-14 [CODE]: Add a real-data held-out evaluation arm from public HTLV-1
  cell-line WGS (PRJDB11215) rather than simulation only.
