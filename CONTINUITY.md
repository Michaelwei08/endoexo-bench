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
- F007 2026-09-14 [CODE]: DESIGN DEFECT found by the trial. Divergence is a
  two-point model (all loci at d, the strain at 0.02), so at d >= 0.10 alignment
  score alone separates the classes and every model reports PR AUC 1.0 / FPR 0.0.
  That is a simulation artefact, not a result. Divergence must become a
  distribution across loci and within the strain.
- F008 2026-09-14 [CODE]: DESIGN GAP found by the trial. no_decoy and full scored
  identically at every divergence, which is correct rather than a bug: the
  endogenous elements sit inside the host reference and the host copy is closer to
  the read than the decoy consensus, so an explicit decoy can never win. The
  simulation therefore had no way to show a decoy helping. The case where it can
  is an insertionally polymorphic, non-reference locus -- carried by the sample,
  absent from the assembly. Added panel modes poly_no_decoy and poly_full.

- F009 2026-09-14 [TOOL]: FIRST CONCLUSIVE RESULT, and it is a conditional. The
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
