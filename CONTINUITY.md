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
- D008 ACTIVE 2026-09-14 [CODE]: The aligner keeps the OPTIMISTIC tie-break -- a tie
  resolves to whichever reference was added first, in practice the exogenous one --
  and exposes the tie instead of hiding it, so the POLICY is the caller's decision
  and its cost is measurable. Every run reports both policies: award_ties_to_exo and
  three_bin_ambiguous. This is how a silently inflated exogenous call count gets
  caught, and it is what made F019 to F022 visible at all. A real aligner breaks the
  tie arbitrarily and reports mapping quality 0, so neither policy is "the" truth;
  reporting one without the other is the mistake.
- D009 ACTIVE 2026-09-14 [CODE]: Never quote a pre-filter false-call fraction on its
  own. It is policy-dependent, by up to the entire effect: the reference-insertion
  case reads 62.8 percent false under one tie policy and 0.0 percent under the
  other. Quote the tie typology or a post-policy number.
- D010 ACTIVE 2026-09-15 [CODE]: Site rates are normalised to mean EXACTLY 1, so the
  `divergence` argument means what its name says. See F028 for why.
- D011 ACTIVE 2026-09-15 [CODE]: The repository carries a stdlib unittest suite,
  `python -m unittest discover -s tests -t .`, and it tests the INVARIANTS the
  findings rest on rather than surface behaviour: that the aligner's vectorised
  max-subarray matches a brute-force reference, that the tie structure counts
  categories and not references, that a locus below the family divergence can exist
  at all (the F011 defect), that the documented point-model ablation actually
  flattens the rates (a defect that shipped in the README), and that the cohort seed
  fixes the virus while the sample seed redraws the host.
- D012 ACTIVE 2026-09-15 [CODE]: The task is evaluated at TWO granularities. Read
  level is bounded by the conserved-window ceiling and is where F026's negative
  result lives. Sample level is what a detection assay actually reports and is the
  only place an aggregate statistic can break that ceiling. Neither replaces the
  other; a result at one granularity may not be quoted as a result at the other.

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

### Test suite, and the defect it found (2026-09-15)

- F028 2026-09-15 [TOOL]: The test suite found a real defect on its first run, and
  it was in the code rather than the test. Gamma site rates have mean 1 in
  EXPECTATION, but the sample mean over the blocks does not: at block_len 400 over
  an 8.5 kb provirus there are only 22 blocks, so the standard error of their mean
  is about 30 percent. Measured over 200 seeds the realized mean rate had sd 0.298
  and ranged 0.345 to 1.912 -- meaning a requested divergence of 0.10 was silently
  delivered as anything from 0.035 to 0.19, by a factor of five, seed to seed.
  Rates are now normalised to mean exactly 1 (D010).
  This does not invalidate any finding: every realized divergence is measured and
  reported (D006), and the findings are comparisons within a run. It does mean the
  realized divergences in results_sweep1 / results_f009_dist / results_typology
  carry that drift, and results_sample_level was launched before the fix.
- F029 2026-09-15 [CODE]: 26 tests, all passing, about 0.12 s. The suite is small on
  purpose: every test corresponds to an invariant a finding depends on or to a
  defect that actually happened here.

### F009 re-run under the distributional divergence model (5 panels x 4 divergences x 3 seeds)

- F013 2026-09-14 [TOOL]: ONE CLAUSE SUPERSEDED by F020 the same day -- the
  "no panel composition eliminates it" sentence counted exact ties as exogenous
  calls and does not survive. The conditional and the redundancy result stand.
  Kept unedited for the audit trail. F009 CORRECTED. What holds, strengthened:
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

### Tie typology (2026-09-14, seed 42; the sweep confirmation is in results_typology.json)

- F019 2026-09-14 [TOOL]: There are TWO failure mechanisms, not one, and they need
  opposite remedies. Measured fraction of false exogenous calls that are EXACT ties
  (AS == XS) at divergence 0.02 / 0.20:
    viral_only          0.0 / 0.0 percent
    host lacks loci     0.0 / 0.0 percent
    host contains loci  100.0 / 100.0 percent
    host lacks + decoy  85.8 / 59.9 percent
  A TIE means the panel holds an equally good explanation. It costs nothing to
  detect, and the unique-best rule removes every one of them -- which is the whole
  reason F015 measured FPR 0.000 there.
  An UNOPPOSED WIN means the read beat everything in the panel with a positive
  margin. No competition-based test can touch it, which is the real reason
  unique-best measures FPR 1.000 on a viral-only panel. F002 called that rule
  vacuous for lack of a competitor; the deeper statement is that its errors there
  are all unopposed wins.
- F020 2026-09-14 [TOOL]: F013 CORRECTED, the second correction of this result
  today. Its claim that no panel composition eliminates cross-mapping rested on
  awarding every tie to the virus. Since 100 percent of the false calls in the
  reference-insertion case are ties, that case IS fully removable: 0 false calls,
  already visible as FPR 0.000 in F015. The honest statement is that cross-mapping
  is eliminable exactly when the panel contains the read's true source, and the
  price is sensitivity rather than specificity.
  The claim does survive for the polymorphic case, where 14.2 percent of false
  calls at divergence 0.02 and 40.1 percent at 0.20 are unopposed wins.
- F021 2026-09-14 [TOOL]: The sensitivity ceiling is a quantity, not a worry. The
  tie fraction among TRUE exogenous reads is a hard upper bound on any
  tie-rejecting rule: 55.5 percent of genuine viral reads are ties at divergence
  0.02 with the host in the panel, so no such rule can exceed 44.5 percent
  sensitivity. At divergence 0.20 the bound is 80.8 percent. The ceiling is set by
  how much of the provirus is conserved relative to the host endogenous copy, and
  it is measurable in advance from the references alone -- no cohort needed.
- F022 2026-09-14 [TOOL]: A CONSENSUS decoy degrades as the locus drifts from the
  consensus, which is a practical design point. It converts unopposed wins into
  ties well at divergence 0.02 (85.8 percent tied) and poorly at 0.20 (59.9
  percent), because a family consensus is a bad proxy for one specific old locus.
  The implication is that the right decoy for polymorphic insertions is a
  LOCUS-RESOLVED reference, not a family consensus.

### Sweep confirmation, 5 panels x 4 divergences x 3 seeds (results_typology.json)

- F023 2026-09-14 [TOOL]: F019 to F022 confirmed. Tie fraction among FALSE calls is
  0.0 +- 0.0 percent on every viral-only and polymorphic-no-decoy cell and 99.8 to
  100.0 percent on every reference-insertion cell, at all four divergences. The
  mechanism split is not a seed accident.
- F024 2026-09-14 [TOOL]: THE SHARPEST RESULT OF THE SESSION. The three-bin rule --
  route a tie to an ambiguous bin -- is a COMPLETE NO-OP when the endogenous locus
  is absent from the panel. False call fraction 82.7 percent before and 82.7 percent
  after, sensitivity 1.000 before and 1.000 after, identical at every divergence,
  because there are no ties to bin. The rule that the field treats as the answer to
  exogenous/endogenous ambiguity does nothing at all in the insertionally
  polymorphic case.
  Where the locus IS in the panel it is near-perfect on specificity and expensive:
  false calls 62.8 to 0.0 percent at divergence 0.02, but sensitivity 0.989 to
  0.343. At divergence 0.20, 31.1 to 0.0 percent for sensitivity 0.999 to 0.780.
  With a consensus decoy it is partial: 69.3 to 38.5 percent at divergence 0.02,
  sensitivity 1.000 to 0.617.
  This is the direct measurement of the three-bin decision rule that the motivating
  prior-art audit could only label a project-specific design choice.
- F025 2026-09-14 [TOOL]: The mate is the remedy for the tie mechanism, and ONLY for
  it. Ablating all seven mate and pair features from the grouped GBM:
    where a competitor is present, FPR at 95 percent sensitivity degrades from
    0.702 to 0.967 at divergence 0.02 and from 0.189 to 0.887 at 0.20
    (reference insertion), and 0.891 to 0.980 / 0.297 to 0.772 with a consensus
    decoy. PR AUC gain from the mate peaks at +0.174.
    where NO competitor is present, the mate buys +0.005 to +0.007 PR AUC at
    divergence 0.02 and FPR is unchanged within noise (0.944 versus 0.958).
  The mechanism is now explicit: the mate helps because it lands a fragment away
  from the read, possibly outside the conserved window, where the COMPETITION is
  decidable. With nothing to compete against, the mate has no more information than
  the read. The +0.121 gain in viral_only at divergence 0.20 comes through a
  different channel -- the mate's absolute score, not its margin.
- F026 2026-09-14 [TOOL]: THE OPEN PROBLEM, stated as a negative result. An
  unopposed win against a young endogenous element is not solvable by anything at
  the read-pair level. Grouped GBM with the full feature set reaches PR AUC 0.138
  +- 0.010 against a 0.173 prevalence floor (lift 0.80x) on viral_only at divergence
  0.02, and 0.134 +- 0.008 against 0.174 (lift 0.77x) on polymorphic-no-decoy --
  BELOW CHANCE, over three seeds, with FPR at 95 percent sensitivity of 0.944 and
  0.962. The three-bin rule is a no-op there (F024) and the mate is worthless there
  (F025). The only routes left are above the read pair: coverage breadth over
  NON-conserved windows of the reference, and host-virus junction evidence. That is
  the next build.
- F027 2026-09-14 [TOOL]: Read-level CV understates FPR at 95 percent sensitivity in
  all sixteen cells, and worst where the task is hardest: 0.646 against a true 0.944
  on viral_only at divergence 0.02, 0.055 against 0.189 on reference-insertion at
  0.20. Quoting a read-level operating point is not a small optimism; it is a
  different answer.

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
