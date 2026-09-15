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

### B: real BWA-MEM reproduction (2026-09-15, scaffolded, BLOCKED on one user action)

- D013 ACTIVE 2026-09-15 [CODE]: The BWA comparison runs the SAME simulated reads
  through both aligners. export_for_bwa.py writes panel.fa, paired FASTQ and a
  separate truth.tsv; run_bwa.sh aligns with `bwa mem -a -k 19`; compare_with_bwa.py
  reproduces the tie typology from the SAM. -a is not optional: a read tied between
  the exogenous reference and the host is reported against one of them arbitrarily,
  with the other present only as a secondary record, so the category tie structure
  cannot be recovered from primary records alone.
- D014 ACTIVE 2026-09-15 [CODE]: Truth lives in a separate file, never in the read
  name. A name that carries its own label is a label leak waiting for the first
  person who parses it.
- F043 2026-09-15 [ASSUMPTION]: PREDICTION ON RECORD, before the numbers exist. The
  two aligners break ties differently and both are defensible: the first-party one
  sorts by score alone so a tie resolves to whichever reference entered the panel
  first, in practice the exogenous one, while BWA picks a primary essentially
  arbitrarily and will therefore report roughly half the tied reads against the host.
  So the award_ties_to_exo call count MUST differ, with the first-party number the
  pessimistic bound, and the three_bin_ambiguous numbers SHOULD agree because that
  policy excludes ties under either tie-break. If they agree, that is a further
  argument for the three-bin policy: it is the only policy whose output does not
  depend on an arbitrary implementation choice inside the aligner.
- F044 2026-09-15 [TOOL]: The decisive number is not the tied fraction, it is the
  distribution of AS - XS among false calls. Under the ungapped aligner it is a spike
  at exactly 0. If gapped alignment spreads it to 1 or 2, "free to detect" becomes
  "needs a threshold" and the binary decomposition softens into a tuning problem.
  compare_with_bwa.py prints that distribution explicitly for that reason.
- 2026-09-15 [TOOL]: Export verified: 8,623 read pairs, 2 references, 310,622 bp,
  5.8 MB. Seven tests added for the SAM parser, because handing an untested parser to
  someone whose sudo is needed to produce its input wastes their time. 33 tests pass.
- 2026-09-15 [USER-ACTION-REQUIRED]: bwa 0.7.17-7 and samtools are in the WSL2 Ubuntu
  24.04 apt repository but not installed, and installing needs a password this session
  cannot supply:
      wsl -e bash -lc 'sudo apt-get update && sudo apt-get install -y bwa samtools'
  After that, everything else is two commands and needs no further decisions.

### Sample-level detection task (2026-09-15, 60 samples, cohort seed 42, load 0.56-7.48x)

- F030 2026-09-15 [TOOL]: THE REFERENCE-INSERTION CASE IS SOLVED, and the fix is
  free. Applying the three-bin rule BEFORE counting, rather than counting every
  read assigned to the exogenous reference, moves sample-level FPR at 95 percent
  sensitivity from 0.633 to 0.000 at divergence 0.02 and from 0.833 to 0.000 at
  0.10. ROC AUC 0.768 to 0.999 and 0.801 to 1.000. This is the prescription the
  benchmark exists to produce: the statistic in production use is a read count,
  and the same count restricted to unopposed wins is near-perfect.
- F031 2026-09-15 [TOOL]: The learned model does NOT beat the right single
  statistic, and saying so is the point. Logistic regression over all 14 sample
  features reaches FPR 0.067 and 0.033 where the unopposed-call count alone reaches
  0.000 and 0.000. With 14 features and 60 samples the model is a worse estimator
  of a quantity one column already carries. A benchmark whose recommendation is
  "use this column" is a better result than one whose recommendation is "fit this".
- F032 2026-09-15 [TOOL]: MY PREDICTION FAILED, and the theory that predicted the
  failure was already in these notes. I expected coverage breadth to break the
  read-level ceiling in the polymorphic case. It does not: breadth ROC AUC is 0.396
  at divergence 0.02 and 0.611 at 0.10, both within 1.5 standard errors of chance
  (null SE 0.075 at 30 versus 30). The reason is in F016, which I wrote and then
  failed to apply: the conserved-window argument belongs to the TIE mechanism only.
  When the true source is absent from the panel, reads from EVERY window of the
  element win unopposed, not just conserved ones, so cross-mapping coverage is
  broad rather than blocky and breadth has nothing to separate.
  The contrast confirms the mechanism rather than the hypothesis: breadth reaches
  AUC 0.953 and 0.987 in the reference-insertion case, where cross-mapping IS
  confined to conserved windows.
  Uniformity failed everywhere it was tested -- max_bin_frac AUC 0.596 and 0.496 in
  the polymorphic case -- so "cross-mapping looks blocky" was the wrong intuition
  for this mechanism in both of its forms.
- F033 2026-09-15 [TOOL]: The polymorphic case is now unsolved at ALL THREE
  granularities: read level at or below chance (F026), read pair no better (F025),
  and sample level topping out at ROC AUC 0.899 / 0.720 with FPR at 95 percent
  sensitivity of 0.467 / 0.833 -- from the raw call count, with nothing beating it.
  The three-bin count is identical to the raw count there, exactly as F024 predicts.
- F034 2026-09-15 [ASSUMPTION]: Unexplained and worth a second cohort seed before it
  is quoted -- the raw count discriminates BETTER at divergence 0.02 (AUC 0.899)
  than at 0.10 (0.720), which is the opposite of the read-level ordering. The
  hypothesis is that a count's discriminating power depends on the VARIANCE of the
  cross-mapping background across samples rather than its magnitude, and per-locus
  divergence spreads more in absolute terms at higher median divergence. One cohort
  seed cannot support that.

### Complete locus catalogue: the polymorphic case is a PANEL problem (2026-09-15)

- F035 2026-09-15 [TOOL]: THE PRESCRIPTION IS COMPLETE. Putting the actual
  polymorphic locus sequences in the panel -- a complete population catalogue,
  which is the upper bound on what panel design can buy -- converts 99.8 to 100.0
  percent of the unopposed wins into detectable exact ties, at divergence 0.02 /
  0.10 / 0.20 over three seeds. Against the same panel without the catalogue the
  tied fraction is 0.0 percent. The three-bin rule, a complete no-op there (F024),
  then removes essentially all of them: false call fraction 0.8 / 0.2 / 0.0 percent
  against 82.6 percent before.
  At sample level the effect is total. The unopposed-call count goes from ROC AUC
  0.899 with FPR at 95 percent sensitivity 0.467 (no catalogue, divergence 0.02) to
  ROC AUC 1.000 and FPR 0.000, and the same at divergence 0.10.
  So the unopposed-win mechanism is not an analysis problem. It is a panel
  COMPLETENESS problem, and completing the panel reduces it to the mechanism the
  field already knows how to handle.
- F036 2026-09-15 [TOOL]: The read-level sensitivity ceiling (F021) does NOT
  propagate to the sample-level decision at these viral loads, and that materially
  softens F015. With the catalogue at divergence 0.02 the three-bin rule retains
  only 0.226 of true exogenous reads -- and sample-level detection is still
  PERFECT, ROC AUC 1.000 and FPR 0.000. Detection needs presence, not completeness;
  throwing away three quarters of the true reads costs nothing when the question is
  whether the sample is infected at all.
  The ceiling should therefore be quoted as a constraint on QUANTIFICATION -- viral
  load, clonality, burden -- and not on detection. Load range here was 0.56 to 7.48
  x; a downward load sweep is running to find where it starts to bind.
- F037 2026-09-15 [TOOL]: A bigger panel costs more sensitivity, as it must. The
  catalogue adds competitors, so more TRUE reads become ties: 76.9 percent of true
  reads are tied at divergence 0.02 with the catalogue. Panel completeness and read
  retention trade against each other directly, and the trade is worth making only
  because of F036.
- F038 2026-09-15 [CODE]: Comparability caveat. results_catalogue_* were produced
  AFTER the rate normalisation of F028 and results_typology BEFORE it, so the
  sensitivity figures are not directly comparable across those files. The tied
  FRACTION comparison is safe because 0.0 percent tied without the locus in the
  panel is structural rather than seed-dependent -- a tie is impossible when no
  competitor exists. A post-fix no_decoy re-run is in flight to close the gap.

### Load sweep and the post-fix comparison (2026-09-15, closes F038)

- F039 2026-09-15 [TOOL]: A COMPLETE CATALOGUE IS EQUIVALENT TO HAVING THE LOCI IN
  THE ASSEMBLY, to three digits. Post-normalisation no_decoy against poly_catalogue
  at divergence 0.02 / 0.10 / 0.20: true tied fraction 76.9 / 40.3 / 19.1 percent
  versus 76.9 / 40.3 / 19.2, three-bin retained sensitivity 0.226 / 0.594 / 0.804
  versus 0.226 / 0.594 / 0.804. The panel does not care HOW the read's true source
  got into it, which is the cleanest possible confirmation that the two-mechanism
  framing is the right one -- everything reduces to presence or absence of the true
  source, and nothing else about the panel matters.
- F040 2026-09-15 [TOOL]: F021's ceiling figure is RESTATED and it moved against us.
  Post-normalisation, 76.9 percent of true exogenous reads are ties at divergence
  0.02, not 65.3, so the ceiling on any tie-rejecting rule is 23.1 percent and not
  34.7. The rate drift of F028 had deflated it. The effect is confined to low
  divergence: at 0.10 the tied fraction moved 39.4 to 40.3 percent and at 0.20 21.9
  to 19.1. F038's comparability caveat is closed -- all read-level tie and
  sensitivity figures should now be quoted from results_nodecoy_postfix and
  results_catalogue_read.
- F041 2026-09-15 [TOOL]: THE CEILING BINDS ONLY WHERE THE ASSAY HAS NO MATERIAL.
  Sweeping viral load downward at divergence 0.10, the three-bin count holds ROC AUC
  1.000 and FPR at 95 percent sensitivity 0.000 over a 0.107 to 0.481x load band,
  and breaks at 0.021 to 0.096x -- AUC 0.803, FPR 1.000.
  The break is not primarily the filter. At 0.02 to 0.1x coverage an infected sample
  contains 0.6 to 2.7 exogenous read PAIRS in expectation, so many infected samples
  carry zero viral reads before any filter runs. The dominant limit there is Poisson
  sampling of a handful of molecules; the ceiling is what turns one or two present
  reads into zero retained. This strengthens F036 rather than qualifying it: the
  sensitivity cost is invisible to the detection call everywhere the assay has
  anything to detect.
- F042 2026-09-15 [TOOL]: THE RECOMMENDATION IS WORTH MOST EXACTLY WHERE DETECTION
  MATTERS MOST. The gap between the production statistic and the three-bin count
  widens as load falls. At 0.107 to 0.481x load, the raw call count reaches FPR at
  95 percent sensitivity of 0.833 while the unopposed count reaches 0.000 -- the raw
  count is near useless (ROC AUC 0.643) in the regime where the three-bin count is
  still perfect. At the original 0.56 to 7.48x band the same contrast was 0.833
  against 0.000 on AUC 0.801 versus 1.000, so the count degrades with load and the
  filtered count does not, until the material runs out.

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
