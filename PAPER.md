# What the rules actually do: a ground-truth benchmark for exogenous versus endogenous retroviral read assignment

**Draft, 2026-09-17.** Not submitted. Every number traces to a file in
`results/`; every limitation is in section 6; the prior-art position is in
`PRIOR_ART.md` and is more conservative than this draft's framing might suggest
to a reader who skips it.

---

## Abstract

When a short read aligns to an exogenous retrovirus reference, it may have come
from the virus or from an endogenous retroelement in the host genome. Detection
pipelines resolve this with a small set of conventions -- competitive alignment
against a host reference, a mapping-quality floor, a unique-best test, an
aligned-length floor, and an ambiguous bin -- and those conventions are
reasonable, widely used, and to our knowledge have not been measured against
ground truth for this specific confusion. We built a simulation in which the
homology is constructed, so every read carries its true source, its true locus,
and the realized divergence of the window it came from, and measured what the
conventions buy and what they cost at three granularities: the read, the sample,
and the cohort.

Three results. First, mis-assignment separates into two mechanisms with opposite
remedies, and which one occurs is decided entirely by whether the read's true
source is present in the alignment panel: 99.8-100% of false calls are exact
score ties when it is, and 0.0% when it is not. Real BWA-MEM reproduces the tie
half exactly -- 100% of its false calls have `AS == XS` and MAPQ 0 -- so the
practical rule is to discard MAPQ 0, which gives 0.0% false calls at 80.2%
sensitivity against 4.8% at 90.5% for counting everything. Second, when the true
source is absent from the panel, no read-level, read-pair-level or sample-level
statistic we tested beats a raw read count, including the intra-host diversity
measure used for the neighbouring question of whether an element is endogenous at
all. Third, the discriminator that does work in that regime is not a per-sample
statistic: sorting host-side junction coordinates by how many individuals carry
them gives 96% specificity, because an endogenous insertion is shared by descent
while an exogenous integration is private to the individual.

We also report that evaluating any of this with a read-level cross-validation
split inflates a gradient-boosted model by up to +0.41 PR AUC, understates the
false-positive rate at fixed sensitivity by a factor of 2.8, and in the hardest
condition reports 2.4x chance where the honest locus-grouped answer is *below*
chance.

We claim no new idea. Each mechanism here has prior art, generally in a
neighbouring subfield; what we contribute is the measurement.

---

## 1. The question, and what is already known about it

An exogenous retrovirus and the endogenous retroelements of its host share
ancestry, and therefore sequence. A read from an endogenous locus can align to an
exogenous reference, and a pipeline counting reads on that reference will count
it. This is not a hypothetical: the conventions used to prevent it exist
precisely because it happens.

What is already established, and what this work does not claim:

- **The question of telling exogenous from endogenous retroviral sequence has
  been asked and answered by other routes.** Hayward et al. separate them by
  intra-host genetic variation -- an exogenous infection's sequences coalesce
  recently, an endogenous element is ancient -- and report the limitation that
  very recent endogenous elements have not accumulated enough change to be
  identified as endogenous. A more recent line separates them by epigenetic
  signature, using CpG depletion and altered trinucleotide frequencies.
- **Exact alignment-score ties, and the arbitrariness of resolving them, are
  documented.** MGmapper states that short reads may map equally well to more
  than one reference with identical alignment scores, and that the assignment is
  then arbitrary because no procedure is available to unravel the ambiguity.
- **The irreducible resolution limit that homology imposes has been computed for
  retroelements.** Umap/Bismap quantify single-read mappability as a function of
  read length, and ERVmancer computes the limit for HERVs specifically, reporting
  that at 75 bp they resolve only to an average of 12.7 internal clades above the
  leaf node and that 634 HERVH-LTR7 elements share on average 2.8% of their reads
  -- framed by those authors as ambiguity no alignment strategy can remove.
- **Adding absent sequence to a panel to soak up spurious mappings is standard.**
  This is the decoy-genome rationale, in use since the 1000 Genomes hs37d5 decoy.
- **Cross-patient recurrence as an artefact heuristic is standard in
  integration-site analysis.** That literature notes both that mispriming sites
  are sometimes reported at identical locations in samples from different
  patients, and that independent integration into an identical nucleotide
  position is statistically near-impossible.
- **Simulating an endogenous locus as absent from the assembly is established.**
  The 2022 Frontiers assessment of HERV-K insertion callers removes the
  proviruses from hg19 by masking and deletion and aligns simulated reads to the
  edited reference.

What we did not find, after two adversarial passes recorded in `PRIOR_ART.md`, is
any source that measures these quantities for this confusion against ground
truth. That gap is the paper. It is a modest one and we state it as such.

## 2. Why this can be given ground truth at all

The obstacle to benchmarking read assignment on real data is that the truth
requires knowing which genomic copy a read came from, which is the thing under
dispute. Adjudicating candidates by review reintroduces the judgement the
benchmark is supposed to test.

We therefore construct the homology. An ancestral provirus is built with HTLV-1
geometry -- two identical 755 bp long terminal repeats flanking a 7 kb internal
region. The exogenous reference, the exogenous strain the reads actually come
from, and an endogenous family inserted at host loci are all derived from that
ancestor at controlled substitution divergence. Every simulated read therefore
carries its true source, its true locus, and the realized divergence of the
window it came from, with no human in the loop.

**Divergence is a distribution, not a point, and the third component of that is
what makes the benchmark non-trivial.** Per-locus divergence is drawn
lognormally, because insertions in a genome are of different ages; a locus below
the family divergence reverts a random subset of the family-specific
substitutions, which is what a young element is. The strain is drawn from a
distribution, because a sample's virus is never the reference. And relative
substitution rates across sites are gamma distributed **at block resolution and
shared between lineages**, because functional constraint applies to regions
rather than to independent sites.

That last point is load-bearing and we found it the hard way. With independent
site rates, the mean of 150 of them concentrates on 1 -- coefficient of variation
1/sqrt(150 x shape), about 0.115 at shape 0.5 -- so every 150 bp read window
carries nearly the family average divergence and alignment score alone separates
the classes. Under that model the benchmark reported perfect separation at every
divergence above 0.10, which is an artefact and not a result. With regional rates
a read can land wholly inside a conserved window: at a median locus divergence of
0.20 the 5th percentile of endogenous read-window divergence is **0.000**, so a
read from a 27% diverged locus can be identical to the exogenous reference over
its whole length. The ablation is stark and is reported rather than asserted:
at median divergence 0.10, endogenous read windows falling inside the exogenous
range go from **76.9%** under the distributional model to **0.3%** under the
point model, and per-locus divergence collapses from a 0.038-0.119 range onto
0.099-0.101. The 5th percentile of endogenous window divergence is 0.0000 under
the distributional model at every median divergence up to 0.30, and 0.0200 to
0.2400 under the point model.

Realized divergences are measured from the generated sequences and reported. They
are never assumed from the parameter that produced them, because they differ from
it -- a defect we found and fixed is described in section 6.

**Inputs are public reference geometry and simulation only.** No patient or
controlled-access data is involved at any stage.

## 3. The task, and how it is scored

A false positive is born when a pipeline **counts** a read. So the task is scoped
to the reads whose best competitive hit is the exogenous reference, and asks
which of those actually came from the virus.

Baselines are the rules in production use, not strawmen: unique-best (`AS > XS`),
an `AS - XS` gap floor, an aligned-length floor, and no-close-alternative.

Two scoring decisions matter.

**The headline metric is false-positive rate at a fixed sensitivity**, not ROC
AUC, because a detection assay runs at one operating point and AUC averages over
thresholds nobody uses. Where a score is a small integer count, that metric
degenerates and we say so rather than reporting it (section 6).

**Cross-validation is grouped by locus.** A read-level split lets a model
memorise where a given endogenous locus lands on the viral reference, and the
size of that effect is one of our results.

We evaluate at three granularities and do not mix them: the read, the sample
(what an assay actually reports), and the cohort.

## 4. Results

### 4.1 Two mechanisms, decided by one binary property of the panel

Mis-assignment does not come in one flavour. It splits by whether the read's true
source is present in the alignment panel, and the two halves need opposite
remedies.

| | Exact tie (`AS == XS`) | Unopposed win (`AS > XS`, wrongly) |
|---|---|---|
| Cause | the panel holds an equally good explanation | the true source is absent from the panel |
| Occurs when | the endogenous locus is in the assembly | the locus is insertionally polymorphic |
| Detectable | free -- it is an equality test | no competition-based test can see it |
| Remedy | route to an ambiguous bin | see 4.4 |

The split is essentially total, and stable across three seeds and four
divergences. The fraction of false exogenous calls that are exact ties is
**99.5-100.0%** across individual runs on every panel containing the locus --
99.9-100.0% as seed means -- and **0.0 +- 0.0%** on every panel lacking it.

Two earlier observations re-attribute cleanly under this split. The unique-best
rule measuring FPR 0.000 in one condition is it rejecting every tie. The same rule
measuring FPR 1.000 on a viral-only panel is not a weak filter -- it is a filter
whose every error is an unopposed win, which no competition test can see.

**The mechanism of the tie half is that the two sequences are identical.** In the
reference-insertion case the surviving false positives are entirely
conserved-window reads: median realized window divergence 0.0000, maximum 0.0067,
100% below 0.05, at both divergence extremes. A read from a window where the
endogenous copy and the exogenous reference are the same sequence carries no
information about which it came from, so the tie is not a defect of the rule. It
is the correct report.

### 4.2 Real BWA-MEM reproduces the tie half exactly, and simplifies the rule

The above was measured with a first-party ungapped aligner, which left an obvious
way to be wrong: a gap buying a base somewhere would turn an exact tie into a
near-tie, and "free to detect" would become "needs a threshold".

Real BWA-MEM 0.7.19 on the same 8,623 read pairs does not do that.

- **100.0%** of false exogenous calls have `AS == XS` **exactly**, using BWA's own
  XS tag. The near-tie band 1..10 contains **zero** of them.
- **100.0%** of them also carry **MAPQ 0**, against 11.3% of true calls.
- All 26 have realized local divergence **0.0000**.

But the rule as we defined it is **not implementable from standard BWA output**.
Routing ties *between categories* needs the competitor's identity, and BWA's XS
tag gives the second-best score without naming the reference it was on. There
were zero `XA:Z` tags in the entire SAM even with `-h 200`. Standard output says a
read is ambiguous; it does not say what it is ambiguous with.

The implementable rule is category-agnostic and already present in every BAM:

| policy | calls | false | sensitivity |
|---|---|---|---|
| count everything | 539 | 4.8% | 0.905 |
| **discard MAPQ 0** | 455 | **0.0%** | **0.802** |
| three-bin by category | -- | not computable | -- |

Sensitivity is against all simulated exogenous reads, not against those that
happened to be called, because a true read lost to host assignment is a loss the
pipeline never sees. The category framing remains the *explanation* -- it is what
makes those reads interpretable -- but the rule a reader should take is one line.

Gapped alignment also costs less than the ungapped harness suggested: among true
calls the tied fraction is 28.3% under BWA against 40.3% under our aligner at the
same divergence, so the sensitivity price is 10 points rather than 40.

### 4.3 The sensitivity ceiling, and why it does not reach the detection call

The tie fraction among *true* exogenous reads is a hard bound on any
tie-rejecting rule, and it is computable in advance from the references alone. At
median divergence 0.02 with the host in the panel, **76.9 +- 5.5%** of genuine
viral reads are themselves ties over three seeds, so no such rule can exceed
**23.1%** read retention; at divergence 0.10 the bound is 59.7% and at 0.20 it is
**80.9%**. The measured three-bin retention tracks it: 0.226 +- 0.055, 0.594 and
0.804.

That sounds fatal and is not, because detection needs presence rather than
completeness. With a full locus catalogue at divergence 0.02 the rule retains
only 22.6% of true reads and **sample-level detection is still perfect**, ROC AUC
1.000 at FPR 0.000. A downward viral-load sweep locates where it starts to bind:
the filtered count holds AUC 1.000 and FPR 0.000 down to a **0.11-0.48x** viral
coverage band and breaks below 0.1x -- but at 0.02-0.1x an infected sample carries
0.6-2.7 viral read *pairs* in expectation, so the binding constraint there is
Poisson sampling of a handful of molecules, not the filter.

The ceiling therefore constrains **quantification** -- load, burden, clonality --
and not the detection call.

And the recommendation is worth most where detection matters most: at 0.11-0.48x
coverage the raw call count is near useless (AUC 0.643, FPR at 95% sensitivity
0.833) in the same regime where the filtered count is still perfect.

### 4.4 When the true source is absent from the panel, nothing per-sample works

This is the harder half and we report it as a negative result.

At the read level the learned model is **at or below chance**: grouped PR AUC
0.138 +- 0.010 against a prevalence floor of 0.173 (lift 0.80x) on a viral-only
panel at divergence 0.02, and 0.134 +- 0.008 against 0.174 on a
polymorphic-no-decoy panel. Ablating all seven mate and pair features changes
nothing there (FPR 0.944 versus 0.958), although the same ablation is worth a lot
where a competitor does exist (0.702 versus 0.967). The mate helps because it
lands a fragment away, possibly outside the conserved window, where the
competition is decidable; with nothing to compete against it carries no more
information than the read.

**Coverage breadth fails, and the mechanism predicted the failure.** Breadth
reaches AUC 0.953-0.987 where cross-mapping is confined to conserved windows, and
0.396-0.611 -- chance, at a null standard error of 0.075 -- where it is not. With
the true source absent, reads from *every* window of the element win unopposed, so
the coverage is broad rather than blocky. Uniformity statistics failed in both
forms of the mechanism, so "cross-mapping looks blocky" is simply the wrong
intuition for it.

**The intra-host diversity statistic does not rescue it either.** This was the one
published discriminator that needs no competitor in the panel, so it was the
candidate. Implemented as three pileup statistics over 60-sample cohorts, it
reaches AUC 0.619-0.774 on the panel-incomplete cells against a raw read count at
0.851-0.902. Its direction also reverses against the literature, and the reversal
is a finding about task identity rather than an error: Hayward et al. ask whether
a set of reads that all come from one candidate element is endogenous or
exogenous, so their pileup reflects coalescence age, whereas here the pileup is a
*mixture* -- an infected sample contributes a viral strain differing
systematically from the reference at its own sites plus endogenous cross-mappings
differing at theirs -- and two populations carry more variation than one.

It does earn a place as a model *feature* in that regime, worth +0.030 to +0.064
AUC and moving FPR at 95% sensitivity from 0.700 to 0.433, and nothing at all
(+0.000) where the panel is complete. Which gives a two-sided prescription: where
the panel is complete, use the single column, because the model loses to it (FPR
0.033 against 0.000); where it is not, use the model, because no single statistic
works but a combination does.

At low load *and* an incomplete panel, everything is at chance: raw count 0.641,
diversity 0.488-0.542, breadth 0.529, models 0.547 and 0.628. That corner is
unsolved by every statistic we tested at every granularity, and it coincides with
the limitation Hayward et al. state for their own method -- independent
corroboration by a different route, which we report as agreement.

### 4.5 What does work is in the cohort, not the sample

An endogenous insertion sits at the same host coordinate in everyone who
inherited it. An exogenous integration sits at a coordinate private to the
individual. Sorting host-side junction coordinates by how many individuals carry
them therefore separates the two **without telling any sequence apart**, which is
why it can work where every competition test is a no-op. The heuristic is
established in integration-site analysis; what we add is the measurement against
ground truth.

Over 50 samples at four clonal fractions:

- **The naive rule is not merely at chance, it is useless.** At a threshold of one
  junction bin it has sensitivity 1.000, FPR 1.000, specificity **0.00** -- it
  calls every sample positive, because every sample carries endogenous junctions.
  ROC AUC 0.507-0.594.
- **Recurrence filtering gives specificity 0.96, constant in clonal fraction** --
  exactly one uninfected sample in twenty-five, at every fraction tested.
- **Sensitivity is clonal-fraction limited**, not discrimination limited:
  0.040 / 0.360 / 0.600 / 0.680 at clonal fraction 0.1 / 0.3 / 0.6 / 1.0.

The observed 0.680 understates the method. Decomposing the gap at full clonality:
all 25 integrations are *detected*, 16 are private, **9 of 9** non-private ones
collide with *another sample's integration*, **0** with an endogenous locus, and
**0** are spurious. The expected count is 4.5 colliding pairs -- C(25,2) x 3/200,
since a plus-or-minus-one-bin tolerance makes each integration occupy three of the
200 bins in a 100 kb backbone -- which is about 9 samples, matching exactly. A 3 Gb
genome at this resolution has about 6 million bins and 0.00015 expected colliding
pairs, so the real sensitivity at full clonality is **25 of 25**.

The cohort machinery is validated independently: endogenous loci recur at their
population frequency, all ten within 0.08 of truth. The detector recovers each
locus's frequency from junction evidence alone.

An obvious refinement we have not implemented: a provirus has *two* boundaries
about 8.5 kb apart, so an integration should produce a **pair** of private
junctions at a known separation, and requiring the pair would use evidence this
rule discards.

### 4.6 Read-level cross-validation gives a different answer, not an optimistic one

This is methodological and it is not specific to this problem, but the size of
the effect here is worth reporting.

A read-level split inflates a gradient-boosted model by **+0.309 +- 0.009** PR AUC
over three seeds under the point-divergence model, and by up to **+0.410** under
the distributional one. It understates false-positive rate at 95% sensitivity in
**all sixteen** read-level cells, worst where the task is hardest (0.646 against a
true 0.944). And in the hardest cell it reports **2.4x chance** where the honest
locus-grouped answer is **below** chance -- so a leaky split there does not inflate
a real signal, it manufactures one.

Ablating two position features removes **92%** of the inflation, which identifies
the mechanism as positional memorisation: the model learns where each endogenous
locus lands on the viral reference. Those features also carry real signal under
grouped CV, so they cannot simply be dropped; they have to be evaluated grouped.

## 5. What this adds up to

The problem is a **reference-completeness** problem before it is a
discrimination problem. Which failure mode a pipeline suffers, which rules can
help, and whether any per-sample statistic can help at all are all determined by
one binary property: is the read's true source in the panel?

If it is, the confusion appears only as exact ties, the aligner already flags
them with MAPQ 0, discarding those reads removes every false call, and the price
is a read-retention cost computable in advance from the references and invisible
to the detection call at any load the assay can work at.

If it is not, no rule that compares alignments can see the problem, and neither
can any within-sample statistic we tested. The remedies are to complete the panel
-- which makes the case identical to the first, to three digits -- or to leave the
sample entirely and use cross-individual recurrence.

## 6. Limitations

Stated at length because several of them are load-bearing.

1. **Simulation only.** No claim is made about any real cohort. The BWA-MEM arm
   validates the harness against a production aligner; it does not validate the
   simulation against biology.
2. **The divergence model is substitution-only.** Indels, recombination and gene
   conversion between elements are not modelled at all. Real endogenous elements
   are full of indels, and indels break the window-alignment logic that the
   conserved-window mechanism rests on. This is the limitation most likely to
   change a number.
3. **The first-party aligner is ungapped**, and its numbers are harness
   validation rather than results. The claim the paper rests on no longer depends
   on it -- BWA reproduces the tie half exactly -- but every other read-level
   figure does.
4. **A simulation parameter did not deliver what it claimed, for a while.** Gamma
   site rates have mean 1 in expectation but not in sample: with 22 blocks across
   an 8.5 kb provirus the standard error of their mean is about 30%, and over 200
   seeds the realized mean ranged 0.345 to 1.912, so a requested divergence of
   0.10 was delivered as anything from 0.035 to 0.19. Rates are now normalised to
   mean exactly 1. The finding this affects is none -- realized divergences are
   measured and reported, and the results are comparisons within a run -- but the
   pre-fix result files carry the drift and are marked.
5. **One metric was wrong and is reported under a name that says so.** False-
   positive rate at 95% sensitivity is meaningless for the recurrence rule,
   because the score is a small integer count and demanding 95% sensitivity forces
   a threshold that admits everything.
6. **The junction-collision loss is a simulation-scale artefact** and the
   corrected figure is given alongside the observed one.
7. **Cohort sizes are small**: 50-60 samples, one cohort seed for the sample-level
   and cohort-level work, three seeds for the read-level work. The 12-sample smoke
   result for the diversity statistic did not survive 60 samples, which is the
   second time in this project a small-n result pointed the wrong way.
8. **No novelty claim.** `PRIOR_ART.md` records two adversarial passes in which
   five positioning claims were withdrawn. Every mechanism here has prior art. The
   one claim we still rest on -- the quantitative all-ties/no-ties dichotomy by
   panel completeness -- survived three targeted searches, but it is a
   search-based absence claim and is weaker than a positive result for that
   reason.
9. **The polymorphic, low-load corner is unsolved**, by every statistic tested at
   every granularity.

## 7. Code and data

`https://github.com/Michaelwei08/endoexo-bench`. Simulation and public reference
geometry only; no data use agreement is required to reproduce anything here.
Every quantitative claim traces to a file in `results/`, mapped in
`results/README.md`, which also separates the two incompatible generations of
result file described in limitation 4. Fifty tests over the invariants the
findings rest on run in under a second.
