# endoexo-bench

A ground-truth benchmark for a question that viral-detection pipelines answer by
convention rather than by measurement: **when a read aligns to an exogenous
retrovirus reference, did it come from the virus, or from an endogenous
retroelement in the host genome?**

Prior work attacks this question with discriminators orthogonal to alignment:
intra-host genetic diversity, which separates a recently coalescing exogenous
infection from an ancient endogenous element, and epigenetic signature. At the
alignment level, ERVmancer resolves read ambiguity *between* HERV loci by
phylogenetic placement, and quantifies the resolution limit that sequence
homology imposes.

What has not been measured, as far as the prior-art check in
[PRIOR_ART.md](PRIOR_ART.md) establishes, is what the rules actually running in
production do: competitive alignment against a host reference, a mapping-quality
floor, a unique-best `AS > XS` test, an aligned-length floor, and an ambiguous
bin. **No claim of novelty is made here** -- four of the five things this
project was originally positioned on turned out to be published, and the claim
matrix records which.

## Why this can be benchmarked at all

Ground truth comes from construction. An ancestral provirus is built with
HTLV-1 geometry (two identical 755 bp LTRs flanking an internal region). The
exogenous reference, the exogenous strain the reads actually come from, and an
endogenous family inserted at many host loci are all derived from it at
controlled substitution divergence. Every simulated read therefore carries its
true source, its true locus, and the realized divergence of the window it came
from.

**Divergence is a distribution, not a point**, in three respects, and the third
is the one that makes the benchmark hard:

1. **Across loci.** Each insertion draws its own divergence lognormally, because
   insertions in a real genome are of different ages. A locus below the family
   divergence reverts a random subset of the family-specific substitutions,
   which is what a young element looks like.
2. **Within the strain.** The sample virus is drawn from a distribution, because
   it is never the reference.
3. **Across sites, with autocorrelation.** Relative substitution rates are gamma
   distributed *at block resolution* and shared between lineages, because
   functional constraint applies to regions -- a gene, a domain, an LTR
   subregion -- not to independent sites. This matters more than it sounds. With
   iid site rates the mean of 150 of them concentrates on 1 (CV about 0.115 at
   shape 0.5), so every read window carries nearly the family average and
   alignment score alone separates the classes. With regional rates a read can
   land wholly inside a conserved window: at a median locus divergence of 0.20
   the 5th percentile of endogenous read-window divergence is **0.000**, so a
   read off a 27 percent diverged locus can be identical to the exogenous
   reference over its whole length.

Realized divergences are measured from the generated sequences and reported;
they are never assumed from the parameter that produced them.

**Inputs are public reference sequence and simulation only.** Nothing here needs
a data use agreement, and no patient data is involved at any stage.

## The task

Take the reads a detection pipeline **would count** -- those whose best
competitive hit is the exogenous reference -- and separate true exogenous reads
from endogenous cross-mappings.

- **Baselines are the rules in production use**, not strawmen.
- **Headline metric is FPR at a fixed sensitivity**, because a detection assay
  runs at one operating point. ROC AUC averages over thresholds nobody uses.
- **Cross-validation is grouped by locus.** A read-level split lets a model
  memorise where a given endogenous locus lands on the viral reference; the
  read-level number is reported only as the leakage contrast.

### Two failure mechanisms, which need opposite remedies

Mis-assignment does not come in one flavour, and separating the two is most of
what this benchmark measures.

| | **Exact tie** (`AS == XS`) | **Unopposed win** (`AS > XS`, wrongly) |
|---|---|---|
| Cause | the panel holds an equally good explanation | the read's true source is absent from the panel |
| When | the endogenous locus is in the assembly | the locus is insertionally polymorphic |
| Detectable? | free -- it is an equality test | no competition-based test can see it |
| Remedy | route to an ambiguous bin | breadth, mate, junction evidence |
| Cost | sensitivity, bounded by the tie fraction among *true* reads | -- |

A read inside a window where the endogenous copy and the exogenous reference
are the *same sequence* is information-free on its own, so the tie is not a
defect of the rule -- it is the correct report. What it costs is quantifiable in
advance from the references alone: the fraction of true exogenous reads that are
also ties is a hard ceiling on any tie-rejecting rule.

### Validated against real BWA-MEM -- and the rule gets simpler

The whole tie story was measured with a first-party **ungapped** aligner, which
left one obvious way for it to be wrong: a gap buying a base somewhere would turn
an exact tie into a near-tie, and "free to detect" would become "needs a
threshold". Real BWA-MEM 0.7.19 on the same 8,623 read pairs says otherwise.

- **100.0%** of false exogenous calls have `AS == XS` **exactly**, using BWA's own
  XS tag. The near-tie band 1..10 contains **zero** of them. The exactness
  survives gapped alignment.
- **100.0%** of them also carry **MAPQ 0** -- BWA flags every one itself.
- All 26 have realized local divergence **0.0000**: every one is a read from a
  window where the endogenous copy is byte-identical to the reference.

But the rule as stated above is **not implementable from standard BWA output**.
It routes ties *between categories*, which needs the competitor's identity, and
BWA's XS tag gives the second-best *score* without saying which reference it was
on. There were **zero `XA:Z` tags** in the entire SAM even with `-h 200`. Standard
output says a read is ambiguous; it does not say what it is ambiguous with.

The implementable rule is category-agnostic and already sitting in every BAM:

| policy | calls | false | sensitivity |
|---|---|---|---|
| count everything | 539 | 4.8% | 0.905 |
| **discard MAPQ 0** | 455 | **0.0%** | **0.802** |
| three-bin by category | -- | *not computable* | -- |

So **discard MAPQ 0**: one field, no panel bookkeeping, zero false calls, ten
points of sensitivity. The category framing stays as the *explanation* -- it is
what makes those reads interpretable and it is what the typology measured -- but
the rule a reader should take away is one line. It also beats the ungapped
harness's own version (0.1% false at sensitivity 0.594 in the comparable cell),
because gapped alignment places reads the ungapped one could not.

### Tie policy is an axis, not a setting

Every run reports both, because the difference is large enough to reverse a
conclusion:

- `award_ties_to_exo` -- the optimistic policy, and the one an inflated call
  count comes from.
- `three_bin_ambiguous` -- ties go to an ambiguous bin, which is the three-bin
  rule that the motivating prior-art audit could only label a project-specific
  design choice. The benchmark measures what it buys and what it costs.

Sensitivity is reported against **all simulated exogenous reads**, not against
the ones that happened to be called, because a true read lost to host assignment
is a loss the pipeline never sees.

## Experimental axes

| Axis | Values |
|---|---|
| Median endogenous divergence | swept |
| Spread across loci / within strain / across sites | `--sigma-locus`, `--sigma-strain`, `--rate-shape` |
| Panel composition, endogenous loci **in** the assembly | `viral_only`, `no_decoy` (+host), `full` (+host +decoy) |
| Panel composition, endogenous loci **absent** from the assembly | `poly_no_decoy`, `poly_full` |
| Sequencing depth, read length, fragment length, error rate | configurable |

`--rate-shape 1e6 --rate-jitter 0 --sigma-locus 0 --sigma-strain 0` recovers the
original point-divergence model, kept reachable so the effect of the
distributional model is measured rather than asserted. At a median divergence of
0.10 that ablation is stark: endogenous read windows falling inside the
exogenous range go from **57.8 percent** under the distributional model to
**0.6 percent** under the point model, and per-locus divergence collapses from a
0.024 to 0.154 range onto 0.099 to 0.101. The point model is why the first
version of this harness reported perfect separation at every divergence above
0.10.

## Layout

```
endoexo/simulate.py   panel + read simulation, truth attached to every read
endoexo/align.py      competitive ungapped local alignment, per-read features
endoexo/evaluate.py   production-rule baselines, learned models, CV schemes
run_slice.py          one divergence x panel point, READ level, end to end
run_sample_level.py   the SAMPLE-level detection task over a simulated cohort
run_leakage_probe.py  isolates the read-level-CV inflation and tests its mechanism
endoexo/cohort.py     a COHORT of genomes sharing one host backbone
run_junction_recurrence.py
                      cross-individual junction recurrence: the cohort-level rule
export_for_bwa.py     writes panel.fa, paired FASTQ and a separate truth.tsv
run_bwa.sh            aligns the exported reads with real BWA-MEM (run in WSL)
compare_with_bwa.py   reproduces the tie typology from the SAM, for comparison
tests/                46 tests over the invariants the findings rest on
results/              committed outputs; every README number traces to one
PRIOR_ART.md          adversarial claim matrix. Read it before citing anything.
```

```sh
python run_slice.py        --divergence 0.02 0.05 0.10 0.20 \
                          --panel-mode viral_only no_decoy full poly_no_decoy poly_catalogue \
                          --seed 42 142 242 --out results/results_typology.json
python run_sample_level.py --divergence 0.02 0.10 --panel-mode no_decoy poly_catalogue \
                          --n-samples 60 --out results/results_catalogue_sample.json
python run_junction_recurrence.py --n-samples 50 --clonal-fraction 0.1 0.3 0.6 1.0 \
                          --out results/results_d_junction.json
python run_leakage_probe.py --seeds 42 142 242
python -m unittest discover -s tests -t .
```

The BWA-MEM comparison needs `bwa` and `samtools`, which on this machine means
WSL and one privileged install:

```sh
sudo apt-get update && sudo apt-get install -y bwa samtools   # once
python export_for_bwa.py --divergence 0.10 --panel-mode no_decoy --outdir bwa_compare
bash run_bwa.sh bwa_compare
python compare_with_bwa.py bwa_compare
```

## What it has found so far

Simulation only, with the pure-Python ungapped aligner, over 6 panel
compositions x 4 median divergences x 3 seeds, plus 60-sample cohorts.
Directional, not publishable -- see the limits below.

### The result

**Exogenous/endogenous confusion is a panel-completeness problem, and it has one
remedy.** Mis-assignment comes in two mechanisms (see the table above), and which
one you get is decided entirely by whether the read's true source is in the
panel: 99.8-100% of false calls are exact ties when it is, and 0.0% are when it
is not.

That makes the whole problem one question -- *is the panel complete?* -- with a
single fix on each side of it:

| panel | mechanism | three-bin rule | sample-level detection |
|---|---|---|---|
| contains the locus | tie, 100% | removes it: 62.8% -> **0.0%** false | AUC **0.999**, FPR@95 **0.000** |
| lacks the locus | unopposed, 0% tied | **no-op**: 82.7% -> 82.7% | AUC 0.899, FPR@95 0.467 |
| lacks it, consensus decoy | mixed, 82% tied | partial: 69.3% -> 38.5% | -- |
| lacks it, **locus catalogue** | tie, **99.9%** | removes it: **0.8%** false | AUC **1.000**, FPR@95 **0.000** |

Reading the last two rows against each other is the finding. A family *consensus*
decoy is a bad proxy for one specific old locus and degrades as the locus ages
(82% of false calls converted at divergence 0.02, 59% at 0.20). The actual locus
sequences convert essentially all of them, and the standard rule -- which was
doing nothing at all one row up -- then becomes sufficient.

**And the statistic in production use is the wrong one.** Pipelines report a read
count. The same count restricted to unopposed wins moves sample-level FPR at 95%
sensitivity from 0.633-0.833 to **0.000**. That is the practical recommendation
and it costs one equality test.

### Three things that cut against the obvious reading

**The learned model does not beat the right single column.** Logistic regression
over all 14 sample-level features reaches FPR 0.033; the unopposed-call count
alone reaches 0.000. With 14 features and 60 samples the model is a worse
estimator of something one column already carries.

**The read-level sensitivity ceiling does not propagate to detection.** The tie
fraction among *true* reads is a hard bound on any tie-rejecting rule -- 76.9% of
genuine viral reads are ties at divergence 0.02, leaving **23.1%** retained.
Sample-level detection is nonetheless perfect. Detection needs presence, not
completeness, so the ceiling constrains *quantification* -- load, burden,
clonality -- and not the detection call.

A downward load sweep puts a number on that. The three-bin count holds AUC 1.000
and FPR@95 0.000 down to a **0.11-0.48x** viral coverage band and breaks below
0.1x. But at 0.02-0.1x an infected sample carries 0.6-2.7 viral read *pairs* in
expectation, so the dominant limit there is Poisson sampling of a handful of
molecules, not the filter. And the recommendation is worth **most** at low load:
at 0.11-0.48x the raw call count is near useless (AUC 0.643, FPR@95 0.833) in
the same regime where the filtered count is still perfect.

**A complete catalogue is equivalent to having the loci in the assembly**, to
three digits -- tied fractions 76.9/40.3/19.1% against 76.9/40.3/19.2%, retained
sensitivity 0.226/0.594/0.804 against 0.226/0.594/0.804. The panel does not care
how the read's true source got into it, which is the sharpest confirmation that
presence or absence of the true source is the only thing that matters.

**Coverage breadth fails against the unopposed mechanism, and the mechanism
predicted it.** Breadth reaches AUC 0.953-0.987 where cross-mapping is confined
to conserved windows, and 0.396-0.611 -- chance, at a null SE of 0.075 -- where
it is not. With the true source absent from the panel, reads from *every* window
of the element win unopposed, so the coverage is broad rather than blocky.
Uniformity statistics failed in both forms of the mechanism.

### The one discriminator that works where the panel is incomplete

It is not a per-sample statistic at all. Every rule above is a competition test,
and a competition test is worthless when the read's true source is absent from
the panel. So is every within-sample statistic tried after it, including the
intra-host diversity measure the literature uses for the neighbouring question.

What works lives in the **cohort**. An endogenous insertion sits at the same
host coordinate in everyone who inherited it; an exogenous integration sits at a
coordinate private to the individual. Sort host-side junction coordinates by how
many samples carry them and the two separate -- **without telling any sequence
apart**. It sidesteps the tie / unopposed-win dichotomy rather than solving it.

Measured over 50 samples at four clonal fractions:

- **"Has a junction" is not merely at chance, it is useless.** At a threshold of
  one junction bin it has sensitivity 1.000, FPR 1.000, specificity **0.00** --
  it calls every sample positive, infected or not, because every sample carries
  endogenous junctions. ROC AUC 0.507-0.594.
- **Recurrence filtering gives specificity 0.96, constant in clonal fraction**
  -- exactly one uninfected sample in twenty-five, at every clonal fraction
  tested.
- **Sensitivity is clonal-fraction limited**, not discrimination-limited:
  0.040 / 0.360 / 0.600 / 0.680 at clonal fraction 0.1 / 0.3 / 0.6 / 1.0, and
  the true integration coordinate is recovered in 1 / 9 / 22 / 25 of 25 infected
  samples. Junction-spanning fragments are intrinsically rare -- only those
  straddling a boundary.
- **Validated independently:** endogenous loci recur at their population
  frequency, all ten within 0.08 of truth. The detector recovers each locus's
  frequency from junction evidence alone.

Two honest caveats. FPR at 95% sensitivity is the *wrong* metric here -- the
score is a small integer count, so demanding 95% sensitivity forces a degenerate
threshold; the meaningful operating point is a threshold of one. And the gap
between 25 integrations detected and 22 counted private is coordinate collision
in a 100 kb backbone (200 bins, 25 integrations, birthday rate about 1.5 pairs); a
real 3 Gb genome has ~6 million bins at the same resolution, so the method is
better than these numbers show.

### Methodology

**Read-level CV gives a different answer, not an optimistic one.** It understates
FPR at 95% sensitivity in all sixteen read-level cells, by the most where the
task is hardest (0.646 against a true 0.944), and in the hardest cell reports
2.4x chance where the honest grouped answer is *below* chance. The inflation is
92% attributable to two position features -- the model memorises where each
endogenous locus lands on the viral reference.

## Status and honest limits

This is a **design-validation slice**, not a result set.

- **The aligner is ungapped and written in pure Python.** Scoring follows
  BWA-MEM (match +1, mismatch -4, local, seed length 19), so the seed
  sensitivity limit is the same, but indels between an endogenous copy and the
  exogenous reference are not modelled and `AS` is a lower bound wherever a gap
  would be needed. No number from this harness is publishable until it is
  reproduced with real BWA-MEM on the same simulated FASTQs.
- **The divergence model is substitution-only.** It now has across-loci,
  within-strain and autocorrelated across-site spread, which removed the
  separable-by-alignment-score artefact of the first version, but indels,
  recombination and gene conversion between elements are not modelled at all.
- **Simulation only, so far.** A held-out real-data arm from public HTLV-1
  cell-line WGS is needed before any claim about real cohorts.
- **No novelty claim is made here.** The prior-art check that motivates the
  project is five months old and must be redone before any such wording.

## Provenance

The method lineage -- competitive alignment with host and retroelement decoys,
unique-best filtering, k-mer masking cost, coverage-breadth versus pile-up
discrimination -- comes from viral-sequencing analysis work in a Stanford lab.
**This repository inherits method only.** It carries no cohort, no sample, and
no unpublished number from that work, which is what makes it independently
owned.
