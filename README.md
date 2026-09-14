# endoexo-bench

A ground-truth benchmark for a question that viral-detection pipelines answer by
convention rather than by measurement: **when a read aligns to an exogenous
retrovirus reference, did it come from the virus, or from an endogenous
retroelement in the host genome?**

Existing HERV tooling splits into DNA insertion callers (ERVcaller, MELT,
RetroSeq, STEAK, xTea) and RNA expression quantifiers (Telescope, ERVmap).
Neither addresses the exogenous/endogenous confusion at the alignment level.
The rules used in practice to resolve it -- competitive alignment against a
host reference, a mapping-quality floor, a unique-best `AS > XS` test, an
aligned-length floor -- are reasonable, widely used, and as far as this
project's prior-art check has established, **never benchmarked against ground
truth**.

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
run_slice.py          one divergence x panel point, end to end
run_leakage_probe.py  isolates the read-level-CV inflation and tests its mechanism
```

```sh
python run_slice.py --divergence 0.02 0.05 0.10 0.20 \
                    --panel-mode viral_only no_decoy full --out results.json
python run_leakage_probe.py --seeds 42 142 242
```

## What it has found so far

Simulation only, with the pure-Python ungapped aligner, over 5 panel compositions
x 4 median divergences x 3 seeds. Directional, not publishable -- see the limits
below.

**1. The three-bin rule is a no-op exactly where the field needs it.** Routing a
tie to an ambiguous bin is the standard answer to exogenous/endogenous
ambiguity. When the endogenous locus sits in the reference assembly it works
almost perfectly on specificity and is expensive:

| endogenous locus | false calls before -> after | sensitivity before -> after |
|---|---|---|
| in the assembly, divergence 0.02 | 62.8% -> **0.0%** | 0.989 -> **0.343** |
| in the assembly, divergence 0.20 | 31.1% -> **0.0%** | 0.999 -> 0.780 |
| absent, consensus decoy, div 0.02 | 69.3% -> 38.5% | 1.000 -> 0.617 |
| **absent, no decoy** | **82.7% -> 82.7%** | **1.000 -> 1.000** |

The last row is the finding. There are no ties to bin, so the rule changes
nothing at any divergence.

**2. The sensitivity ceiling is computable in advance, from references alone.**
The fraction of *true* exogenous reads that are themselves ties is a hard bound
on any tie-rejecting rule: 65.3% of genuine viral reads are ties at divergence
0.02 with the host in the panel, so no such rule can exceed **34.7%**
sensitivity. No cohort is needed to know this before running an assay.

**3. The mate rescues the tie mechanism, and only that one.** Ablating all seven
mate and pair features from the grouped model degrades FPR at 95% sensitivity
from 0.702 to 0.967 (divergence 0.02) and 0.189 to 0.887 (0.20) where a
competitor is present -- and changes nothing (0.944 vs 0.958) where none is. The
mate helps because it lands a fragment away, possibly outside the conserved
window, where the competition is decidable; with nothing to compete against it
carries no more information than the read.

**4. A consensus decoy degrades as the locus ages.** It converts unopposed wins
into detectable ties for 82.1% of false calls at divergence 0.02 but only 59.0%
at 0.20. The right decoy for polymorphic insertions is locus-resolved, not a
family consensus.

**5. The open problem, as a negative result.** An unopposed win against a young
endogenous element is not solvable at the read-pair level. The grouped model
reaches PR AUC 0.138 against a 0.173 prevalence floor -- **below chance**, over
three seeds. The three-bin rule is a no-op there and the mate is worthless
there. The remaining routes are above the read pair: coverage breadth over
non-conserved windows, and host-virus junction evidence.

**6. Read-level CV gives a different answer, not an optimistic one.** It
understates FPR at 95% sensitivity in all sixteen cells and by the most where
the task is hardest (0.646 against a true 0.944), and in the hardest cell it
reports 2.4x chance where the honest grouped answer is below chance.

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
