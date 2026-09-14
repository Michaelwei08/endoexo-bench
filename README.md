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
