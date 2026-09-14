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
true source and its true locus, with no human adjudication anywhere.

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
| Endogenous divergence from the ancestral provirus | swept |
| Panel composition | `viral_only`, `no_decoy` (+host), `full` (+host +endogenous decoy) |
| Sequencing depth, read length, fragment length, error rate | configurable |

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
- **Divergence is currently two-point** -- one value for the endogenous family,
  one for the exogenous strain. That makes the high-divergence regime separable
  by alignment score alone, which is an artefact of the simulation and not a
  finding. A divergence distribution across loci is required.
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
