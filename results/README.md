# Committed results

Every quantitative claim in the top-level `README.md` traces to a file here, so a
reader can check it without re-running anything. All of it is simulation output:
no patient data, no controlled-access data, nothing that needs a use agreement.

Each file is regenerable from its recorded seeds. Runtimes are minutes.

## READ THIS FIRST: two incompatible generations of file

A defect in the site-rate model was found and fixed on 2026-09-15 (`F028` in
`../CONTINUITY.md`). Gamma site rates have mean 1 in *expectation* but not in
sample: with only 22 blocks across an 8.5 kb provirus the standard error of their
mean is about 30%, so a requested divergence of 0.10 was being delivered as
anything from 0.035 to 0.19. Rates are now normalised to mean exactly 1.

**Read-level tie fractions and sensitivity figures must be quoted from the
post-fix files.** The pre-fix files remain because the findings they back are
comparisons *within* a run, which the drift does not invalidate, and because
removing them would remove the audit trail.

| | files |
|---|---|
| **post-fix** (quote these) | `results_catalogue_read`, `results_catalogue_sample`, `results_nodecoy_postfix`, `results_load_*`, `results_f049_*`, `results_d_*`, `results_bwa_comparison`, `results_divergence_ablation` |
| **pre-fix** (directional only) | `results_f009_dist`, `results_typology`, `results_sample_level` |
| **point-divergence model** (superseded, kept for the trail) | `results_sweep1`, `results_poly`, `results_leakage_probe` |

## What backs what

| file | findings | what it is |
|---|---|---|
| `results_sweep1.json` | F001, F002 | first divergence x panel sweep. Point-divergence model, so every cell at divergence >= 0.10 is uninformative (F007). |
| `results_poly.json` | F009, F010 | the polymorphic panel modes, still on the point model. F009's zeros were artefacts of it. |
| `results_leakage_probe.json` | F003-F005 | three-seed read-level-CV leakage probe and the position-feature ablation that identifies its mechanism. |
| `results_f009_dist.json` | F013-F018 | F009 re-run under the distributional divergence model, 5 panels x 4 divergences x 3 seeds. |
| `results_typology.json` | F023-F027 | the tie typology and the mate-feature ablation, same grid. |
| `results_sample_level.json` | F030-F034 | the sample-level detection task, 60-sample cohorts. Includes the failed breadth hypothesis (F032). |
| `results_catalogue_read.json` | F035 | complete locus catalogue, read level, 3 divergences x 3 seeds. |
| `results_catalogue_sample.json` | F035 | the same at sample level. |
| `results_nodecoy_postfix.json` | F039, F040 | post-fix `no_decoy`, which is what shows a complete catalogue is *equivalent* to having the loci in the assembly, and what restates the sensitivity ceiling. |
| `results_load_0.02_0.10.json` | F041, F042 | viral-load sweep, lowest band. This is where the three-bin count breaks. |
| `results_load_0.10_0.50.json` | F041, F042 | viral-load sweep, 0.11-0.48x band. This is where it still holds perfectly. |
| `results_f049_diversity.json` | F050-F055 | the intra-host diversity statistic as a standalone score. It loses to a raw read count where the panel is incomplete. |
| `results_f049_lowload.json` | F054 | the same at low viral load, where everything is at chance. |
| `results_f049_ablation.json` | F056, F057 | diversity features ablated from the sample-level model, for a clean attribution. |
| `results_d_junction.json` | F058-F063 | cross-individual junction recurrence over 50 samples at four clonal fractions. |
| `results_d_decomposition.json` | F075 | why a detected integration is not always private. It is integration-integration collision, 9 of 9. |
| `results_bwa_comparison.json` | F066-F073 | REAL BWA-MEM 0.7.19 on the exported reads. The only file here not produced by the first-party aligner. |
| `results_divergence_ablation.json` | F012, and section 2 of PAPER.md | what the distributional divergence model buys, measured. Existed as an ad-hoc check first; promoted because the paper quotes it. |

## Regenerating

```sh
python run_slice.py        --divergence 0.02 0.05 0.10 0.20 \
                           --panel-mode viral_only no_decoy full poly_no_decoy poly_full \
                           --seed 42 142 242 --host-filler-bp 200000 --n-herv-loci 12 \
                           --exo-depth 20 --host-depth 8 --out results/results_typology.json
python run_sample_level.py --divergence 0.02 0.10 --panel-mode no_decoy poly_catalogue \
                           --n-samples 60 --n-herv-loci 6 --host-filler-bp 100000 \
                           --cohort-seed 42 --out results/results_catalogue_sample.json
python run_leakage_probe.py --seeds 42 142 242
```

Re-running the pre-fix files against current code will NOT reproduce their
numbers, by design -- the rate normalisation changed them. That is what the two
generations above are about.
