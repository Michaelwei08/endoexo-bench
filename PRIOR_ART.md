# Prior art check

Date: 2026-09-15. Adversarial: the goal of this pass was to find work that
PRE-EMPTS each claim, not work that supports the project.

This replaces the motivating audit the project was founded on, which was dated
2026-04-15, belonged to a different repository, and scoped the question to one
lab's rubric rather than to the literature.

**Verdict up front, revised after the second pass at the end of this file: every
IDEA in this project has prior art. What does not is the MEASUREMENT. The framing
has to be "here is what these rules actually do, against ground truth, in one
place" -- never "here is a new idea".**

## Claim matrix

| ID | Claim as previously positioned | Verdict | What pre-empts it |
|---|---|---|---|
| P1 | Nobody has asked how to separate exogenous retroviral sequence from endogenous retroelements. | **DEAD** | Hayward et al. 2015 ask exactly this and answer it with a different statistic: intra-host genetic variation. An exogenous infection's sequences coalesce recently, an endogenous element is ancient, so nucleotide diversity, minor allele frequency and proportion of variable sites separate them. They ran simulations over generations, substitution rates, sequence lengths and read depths. Separately, Kaplan et al. 2025 separate the two by epigenetic signature -- CpG depletion and altered trinucleotide frequencies. |
| P2 | Simulating an endogenous locus as ABSENT from the reference assembly was a design gap this project found (recorded as F008). | **DEAD** | The 2022 Frontiers assessment already does it: HERV-K sites were simulated as novel insertions by removing the proviruses from the reference genome, and the simulated FASTQs were aligned to the edited reference with BWA-MEM. That is the `poly_no_decoy` construction, published four years earlier. |
| P3 | Exact alignment-score ties, and the arbitrariness of resolving them, are an unrecognised failure mode. | **DEAD** | MGmapper (2017) states it in as many words: short reads may map equally well to more than one reference with identical alignment scores, and when that happens the assignment is arbitrary because no procedure is available to unravel the ambiguity. |
| P4 | A resolution ceiling computable from reference sequence alone, as a function of read length, is a new idea. | **DEAD** | Umap/Bismap (2018) is exactly this machinery -- single-read mappability, the fraction of a region overlapping at least one uniquely mappable k-mer, parameterised by read length. And ERVmancer (2026) computes the ceiling for retroelements specifically: at a 75 bp read length HERVs resolve only to an average of 12.7 internal clades above the leaf node, and 634 HERVH-LTR7 elements share on average only 2.8% of their reads, which the authors frame as inherent ambiguity that no alignment strategy can remove. Same idea, different reference pair, RNA rather than DNA. |
| P5 | Existing tools address neither case at the alignment level. | **DEAD, and it shipped** | ERVmancer resolves read-mapping ambiguity at the alignment level using phylogenetic placement at the lowest common ancestor. It is HERV-versus-HERV rather than exogenous-versus-endogenous, and RNA rather than DNA, but the sentence as written is false. It was on a resume entry for about an hour before this check caught it. |
| P6 | The error decomposes quantitatively by panel completeness: ~100% exact ties when the read's true source is in the panel, ~0% when it is not. | **PLAUSIBLY SURVIVES** | The tie phenomenon is documented (P3) and irreducible ambiguity is quantified (P4), but I found no source measuring the dichotomy against ground truth, nor stating that which mechanism you get is decided by panel completeness. Needs a targeted second pass before any claim. |
| P7 | The measured cost and benefit of an ambiguous-bin policy, at read and sample level, including that it is a complete no-op when the true source is absent. | **PLAUSIBLY SURVIVES** | The policy is standard and its rationale is stated in the literature. Its measured cost against ground truth -- sample-level FPR at 95% sensitivity 0.467 to 0.000, read retention 0.989 to 0.343 -- I did not find. |
| P8 | Locus-grouped versus read-level cross-validation changes the answer on this task. | **PLAUSIBLY SURVIVES** | Nothing found. Leakage-aware evaluation is of course not novel in general; the specific measurement on this task appears unreported. |
| P9 | A complete locus catalogue in the panel is numerically equivalent to the loci being present in the assembly. | **PLAUSIBLY SURVIVES** | Nothing found. |
| P10 | Young endogenous elements are the unsolvable regime. | **NOT NOVEL -- AND THAT IS THE POINT** | Hayward et al. 2015 state the same limitation from a completely different method: very recent endogenous retroviruses will not have accumulated enough genetic changes to be identified as endogenous. My hardest cell reproduces their stated limitation by a different route, which is independent corroboration and should be reported as agreement rather than as a finding. |

## What this does to the project

**Kept, with the framing changed.** The contribution is not that the question is
unasked. It is a measurement: given the rules that detection pipelines actually
run, what do they buy and what do they cost, against ground truth, and which
property of the setup decides the answer. P6 to P9 are the candidate content and
P10 is a consistency check that strengthens rather than weakens them.

**Killed.** Any wording of the form "first", "novel", "unaddressed", "no
benchmark exists", or "existing tools do not address". P5 in particular was live
on a resume entry and is now corrected.

**Owed before any submission.**

1. **DONE 2026-09-17.** A targeted second pass on P6 to P9 specifically. It moved
   P7 to weakened, P9 to dead, and D to severely weakened.
2. **DONE 2026-09-17.** Read Vy-PER in full. It does not pre-empt the three things
   this project measures; see below.
3. **DONE 2026-09-17.** Read the 2022 Frontiers assessment in full. It confirms P2
   exactly and supplies the direct quantitative comparators.
4. **DONE 2026-09-15 -- see the update at the end of this file.** Decide whether
   the intra-host diversity statistic of Hayward et al. belongs in the benchmark
   as an additional baseline. It was the one item that could change the
   project's conclusion rather than its wording, because it is a genuinely
   different discriminator and it does not need the alignment panel to be
   complete. It was implemented, measured, and does not rescue the hard case.

All four are done. Two new items are owed; see the end of this file.

## Sources

- Hayward, Grabherr, Jern. Characterizing novel endogenous retroviruses from genetic variation inferred from short sequence reads. https://pmc.ncbi.nlm.nih.gov/articles/PMC4616055/
- An assessment of bioinformatics tools for the detection of human endogenous retroviral insertions in short-read genome sequencing data. Frontiers in Bioinformatics, 2022. https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.1062328/full
- A phylogeny-guided framework for decoding mechanisms of human endogenous retrovirus regulation in health and disease (ERVmancer). https://pmc.ncbi.nlm.nih.gov/articles/PMC13228468/
- Karimzadeh et al. Umap and Bismap: quantifying genome and methylome mappability. NAR 2018. https://academic.oup.com/nar/article/46/20/e120/5086676
- MGmapper: Reference based mapping and taxonomy annotation of metagenomics sequence reads. https://pmc.ncbi.nlm.nih.gov/articles/PMC5415185/
- Epigenetic motifs distinguishing endogenous from exogenous retroviral integrants. J Virol 2025. https://journals.asm.org/doi/10.1128/jvi.00775-25
- Vy-PER: eliminating false positive detection of virus integration events in next generation sequencing data. Sci Rep 2015. https://www.nature.com/articles/srep11534 (surfaced, NOT yet read)

## Update 2026-09-15: owed item 4 is done

The intra-host diversity statistic was implemented and measured. It does NOT
rescue the panel-incomplete regime: ROC AUC 0.619 to 0.774 there against a raw
read count at 0.851 to 0.902, and at low load everything including it sits at
chance. Its direction reverses against the literature because Hayward's question
is about the coalescence age of reads that all come from one element, while this
question is about whether a mixture is present. See F050 to F055 in
CONTINUITY.md. The project's conclusion survives the challenge.


## Second pass, 2026-09-17: the targeted one, aimed at P6 to P9 and at D

Owed items 1 to 3 are now done. Vy-PER and the 2022 Frontiers assessment were
read in full, and the P6-to-P9 queries were aimed at the mechanisms rather than
at the general question. The verdicts move AGAINST the project again.

| ID | Verdict after the second pass | What moved it |
|---|---|---|
| P6 | **STILL THE MOST DEFENSIBLE** | Related to the decoy-genome rationale below, but I found no source stating the QUANTITATIVE dichotomy -- that error is essentially all exact ties when the true source is in the panel and essentially none when it is not. |
| P7 | **WEAKENED** | Quantifying what a filter costs against ground truth is standard practice in this literature: the 2024 BMC pathogen-detection paper does it with precision-recall curves and an LOD50 of 10.2 reads. And after B, this project's rule reduces to discarding MAPQ 0, which is universal. What is left is the specific numbers for this specific rule, not the act of measuring. |
| P8 | **UNAFFECTED** | Nothing found. |
| P9 | **DEAD** | This is the decoy-genome rationale, standard in human genomics since the 1000 Genomes hs37d5 decoy: "The decoy contains known true human genome sequence that is not in the reference genome and will suck up reads that would otherwise map with low quality in the reference." My catalogue-equals-assembly result is that idea applied to endogenous retroelements. |
| D | **SEVERELY WEAKENED** | The HIV integration-site literature already uses cross-patient recurrence as an artefact heuristic, and states both halves of my argument: mispriming sites "are sometimes reported to be present at identical locations in samples from different patients", and "the possibility of independent integrations into the identical nucleotide position in the human genome is highly statistically unfavorable". What remains is that their recurrent class is PCR mispriming in linker-mediated PCR assays while mine is endogenous insertions in WGS, and that the measurement -- naive junction rule at specificity 0.00 against recurrence filtering at 0.96, with ground truth -- I did not find. |

### Vy-PER, read in full

Squarely in this territory and attacks it differently. It filters by short tandem
repeat and homopolymer content, re-aligns candidates to the host genome with BLAT,
and requires a minimum number of supporting chimeras. It does **not** discuss MAPQ
filtering, does **not** identify endogenous retroviruses as a false-positive
source, does **not** quantify the sensitivity cost of its filters against
simulated ground truth, and does **not** use cross-individual recurrence. So it
is a different attack on a neighbouring problem, and it does not pre-empt the
three things this project measures.

### The 2022 Frontiers assessment, read in full

Confirms P2 exactly. Non-reference HERV-K insertions were simulated with DWGSIM
after "we removed these proviruses from the hg19 reference using Bedtools masking
followed by deletion" -- 15 LTR3A target loci plus 4 LTR3B specificity controls,
at 150 and 100 bp and 32X / 10.5X / 7X. Tool sensitivity 0.20 to 0.80 and
precision 0.56 to 0.92 at 32X and 150 bp, with a false discovery rate between 8
and 55 percent across tools and datasets. Those are the direct comparators for
anything quantitative here.

It does **not** address exogenous retrovirus contamination, cross-mapping between
exogenous viruses and HERVs, mapping quality or ambiguous read assignment, or
cross-individual recurrence.

## What the project's contribution actually is, after two passes

Every IDEA in this project turns out to have prior art. The tie phenomenon, the
decoy rationale, the resolution ceiling, the recurrence heuristic, the practice of
quantifying a filter's cost -- each is published, mostly in a neighbouring
subfield, and mostly not cited by the others.

What is not published is the MEASUREMENT: these quantities, for this problem,
against ground truth, in one place. That is a legitimate and modest contribution
and it is the only one this project should claim. The writing should be framed as
"here is what these rules actually do", never as "here is a new idea".

Five claims have now been walked back across two passes. That is the check
working rather than failing, but the cumulative lesson is worth stating plainly:
in a problem this well-trodden, the prior is that any mechanism one notices has
been noticed. The numbers are what is new.

### Still owed

- A third pass specifically on P6, which is now the load-bearing claim. The right
  query is about the DICHOTOMY -- whether anyone has reported that mis-assignment
  is all-ties or no-ties depending on reference completeness -- and not about ties
  or decoys separately, both of which are known.
- Read the 1000 Genomes / hs37d5 decoy documentation properly, to cite P9's
  pre-emption precisely rather than from a secondary description.
