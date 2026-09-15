#!/usr/bin/env python3
"""Write the simulated panel and reads out as real FASTA and FASTQ, so the same
reads can be pushed through BWA-MEM and through this repository's own aligner.

WHY THIS EXISTS
  Every result here rests on a first-party UNGAPPED aligner, and the sharpest
  claim in the project rests on exact score ties: 99.8-100 percent of false
  calls are `AS == XS` when the read's true source is in the panel, which makes
  them free to detect. Gapped local alignment might not preserve that. If BWA
  turns an exact tie into a near-tie -- scores differing by one or two because a
  gap bought a base somewhere -- then "free to detect" becomes "needs a
  threshold", and the binary decomposition softens into a tuning problem.

  That is a falsifiable prediction and this file is how it gets tested. The
  comparison is the SAME READS through both aligners, so any difference is the
  aligner and not the simulation.

WHAT IT WRITES
  panel.fa        every panel reference, one record each
  reads_1.fastq   first end
  reads_2.fastq   second end, REVERSE COMPLEMENTED
  truth.tsv       read name -> true source, locus, realized window divergence

  Read names are opaque (`r000001`). The truth lives in a separate file on
  purpose: a name that carries its own label is a label leak waiting to happen
  the first time someone parses it.

TWO CONVENTIONS THAT MATTER FOR BWA AND NOT FOR OUR ALIGNER
  1 The second end is reverse complemented. Our aligner tries both strands per
    read, so orientation was irrelevant to it; BWA's pairing, proper-pair flag
    and insert-size estimate all depend on the forward-reverse convention.
  2 Base qualities are a constant Q40. BWA-MEM does not use base quality in its
    match/mismatch scoring, and the simulator has already injected errors into
    the bases themselves, so a quality track would be decoration. Stated here so
    nobody later reads the constant as an accident.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from endoexo.align import revcomp
from endoexo.simulate import simulate_reads, simulate_world
from run_slice import PANEL_MODES, restrict_panel

BASES = "ACGT"
CONST_QUAL = "I"          # Phred+33 'I' = Q40


def to_string(seq: np.ndarray) -> str:
    return "".join(BASES[b] for b in seq.tolist())


def write_fasta(path: Path, refs: dict[str, np.ndarray], width: int = 70) -> None:
    with path.open("w", encoding="ascii", newline="\n") as fh:
        for ref_id, seq in refs.items():
            fh.write(f">{ref_id}\n")
            s = to_string(seq)
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def write_fastq_pair(path1: Path, path2: Path, reads, read_len: int) -> list[tuple]:
    """Write both ends and return the truth rows in the same order."""
    truth = []
    with path1.open("w", encoding="ascii", newline="\n") as f1, \
         path2.open("w", encoding="ascii", newline="\n") as f2:
        for i, t in enumerate(reads):
            name = f"r{i:06d}"
            q = CONST_QUAL * read_len
            f1.write(f"@{name}/1\n{to_string(t.read)}\n+\n{q}\n")
            # second end reverse complemented, per the forward-reverse convention
            f2.write(f"@{name}/2\n{to_string(revcomp(t.mate))}\n+\n{q}\n")
            truth.append((name, t.source, t.locus, t.frag_start, t.local_div))
    return truth


def write_truth(path: Path, truth: list[tuple], meta: dict, panel_mode: str,
                divergence: float, seed: int) -> None:
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write(f"# panel_mode={panel_mode} divergence={divergence} seed={seed}\n")
        fh.write(f"# realized_strain_divergence={meta['realized_strain_divergence']:.6f}\n")
        fh.write(f"# locus_divergence_min={meta['locus_divergence_min']:.6f} "
                 f"median={meta['locus_divergence_median']:.6f} "
                 f"max={meta['locus_divergence_max']:.6f}\n")
        fh.write("# local_div is a DIAGNOSTIC, never a feature\n")
        fh.write("read_name\ttrue_source\ttrue_locus\tfrag_start\tlocal_divergence\n")
        for name, source, locus, start, ldiv in truth:
            ld = "" if ldiv != ldiv else f"{ldiv:.6f}"      # NaN -> empty
            fh.write(f"{name}\t{source}\t{locus}\t{start}\t{ld}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", type=float, default=0.10)
    ap.add_argument("--panel-mode", default="no_decoy", choices=list(PANEL_MODES))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--exo-depth", type=float, default=20.0)
    ap.add_argument("--host-depth", type=float, default=8.0)
    ap.add_argument("--n-herv-loci", type=int, default=12)
    ap.add_argument("--host-filler-bp", type=int, default=200_000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--outdir", default="bwa_compare")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    panel, host, herv_spans, exo_strain, meta = simulate_world(
        args.divergence, n_herv_loci=args.n_herv_loci,
        host_filler_bp=args.host_filler_bp, seed=args.seed)
    restrict_panel(panel, args.panel_mode)

    reads = simulate_reads(exo_strain, host, herv_spans, meta,
                           exo_depth=args.exo_depth, host_depth=args.host_depth,
                           read_len=args.read_len, seed=args.seed + 1)

    write_fasta(out / "panel.fa", panel.refs)
    truth = write_fastq_pair(out / "reads_1.fastq", out / "reads_2.fastq",
                             reads, args.read_len)
    write_truth(out / "truth.tsv", truth, meta, args.panel_mode,
                args.divergence, args.seed)

    # The category of each reference, which BWA has no idea about and the
    # comparison needs in order to reproduce the tie typology.
    with (out / "reference_map.tsv").open("w", encoding="ascii", newline="\n") as fh:
        fh.write("reference_id\tcategory\tlength\n")
        for ref_id, seq in panel.refs.items():
            fh.write(f"{ref_id}\t{panel.categories[ref_id]}\t{len(seq)}\n")

    print(f"wrote {out}/panel.fa ({len(panel.refs)} references, "
          f"{sum(len(s) for s in panel.refs.values()):,} bp)")
    print(f"wrote {out}/reads_1.fastq and reads_2.fastq ({len(reads):,} pairs)")
    print(f"wrote {out}/truth.tsv and reference_map.tsv")
    print(f"realized locus divergence "
          f"{meta['locus_divergence_min']:.3f}/{meta['locus_divergence_median']:.3f}/"
          f"{meta['locus_divergence_max']:.3f}, strain "
          f"{meta['realized_strain_divergence']:.4f}")


if __name__ == "__main__":
    main()
