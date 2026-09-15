#!/usr/bin/env bash
# Align the exported reads with real BWA-MEM. Run inside WSL, where bwa lives.
#
#   wsl -e bash -lc 'cd /mnt/d/Stanford/research/own/endoexo_bench && bash run_bwa.sh bwa_compare'
#
# One user action is required first, because it needs sudo:
#   sudo apt-get update && sudo apt-get install -y bwa samtools
#
# -a is not optional. It emits every alignment rather than only the primary one,
# and the category tie structure cannot be recovered from primary records alone:
# a read tied between the exogenous reference and the host is reported against
# one of them arbitrarily, with the other appearing only as a secondary record.
#
# -k 19 matches the first-party aligner's seed length, which is also BWA-MEM's
# own default. It is written out so the two are visibly the same choice rather
# than coincidentally the same.
set -euo pipefail

DIR="${1:-bwa_compare}"
THREADS="${THREADS:-4}"

for tool in bwa samtools; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "ERROR: $tool not found." >&2
    echo "  sudo apt-get update && sudo apt-get install -y bwa samtools" >&2
    exit 127
  }
done

for f in panel.fa reads_1.fastq reads_2.fastq truth.tsv reference_map.tsv; do
  [ -f "$DIR/$f" ] || { echo "ERROR: $DIR/$f missing. Run export_for_bwa.py first." >&2; exit 1; }
done

echo "== bwa index =="
bwa index "$DIR/panel.fa" 2>&1 | tail -2

echo "== bwa mem -a -k 19 =="
bwa mem -a -k 19 -t "$THREADS" "$DIR/panel.fa" \
    "$DIR/reads_1.fastq" "$DIR/reads_2.fastq" > "$DIR/aln.sam" 2> "$DIR/bwa.log"
tail -3 "$DIR/bwa.log"

echo "== sorted BAM, for anything that wants one =="
samtools sort -@ "$THREADS" -o "$DIR/aln.bam" "$DIR/aln.sam"
samtools index "$DIR/aln.bam"

echo "== record the versions, because the comparison is only meaningful with them =="
{
  echo "bwa: $(bwa 2>&1 | grep -i '^Version' | head -1)"
  echo "samtools: $(samtools --version | head -1)"
  echo "flags: bwa mem -a -k 19 -t $THREADS"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$DIR/aligner_versions.txt"
cat "$DIR/aligner_versions.txt"

echo
echo "records in SAM: $(grep -vc '^@' "$DIR/aln.sam")"
echo "next: python compare_with_bwa.py $DIR"
