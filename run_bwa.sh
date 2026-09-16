#!/usr/bin/env bash
# Align the exported reads with real BWA-MEM. Run inside WSL, where bwa lives.
#
#   wsl -e bash -lc 'cd /mnt/d/Stanford/research/own/endoexo_bench && bash run_bwa.sh bwa_compare'
#
# NO SUDO IS NEEDED. bwa builds from source against the zlib headers Ubuntu
# already ships, so the apt route is unnecessary:
#   mkdir -p ~/tools && cd ~/tools && git clone --depth 1 https://github.com/lh3/bwa.git
#   cd bwa && make                       # then BWA=~/tools/bwa/bwa
#
# samtools is OPTIONAL. compare_with_bwa.py parses the SAM directly; samtools is
# used here only to produce a sorted BAM for anything that wants one, and the
# script skips that step rather than failing when samtools is absent. Requiring
# it was the only reason this step looked like it needed a privileged install.
#
# -a is kept but is NOT sufficient: it does not emit secondary alignments for
# properly paired reads, which the first real run demonstrated -- one record per
# read end, zero secondary flags. -h 200 is what actually matters. It emits XA:Z
# listing alternative hits BY REFERENCE NAME, which is the only place default BWA
# output names the competing reference, and the category of the competitor is
# what the tie typology is defined on.
#
# -k 19 matches the first-party aligner's seed length, which is also BWA-MEM's
# own default. It is written out so the two are visibly the same choice rather
# than coincidentally the same.
set -euo pipefail

DIR="${1:-bwa_compare}"
THREADS="${THREADS:-4}"

BWA="${BWA:-$(command -v bwa || echo "$HOME/tools/bwa/bwa")}"
[ -x "$BWA" ] || {
  echo "ERROR: no bwa at '$BWA'." >&2
  echo "  build it without sudo:" >&2
  echo "  mkdir -p ~/tools && cd ~/tools && git clone --depth 1 https://github.com/lh3/bwa.git && cd bwa && make" >&2
  exit 127
}
SAMTOOLS="${SAMTOOLS:-$(command -v samtools || true)}"

for f in panel.fa reads_1.fastq reads_2.fastq truth.tsv reference_map.tsv; do
  [ -f "$DIR/$f" ] || { echo "ERROR: $DIR/$f missing. Run export_for_bwa.py first." >&2; exit 1; }
done

echo "== bwa index =="
"$BWA" index "$DIR/panel.fa" 2>&1 | tail -2

echo "== bwa mem -a -h 200 -k 19 =="
"$BWA" mem -a -h 200 -k 19 -t "$THREADS" "$DIR/panel.fa" "$DIR/reads_1.fastq" "$DIR/reads_2.fastq" > "$DIR/aln.sam" 2> "$DIR/bwa.log"
tail -3 "$DIR/bwa.log"

if [ -n "$SAMTOOLS" ]; then
  echo "== sorted BAM, for anything that wants one =="
  "$SAMTOOLS" sort -@ "$THREADS" -o "$DIR/aln.bam" "$DIR/aln.sam"
  "$SAMTOOLS" index "$DIR/aln.bam"
else
  echo "== samtools absent, skipping the BAM. The comparison does not need it. =="
fi

echo "== record the versions, because the comparison is only meaningful with them =="
{
  echo "bwa: $("$BWA" 2>&1 | grep -i '^Version' | head -1)"
  echo "bwa path: $BWA"
  if [ -n "$SAMTOOLS" ]; then echo "samtools: $("$SAMTOOLS" --version | head -1)"
  else echo "samtools: absent (not required)"; fi
  echo "flags: bwa mem -a -h 200 -k 19 -t $THREADS"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$DIR/aligner_versions.txt"
cat "$DIR/aligner_versions.txt"

echo
echo "records in SAM: $(grep -vc '^@' "$DIR/aln.sam")"
echo "next: python compare_with_bwa.py $DIR"
