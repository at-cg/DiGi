#!/bin/bash

# Count the FracMinHash-sampled k-mers of MHC-filtered paired-end reads.
# Usage: bash process_reads.sh <R1.fastq[.gz]> <R2.fastq[.gz]> <sample> <out_dir>
#   THREADS (default 1), KMC_MEM in GB (default 16) and KMC_TMP (dir for KMC's bins/db,
#   default <out_dir>) can be set as environment variables.
# -> <out_dir>/<sample>_read_kmers_0.05_k31.{tsv,stats,params}
#    .stats holds the read and base totals of R1 + R2; run_iqp.sh derives coverage,
#    read length and the costs C1 and C2 from it.

set -euo pipefail

R1="$1"
R2="$2"
SAMPLE="$3"
OUTDIR="$4"
fraction=0.05  # fixed by the node k-mers
k=31
threads=${THREADS:-1}
kmc_mem=${KMC_MEM:-16}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/kmer_set_sampling"
[ -x "$BIN" ] || { echo "missing $BIN, run make in $ROOT first" >&2; exit 1; }
for tool in kmc kmc_tools; do
  command -v "$tool" >/dev/null || { echo "missing $tool, activate the digi env" >&2; exit 1; }
done
for f in "$R1" "$R2"; do
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

# Skip if this exact run already finished. The stamp is removed first and written
# last, so an interrupted run never counts as done.
out="$OUTDIR/${SAMPLE}_read_kmers_${fraction}_k${k}.tsv"
stamp="${out%.tsv}.params"
stats="${out%.tsv}.stats"
params=$(
  for f in "$R1" "$R2"; do realpath "$f"; stat -L -c '%s %Y' "$f"; done
  echo "fraction=$fraction k=$k"
  cat "$0" "$BIN" | md5sum
)
if [ -s "$out" ] && [ -s "$stats" ] && [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$params" ]; then
  echo "Reads and read settings unchanged, skipping read processing"
  echo "Read k-mers: $out"
  exit 0
fi
rm -f "$stamp"

# Read and base totals (gzip -f passes plain FASTQ through)
gzip -cdf "$R1" "$R2" \
  | awk 'BEGIN { print "reads\tbases" } NR % 4 == 2 { n++; b += length($0) } END { printf "%.0f\t%.0f\n", n, b }' \
  > "$stats"

# Count canonical k-mers of R1 + R2 with KMC, then sample them. KMC drops k-mers
# containing N and canonicalises as the lexicographic min, like node_kmers_gen.
# Sampling depends on the k-mer alone, so counting everything first gives the
# same counts as sampling while counting.
work=$(mktemp -d "${KMC_TMP:-$OUTDIR}/kmc.${SAMPLE}.XXXXXX")
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/bins"
readlink -f "$R1" "$R2" > "$work/reads.lst"

echo "[1/2] Counting k-mers with KMC (k=$k)..."
# -ci1: keep singletons (KMC default -ci2 drops them); -cs: no 255 cap on counts
kmc -k"$k" -ci1 -cs4294967295 -t"$threads" -m"$kmc_mem" -fq @"$work/reads.lst" "$work/db" "$work/bins" \
  > "$work/kmc.log" 2>&1 || { cat "$work/kmc.log" >&2; echo "kmc failed for $SAMPLE" >&2; exit 1; }
grep -E "No. of unique counted k-mers|Total no. of k-mers|Total no. of reads" "$work/kmc.log"

echo "[2/2] Sampling read k-mers (fraction $fraction)..."
# Stream the dump straight into the sampler instead of writing it out
kmc_tools -t"$threads" transform "$work/db" dump /dev/stdout 2> "$work/dump.log" \
  | "$BIN" "$fraction" /dev/stdin "$out.tmp" 2>&1 \
  && [ -s "$out.tmp" ] && mv "$out.tmp" "$out" \
  || { cat "$work/dump.log" >&2; echo "read k-mer sampling failed for $SAMPLE" >&2; rm -f "$out.tmp"; exit 1; }

printf '%s\n' "$params" > "$stamp"
echo "Read k-mers: $out"
