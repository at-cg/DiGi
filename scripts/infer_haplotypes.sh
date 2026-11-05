#!/bin/bash

# Usage: bash infer_haplotypes.sh -r1 <R1.fastq[.gz]> -r2 <R2.fastq[.gz]> -g <graph_dir> -o <out_dir>
#                                 [-s <sample>] [-R <recombination_per_haplotype>]
#   r1, r2: paired-end reads filtered to the MHC region
#   g: graph folder holding panel/, bed/ and node_kmers_0.05_k31.tsv
#   s: sample name used in the output files (default: sample)
#   R: maximum switches between reference haplotypes per inferred haplotype (default 10)
#   FASTP (1 to run fastp on uncleaned reads), THREADS, C1 and C2 can be set as environment
#   variables. C1 and C2 default to cost parameters estimated from the reads' coverage and
#   length (see run_iqp.sh).
#   Needs bin/dip_iqp (make iqp) and a Gurobi license.
# -> <out_dir>/<sample>_{1,2}_*.fa, <out_dir>/reads (read k-mers)

set -euo pipefail

usage() {
    echo "usage: $0 -r1 <R1> -r2 <R2> -g <graph_dir> -o <out_dir> [-s <sample>] [-R <R>]" >&2
    exit 1
}

R1=""
R2=""
GRAPH_DIR=""
OUTDIR=""
SAMPLE=sample
while [ $# -gt 0 ]; do
    [ $# -ge 2 ] || usage
    case "$1" in
        -r1) R1="$2" ;;
        -r2) R2="$2" ;;
        -g) GRAPH_DIR="$2" ;;
        -o) OUTDIR="$2" ;;
        -s) SAMPLE="$2" ;;
        -R) export R="$2" ;;
        *) usage ;;
    esac
    shift 2
done
[ -n "$R1" ] && [ -n "$R2" ] && [ -n "$GRAPH_DIR" ] && [ -n "$OUTDIR" ] || usage

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Fail before read processing if the graph is not there
[ -d "$GRAPH_DIR/panel" ] || { echo "missing $GRAPH_DIR/panel" >&2; exit 1; }
[ -d "$GRAPH_DIR/bed" ] || { echo "missing $GRAPH_DIR/bed" >&2; exit 1; }
[ -s "$GRAPH_DIR/node_kmers_0.05_k31.tsv" ] || { echo "missing $GRAPH_DIR/node_kmers_0.05_k31.tsv" >&2; exit 1; }

echo "Read input files: $R1, $R2"
echo "=== Step 1/2: Read processing ==="
bash "$ROOT/scripts/process_reads.sh" "$R1" "$R2" "$SAMPLE" "$OUTDIR/reads"
echo "=== Step 2/2: Diploid inference ==="
bash "$ROOT/scripts/run_iqp.sh" "$GRAPH_DIR" "$OUTDIR/reads" "$SAMPLE" "$OUTDIR"

echo "=== Done in ${SECONDS}s ==="
