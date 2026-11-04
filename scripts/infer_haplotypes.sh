#!/bin/bash

# Infer a sample's two haplotypes: quality-filter the reads, sample their k-mers
# and solve the IQP against the graph.
# Usage: bash infer_haplotypes.sh <R1.fastq[.gz]> <R2.fastq[.gz]> <sample> <graph_dir> <out_dir> [R]
#   R1, R2: paired-end reads filtered to the region of interest
#   graph_dir: output folder of preprocess_graph.sh
#   R: maximum switches between reference haplotypes per inferred haplotype (default 10)
#   FRACTION, K, THREADS, C1 and C2 can be set as environment variables. Use the
#   same FRACTION and K as preprocess_graph.sh.
#   Needs bin/dip_iqp (make iqp) and a Gurobi license.
# -> <out_dir>/<sample>_{1,2}_*.fa, <out_dir>/<sample>_{1,2}_*_haplotype_regions.tsv,
#    <out_dir>/reads (filtered reads and read k-mers)

set -euo pipefail

[ $# -eq 5 ] || [ $# -eq 6 ] \
    || { echo "usage: $0 <R1> <R2> <sample> <graph_dir> <out_dir> [R]" >&2; exit 1; }

R1="$1"
R2="$2"
SAMPLE="$3"
GRAPH_DIR="$4"
OUTDIR="$5"
export R="${6:-${R:-10}}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Fail before read processing if the graph is not there
[ -s "$GRAPH_DIR/ref_graph.txt" ] || { echo "missing $GRAPH_DIR/ref_graph.txt" >&2; exit 1; }

bash "$ROOT/scripts/process_reads.sh" "$R1" "$R2" "$SAMPLE" "$OUTDIR/reads"
bash "$ROOT/scripts/run_iqp.sh" "$GRAPH_DIR/ref_graph.txt" "$OUTDIR/reads" "$SAMPLE" "$OUTDIR"
