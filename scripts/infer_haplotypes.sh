#!/bin/bash

# Infer a sample's two haplotypes: quality-filter the reads, sample their k-mers
# and solve the IQP against the graph. Read processing is skipped if its reads and
# settings are unchanged, so rerunning with new C1, C2 or R repeats only the IQP.
# Usage: bash infer_haplotypes.sh <R1.fastq[.gz]> <R2.fastq[.gz]> <sample> <graph_dir> <out_dir> [R]
#   R1, R2: paired-end reads filtered to the region of interest
#   graph_dir: output folder of preprocess_graph.sh; its panel must be unchanged
#   R: maximum switches between reference haplotypes per inferred haplotype (default 10)
#   FRACTION, K, THREADS, C1 and C2 can be set as environment variables. Use the
#   same FRACTION and K as preprocess_graph.sh.
#   Needs bin/dip_iqp (make iqp) and a Gurobi license.
# -> <out_dir>/<sample>_{1,2}_*.fa, <out_dir>/reads (filtered reads and read k-mers)

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

echo "=== Step 1/2: Read processing ==="
bash "$ROOT/scripts/process_reads.sh" "$R1" "$R2" "$SAMPLE" "$OUTDIR/reads"
echo "=== Step 2/2: Diploid inference ==="
bash "$ROOT/scripts/run_iqp.sh" "$GRAPH_DIR/ref_graph.txt" "$OUTDIR/reads" "$SAMPLE" "$OUTDIR"

echo "=== Done in ${SECONDS}s ==="
