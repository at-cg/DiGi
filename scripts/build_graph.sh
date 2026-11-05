#!/bin/bash

# Build the junction-linked haplotype graph (node boundaries on each assembly) and its
# FracMinHash-sampled node k-mers.
# Usage: bash build_graph.sh <graph.yml> [samples_included.txt]
#   FRACTION (default 0.05) and K (default 31) can be set as environment variables.
# -> <graph_file>, node_kmers_<fraction>_k<k>.tsv (next to graph_file)

set -euo pipefail

CONFIG="$1"
SAMPLES="${2:-}"
fraction=${FRACTION:-0.05}
k=${K:-31}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/node_kmers_gen"
[ -x "$BIN" ] || { echo "missing $BIN, run make in $ROOT first" >&2; exit 1; }

# Graph output path, read from the config
GRAPH_FILE=$(python3 -c "import sys, yaml; c = yaml.safe_load(open(sys.argv[1])); print(c['graph']['graph_file'].format(base_dir=c['paths']['base_dir']))" "$CONFIG")
OUTDIR=$(dirname "$GRAPH_FILE")
KMERS="$OUTDIR/node_kmers_${fraction}_k${k}.tsv"
mkdir -p "$OUTDIR"

# Build graph nodes
echo "[1/2] Building graph nodes..."
python3 "$ROOT/src/graph/graph_builder.py" --config "$CONFIG" ${SAMPLES:+--samples_included "$SAMPLES"}

# Sample node k-mers
echo "[2/2] Sampling node k-mers (fraction $fraction, k=$k)..."
"$BIN" --fraction "$fraction" --k "$k" --input "$GRAPH_FILE" --output "$KMERS.tmp" \
  && [ -s "$KMERS.tmp" ] && mv "$KMERS.tmp" "$KMERS" \
  || { echo "node k-mer sampling failed" >&2; exit 1; }

echo "Graph:      $GRAPH_FILE"
echo "Node k-mers: $KMERS"
