#!/bin/bash

# Build the junction-linked haplotype graph (node boundaries on each assembly) and its
# FracMinHash-sampled node k-mers.
# Usage: bash build_graph.sh <graph.yml> [samples_included.txt]
# -> <graph_file>, node_kmers_<fraction>_k<k>.{tsv,params} (next to graph_file)

set -euo pipefail

CONFIG="$1"
SAMPLES="${2:-}"
fraction=${FRACTION:-0.05}
k=${K:-31}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/node_kmers_gen"
[ -x "$BIN" ] || { echo "missing $BIN, run make in $ROOT first" >&2; exit 1; }

# Paths from the config
mapfile -t paths < <(python3 -c "
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
fmt = {'base_dir': c['paths']['base_dir']}
for key in ('graph_file', 'haplotype_fasta_dir', 'rh_bed_file'):
    print(c['graph'][key].format(**fmt))
" "$CONFIG")
GRAPH_FILE="${paths[0]}"
PANEL_DIR="${paths[1]}"
BED_DIR="${paths[2]}"
OUTDIR=$(dirname "$GRAPH_FILE")
KMERS="$OUTDIR/node_kmers_${fraction}_k${k}.tsv"
mkdir -p "$OUTDIR"

# Skip if this exact graph is already built. The stamp is removed first and written
# last, so an interrupted run never counts as done.
stamp="${KMERS%.tsv}.params"
params=$(
  find -L "$PANEL_DIR" "$BED_DIR" -type f \( -name '*.fa' -o -name '*_junctions.bed' \) \
    -printf '%p %s %T@\n' | sort
  echo "fraction=$fraction k=$k"
  cat "$0" "$ROOT/src/graph/graph_builder.py" "$BIN" ${SAMPLES:+"$SAMPLES"} | md5sum
)
if [ -s "$GRAPH_FILE" ] && [ -s "$KMERS" ] && [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$params" ]; then
  echo "Graph inputs and settings unchanged, skipping graph building"
  echo "Graph:      $GRAPH_FILE"
  echo "Node k-mers: $KMERS"
  exit 0
fi
rm -f "$stamp"

# Build graph nodes
echo "[1/2] Building graph nodes..."
python3 "$ROOT/src/graph/graph_builder.py" --config "$CONFIG" ${SAMPLES:+--samples_included "$SAMPLES"}

# Sample node k-mers
echo "[2/2] Sampling node k-mers (fraction $fraction, k=$k)..."
"$BIN" --fraction "$fraction" --k "$k" --input "$GRAPH_FILE" --output "$KMERS.tmp" \
  && [ -s "$KMERS.tmp" ] && mv "$KMERS.tmp" "$KMERS" \
  || { echo "node k-mer sampling failed" >&2; exit 1; }

printf '%s\n' "$params" > "$stamp"
echo "Graph:      $GRAPH_FILE"
echo "Node k-mers: $KMERS"
