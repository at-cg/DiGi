#!/bin/bash

# Preprocess the haplotype graph: cut each assembly at its junctions into nodes
# and sample the node k-mers.
# Usage: bash preprocess_graph.sh <panel_dir> <bed_dir> <graph_dir>
#   panel_dir: one <hap>.fa per assembly, holding only the region of interest
#   bed_dir: <hap>/<hap>_junctions.bed per assembly, with the junction position on
#   that assembly in column 7. All BEDs list the same junctions in the same order.
#   FRACTION (default 0.05) and K (default 31) can be set as environment variables.
# -> <graph_dir>/ref_graph.txt, debug_ref_graph.txt, node_kmers_<fraction>_k<k>.tsv

set -euo pipefail

[ $# -eq 3 ] || { echo "usage: $0 <panel_dir> <bed_dir> <graph_dir>" >&2; exit 1; }

PANEL_DIR="$(realpath "$1")"
BED_DIR="$(realpath "$2")"
mkdir -p "$3"
GRAPH_DIR="$(realpath "$3")"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# The graph builder reads <hap>.fa only
shopt -s nullglob
panel=("$PANEL_DIR"/*.fa)
others=("$PANEL_DIR"/*.fasta "$PANEL_DIR"/*.fna "$PANEL_DIR"/*.fas)
shopt -u nullglob
[ ${#panel[@]} -gt 0 ] || { echo "no .fa files in $PANEL_DIR" >&2; exit 1; }
[ ${#others[@]} -eq 0 ] || { echo "rename panel FASTAs to .fa: ${others[*]}" >&2; exit 1; }

# The graph builder skips lines without column 7, which would lose junctions silently
for fa in "${panel[@]}"; do
    hap=$(basename "$fa" .fa)
    bed="$BED_DIR/$hap/${hap}_junctions.bed"
    [ -s "$bed" ] || { echo "missing $bed" >&2; exit 1; }
    awk -F'\t' 'NF < 7 { exit 1 }' "$bed" || { echo "fewer than 7 columns in $bed" >&2; exit 1; }
done

CONFIG="$GRAPH_DIR/graph.yml"
cat > "$CONFIG" <<EOF
experiment:
  name: DiGi
paths:
  base_dir: "$GRAPH_DIR"

graph:
  haplotype_fasta_dir: "$PANEL_DIR"
  rh_bed_file: "$BED_DIR"
  graph_file: "{base_dir}/ref_graph.txt"
EOF

bash "$ROOT/scripts/build_graph.sh" "$CONFIG"
