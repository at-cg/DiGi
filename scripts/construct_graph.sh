#!/bin/bash

# Construct the haplotype graph: lift the junctions over to each panel assembly,
# drop junctions that fail in any assembly, build the nodes and sample their k-mers.
# Usage: bash construct_graph.sh <ref.fa> <junctions.bed> <panel_dir> <graph_dir>
#   ref.fa: linear reference the junctions are defined on (e.g. GRCh38)
#   junctions.bed: junctions on ref.fa
#   panel_dir: one <hap>.fa per assembly, holding only the region of interest
#   FRACTION (default 0.05) and K (default 31) can be set as environment variables.
# -> <graph_dir>/ref_graph.txt, debug_ref_graph.txt, node_kmers_<fraction>_k<k>.tsv

set -euo pipefail

[ $# -eq 4 ] || { echo "usage: $0 <ref.fa> <junctions.bed> <panel_dir> <graph_dir>" >&2; exit 1; }

REF="$(realpath "$1")"
JUNCTIONS="$(realpath "$2")"
PANEL_DIR="$(realpath "$3")"
mkdir -p "$4"
GRAPH_DIR="$(realpath "$4")"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

for f in "$REF" "$JUNCTIONS"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done

# The graph builder reads <hap>.fa only, so reject anything else before liftover
shopt -s nullglob
panel=("$PANEL_DIR"/*.fa)
others=("$PANEL_DIR"/*.fasta "$PANEL_DIR"/*.fna "$PANEL_DIR"/*.fas)
shopt -u nullglob
[ ${#panel[@]} -gt 0 ] || { echo "no .fa files in $PANEL_DIR" >&2; exit 1; }
[ ${#others[@]} -eq 0 ] || { echo "rename panel FASTAs to .fa: ${others[*]}" >&2; exit 1; }

# Sorted 4-column BED: the builder reads the lifted position by column number
SORTED_BED="$GRAPH_DIR/junctions.bed"
awk -v OFS='\t' '!/^(#|track|browser)/ && NF >= 3 { print $1, $2, $3, (NF >= 4 ? $4 : "junction_" NR) }' "$JUNCTIONS" \
    | sort -k1,1 -k2,2n > "$SORTED_BED"
[ -s "$SORTED_BED" ] || { echo "no junctions in $JUNCTIONS" >&2; exit 1; }

# Run from graph_dir so Nextflow's work/ and logs land there
LIFTOVER_DIR="$GRAPH_DIR/liftover"
(cd "$GRAPH_DIR" && bash "$ROOT/scripts/liftover.sh" "$REF" "$SORTED_BED" "$PANEL_DIR" "$LIFTOVER_DIR")

python3 "$ROOT/src/graph/filter_junctions.py" "$SORTED_BED" "$LIFTOVER_DIR"

CONFIG="$GRAPH_DIR/graph.yml"
cat > "$CONFIG" <<EOF
experiment:
  name: DiGi
paths:
  base_dir: "$GRAPH_DIR"

graph:
  haplotype_fasta_dir: "$PANEL_DIR"
  rh_bed_file: "{base_dir}/liftover"
  graph_file: "{base_dir}/ref_graph.txt"
EOF

bash "$ROOT/scripts/build_graph.sh" "$CONFIG"
