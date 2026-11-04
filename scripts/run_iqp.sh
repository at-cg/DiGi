#!/bin/bash

# Infer a sample's two MHC haplotypes with the recombination-capped IQP
# (src/iqp/main.cpp) from its read k-mers and the graph node k-mers.
# Usage: bash run_iqp.sh <graph_file> <reads_dir> <sample> <out_dir>
#   graph_file comes from build_graph.sh; debug_<graph_file> and
#   node_kmers_<fraction>_k<k>.tsv must sit next to it.
#   reads_dir is the out_dir of process_reads.sh, holding
#   <sample>_read_kmers_<fraction>_k<k>.tsv.
#   FRACTION (default 0.05) and K (default 31) select these k-mer files.
#   C1 (missing k-mer penalty, default 5.2), C2 (novel k-mer penalty, default 3.5) and
#   R (max recombinations per haplotype, default 10) can be set as environment variables.
#   Needs bin/dip_iqp (make iqp) and a Gurobi license (e.g. GRB_LICENSE_FILE or WLS variables).
# -> <out_dir>/<sample>_{1,2}_<weights>_<read_kmers name>.fa, <out_dir>/<sample>_<weights>.log

set -euo pipefail

LABELS="$(realpath "$1")"
READS_DIR="$(realpath "$2")"
SAMPLE="$3"
OUTDIR="$4"
fraction=${FRACTION:-0.05}
k=${K:-31}
c1=${C1:-5.2}
c2=${C2:-3.5}
R=${R:-10}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/dip_iqp"
[ -x "$BIN" ] || { echo "missing $BIN, run make iqp in $ROOT first" >&2; exit 1; }

COORDS="$(dirname "$LABELS")/debug_$(basename "$LABELS")"
GK="$(dirname "$LABELS")/node_kmers_${fraction}_k${k}.tsv"
SK="$READS_DIR/${SAMPLE}_read_kmers_${fraction}_k${k}.tsv"

for f in "$LABELS" "$COORDS" "$GK" "$SK"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done

# H = haplotypes in the graph, r = junctions per haplotype = nodes per haplotype - 1
H=$(grep -c '^Processing ' "$COORDS")
nodes=$(wc -l < "$LABELS")
[ $((nodes % H)) -eq 0 ] || { echo "$nodes nodes not divisible by H=$H" >&2; exit 1; }
r=$((nodes / H - 1))

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

rel_weights="c1-${c1}_c2-${c2}_R-${R}"
cov=$(basename "$SK" .tsv)
echo "Running IQP on $SAMPLE (H=$H, r=$r, $rel_weights)..."

# Run from the output directory so Gurobi's dip_iqp.log lands there
cd "$OUTDIR"
"$BIN" --graph_kmers "$GK" \
  --sk "$SK" \
  --c1 "$c1" --c2 "$c2" \
  --R "$R" \
  --H "$H" --r "$r" \
  --labels "$LABELS" \
  --node_coords "$COORDS" \
  --sample "$SAMPLE" \
  --out_dir "$OUTDIR" \
  --cov "$cov" \
  --rel_weights "$rel_weights" \
  > "$OUTDIR/${SAMPLE}_${rel_weights}.log" 2>&1 \
  || { echo "IQP failed, see $OUTDIR/${SAMPLE}_${rel_weights}.log" >&2; exit 1; }

echo "Inferred haplotypes in $OUTDIR"
