#!/bin/bash

# Infer a sample's two MHC haplotypes with the recombination-capped IQP
# (src/iqp/main.cpp) from its read k-mers and the graph node k-mers.
# Usage: bash run_iqp.sh <graph_dir> <reads_dir> <sample> <out_dir>
#   graph_dir holds panel/, bed/ and node_kmers_0.05_k31.tsv.
#   reads_dir is the out_dir of process_reads.sh, holding
#   <sample>_read_kmers_0.05_k31.{tsv,stats}.
#   C1 (missing k-mer penalty) and C2 (novel k-mer penalty) default to cost parameters
#   estimated from this sample's coverage and read length, with mu = eps = 0.001. MHC_LEN
#   is the bp the reads were filtered to (default GRCh38 MHC, 4970452).
#   R (max recombinations per haplotype, default 10) can also be set.
#   All of these are environment variables.
#   Needs bin/dip_iqp (make iqp) and a Gurobi license (e.g. GRB_LICENSE_FILE or WLS variables).
# -> <out_dir>/<sample>_{1,2}_<weights>.fa, <out_dir>/<sample>_<weights>.log

set -euo pipefail

GRAPH="$(realpath "$1")"
READS_DIR="$(realpath "$2")"
SAMPLE="$3"
OUTDIR="$4"
fraction=0.05  # fixed by the node k-mers
k=31
R=${R:-10}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/dip_iqp"
[ -x "$BIN" ] || { echo "missing $BIN, run make iqp in $ROOT first" >&2; exit 1; }

GK="$GRAPH/node_kmers_${fraction}_k${k}.tsv"
SK="$READS_DIR/${SAMPLE}_read_kmers_${fraction}_k${k}.tsv"

for f in "$GRAPH" "$GK" "$SK"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done

c1=${C1:-}
c2=${C2:-}
if [ -z "$c1" ] || [ -z "$c2" ]; then
    STATS="${SK%.tsv}.stats"
    [ -s "$STATS" ] || { echo "missing $STATS, rerun process_reads.sh or set C1 and C2" >&2; exit 1; }
    # C = total bases / region length, L = mean read length, lam = per-haplotype k-mer coverage
    costs=$(awk -v k="$k" -v mu=0.001 -v eps=0.001 -v G="${MHC_LEN:-4970452}" '
        NR == 2 {
            if ($1 == 0 || $2 / $1 < k) exit 1
            L = $2 / $1
            C = $2 / G
            lam = C / 2 * (L - k + 1) / L
            hit = 1 - exp(-(1 - eps)^k * lam)
            p = (1 - mu)^k * hit
            q = (1 - mu)^(k - 1) * ((1 - mu) * (1 - exp(-eps / 3 * (1 - eps)^(k - 1) * lam)) + mu / 3 * hit)
            printf "%.4f %.4f %.1f %.1f\n", log(p / q), log((1 - q) / (1 - p)), C, L
        }' "$STATS") || { echo "no reads longer than k=$k in $STATS" >&2; exit 1; }
    read -r est1 est2 C L <<< "$costs"
    echo "Coverage ${C}x, mean read length $L: estimated cost parameters c1=$est1, c2=$est2"
    c1=${c1:-$est1}
    c2=${c2:-$est2}
fi

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

rel_weights="c1-${c1}_c2-${c2}_R-${R}"
echo "Running IQP on $SAMPLE ($rel_weights)..."

# Run from the output directory so Gurobi's dip_iqp.log lands there
cd "$OUTDIR"
"$BIN" --graph_kmers "$GK" \
  --sk "$SK" \
  --c1 "$c1" --c2 "$c2" \
  --R "$R" \
  --graph "$GRAPH" \
  --sample "$SAMPLE" \
  --out_dir "$OUTDIR" \
  --rel_weights "$rel_weights" \
  > "$OUTDIR/${SAMPLE}_${rel_weights}.log" 2>&1 \
  || { echo "IQP failed, see $OUTDIR/${SAMPLE}_${rel_weights}.log" >&2; exit 1; }

echo "Inferred haplotypes in $OUTDIR"
