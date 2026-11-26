#!/bin/bash
# Build a DiGi graph folder from a linear reference and haplotype assemblies.
#
#   1. Junction BED on the reference: one 1-bp junction every <spacing> bp.
#   2. nf-LO (minimap2) builds a ref -> haplotype chain for every <fasta_dir>/*.fa.
#   3. CrossMap lifts the junction BED through each chain.
#   4. junction_filter drops the junctions that do not lift cleanly in EVERY haplotype,
#      then keeps the largest set of the rest that is in reference order in every
#      haplotype, so every haplotype keeps the same junctions (same r, same node count).
#   5. gfa_const writes the graph as GFA and node_kmer_sampler its sampled node
#      k-mers (fraction 0.05, k = 31): the two files the IQP reads.
#
# Usage: liftover.sh <ref.fa> <fasta_dir> <spacing_bp> <outdir>
#   ref.fa      linear reference, first sequence used (e.g. GRCh38_0.fa, the MHC excerpt)
#   fasta_dir   haplotype assemblies, one <hap>.fa per haplotype, first sequence used
#   spacing_bp  distance between junctions on the reference, in bp (any positive integer)
#   outdir      graph folder to create
#
# Output (pass <outdir> to infer_haplotypes.sh -g):
#   <outdir>/graph.gfa                         the graph, with node sequences
#   <outdir>/node_kmers_0.05_k31.tsv           sampled k-mers of each node
#   <outdir>/panel/<hap>.fa                    link to the haplotype FASTA
#   <outdir>/bed/<hap>/<hap>_junctions.bed     lifted junctions, position on <hap> in column 7
#   <outdir>/junctions_ref.bed                 junctions on the reference
#   <outdir>/lift_fails.txt                    one line per failed lift: junction name, haplotype
#   <outdir>/order_drops.txt                   junctions dropped by the order check, one per line
#   <outdir>/liftover/<hap>/                   nf-LO files, raw CrossMap output and log
# A haplotype whose liftover/<hap>/chainnet/liftover.chain already exists is not re-aligned.
#
# NF_LO_VERSION overrides the pinned nf-LO release. Needs awk; Nextflow, its JDK,
# nf-LO's tools and CrossMap come from construct_graph/bin and construct_graph/extra
# (run `make graph` in the DiGi root first).

# Stop at the first failing command (-e), on an unset variable (-u),
# and when any command inside a pipe fails (-o pipefail).
set -euo pipefail

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

# This script is in construct_graph/scripts/, so construct_graph/ is two levels up
SCRIPT_DIR=$(dirname "$(realpath "$0")")
GRAPH_DIR=$(dirname "$SCRIPT_DIR")
TOOLS="$GRAPH_DIR/bin"

# Check that every tool we need is installed
for f in "$TOOLS/CrossMap" "$TOOLS/nextflow" "$GRAPH_DIR/extra/java/bin/java" \
         "$TOOLS/junction_filter" "$TOOLS/gfa_const" "$TOOLS/node_kmer_sampler"; do
    if [ ! -x "$f" ]; then
        echo "Missing $f: run 'make graph' in the DiGi root" >&2
        exit 1
    fi
done

export PATH="$TOOLS:$PATH"
export JAVA_HOME="$GRAPH_DIR/extra/java"
export NXF_HOME="$GRAPH_DIR/extra/nextflow"

# Read the four arguments
if [ $# -ne 4 ]; then
    echo "Usage: $0 <ref.fa> <fasta_dir> <spacing_bp> <outdir>" >&2
    exit 1
fi

REF=$(realpath "$1")
FASTA_DIR=$(realpath "$2")
SPACING=$3
OUTBASE=$(realpath -m "$4")   # -m: the folder does not have to exist yet

# Use NF_LO_VERSION if it is set, otherwise 1.8.8
NF_LO_VERSION=${NF_LO_VERSION:-1.8.8}

# SPACING must be made of digits only, and be more than 0
if ! [[ "$SPACING" =~ ^[0-9]+$ ]] || [ "$SPACING" -eq 0 ]; then
    echo "spacing_bp must be a positive integer: $SPACING" >&2
    exit 1
fi

# List the haplotype FASTA files.
# nullglob: if there is no .fa file, the list is empty (instead of the text "*.fa").
shopt -s nullglob
FASTAS=("$FASTA_DIR"/*.fa)
if [ ${#FASTAS[@]} -eq 0 ]; then
    echo "No .fa files in $FASTA_DIR" >&2
    exit 1
fi

# Name of the first sequence in a FASTA:
# take the first ">" line, keep its first word, remove the ">"
seq_name() {
    grep -m 1 '^>' "$1" | awk '{ print $1 }' | tr -d '>'
}

# Length of the first sequence in a FASTA:
# print its lines (awk stops at the second ">" line), remove line breaks, count characters
seq_length() {
    awk '/^>/ { headers++; next }  headers == 2 { exit }  { print }' "$1" | tr -d '\r\n' | wc -c
}

mkdir -p "$OUTBASE/panel" "$OUTBASE/bed" "$OUTBASE/liftover"
date

# ---------------------------------------------------------------------------
# 1. Junctions on the reference: SPACING, 2*SPACING, ... below the reference length
# ---------------------------------------------------------------------------

REF_NAME=$(seq_name "$REF")
REF_LEN=$(seq_length "$REF")
JUNCTIONS="$OUTBASE/junctions_ref.bed"

# One BED line per junction: <chrom> <pos> <pos+1> grid_<pos>
awk -v chrom="$REF_NAME" -v len="$REF_LEN" -v step="$SPACING" '
    BEGIN {
        for (pos = step; pos < len; pos += step)
            printf "%s\t%d\t%d\tgrid_%d\n", chrom, pos, pos + 1, pos
    }' > "$JUNCTIONS"

N_SRC=$(wc -l < "$JUNCTIONS")
echo "Reference $REF_NAME ($REF_LEN bp), spacing $SPACING bp: $N_SRC junctions"
if [ "$N_SRC" -eq 0 ]; then
    echo "No junctions: spacing is not smaller than the reference" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 2.-3. nf-LO chain and CrossMap for every haplotype
# ---------------------------------------------------------------------------

# Download the nf-LO pipeline once (needs network), then run it offline
nextflow pull evotools/nf-LO -r "$NF_LO_VERSION"
export NXF_OFFLINE=1

HAPS=()   # haplotype names, filled in by the loop below
for fasta in "${FASTAS[@]}"; do
    hap=$(basename "$fasta" .fa)   # e.g. /path/HG002_1.fa -> HG002_1
    HAPS+=("$hap")
    echo "Processing: $hap"

    SAMPLE_OUT="$OUTBASE/liftover/$hap"
    mkdir -p "$SAMPLE_OUT"
    CHAIN="$SAMPLE_OUT/chainnet/liftover.chain"

    # Align only if the chain file is missing or empty (-s: exists and is not empty)
    if [ ! -s "$CHAIN" ]; then
        nextflow run evotools/nf-LO -r "$NF_LO_VERSION" \
            --source "$REF" \
            --target "$fasta" \
            --distance near \
            --aligner minimap2 \
            --no_maf \
            --outdir "$SAMPLE_OUT"
    fi

    # Lift the junctions to this haplotype: results to the .bed, messages to the .log
    CrossMap region "$CHAIN" "$JUNCTIONS" \
        > "$SAMPLE_OUT/${hap}_junctions_raw.bed" \
        2> "$SAMPLE_OUT/crossmap.log"

    # Link the FASTA into panel/ (-f: replace an old link)
    ln -sfn "$fasta" "$OUTBASE/panel/$hap.fa"
    echo "Finished: $hap"
done

# ---------------------------------------------------------------------------
# 4. Keep the junctions that lift cleanly in every haplotype and keep the
#    reference order in all of them (see src/junction_filter.cpp)
# ---------------------------------------------------------------------------

junction_filter "$OUTBASE"

# ---------------------------------------------------------------------------
# 5. GFA and node k-mers (the fraction and k the read k-mers use, see scripts/process_reads.sh)
# ---------------------------------------------------------------------------

gfa_const "$OUTBASE" "$OUTBASE/graph.gfa"
node_kmer_sampler "$OUTBASE" "$OUTBASE/node_kmers_0.05_k31.tsv" 0.05 31
echo "Graph folder: $OUTBASE"
