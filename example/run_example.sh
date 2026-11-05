#!/bin/bash

# Infer HG002's two MHC haplotypes from 10x MHC-filtered Illumina reads with an
# 18-haplotype HGSVC graph cut every 10 kb.
# Usage: [R] bash example/run_example.sh

set -euo pipefail

EX="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$EX")"
cd "$EX"

export FRACTION=${FRACTION:-0.05}
export K=${K:-31}
export C1=${C1:-6.5127} C2=${C2:-4.0225} R=${R:-10}

mkdir -p work
[ -d work/panel ] || tar -C work -xzf data/HGSVC_MHC_panel_18hap_grid10kb.tar.xz

# 1. Graph and node k-mers
echo "=== Step 1/3: Graph and node k-mers ==="
bash "$ROOT/scripts/build_graph.sh" graph.yml

# 2. Read k-mers
echo "=== Step 2/3: Read processing ==="
bash "$ROOT/scripts/process_reads.sh" data/HG002.MHC_R1.10x.fastq.gz data/HG002.MHC_R2.10x.fastq.gz HG002 work/reads

# 3. Diploid inference
echo "=== Step 3/3: Diploid inference ==="
if [ -x "$ROOT/bin/dip_iqp" ]; then
  bash "$ROOT/scripts/run_iqp.sh" work/ref_graph/ref_graph.txt work/reads HG002 work/inferred
else
  echo "bin/dip_iqp not built (make iqp, needs Gurobi), skipping the IQP step" >&2
fi

echo "=== Done in ${SECONDS}s ==="
