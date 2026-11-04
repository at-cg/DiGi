# DiGi

DiGi infers an individual's two haplotypes (their diploid genome sequence) for the MHC region, using only that individual's short reads.

It compares the reads against a **haplotype trellis graph**. To build the graph, a set of reference haplotype assemblies in a common orientation is cut at shared junction points, which splits each assembly into the same number of nodes. DiGi then models each of the individual's two haplotypes as a path through the graph. A path can switch between reference haplotypes at junctions, a limited number of times specified by a recombination budget $R\leq r$ per haplotype, where $r$ is the number of junctions.

## Installation

```bash
conda env create -f environment.yml && conda activate digi
make         # k-mer samplers
make iqp     # inference solver; needs GUROBI_HOME and a Gurobi license
```

## Pipeline
The pipeline has two steps:

1. **Graph construction.** Junctions defined on a linear reference (e.g., GRCh38 no-alt) are lifted over to each haplotype assembly in the reference panel. The assemblies are cut at these junctions into nodes, and the k-mers of each node are sampled as a multiset using FracMinHash (scaling fraction 0.05 by default).
2. **Inference.** The target individual's short reads are reduced to sampled k-mers, which are compared against the graph from step 1 to infer the two haplotypes.

A sample reference graph (108 haplotype assemblies with 438 junctions, spaced roughly every 10 kb) is available on Zenodo.

```bash
# 1a. Lift the junctions over to each panel assembly
#     Output: <liftover_dir>/<hap>/<hap>_junctions.bed
bash tmp/liftover.sh <GRCh38.fa> <junctions.bed> <panel_dir> <liftover_dir>

# 1b. Build the graph and sample its node k-mers
#     Output: <graph_dir>/ref_graph.txt, <graph_dir>/node_kmers_<fraction>_k<k>.tsv
bash scripts/build_graph.sh config/graph.yml

# 2a. Filter the reads and sample their k-mers
#     Output: <reads_dir>/<sample>_read_kmers_<fraction>_k<k>.tsv
bash scripts/process_reads.sh <R1.fastq.gz> <R2.fastq.gz> <sample> <reads_dir>

# 2b. Infer the two haplotypes
#     Output: <out_dir>/<sample>_{1,2}_*.fa, <out_dir>/<sample>_{1,2}_*_haplotype_regions.tsv
bash scripts/run_iqp.sh <graph_dir>/ref_graph.txt <reads_dir> <sample> <out_dir>
```

`<R1.fastq.gz>` and `<R2.fastq.gz>` are paired-end reads filtered to the region. `<graph_dir>` is the folder of `graph_file` in the config. Step 2b finds the k-mer files automatically from `<graph_dir>` and `<reads_dir>`.

**Panel requirements.** Before liftover, each panel FASTA must contain only the MHC region of that assembly, not the whole genome, and all assemblies must share a common orientation.

Set input paths in [`config/graph.yml`](config/graph.yml).

**Settings.** Each one is set with an environment variable:

| Variable | Default | Meaning |
|---|---|---|
| `FRACTION` | 0.05 | FracMinHash sampling fraction, the share of k-mers kept |
| `K` | 31 | k-mer length |
| `THREADS` | 1 | Threads used by fastp and KMC in read preprocessing |
| `C1` | 5.2 | Penalty for a k-mer the path expects but the reads lack |
| `C2` | 3.5 | Penalty for a read k-mer the path lacks |
| `R` | 10 | Maximum switches between reference haplotypes per inferred haplotype |

Steps 1 and 2 must use the same `FRACTION` and `K`.

## Example

Infers HG002's MHC haplotypes from 10x coverage reads, using an 18-haplotype graph (about 5 minutes on one core):

```bash
bash example/run_example.sh
```

See [`example/README.md`](example/README.md) for details.
