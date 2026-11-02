# DiGi

DiGi infers an individual's two haplotypes (their diploid genome sequence) for a region such as the MHC, using only that individual's short reads.

It compares the reads against a **haplotype trellis graph**. To build the graph, a set of reference haplotype assemblies is cut at shared junction points, which splits each assembly into the same number of nodes. DiGi then models each of the individual's two haplotypes as a path through the graph. A path can switch between reference haplotypes at junctions, a limited number of times.

## Installation

```bash
conda env create -f environment.yml && conda activate digi
make         # k-mer samplers
make iqp     # inference solver; needs GUROBI_HOME and a Gurobi license
```

## Pipeline

| Step | Script | Input | Output |
|---|---|---|---|
| 1. Graph construction | `scripts/build_graph.sh <graph.yml>` | One FASTA per reference haplotype, plus a BED file per haplotype with the junction coordinates lifted onto that haplotype | Graph node sequences, and k-mers sampled from each node |
| 2. Read preprocessing + inference | `scripts/process_reads.sh <R1> <R2> <sample> <out_dir>`, then<br>`scripts/run_iqp.sh <graph> <node_kmers> <read_kmers> <sample> <out_dir>` | Paired-end FASTQ reads filtered to the region, plus the graph and node k-mers from step 1 | Quality-filtered reads and their sampled k-mer counts, then two haplotype FASTAs, and for each one the reference haplotype it follows in each region |

Set input paths in [`config/graph.yml`](config/graph.yml). Junction BEDs on GRCh38
can be lifted to each assembly with `tmp/liftover.sh`.

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
