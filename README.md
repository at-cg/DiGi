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
# 1. Graph construction
bash scripts/construct_graph.sh <GRCh38.fa> <junctions.bed> <panel_dir> <graph_dir>

# 2. Inference
bash scripts/infer_haplotypes.sh <R1.fastq.gz> <R2.fastq.gz> <sample> <graph_dir> <out_dir> [R]
```

`<R1.fastq.gz>` and `<R2.fastq.gz>` are paired-end reads filtered to the region. Step 1 removes any junction that fails to lift to exactly one position in every assembly. Step 2 writes `<sample>_{1,2}_*.fa` and `<sample>_{1,2}_*_haplotype_regions.tsv` to `<out_dir>`.

<!-- **Panel requirements.** `<panel_dir>` holds one `<hap>.fa` per assembly. Each FASTA must contain only the MHC region of that assembly, as a single sequence, not the whole genome, and all assemblies must share a common orientation. -->

<!-- **Settings.** Each one is set with an environment variable; `R` can also be given as the last argument of step 2:

| Variable | Default | Meaning |
|---|---|---|
| `FRACTION` | 0.05 | FracMinHash sampling fraction, the share of k-mers kept |
| `K` | 31 | k-mer length |
| `C1` | 5.2 | Penalty for a k-mer the path expects but the reads lack |
| `C2` | 3.5 | Penalty for a read k-mer the path lacks |
| `R` | 10 | Maximum switches between reference haplotypes per inferred haplotype | -->


## Example

Infers HG002's MHC haplotypes from 10x coverage reads, using an 18-haplotype graph (about 5 minutes on one core):

```bash
bash example/run_example.sh
```

See [`example/README.md`](example/README.md) for details.
