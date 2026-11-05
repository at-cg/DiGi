# DiGi
DiGi infers an individual's two haplotypes (their diploid genome sequence) for the MHC region, using only that individual's short reads and a simplified multi-haplotype reference graph we refer to as a **haplotype trellis graph**. DiGi models each of the individual's two haplotypes as a path through the graph. A path can switch between reference haplotypes at junctions, a limited number of times specified by a recombination budget $R\leq r$ per haplotype, where $r$ is the number of junctions.

> DiGi needs a Gurobi license to solve the integer program during inference. Free academic licenses are available at [gurobi.com/academics](https://www.gurobi.com/academics).

## Installation

```bash
git clone https://github.com/at-cg/DiGi.git && cd DiGi
conda env create -f environment.yml && conda activate digi
make
export GRB_LICENSE_FILE=/path/to/gurobi.lic
make iqp
```

## Pipeline
The pipeline has two steps:

1. **Graph preprocessing.** We currently assume that each haplotype assembly comes with a BED file of its junctions, and that all assemblies share the same set of junctions, defined on a linear reference. The assemblies are cut at these junctions into nodes, and the k-mers of each node are sampled as a multiset using FracMinHash (scaling fraction 0.05 by default).
2. **Inference.** The target individual's short reads are reduced to sampled k-mers, which are compared against the graph from step 1 to infer the two haplotypes.

<!-- A sample reference graph (108 haplotype assemblies with 438 junctions, spaced roughly every 10 kb) is available on Zenodo(**!**). -->

```bash
# 1. Graph preprocessing
bash scripts/preprocess_graph.sh <panel_dir> <bed_dir> <graph_dir>

# 2. Inference
bash scripts/infer_haplotypes.sh <R1.fastq.gz> <R2.fastq.gz> <sample> <graph_dir> <out_dir> [R]
```

`<panel_dir>` holds one `<hap>.fa` per assembly and `<bed_dir>` holds its `<hap>/<hap>_junctions.bed`. example/data contains a small dataset in this format. `<R1.fastq.gz>` and `<R2.fastq.gz>` are paired-end reads filtered to the region. Step 2 writes `<sample>_{1,2}_*.fa` to `<out_dir>`. $R$, the recombination budget per haplotype is set to 10, by default. 
<!-- Rerunning step 2 with the same `<out_dir>` skips read processing if the reads, FracMinHash fraction and k are unchanged, so only inference runs again. -->

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

Infers HG002's MHC haplotypes from 10x coverage reads, using an 18-haplotype graph (about 5 minutes on one CPU core):

```bash
bash example/run_example.sh
```

See [`example/README.md`](example/README.md) for details.
