# DiGi
DiGi infers an individual's two haplotypes (their diploid genome sequence) for the MHC region, using only that individual's short reads and a simplified multi-haplotype reference graph. These sequences can be used for genotyping, however, one must note that these are partially phased and end-to-end phasing is not currently guaranteed.

> DiGi needs a Gurobi license to solve the integer program during inference. Free academic licenses are available at [gurobi.com/academics](https://www.gurobi.com/academics).

## Installation

```bash
git clone https://github.com/at-cg/DiGi.git && cd DiGi
conda env create -f environment.yml && conda activate digi
make
export GRB_LICENSE_FILE=/path/to/gurobi.lic
make iqp
```

## Usage
<!-- A sample reference graph (108 haplotype assemblies with 438 junctions, spaced roughly every 10 kb) is available on Zenodo(**!**). -->
The graph is hosted on Zenodo[link], download it into the <graph_dir> folder before running the code.
```bash
<!-- # 1. Graph preprocessing
bash scripts/preprocess_graph.sh <panel_dir> <bed_dir> <graph_dir> -->

# Inference
bash scripts/infer_haplotypes.sh -r1 data/<R1.fastq.gz> data/<R2.fastq.gz> -g <graph_dir> -o <out_dir> [-s <sample>] [-R <recombination_per_haplotype>]
```

example/data contains a small dataset in this format. `<R1.fastq.gz>` and `<R2.fastq.gz>` are paired-end reads filtered to the MHC region stored in data/ folder. As a result, `<sample>_{1,2}_*.fa` are written to `<out_dir>`. $R$, the maximum number of recombinations per haplotype is set to 10, by default. 
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
