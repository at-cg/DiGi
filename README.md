# DiGi
DiGi is a pangenome-based genome inference method. It is useful for genotyping and partial phasing of variants in the [MHC region](https://genome.ucsc.edu/cgi-bin/hgc?db=hg38&c=chr6&l=31944158&r=31944159&o=28510019&t=33480577&g=comments&i=MHC).
Given a pangenome index of MHC reference haplotypes (e.g., from haplotype-resolved long-read assemblies) and short-reads of donor genome, DiGi infers two mosaic paths in the pangenome that are most consistent with donor's genome. It gives two fasta files as output, each corresponding to an inferred haplotype sequence from the pangenome. The user can use these outputs and either run [variant calling](https://github.com/EichlerLab/pav) or [gene annotation tools](https://github.com/YingZhou001/Immuannot) for downstream applications.

> DiGi needs a Gurobi license to solve the integer program during inference. Free academic licenses are available at [gurobi.com/academics](https://www.gurobi.com/academics).

## Installation

```bash
git clone https://github.com/at-cg/DiGi.git && cd DiGi
conda env create -f environment.yml && conda activate digi
make
export GRB_LICENSE_FILE=/path/to/gurobi.lic
make iqp
```
<!-- # 1. Graph preprocessing
bash scripts/preprocess_graph.sh <panel_dir> <bed_dir> <graph_dir> -->
## Usage
<!-- A sample reference graph (108 haplotype assemblies with 438 junctions, spaced roughly every 10 kb) is available on Zenodo(**!**). -->
The graph is hosted on Zenodo[link], download it into the `<graph_dir>` folder before running the code.
```bash
# Inference
bash scripts/infer_haplotypes.sh -r1 R1.fastq -r2 R2.fastq -g <graph_dir> -outdir <out_dir> -outname <prefix> 
```

example/data contains a small dataset in this format. `R1.fastq` and `R2.fastq` are paired-end reads filtered to the MHC region (gzipped files also work). As a result, `<prefix>_{1,2}_*.fa` are written to `<out_dir>`.
<!-- Rerunning step 2 with the same `<out_dir>` skips read processing if the reads, FracMinHash fraction and k are unchanged, so only inference runs again. -->

<!-- **Panel requirements.** `<panel_dir>` holds one `<hap>.fa` per assembly. Each FASTA must contain only the MHC region of that assembly, as a single sequence, not the whole genome, and all assemblies must share a common orientation. -->

<!-- **Settings.** Each one is set with an environment variable:

| Variable | Default | Meaning |
|---|---|---|
| `FRACTION` | 0.05 | FracMinHash sampling fraction, the share of k-mers kept |
| `K` | 31 | k-mer length |
| `C1` | 5.2 | Penalty for a k-mer the path expects but the reads lack |
| `C2` | 3.5 | Penalty for a read k-mer the path lacks | -->


## Example

Infers HG002's MHC haplotypes from 10x coverage reads, using an 18-haplotype graph (about 5 minutes on one CPU core):

```bash
mkdir -p example/work && tar -C example/work -xJf example/data/graph_18_haps.tar.xz
bash scripts/infer_haplotypes.sh -r1 example/data/HG002.MHC_R1.10x.fastq.gz -r2 example/data/HG002.MHC_R2.10x.fastq.gz -g example/work/graph_18_haps -outdir example/work/inferred -outname HG002
```

See [`example/README.md`](example/README.md) for details.
