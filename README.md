# DiGi
DiGi is a pangenome-based method for diploid genome inference. Given a set of reference haplotypes and short reads from a target genome, it infers the two mosaics of reference haplotypes that best explain the reads. The underlying algorithm is not specific to any genomic region: it is applicable wherever complete reference haplotypes can be aligned at a common set of junction points. The current implementation has been primarily developed and tested on the human [major histocompatibility complex (MHC)](https://genome.ucsc.edu/cgi-bin/hgc?db=hg38&c=chr6&l=31944158&r=31944159&o=28510019&t=33480577&g=comments&i=MHC), a highly polymorphic region of about 5 Mbp on chromosome 6. 

DiGi outputs two FASTA files, one per inferred haplotype. These sequences can be passed to variant-calling tools such as [PAV](https://github.com/EichlerLab/pav) or gene-annotation tools such as [Immuannot](https://github.com/YingZhou001/Immuannot) for downstream analysis.

## Installation
Requires Linux (x86-64) and GCC 8.3 or newer.

### Step 1: Download source code and compile
```bash
git clone https://github.com/at-cg/DiGi.git
cd DiGi
make
```

`make` also downloads the external dependencies: [KMC](https://github.com/refresh-bio/KMC) and Gurobi.

### Step 2: Gurobi license setup 

> [!IMPORTANT]
> DiGi requires a Gurobi license to solve its integer program. It reads your Gurobi license from the `GRB_LICENSE_FILE` environment variable. For academic use, a license is easy to obtain.

1. Register or log in at the [Gurobi Academic Program](https://www.gurobi.com/academics) with your institutional email address and request a **WLS Academic** license.
2. Open the [Gurobi Web License Manager](https://license.gurobi.com/manager/licenses).
3. Go to **API Keys** and create an API key for your WLS Academic license.
4. Download the generated `gurobi.lic` file.
5. Point DiGi to it by using the following command:

```bash
export GRB_LICENSE_FILE=/path/to/gurobi.lic
```

Optional: To avoid running this command in every new session, add it to your `~/.bashrc`.


## General usage
DiGi takes paired-end short reads and an MHC pangenome graph as input:

```bash
bash scripts/infer_haplotypes.sh -r1 read1.fq -r2 read2.fq -g <graph_dir> -outdir <out_dir> -outname <prefix> 
```

`read1.fq` and `read2.fq` should contain paired-end reads filtered to the MHC region. After a successful run, DiGi writes `<prefix>_hap1.fa` and `<prefix>_hap2.fa` to `<out_dir>`.

A pre-built MHC pangenome graph of 54 individuals (108 haplotypes) is available on [Zenodo](https://zenodo.org/records/23974110). Download and extract the archive, then pass the extracted folder to `-g` as `<graph_dir>`.

> [!NOTE]
> Phasing is reliable within genes but may contain switch errors over long distances, so the two output haplotypes should not be assumed to be phased across the full MHC.

## Try DiGi on small data

This example infers the MHC haplotypes of HG002 genome from 10x coverage reads, using a small graph of 18 haplotypes. It takes about five minutes.

```bash
# From the DiGi project folder
# Download the small graph
wget https://github.com/at-cg/DiGi/releases/download/v0.1/graph_18_haps.tar.gz
tar -xzvf graph_18_haps.tar.gz

# Run inference
bash scripts/infer_haplotypes.sh \
  -r1 example/data/HG002.MHC_R1.10x.fastq.gz \
  -r2 example/data/HG002.MHC_R2.10x.fastq.gz \
  -g graph_18_haps -outdir DiGi_TestRun -outname HG002
```

When the run finishes, the two inferred haplotypes are in `DiGi_TestRun/HG002_hap1.fa` and `DiGi_TestRun/HG002_hap2.fa`. The folder also holds the solver logs and, under `reads/`, intermediate read k-mers, which you can ignore.

## Algorithm overview

DiGi formulates diploid genome inference as an integer quadratic program (IQP) that jointly reconstructs two haplotypes as mosaic paths through the pangenome graph. The optimization balances agreement between the inferred haplotypes and observed read $k$-mers while limiting recombinations between reference haplotypes.

<p align="center">
  <img src="docs/figures/digi_overview.svg" width="100%" alt="Overview of the DiGi algorithm">
</p>
