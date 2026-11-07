# Example: HG002 MHC at 10x

Infers the two MHC haplotypes of HG002 from 10x Illumina reads, using a small
graph of 18 HGSVC haplotypes with a junction every 10 kb.

```bash
git clone https://github.com/at-cg/DiGi.git && cd DiGi
conda env create -f environment.yml && conda activate digi
make            # k-mer sampler
export GRB_LICENSE_FILE=/path/to/gurobi.lic
make iqp        # IQP solver
mkdir -p example/work
tar -C example/work -xJf example/data/graph_18_haps.tar.xz
bash scripts/infer_haplotypes.sh -r1 example/data/HG002.MHC_R1.10x.fastq.gz -r2 example/data/HG002.MHC_R2.10x.fastq.gz -g example/work/graph_18_haps -outdir example/work/inferred -outname HG002
```

Everything is written to `example/work/`. The inferred haplotypes are
`example/work/inferred/HG002_{1,2}_*.fa`. Expect a runtime of 5 minutes on a single CPU core.

## Data

| File | Contents |
|---|---|
| `data/graph_18_haps.tar.xz` | Supporting files for graph construction (haplotype FASTA, junction BEDs, sampled $k$-mers at each node). |
| `data/HG002.MHC_R{1,2}.10x.fastq.gz` | HG002 NIST Ashkenazi trio Illumina 2x250 bp reads. Reads mapping to the GRCh38 MHC (chr6:28,510,125–33,480,577), downsampled to 10× mean depth.|

## Settings

Cost parameters C1 and C2 are estimated from the reads using the parameter estimation model in the main text.
