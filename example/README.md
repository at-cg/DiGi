# Example: HG002 MHC at 10x

Infers the two MHC haplotypes of HG002 from 10x Illumina reads, using a small
graph of 18 HGSVC haplotypes with a junction every 10 kb.

```bash
git clone https://github.com/at-cg/DiGi.git && cd DiGi
conda env create -f environment.yml && conda activate digi
make            # k-mer sampler
export GRB_LICENSE_FILE=/path/to/gurobi.lic
make iqp        # IQP solver
mkdir -p example/work && tar -C example/work -xJf example/data/graph_18_haps.tar.xz
bash scripts/infer_haplotypes.sh -r1 example/data/HG002.MHC_R1.10x.fastq.gz -r2 example/data/HG002.MHC_R2.10x.fastq.gz -g example/work/graph_18_haps -o example/work/inferred -s HG002
```

Everything is written to `example/work/`. The inferred haplotypes are
`example/work/inferred/HG002_{1,2}_c1-6.3245_c2-3.1333_R-10.fa`. Expect a runtime of 5 minutes on a single CPU core.

## Data

| File | Contents |
|---|---|
| `data/graph_18_haps.tar.xz` | `graph_18_haps/panel/<hap>.fa`: MHC sequences of 18 HGSVC haplotypes. `graph_18_haps/bed/<hap>/<hap>_junctions.bed`: 438 1-bp junctions, one every 10 kb of GRCh38 MHC, lifted to each haplotype. `graph_18_haps/node_kmers_0.05_k31.tsv`: sampled k-mers of each graph node. |
| `data/HG002.MHC_R{1,2}.10x.fastq.gz` | HG002 NIST Ashkenazi trio Illumina 2x250 bp reads. Cleaned with fastp, then pairs with both mates aligned to GRCh38 chr6:28510125-33480577, downsampled to 10x MHC depth. |

## Settings

Cost parameters C1 and C2 are estimated from the reads using the parameter estimation model in the main text.
To change $R$, add `-R 5` to the command above.
