# Example: HG002 MHC at 10x

Infers the two MHC haplotypes of HG002 from 10x Illumina reads, using a small
graph of 18 HGSVC haplotypes with a junction every 10 kb.

```bash
conda env create -f environment.yml && conda activate digi
make            # k-mer samplers
make iqp        # IQP solver, needs GUROBI_HOME and a Gurobi license
bash example/run_example.sh
```

Everything is written to `example/work/`. The inferred haplotypes are
`example/work/inferred/HG002_{1,2}_c1-6.5127_c2-4.0225_R-10_*.fa`. Expect a runtime of 5 minutes on a single CPU core.

## Data

| File | Contents |
|---|---|
| `data/HGSVC_MHC_panel_18hap_grid10kb.tar.xz` | `panel/<hap>.fa`: MHC sequences of 18 HGSVC haplotypes. `liftover/<hap>/<hap>_junctions_lifted_final.bed`: 438 1-bp junctions, one every 10 kb of GRCh38 MHC, lifted to each haplotype. |
| `data/HG002.MHC_R{1,2}.10x.fastq.gz` | HG002 NIST Ashkenazi trio Illumina 2x250 bp reads. Pairs with both mates aligned to GRCh38 chr6:28510125-33480577, downsampled to 10x MHC depth. |

## Settings

FracMinHash fraction 0.05, k = 31, C1 = 6.5127 (missing k-mer penalty),
C2 = 4.0225 (novel k-mer penalty), R = 10 (maximum recombinations per haplotype), where cost parameters C1 and C2 are set using the parameter estimation model in the main text.
The current pieline can only override $R$ with the `R` environment variables.
