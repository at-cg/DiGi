# Graph construction

## System requirements

- Linux x86-64, glibc 2.25 or newer
- GCC 8 or newer
- Python 3.8 or newer, with the `venv` module (`make graph PYTHON=/path/to/python3`)

## Usage

From the DiGi root:

```bash
make graph
bash construct_graph/scripts/liftover.sh <ref.fa> <fasta_dir> <spacing_bp> <outdir>
```

- `<ref.fa>`: linear reference
- `<fasta_dir>`: folder of reference assemblies, one `<hap>.fa` per haplotype
- `<spacing_bp>`: junction spacing on the linear reference, in bp (e.g. `10000`)
- `<outdir>`: output graph folder, to pass to `scripts/infer_haplotypes.sh -g`

**Note:** Only the first sequence of each FASTA file is used.

If `make graph` fails to verify a certificate, install a newer certificate list with pip and pass it as `CA_BUNDLE`:

```bash
pip install --upgrade certifi
make graph CA_BUNDLE=$(python3 -m certifi)
```
