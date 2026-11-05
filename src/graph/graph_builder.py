import os
import sys
import yaml
import argparse

# Usage: python graph_builder.py --config graph.yml [OPTIONAL] --sample <sample_name>
#        [OPTIONAL] --samples_included <path/to/samples_included.txt>
# sample_name is especially crucial if the sample is a part of the reference graph
# samples_included.txt lists (one per line) the ONLY samples that should be
# included when building the graph

def load_config(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    fmt = {
        "base_dir": config["paths"]["base_dir"]
    }

    for key, value in config["graph"].items():
        if isinstance(value, str):
            config["graph"][key] = value.format(**fmt)

    return config

def load_single_fasta(path):
    """Load the first (or only) sequence from a single FASTA file."""
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    break  # only read the first sequence
                name = line[1:].split()[0]
            elif line:
                buf.append(line)
    return name, "".join(buf)

def node_boundaries(seq_len, bed_path):
    """Node boundaries of a haplotype: 0, its junctions (col7 of the BED), seq_len."""
    pos = [0]
    with open(bed_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) > 6:
                pos.append(int(parts[6]))
    if pos[-1] < seq_len:
        pos.append(seq_len)
    if any(b <= a for a, b in zip(pos, pos[1:])) or pos[-1] != seq_len:
        raise ValueError(f"junctions of {bed_path} are not increasing within 0..{seq_len}")
    return pos

def filter_headers(headers, exclude_sample=None):
    if not exclude_sample:
        return headers

    prefix = f"{exclude_sample}"
    return [h for h in headers if not h.startswith(prefix)]

def load_included_samples(path):
    """Read a text file of sample names (one per line) to include.
    Blank lines and lines starting with '#' are ignored."""
    included = set()
    with open(path) as f:
        for line in f:
            name = line.strip()
            if name and not name.startswith("#"):
                included.add(name)
    return included

if __name__ == "__main__":

    # load YAML config
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument("--sample", default=None, help="Sample to be inferred")

    parser.add_argument("--samples_included", default=None,
                         help="Optional path to a text file listing sample names "
                              "(one per line) to include; graph is built only for these")
    args = parser.parse_args()

    config = load_config(args.config)
    exclude_sample = args.sample

    # paths from config
    fasta_folder = config["graph"]["haplotype_fasta_dir"]
    out_path = config["graph"]["graph_file"]
    bed_folder = config["graph"]["rh_bed_file"]

# discover haplotypes from FASTA files
    headers = sorted(
        os.path.splitext(f)[0]
        for f in os.listdir(fasta_folder)
        if f.endswith((".fa", ".fasta"))
    )
    # apply exclusions
    headers = filter_headers(headers, exclude_sample=exclude_sample)

    if args.samples_included:
        included = load_included_samples(args.samples_included)
        before = len(headers)
        headers = [h for h in headers if h in included]
        print(f"Samples_included filter: kept {len(headers)}/{before} haplotypes "
              f"based on {args.samples_included}")

        missing = included - set(headers)
        if missing:
            print(f"!! {len(missing)} sample(s) from samples_included not found "
                  f"in {fasta_folder}: {sorted(missing)}")

    # One record per haplotype: its FASTA name and file, then its first node id and
    # node boundaries. Written only once every haplotype has passed.
    records = []
    next_node = 1
    for h in headers:
        fasta_file = os.path.join(fasta_folder, f"{h}.fa")
        bed_file = os.path.join(bed_folder, f"{h}/{h}_junctions.bed") # *** extract bed files in one folder

        seq_name, seq = load_single_fasta(fasta_file)
        if not seq:
            sys.exit(f"Empty or unreadable FASTA: {fasta_file}")

        pos = node_boundaries(len(seq), bed_file)
        records.append((seq_name, f"{h}.fa", next_node, pos))
        next_node += len(pos) - 1

    # Every haplotype must have the same number of nodes (junctions + 1)
    n_nodes = {len(r[3]) - 1 for r in records}
    if len(n_nodes) != 1:
        counts = sorted((len(r[3]) - 1, r[0]) for r in records)
        sys.exit(f"Unequal node counts: {counts[0]} ... {counts[-1]}")

    panel = os.path.relpath(fasta_folder, os.path.dirname(os.path.abspath(out_path)))
    with open(out_path, "w") as out:
        out.write(f"#panel\t{panel}\n")
        for name, fasta, first_node, pos in records:
            out.write(f">{name}\t{fasta}\n{first_node} {' '.join(map(str, pos))}\n")

    n = n_nodes.pop()
    print(f"      {len(records)} haplotypes, {n} nodes each ({n - 1} junctions)")
