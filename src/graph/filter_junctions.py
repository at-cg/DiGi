import os
import sys

# Keep only the junctions that lift to exactly one position in every haplotype,
# so that all haplotypes are cut at the same junctions.
# Usage: python filter_junctions.py <junctions.bed> <liftover_dir>
# <liftover_dir>/<hap>/<hap>_junctions.bed -> <hap>_junctions_lifted_final.bed


def junction_key(parts):
    return tuple(parts[:3])


def main(junctions_path, liftover_dir):
    with open(junctions_path) as f:
        junctions = [junction_key(line.rstrip("\n").split("\t")) for line in f if line.strip()]

    haps = sorted(
        h for h in os.listdir(liftover_dir)
        if os.path.exists(os.path.join(liftover_dir, h, f"{h}_junctions.bed"))
    )
    if not haps:
        sys.exit(f"no <hap>/<hap>_junctions.bed found in {liftover_dir}")

    lifted = {}
    for hap in haps:
        lines = {}
        with open(os.path.join(liftover_dir, hap, f"{hap}_junctions.bed")) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                # A lifted line has "->" in column 5; CrossMap marks others Fail/Unmap
                if len(parts) > 6 and parts[4] == "->":
                    lines.setdefault(junction_key(parts), []).append(line)
        lifted[hap] = lines

    keep = [j for j in junctions if all(len(lifted[h].get(j, [])) == 1 for h in haps)]
    if not keep:
        sys.exit("no junction lifted uniquely to every haplotype")

    for hap in haps:
        out_path = os.path.join(liftover_dir, hap, f"{hap}_junctions_lifted_final.bed")
        with open(out_path, "w") as out:
            for j in keep:
                out.write(lifted[hap][j][0])

    print(f"Kept {len(keep)}/{len(junctions)} junctions across {len(haps)} haplotypes")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: filter_junctions.py <junctions.bed> <liftover_dir>")
    main(sys.argv[1], sys.argv[2])
