// Input is a DiGi graph folder (panel/<hap>.fa and bed/<hap>/<hap>_junctions.bed),
// read with graph_file.hpp. Node ids follow the IQP lattice:
//   haplotype node (row, col) -> row*(r+1) + col + 1     sequence of [pos[col], pos[col+1])
//   junction node c           -> H*(r+1) + c + 1         no sequence ("*", LN:i:0)
//
//   S <node_id> <seq | *> [LN:i:<len> if *] hp:Z:<hap name> rr:Z:<start>-<end>
//   S <junction_id> * LN:i:0
//   L <from_id> + <to_id> + 0M
// The IQP source (0) and sink are left out, they carry no sequence.
//
// Build:  make graph (into construct_graph/bin/), run by construct_graph/scripts/liftover.sh
// Usage:  gfa_const <graph_dir> <out.gfa> [--no-seq]
//         --no-seq writes "*" instead of node sequences (smaller file, LN:i: kept)

#include <fstream>
#include <iostream>
#include <string>

#include "graph_file.hpp"

// Node id of junction c: numbered after all H*(r+1) haplotype nodes
static long long junctionNode(int c, int H, int r) {
    return 1LL * H * (r + 1) + c + 1;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_dir> <out.gfa> [--no-seq]\n";
        return 1;
    }
    const std::string graph_dir = argv[1];
    const std::string out_path = argv[2];
    const bool write_seq = !(argc > 3 && std::string(argv[3]) == "--no-seq");

    try {
        const HaplotypeGraph g = loadGraph(graph_dir);
        const int H = g.H, r = g.r;

        std::ofstream out(out_path);
        if (!out.is_open()) throw std::runtime_error("Cannot write " + out_path);
        out << "H\tVN:Z:1.0\n";

        // Segments: one per haplotype node, tagged with its haplotype (hp) and its
        // range [start, end) on that haplotype's assembly (rr)
        for (int row = 0; row < H; ++row) {
            const GraphHaplotype& h = g.haps[row];
            const std::string seq = write_seq ? loadHaplotypeSequence(h) : std::string();
            for (int col = 0; col <= r; ++col) {
                const long long start = h.pos[col], end = h.pos[col + 1];
                out << "S\t" << hapNode(row, col, r) << "\t";
                if (write_seq) out << seq.substr(start, end - start);
                else           out << "*\tLN:i:" << (end - start);
                out << "\thp:Z:" << h.name << "\trr:Z:" << start << "-" << end << "\n";
            }
        }

        // Segments: one empty node per junction
        for (int c = 0; c < r; ++c)
            out << "S\t" << junctionNode(c, H, r) << "\t*\tLN:i:0\n";

        // Links along each haplotype row
        for (int row = 0; row < H; ++row)
            for (int col = 0; col < r; ++col)
                out << "L\t" << hapNode(row, col, r) << "\t+\t" << hapNode(row, col + 1, r) << "\t+\t0M\n";

        // Links through each junction: every row's column c -> junction c -> every row's column c+1
        for (int c = 0; c < r; ++c)
            for (int row = 0; row < H; ++row) {
                out << "L\t" << hapNode(row, c, r) << "\t+\t" << junctionNode(c, H, r) << "\t+\t0M\n";
                out << "L\t" << junctionNode(c, H, r) << "\t+\t" << hapNode(row, c + 1, r) << "\t+\t0M\n";
            }

        std::cout << "Wrote " << out_path << ": " << H * (r + 1) + r << " segments ("
                  << H << " haplotypes x " << r + 1 << " nodes + " << r << " junctions), "
                  << H * r + 2LL * H * r << " links\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
