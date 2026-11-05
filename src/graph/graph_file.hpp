// Reader for a DiGi graph folder:
//
//   panel/<hap>.fa                  assembly of one haplotype, first sequence used
//   bed/<hap>/<hap>_junctions.bed   its junctions, position on <hap> in column 7
//   node_kmers_0.05_k31.tsv         node k-mers, read separately
//
// Haplotypes are numbered in sorted <hap> order, which is the node numbering of
// the node k-mers. Nodes are contiguous: node first_id+c spans [pos[c], pos[c+1])
// of the assembly.

#pragma once

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct GraphHaplotype {
    std::string name;            // e.g. "HG00096#1#MHC"
    std::string fasta;           // full path of its panel FASTA
    long long first_node = 0;    // id of its first node
    std::vector<long long> pos;  // r+2 increasing node boundaries
};

struct HaplotypeGraph {
    int H = 0;                         // haplotypes
    int r = 0;                         // junctions per haplotype
    std::vector<GraphHaplotype> haps;  // in row order
};

// Loads the name and sequence of the first record of a FASTA file.
inline void readFirstSequence(const std::string& path, std::string& name, std::string& seq) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot open panel FASTA: " + path);

    name.clear();
    seq.clear();
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) break;  // only the first sequence is used
            name = line.substr(1, line.find_first_of(" \t") - 1);
        } else {
            seq += line;
        }
    }
}

// Node boundaries of a haplotype: 0, its junctions, its length.
inline std::vector<long long> nodeBoundaries(const std::string& bed, long long length) {
    std::ifstream fin(bed);
    if (!fin.is_open()) throw std::runtime_error("Cannot open junction BED: " + bed);

    std::vector<long long> pos{0};
    std::string line, field;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        for (int col = 0; col < 7; ++col)
            if (!std::getline(ss, field, '\t'))
                throw std::runtime_error("Fewer than 7 columns in " + bed);
        pos.push_back(std::stoll(field));
    }
    if (pos.back() < length) pos.push_back(length);

    for (size_t i = 1; i < pos.size(); ++i)
        if (pos[i] <= pos[i - 1])
            throw std::runtime_error("Junctions not increasing in " + bed);
    if (pos.back() != length)
        throw std::runtime_error("Junctions past the end of the assembly in " + bed);
    return pos;
}

inline HaplotypeGraph loadGraph(const std::string& graph_dir) {
    namespace fs = std::filesystem;
    const fs::path panel_dir = fs::path(graph_dir) / "panel";
    const fs::path bed_dir = fs::path(graph_dir) / "bed";
    if (!fs::is_directory(panel_dir)) throw std::runtime_error("Missing folder: " + panel_dir.string());

    std::vector<std::string> names;
    for (const auto& entry : fs::directory_iterator(panel_dir))
        if (entry.path().extension() == ".fa") names.push_back(entry.path().stem().string());
    if (names.empty()) throw std::runtime_error("No .fa files in " + panel_dir.string());
    std::sort(names.begin(), names.end());

    HaplotypeGraph g;
    long long next_node = 1;
    std::string seq;
    for (const std::string& hap : names) {
        GraphHaplotype h;
        h.fasta = (panel_dir / (hap + ".fa")).string();
        readFirstSequence(h.fasta, h.name, seq);
        if (seq.empty()) throw std::runtime_error("Empty FASTA: " + h.fasta);

        const fs::path bed = bed_dir / hap / (hap + "_junctions.bed");
        h.pos = nodeBoundaries(bed.string(), static_cast<long long>(seq.size()));
        h.first_node = next_node;
        next_node += static_cast<long long>(h.pos.size()) - 1;

        // The IQP finds nodes by id, so every row (hap) needs the same number of nodes
        if (!g.haps.empty() && h.pos.size() != g.haps[0].pos.size())
            throw std::runtime_error("Unequal node counts at " + h.name);
        g.haps.push_back(std::move(h));
    }

    g.H = static_cast<int>(g.haps.size());
    g.r = static_cast<int>(g.haps[0].pos.size()) - 2;
    if (g.r < 0) throw std::runtime_error("Too few node boundaries in " + graph_dir);
    return g;
}

// Loads the first sequence of a haplotype's panel FASTA.
inline std::string loadHaplotypeSequence(const GraphHaplotype& h) {
    std::string name, seq;
    readFirstSequence(h.fasta, name, seq);
    return seq;
}
