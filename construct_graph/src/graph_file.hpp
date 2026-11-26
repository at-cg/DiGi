// Input: a folder containing a panel of haplotype assemblies and their junctions.
//   panel/<hap>.fa                  assembly of one haplotype, first sequence used
//   bed/<hap>/<hap>_junctions.bed   its junctions, position on <hap> in column 7
//
// Haplotypes are numbered in sorted <hap> order. Node c of a haplotype spans
// [pos[c], pos[c+1]) of its assembly.

#pragma once

#include <dirent.h>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

struct GraphHaplotype {
    std::string name;            // FASTA record name, e.g. "HG00096#1#MHC"
    std::string fasta;           // path of its panel FASTA
    std::vector<long long> pos;  // r+2 increasing node boundaries: 0, junctions, length
};

struct HaplotypeGraph {
    int H = 0;                         // haplotypes
    int r = 0;                         // junctions per haplotype
    std::vector<GraphHaplotype> haps;  // in row order
};

// Node id of haplotype node (row, col) in a graph with r junctions: rows are
// numbered one after another, from 1.
inline long long hapNode(int row, int col, int r){
    return 1LL * row * (r + 1) + col + 1;
}

// Header name and sequence of the first record of a FASTA file.
inline void readFirstSequence(const std::string& path, std::string& name, std::string& seq){
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

// Node boundaries of a haplotype: 0, its junctions (BED column 7), its length.
inline std::vector<long long> nodeBoundaries(const std::string& bed, long long length) {
    std::ifstream fin(bed);
    if (!fin.is_open()) throw std::runtime_error("Cannot open junction BED: " + bed);

    std::vector<long long> pos{0};
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        size_t begin = 0;
        for (int col = 1; col < 7; ++col) {  // skip to the start of column 7
            begin = line.find('\t', begin);
            if (begin == std::string::npos) throw std::runtime_error("Fewer than 7 columns in " + bed);
            ++begin;
        }
        pos.push_back(std::stoll(line.substr(begin)));
    }
    if (pos.back() < length) pos.push_back(length);

    for (size_t i = 1; i < pos.size(); ++i)
        if (pos[i] <= pos[i - 1]) throw std::runtime_error("Junctions not increasing in " + bed);
    return pos;
}

inline HaplotypeGraph loadGraph(const std::string& graph_dir) {
    const std::string panel_dir = graph_dir + "/panel";
    const std::string bed_dir = graph_dir + "/bed";

    DIR* dir = opendir(panel_dir.c_str());
    if (!dir) throw std::runtime_error("Missing folder: " + panel_dir);
    std::vector<std::string> names;
    while (const dirent* e = readdir(dir)) {
        const std::string f = e->d_name;
        if (f.size() > 3 && f.compare(f.size() - 3, 3, ".fa") == 0) names.push_back(f.substr(0, f.size() - 3));
    }
    closedir(dir);
    if (names.empty()) throw std::runtime_error("No .fa files in " + panel_dir);
    std::sort(names.begin(), names.end());

    HaplotypeGraph g;
    std::string seq;
    for (const std::string& hap : names) {
        GraphHaplotype h;
        h.fasta = panel_dir + "/" + hap + ".fa";
        readFirstSequence(h.fasta, h.name, seq);
        if (seq.empty()) throw std::runtime_error("Empty FASTA: " + h.fasta);
        h.pos = nodeBoundaries(bed_dir + "/" + hap + "/" + hap + "_junctions.bed", static_cast<long long>(seq.size()));

        // Node ids are row*(r+1)+col+1, so every haplotype needs the same node count
        if (!g.haps.empty() && h.pos.size() != g.haps[0].pos.size())
            throw std::runtime_error("Unequal node counts at " + h.name);
        g.haps.push_back(std::move(h));
    }

    g.H = static_cast<int>(g.haps.size());
    g.r = static_cast<int>(g.haps[0].pos.size()) - 2;
    if (g.r < 0) throw std::runtime_error("Too few node boundaries in " + graph_dir);
    return g;
}

// Sequence of a haplotype's panel FASTA (first record).
inline std::string loadHaplotypeSequence(const GraphHaplotype& h) {
    std::string name, seq;
    readFirstSequence(h.fasta, name, seq);
    return seq;
}
