// Reader for the DiGi graph file written by graph_builder.py.
//
//   #panel <dir>                     relative paths resolve against the graph file
//   >HG00268#1#MHC<TAB>HG00268_1.fa  FASTA name, then file inside the panel dir
//   1 0 10000 19968 ... 4947040      first node id, then r+2 node boundaries
//
// Nodes are contiguous: node first_id+c spans [pos[c], pos[c+1]) of the assembly.

#pragma once

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

inline HaplotypeGraph loadGraph(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) throw std::runtime_error("Cannot open graph file: " + filename);

    namespace fs = std::filesystem;
    fs::path panel_dir;
    bool have_panel = false;
    HaplotypeGraph g;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line.rfind("#panel", 0) == 0) {
            std::istringstream ss(line.substr(6));
            std::string dir;
            std::getline(ss >> std::ws, dir);
            panel_dir = fs::path(filename).parent_path() / dir;  // keeps an absolute dir as is
            have_panel = true;
        } else if (line[0] == '>') {
            size_t tab = line.find('\t');
            if (tab == std::string::npos || !have_panel)
                throw std::runtime_error("Bad graph header: " + line);
            GraphHaplotype h;
            h.name = line.substr(1, tab - 1);
            h.fasta = (panel_dir / line.substr(tab + 1)).string();

            if (!std::getline(fin, line))
                throw std::runtime_error("Graph file ends after header of " + h.name);
            std::istringstream ss(line);
            ss >> h.first_node;
            for (long long p; ss >> p;) h.pos.push_back(p);
            g.haps.push_back(std::move(h));
        } else {
            throw std::runtime_error("Unexpected line in " + filename + ", rebuild it with build_graph.sh: " + line.substr(0, 40));
        }
    }

    if (g.haps.empty()) throw std::runtime_error("No haplotypes in graph file: " + filename);
    g.H = static_cast<int>(g.haps.size());
    g.r = static_cast<int>(g.haps[0].pos.size()) - 2;
    if (g.r < 0) throw std::runtime_error("Too few node boundaries in graph file: " + filename);

    // The IQP finds nodes by id, so ids must follow the row and column layout
    for (int row = 0; row < g.H; ++row) {
        const GraphHaplotype& h = g.haps[row];
        if (static_cast<int>(h.pos.size()) != g.r + 2)
            throw std::runtime_error("Unequal node counts in graph file at " + h.name);
        if (h.first_node != static_cast<long long>(row) * (g.r + 1) + 1)
            throw std::runtime_error("Unexpected first node id for " + h.name);
        for (size_t i = 1; i < h.pos.size(); ++i)
            if (h.pos[i] <= h.pos[i - 1])
                throw std::runtime_error("Node boundaries not increasing for " + h.name);
    }
    return g;
}

// Loads the first sequence of a haplotype's panel FASTA. The graph holds only
// coordinates, so a panel that changed after the graph was built must fail here
// rather than yield wrong nodes.
inline std::string loadHaplotypeSequence(const GraphHaplotype& h) {
    std::ifstream fin(h.fasta);
    if (!fin.is_open()) throw std::runtime_error("Cannot open panel FASTA: " + h.fasta);

    std::string seq, line, name;
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

    if (name != h.name)
        throw std::runtime_error(h.fasta + ": sequence " + name + " does not match graph name " + h.name);
    if (static_cast<long long>(seq.size()) != h.pos.back())
        throw std::runtime_error(h.fasta + ": length " + std::to_string(seq.size()) +
                                 " does not match graph length " + std::to_string(h.pos.back()));
    return seq;
}
