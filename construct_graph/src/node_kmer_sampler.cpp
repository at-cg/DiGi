// Node k-mers for src/iqp/main.cpp (--graph_kmers), from a DiGi graph folder
// (panel/<hap>.fa and bed/<hap>/<hap>_junctions.bed, read with graph_file.hpp).
// - For each non-empty node label: enumerate k-mers (k = 31), canonicalize, keep if FracMinHash hash < scale x 2^32 (scale = 0.05).
// - Output canonical k-mer **strings**, repeats kept (novel penalty counts occurrences)
//
// Node ids follow the IQP lattice (same as gfa_const.cpp):
//   haplotype node (row, col) -> row*(r+1) + col + 1     sequence of [pos[col], pos[col+1])
// Junction nodes carry no sequence and are not written.
//
// Output, one line per node with at least one sampled k-mer:
//   <node_id>\t<kmer>,<kmer>,...
//
// The hash is MurmurHash3_32, seed 0, on the canonical k-mer, with threshold
// fraction * 2^32, exactly as src/kmer_sampling/kmer_set_sampling.cpp, so the
// k-mers sampled from reads and from the graph agree (a lower fraction gives a
// subset of a higher one).
//
// Build:  make graph (into construct_graph/bin/), run by construct_graph/scripts/liftover.sh
// Usage:  node_kmer_sampler <graph_dir> <out.tsv> [fraction=0.05] [k=31]

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "graph_file.hpp"

// MurmurHash3 by Austin Appleby (public domain): https://github.com/aappleby/smhasher
static inline uint32_t rotl32(uint32_t x, int8_t r) { return (x << r) | (x >> (32 - r)); }

static inline uint32_t fmix32(uint32_t h) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

static uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    const uint32_t c1 = 0xcc9e2d51, c2 = 0x1b873593, n = 0xe6546b64;
    uint32_t hash = seed;
    const size_t nblocks = len / 4;

    for (size_t i = 0; i < nblocks; ++i) {
        uint32_t k;
        std::memcpy(&k, key + i * 4, sizeof(uint32_t));
        k *= c1;
        k = rotl32(k, 15);
        k *= c2;
        hash ^= k;
        hash = rotl32(hash, 13);
        hash = hash * 5 + n;
    }

    const uint8_t* tail = key + nblocks * 4;
    uint32_t k1 = 0;
    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16; [[fallthrough]];
        case 2: k1 ^= tail[1] << 8;  [[fallthrough]];
        case 1: k1 ^= tail[0];
                k1 *= c1;
                k1 = rotl32(k1, 15);
                k1 *= c2;
                hash ^= k1;
    }
    hash ^= static_cast<uint32_t>(len);
    return fmix32(hash);
}

static char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'A';  // 'T'
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_dir> <out.tsv> [fraction=0.05] [k=31]\n";
        return 1;
    }
    const std::string graph_dir = argv[1];
    const std::string out_path = argv[2];

    try {
        const double fraction = argc > 3 ? std::stod(argv[3]) : 0.05;
        const int k = argc > 4 ? std::stoi(argv[4]) : 31;
        if (fraction <= 0.0 || fraction > 1.0) throw std::runtime_error("fraction must be in (0, 1]");
        if (k < 1) throw std::runtime_error("k must be positive");
        const double threshold = fraction * (static_cast<double>(0xFFFFFFFFUL) + 1.0);

        const HaplotypeGraph g = loadGraph(graph_dir);
        std::ofstream out(out_path);
        if (!out.is_open()) throw std::runtime_error("Cannot write " + out_path);

        std::string kmer(k, 'A'), rev(k, 'A');
        long long nodes_written = 0, kmers_written = 0;

        for (int row = 0; row < g.H; ++row) {
            const GraphHaplotype& h = g.haps[row];
            std::string asm_seq = loadHaplotypeSequence(h);
            for (char& c : asm_seq) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

            for (int col = 0; col <= g.r; ++col) {
                const long long start = h.pos[col], end = h.pos[col + 1];
                if (end - start < k) continue;

                std::string line;
                long long count = 0;
                for (long long i = start; i + k <= end; ++i) {
                    bool ok = true;  // only A, C, G, T k-mers
                    for (int j = 0; j < k; ++j) {
                        const char c = asm_seq[i + j];
                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { ok = false; break; }
                        kmer[j] = c;
                        rev[k - 1 - j] = complement(c);
                    }
                    if (!ok) continue;

                    const std::string& canon = kmer < rev ? kmer : rev;
                    const uint32_t hv = murmur3_32(reinterpret_cast<const uint8_t*>(canon.data()), k, 0);
                    if (static_cast<double>(hv) >= threshold) continue;

                    if (count++) line += ',';
                    line += canon;
                }
                if (!count) continue;

                out << hapNode(row, col, g.r) << '\t' << line << '\n';
                ++nodes_written;
                kmers_written += count;
            }
        }

        std::cerr << "Wrote " << out_path << ": " << kmers_written << " k-mers (k=" << k
                  << ", fraction=" << fraction << ") on " << nodes_written << " nodes\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
