#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdint>

// Platform-specific
#if defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#define ROTL32(x,y) _rotl(x,y)
#else
#define FORCE_INLINE inline __attribute__((always_inline))
inline uint32_t rotl32(uint32_t x, int8_t r) {
    return (x << r) | (x >> (32 - r));
}
#define ROTL32(x,y) rotl32(x,y)
#endif

FORCE_INLINE uint32_t fmix32(uint32_t h) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const uint32_t r1 = 15;
    const uint32_t r2 = 13;
    const uint32_t m  = 5;
    const uint32_t n  = 0xe6546b64;

    uint32_t hash = seed;
    const int nblocks = len / 4;

    for (int i = 0; i < nblocks; i++) {
        uint32_t k;
        std::memcpy(&k, key + (i * 4), sizeof(uint32_t));
        k *= c1;
        k = ROTL32(k, r1);
        k *= c2;
        hash ^= k;
        hash = ROTL32(hash, r2);
        hash = hash * m + n;
    }

    const uint8_t* tail = key + nblocks * 4;
    uint32_t k1 = 0;
    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = ROTL32(k1, r1); k1 *= c2; hash ^= k1;
    };

    hash ^= len;
    hash = fmix32(hash);
    return hash;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <fraction> <kmer_tsv> [output_tsv]\n"
                  << "  fraction : e.g. 0.001 keeps ~1/1000 of k-mers\n"
                  << "  kmer_tsv : input file with <kmer>\\t<count> per line "
                     "(e.g. from `kmc_tools transform <db> dump`)\n"
                  << "             k-mers must already be canonical (KMC default, i.e. no -b)\n";
        return 1;
    }

    double fraction = std::stod(argv[1]);
    if (fraction <= 0.0 || fraction > 1.0) {
        std::cerr << "Error: fraction must be in (0, 1]\n";
        return 1;
    }

    std::ifstream in_file(argv[2]);
    if (!in_file.is_open()) {
        std::cerr << "Error: cannot open input file " << argv[2] << "\n";
        return 1;
    }

    std::ostream* out = &std::cout;
    std::ofstream file_out;
    if (argc >= 4) {
        file_out.open(argv[3]);
        if (!file_out.is_open()) {
            std::cerr << "Error: cannot open output file " << argv[3] << "\n";
            return 1;
        }
        out = &file_out;
    }

    // Threshold on the 32-bit hash space. A hash h passes if h < threshold,
    // which happens with probability ~= fraction.
    const double threshold = fraction * (static_cast<double>(0xFFFFFFFFUL) + 1.0);
    const uint32_t seed = 0;    // matches node_kmers_gen.cpp, so read and node k-mer samples agree

    std::string line;
    uint64_t total = 0, kept = 0;
    int expected_k = -1;

    while (std::getline(in_file, line)) {
        if (line.empty()) continue;

        auto tab = line.find('\t');
        if (tab == std::string::npos) {
            std::cerr << "Warning: skipping malformed line: " << line << "\n";
            continue;
        }

        // Zero-copy view of the k-mer portion
        const char* kmer_ptr = line.data();
        size_t      kmer_len = tab;
        uint64_t    count    = std::stoull(line.c_str() + tab + 1);

        if (expected_k < 0) expected_k = static_cast<int>(kmer_len);
        else if (static_cast<int>(kmer_len) != expected_k) {
            std::cerr << "Warning: inconsistent k-mer length ("
                      << kmer_len << " vs " << expected_k << ")\n";
        }

        ++total;

        // KMC (without -b) already emits canonical k-mers: the lexicographic min
        // of the k-mer and its reverse complement under A<C<G<T, which is the
        // same rule as get_canonical() in node_kmers_gen.cpp. Hash as-is.
        uint32_t h = murmur3_32(reinterpret_cast<const uint8_t*>(kmer_ptr),
                                kmer_len, seed);
        if (static_cast<double>(h) < threshold) {
            out->write(kmer_ptr, kmer_len);
            *out << '\t' << count << '\n';
            ++kept;
        }
    }

    std::cerr << "Input k-mers:    " << total << "\n"
              << "Sampled k-mers:  " << kept
              << "  (expected ~" << static_cast<uint64_t>(total * fraction) << ")\n"
              << "Effective scaled: "
              << (kept ? static_cast<double>(total) / kept : 0.0) << "\n";
    return 0;
}