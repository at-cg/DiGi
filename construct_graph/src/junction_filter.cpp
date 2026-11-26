// Junction filter for construct_graph/scripts/liftover.sh: picks the junctions that
// every haplotype keeps, in reference order, and writes each haplotype's junction BED.
//
// Input (inside <outdir>, written by liftover.sh):
//   junctions_ref.bed                        junctions on the reference, in reference order
//   panel/<hap>.fa                           haplotype assembly, first sequence used
//   liftover/<hap>/<hap>_junctions_raw.bed   CrossMap output for that haplotype
//
// Pass 1, failed lifts. A junction lifts cleanly onto a haplotype when it has exactly
// one CrossMap row that
//   - has 10 columns, "->" in column 5 and "map_ratio=1.0000" in column 10,
//   - lands on the first sequence of the haplotype FASTA (column 6, maybe with "chr" in front),
//   - lands inside the assembly: 0 < position (column 7) < length.
// A junction that does not lift cleanly onto every haplotype is dropped from all of them.
//
// Pass 2, junction order. The m junctions left are j = 0..m-1 in reference order, and
// pos[k][j] is the position of junction j on haplotype k (k = 0..n-1). We want the
// largest set S of junctions that keeps the reference order on every haplotype:
// for each pair i < j in S and each k, pos[k][i] < pos[k][j].
// Longest chain DP:
//   dp[j]   = 1 + max { dp[i] : i < j and pos[k][i] < pos[k][j] for all k }   (1 if no such i)
//   prev[j] = the i that gives the max (-1 if none)
// S is the chain that ends at the j with the largest dp[j], read back through prev[].
// Time O(n * m^2).
//
// Output (inside <outdir>):
//   bed/<hap>/<hap>_junctions.bed   CrossMap rows of the junctions in S, reference order
//   lift_fails.txt                  one line per failed lift: <junction>\t<hap>
//   order_drops.txt                 junctions removed by the order check, one per line
//
// Build:  make graph (into construct_graph/bin/), run by construct_graph/scripts/liftover.sh
// Usage:  junction_filter <outdir>

#include <sys/stat.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph_file.hpp"

// Split a line at tabs.
static std::vector<std::string> splitTabs(const std::string& line) {
    std::vector<std::string> cols;
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            cols.push_back(line.substr(start));
            return cols;
        }
        cols.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
}

// Haplotype names: <outdir>/panel/<hap>.fa, sorted (the order graph_file.hpp uses).
static std::vector<std::string> listHaplotypes(const std::string& panel_dir) {
    DIR* dir = opendir(panel_dir.c_str());
    if (!dir) throw std::runtime_error("Missing folder: " + panel_dir);
    std::vector<std::string> haps;
    while (const dirent* e = readdir(dir)) {
        const std::string f = e->d_name;
        if (f.size() > 3 && f.compare(f.size() - 3, 3, ".fa") == 0) haps.push_back(f.substr(0, f.size() - 3));
    }
    closedir(dir);
    if (haps.empty()) throw std::runtime_error("No .fa files in " + panel_dir);
    std::sort(haps.begin(), haps.end());
    return haps;
}

// What one haplotype's CrossMap output says about each junction.
struct HapLift {
    std::vector<int> clean_rows;       // clean rows per junction (we want exactly 1)
    std::vector<long long> pos;        // position of the clean row (column 7)
    std::vector<std::string> row;      // the clean row itself
};

// Read <hap>_junctions_raw.bed and count the clean rows of each junction.
static HapLift readLift(const std::string& raw_bed, const std::string& contig, long long length,
                        const std::unordered_map<std::string, int>& junction_index) {
    const size_t m = junction_index.size();
    HapLift lift{std::vector<int>(m, 0), std::vector<long long>(m, 0), std::vector<std::string>(m)};

    std::ifstream fin(raw_bed);
    if (!fin.is_open()) throw std::runtime_error("Cannot open CrossMap output: " + raw_bed);
    std::string line;
    while (std::getline(fin, line)) {
        std::vector<std::string> c = splitTabs(line);
        if (c.size() != 10 || c[4] != "->" || c[9] != "map_ratio=1.0000") continue;
        if (c[5] != contig && c[5] != "chr" + contig) continue;
        long long p = std::stoll(c[6]);
        if (p <= 0 || p >= length) continue;

        auto it = junction_index.find(c[3]);
        if (it == junction_index.end()) continue;
        int j = it->second;
        lift.clean_rows[j] += 1;
        lift.pos[j] = p;
        lift.row[j] = line;
    }
    return lift;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <outdir>\n";
        return 1;
    }
    const std::string out = argv[1];

    try {
        // Reference junctions, in reference order
        std::vector<std::string> junctions;
        std::unordered_map<std::string, int> junction_index;
        {
            std::ifstream fin(out + "/junctions_ref.bed");
            if (!fin.is_open()) throw std::runtime_error("Cannot open " + out + "/junctions_ref.bed");
            std::string line;
            while (std::getline(fin, line)) {
                if (line.empty()) continue;
                std::vector<std::string> c = splitTabs(line);
                if (c.size() < 4) throw std::runtime_error("Fewer than 4 columns in junctions_ref.bed");
                junction_index[c[3]] = static_cast<int>(junctions.size());
                junctions.push_back(c[3]);
            }
        }
        const int total = static_cast<int>(junctions.size());

        // CrossMap results of every haplotype
        std::vector<std::string> haps = listHaplotypes(out + "/panel");
        const int n = static_cast<int>(haps.size());
        std::vector<HapLift> lifts;
        for (const std::string& hap : haps) {
            std::string contig, seq;
            readFirstSequence(out + "/panel/" + hap + ".fa", contig, seq);
            std::string raw = out + "/liftover/" + hap + "/" + hap + "_junctions_raw.bed";
            lifts.push_back(readLift(raw, contig, static_cast<long long>(seq.size()), junction_index));
        }

        // ---------------------------------------------------------------
        // Pass 1: keep the junctions that lift cleanly onto every haplotype
        // ---------------------------------------------------------------
        std::ofstream fails(out + "/lift_fails.txt");
        std::vector<int> lifted;  // indices into junctions, reference order
        for (int j = 0; j < total; ++j) {
            bool ok = true;
            for (int k = 0; k < n; ++k) {
                if (lifts[k].clean_rows[j] != 1) {
                    fails << junctions[j] << '\t' << haps[k] << '\n';
                    ok = false;
                }
            }
            if (ok) lifted.push_back(j);
        }

        // ---------------------------------------------------------------
        // Pass 2: longest chain that keeps the reference order on every haplotype
        // ---------------------------------------------------------------
        const int m = static_cast<int>(lifted.size());
        // pos[k][a] = position of the a-th lifted junction on haplotype k
        std::vector<std::vector<long long>> pos(n, std::vector<long long>(m));
        for (int k = 0; k < n; ++k)
            for (int a = 0; a < m; ++a) pos[k][a] = lifts[k].pos[lifted[a]];

        std::vector<int> dp(m, 1), prev(m, -1);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < j; ++i) {
                if (dp[i] + 1 <= dp[j]) continue;  // cannot improve dp[j]
                bool before_everywhere = true;
                for (int k = 0; k < n; ++k) {
                    if (pos[k][i] >= pos[k][j]) {
                        before_everywhere = false;
                        break;
                    }
                }
                if (before_everywhere) {
                    dp[j] = dp[i] + 1;
                    prev[j] = i;
                }
            }
        }

        // Read the best chain back from its last junction
        std::vector<bool> in_chain(m, false);
        if (m > 0) {
            int last = static_cast<int>(std::max_element(dp.begin(), dp.end()) - dp.begin());
            for (int a = last; a != -1; a = prev[a]) in_chain[a] = true;
        }

        std::ofstream drops(out + "/order_drops.txt");
        std::vector<int> kept;  // indices into junctions, reference order
        for (int a = 0; a < m; ++a) {
            if (in_chain[a]) kept.push_back(lifted[a]);
            else drops << junctions[lifted[a]] << '\n';
        }

        // ---------------------------------------------------------------
        // Output: one junction BED per haplotype
        // ---------------------------------------------------------------
        mkdir((out + "/bed").c_str(), 0755);
        for (int k = 0; k < n; ++k) {
            std::string dir = out + "/bed/" + haps[k];
            mkdir(dir.c_str(), 0755);
            std::ofstream bed(dir + "/" + haps[k] + "_junctions.bed");
            if (!bed.is_open()) throw std::runtime_error("Cannot write " + dir + "/" + haps[k] + "_junctions.bed");
            for (int j : kept) bed << lifts[k].row[j] << '\n';
        }

        const int r = static_cast<int>(kept.size());
        std::cout << total << " junctions, " << (total - m) << " fail to lift in >=1 of " << n
                  << " haplotypes, " << (m - r) << " out of order, " << r << " kept (r = " << r << ", "
                  << (r + 1) << " nodes per haplotype)\n";
    } catch (const std::exception& e) {
        std::cerr << "junction_filter: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
