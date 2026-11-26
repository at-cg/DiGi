#include <fstream>
#include <map>
#include <unordered_set>
#include <set>
#include <vector>
#include <iostream>
#include <chrono>
#include <cstdint>
#include <string>
#include <filesystem>
#include <thread>
#include <unordered_map>

#include "gurobi_c++.h"

// --- Graph Structure ---
class Graph{
    public:
        int n;
        std::vector<std::vector<int>> adj;
        Graph(int n);
        void addEdge(int u, int v);
};

// --- Graph Structure Methods ---
Graph::Graph(int n): n(n), adj(n) {}
void Graph::addEdge(int u, int v){ adj[u].push_back(v); }

// --- GFA graph (construct_graph/src/gfa_const.cpp) ---
// Segments get indices 1..N in file order; 0 is the source and N+1 the sink, added
// here: the source links to every segment without incoming links, and every segment
// without outgoing links links to the sink. A segment with an hp:Z:<haplotype> tag
// is a haplotype node (rr:Z:<start>-<end> is its range on that haplotype); one
// without is a junction, and passing through it is a recombination.
// Only the graph is kept in memory; segment sequences are read after solving, for
// the nodes on the two paths (loadSegmentSequences).
struct GfaGraph {
    int N = 0;                                     // segments
    int sink = 0;                                  // N + 1
    std::vector<std::string> name;                 // by index, "" for source and sink
    std::unordered_map<std::string, int> index;    // segment name -> index
    std::vector<int> hap;                          // haplotype number, -1 for junctions
    std::vector<long long> start, end;             // rr range, -1 if absent
    std::vector<std::string> hap_names;            // in order of first appearance
    std::vector<std::vector<int>> adj;             // links, plus source and sink edges
};

// Tab-separated fields of a GFA line, as [begin, end) offsets
static std::vector<std::pair<size_t, size_t>> gfaFields(const std::string& line) {
    std::vector<std::pair<size_t, size_t>> f;
    size_t b = 0;
    while (true) {
        size_t e = line.find('\t', b);
        if (e == std::string::npos) { f.emplace_back(b, line.size()); break; }
        f.emplace_back(b, e);
        b = e + 1;
    }
    return f;
}

GfaGraph loadGfa(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot open GFA: " + path);

    GfaGraph g;
    g.name.push_back("");  // source
    g.hap.push_back(-1); g.start.push_back(-1); g.end.push_back(-1);
    std::unordered_map<std::string, int> hap_index;
    std::vector<std::pair<std::string, std::string>> links;

    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.size() < 2 || line[1] != '\t') continue;
        const auto f = gfaFields(line);
        auto field = [&](size_t i) { return line.substr(f[i].first, f[i].second - f[i].first); };

        if (line[0] == 'S') {
            if (f.size() < 3) throw std::runtime_error("Short S line in " + path);
            const std::string seg = field(1);
            if (!g.index.emplace(seg, static_cast<int>(g.name.size())).second)
                throw std::runtime_error("Duplicate segment " + seg + " in " + path);
            g.name.push_back(seg);
            int hp = -1;
            long long s = -1, e = -1;
            for (size_t i = 3; i < f.size(); ++i) {
                const std::string tag = field(i);
                if (tag.rfind("hp:Z:", 0) == 0) {
                    auto it = hap_index.emplace(tag.substr(5), static_cast<int>(g.hap_names.size())).first;
                    if (it->second == static_cast<int>(g.hap_names.size())) g.hap_names.push_back(tag.substr(5));
                    hp = it->second;
                } else if (tag.rfind("rr:Z:", 0) == 0) {
                    const size_t dash = tag.find('-', 5);
                    if (dash == std::string::npos) throw std::runtime_error("Bad rr tag on segment " + seg);
                    s = std::stoll(tag.substr(5, dash - 5));
                    e = std::stoll(tag.substr(dash + 1));
                }
            }
            g.hap.push_back(hp); g.start.push_back(s); g.end.push_back(e);
        } else if (line[0] == 'L') {
            if (f.size() < 5) throw std::runtime_error("Short L line in " + path);
            if (field(2) != "+" || field(4) != "+")
                throw std::runtime_error("Only + to + links are supported: " + field(1) + " -> " + field(3));
            links.emplace_back(field(1), field(3));
        }
    }
    g.N = static_cast<int>(g.name.size()) - 1;
    if (g.N == 0) throw std::runtime_error("No segments in " + path);
    g.sink = g.N + 1;
    g.name.push_back("");  // sink
    g.hap.push_back(-1); g.start.push_back(-1); g.end.push_back(-1);

    g.adj.assign(g.N + 2, {});
    std::vector<char> has_in(g.N + 2, 0), has_out(g.N + 2, 0);
    for (const auto& [from, to] : links) {
        auto u = g.index.find(from), v = g.index.find(to);
        if (u == g.index.end() || v == g.index.end())
            throw std::runtime_error("Link to an unknown segment: " + from + " -> " + to);
        g.adj[u->second].push_back(v->second);
        has_out[u->second] = has_in[v->second] = 1;
    }
    for (int i = 1; i <= g.N; ++i) {
        if (!has_in[i]) g.adj[0].push_back(i);
        if (!has_out[i]) g.adj[i].push_back(g.sink);
    }
    if (g.adj[0].empty()) throw std::runtime_error("Every segment has an incoming link (cycle?) in " + path);
    return g;
}

// Sequences of the given segments (by index), read from the GFA
std::map<int, std::string> loadSegmentSequences(const std::string& path, const GfaGraph& g,
                                                const std::set<int>& wanted) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot open GFA: " + path);
    std::map<int, std::string> seqs;
    std::string line;
    while (seqs.size() < wanted.size() && std::getline(fin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.size() < 2 || line[0] != 'S' || line[1] != '\t') continue;
        const size_t t1 = line.find('\t', 2);
        if (t1 == std::string::npos) continue;
        const int idx = g.index.at(line.substr(2, t1 - 2));
        if (!wanted.count(idx)) continue;
        size_t t2 = line.find('\t', t1 + 1);
        if (t2 == std::string::npos) t2 = line.size();
        std::string seq = line.substr(t1 + 1, t2 - t1 - 1);
        if (seq == "*") {
            if (g.end[idx] > g.start[idx])
                throw std::runtime_error("Segment " + g.name[idx] + " has no sequence in " + path +
                                         " (written with gfa_const --no-seq?)");
            seq.clear();
        }
        seqs.emplace(idx, std::move(seq));
    }
    return seqs;
}

// 1. Load Target Read K-mers (Sk) as an unordered_set (ignoring frequencies)
std::unordered_set<std::string> loadSk(const std::string& sk_file) {
    std::unordered_set<std::string> Sk;
    std::ifstream file(sk_file);
    if (!file.is_open()) throw std::runtime_error("Cannot open Sk file: " + sk_file);

    std::string line;
    while(std::getline(file, line)){
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();

        // Target kmers are stored as <kmer>\t<number_of_occurrences>
        // We only extract the kmer and ignore the occurrences
        size_t tab = line.find('\t');
        if (tab != std::string::npos) {
            line = line.substr(0, tab);
        }
        Sk.insert(line);
    }
    return Sk;
}

// 2. Process Graph K-mer file in format: <segment name>\t<multiset of kmers>
//    (construct_graph/src/node_kmer_sampler.cpp); every name must be a GFA segment
void processGraphKmers(const std::string& graph_file,
                       const std::unordered_set<std::string>& Sk,
                       std::map<std::string, std::vector<int>>& targetKmerToNodes,
                       std::vector<uint32_t>& novel_kmer_count,
                       const std::unordered_map<std::string, int>& segment_index)
{
    std::ifstream file(graph_file);
    if (!file.is_open()) throw std::runtime_error("Cannot open Graph Kmers: " + graph_file);

    // Temporary map to collect unique node IDs for each target k-mer across different lines
    std::map<std::string, std::set<int>> temp_target_nodes;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();

        size_t tab_pos = line.find('\t');
        if (tab_pos == std::string::npos) {
            tab_pos = line.find(' ');
            if (tab_pos == std::string::npos) continue;
        }

        auto seg = segment_index.find(line.substr(0, tab_pos));
        if (seg == segment_index.end())
            throw std::runtime_error("Node k-mers name a segment that is not in the GFA: " + line.substr(0, tab_pos));
        const int node_id = seg->second;

        std::string kmers_str = line.substr(tab_pos + 1);

        // Optional: remove surrounding brackets if present in file format
        if (!kmers_str.empty() && kmers_str.front() == '[') kmers_str.erase(0, 1);
        if (!kmers_str.empty() && kmers_str.back() == ']') kmers_str.pop_back();

        size_t start = 0;
        while (start < kmers_str.size()) {
            size_t comma_pos = kmers_str.find(',', start);
            if (comma_pos == std::string::npos) comma_pos = kmers_str.size();

            if (comma_pos > start) {
                std::string kmer = kmers_str.substr(start, comma_pos - start);

                // Trim surrounding whitespace (e.g. "[AAA, CCC]") so it doesn't count as novel
                size_t first = kmer.find_first_not_of(" \t");
                if (first == std::string::npos) kmer.clear();
                else kmer = kmer.substr(first, kmer.find_last_not_of(" \t") - first + 1);

                if (!kmer.empty()) {
                    if (Sk.find(kmer) != Sk.end()) {
                        temp_target_nodes[kmer].insert(node_id);
                    } else {
                        novel_kmer_count[node_id]++; // =1 if no frqeuency needed
                    }
                }
            }
            start = comma_pos + 1;
        }
    }

    // Convert sets to vectors for the final targetKmerToNodes map
    for (const auto& pair : temp_target_nodes) {
        targetKmerToNodes[pair.first].assign(pair.second.begin(), pair.second.end());
    }
}

int run(int argc, char *argv[]){
    auto start = std::chrono::high_resolution_clock::now();
    std::string GRAPH_KMERS_FILE = "";
    std::string SK_FILE          = "";
    // c1 = missing target k-mer penalty, c2 = novel k-mer penalty.
    // Recombinations carry no objective cost; they are only limited by the --R budget.
    double c1 = 1.0, c2 = 1.0;
    int R = -1;
    std::string GFA_FILE = "";
    std::string sample_name      = "Sample";
    bool write_regions = false; // off for now; --write_regions writes *_haplotype_regions.tsv

    std::string out_dir      = ".";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph_kmers" && i + 1 < argc) GRAPH_KMERS_FILE = argv[++i];
        else if (arg == "--sk" && i + 1 < argc)     SK_FILE = argv[++i];
        else if (arg == "--c1" && i + 1 < argc)     c1 = std::stod(argv[++i]);
        else if (arg == "--c2" && i + 1 < argc)     c2 = std::stod(argv[++i]);
        else if (arg == "--R" && i + 1 < argc)      R = std::stoi(argv[++i]);  
        else if (arg == "--gfa" && i + 1 < argc)    GFA_FILE = argv[++i];
        else if (arg == "--sample" && i + 1 < argc) sample_name = argv[++i];
        else if (arg == "--write_regions")           write_regions = true;
        else if (arg == "--out_dir" && i + 1 < argc) out_dir = argv[++i];
    }

    if (GRAPH_KMERS_FILE.empty() || SK_FILE.empty() || GFA_FILE.empty()) {
        std::cerr << "Usage: " << argv[0] << " --graph_kmers <file> --sk <file> --gfa <graph.gfa> [--R <int>] [options]\n";
        return 1;
    }

    const GfaGraph gfa = loadGfa(GFA_FILE);
    const int sink = gfa.sink;
    int n_junctions = 0;
    for (int i = 1; i <= gfa.N; ++i) n_junctions += gfa.hap[i] < 0;
    std::cout << "Graph read: (haplotypes: " << gfa.hap_names.size() << ", junctions: " << n_junctions
              << ", segments: " << gfa.N << ")\n";

    // Build Graph: GFA links plus source and sink edges
    Graph g(gfa.N + 2);
    for (int u = 0; u < g.n; ++u)
        for (int v : gfa.adj[u]) g.addEdge(u, v);

    // --- PIPELINE EXECUTION ---
    std::cout << "Loading Target Kmers...\n";
    auto Sk = loadSk(SK_FILE);

    std::cout << "Processing Graph Kmers in memory...\n";
    std::map<std::string, std::vector<int>> kmerToNodes;
    std::vector<uint32_t> novel_kmer_count(g.n, 0);

    processGraphKmers(GRAPH_KMERS_FILE, Sk, kmerToNodes, novel_kmer_count, gfa.index);

    // --- GUROBI MODEL ---
    int exit_code = 0;
    try {
        GRBEnv env(true);
        env.set("LogFile", "dip_iqp.log");

        const int max_attempts = 3;
        for (int attempt = 1; attempt <= max_attempts; ++attempt) {
            try {
                env.start();
                break;
            } catch (GRBException& e) {
                if (attempt == max_attempts) throw;
                std::cerr << "WLS checkout failed (attempt " << attempt << "/"
                          << max_attempts << "): " << e.getMessage()
                          << " -- retrying...\n";
                std::this_thread::sleep_for(std::chrono::seconds(5 * attempt * attempt));
            }
        }

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "Diploid Path Inference");
        model.set(GRB_IntParam_Threads, 1); // serial solve per run - remove this to allow Gurobi to use multiple threads for solving
        
        std::vector<GRBVar> x_1_vars, x_2_vars;
        std::vector<std::pair<int,int>> edge_list;
        std::vector<std::vector<int>> out_edge_idx1(g.n), in_edge_idx1(g.n);
        std::vector<std::vector<int>> out_edge_idx2(g.n), in_edge_idx2(g.n);
        auto is_junction = [&](int node) { return node != 0 && node != sink && gfa.hap[node] < 0; };
        auto recomb_weight = [&](int u, int v) { return (is_junction(u) || is_junction(v)) ? 0.5 : 0.0; };    // assign weight 1/2 to each recombination edge, so that the total recombination count is the sum of weights along the path
        std::vector<double> edge_w;

        GRBQuadExpr obj = 0;

        for(int u = 0; u < g.n; ++u){
            for(int v : g.adj[u]){
                GRBVar var1 = model.addVar(0, 1, 0, GRB_BINARY, "x1_" + std::to_string(u) + "_" + std::to_string(v));
                GRBVar var2 = model.addVar(0, 1, 0, GRB_BINARY, "x2_" + std::to_string(u) + "_" + std::to_string(v));
                x_1_vars.push_back(var1); out_edge_idx1[u].push_back(x_1_vars.size()-1); in_edge_idx1[v].push_back(x_1_vars.size()-1);
                x_2_vars.push_back(var2); out_edge_idx2[u].push_back(x_2_vars.size()-1); in_edge_idx2[v].push_back(x_2_vars.size()-1);
                edge_list.emplace_back(u,v);
                edge_w.push_back(recomb_weight(u, v));

                double cost = c2 * novel_kmer_count[v];

                obj += cost * (var1 + var2);    // Add the novel k-mer penalty for each edge leading to a node with novel k-mers
            }
        }

        std::vector<std::string> kmer_list;
        for (const auto& pair : kmerToNodes) kmer_list.push_back(pair.first);

        std::vector<GRBVar> y_1_vars(kmer_list.size()), y_2_vars(kmer_list.size());
        for(size_t i = 0; i < kmer_list.size(); ++i){
            y_1_vars[i] = model.addVar(0, 1, 0, GRB_BINARY, "y1_" + kmer_list[i]);
            y_2_vars[i] = model.addVar(0, 1, 0, GRB_BINARY, "y2_" + kmer_list[i]);

            obj += c1 * (1.0 - y_1_vars[i]) * (1.0 - y_2_vars[i]); // Add the missing target k-mer penalty for each k-mer not covered by either path
            GRBLinExpr cover1 = 0, cover2 = 0;
            for(int v : kmerToNodes[kmer_list[i]]){
                for(int idx1 : in_edge_idx1[v]) cover1 += x_1_vars[idx1];
                for(int idx2 : in_edge_idx2[v]) cover2 += x_2_vars[idx2];
            }
            model.addConstr(cover1 >= y_1_vars[i]);     // Ensure that if a target k-mer is covered by path 1, then y_1_vars[i] must be 1
            model.addConstr(cover2 >= y_2_vars[i]);     // Ensure that if a target k-mer is covered by path 2, then y_2_vars[i] must be 1
        }

        model.setObjective(obj, GRB_MINIMIZE);

        for(int i = 0; i < g.n; i++){
            GRBLinExpr in1 = 0, out1 = 0, in2 = 0, out2 = 0;
            for(int idx : in_edge_idx1[i]) in1 += x_1_vars[idx];    // Sum of incoming edges to node i for path 1
            for(int idx : out_edge_idx1[i]) out1 += x_1_vars[idx];  // Sum of outgoing edges from node i for path 1
            for(int idx : in_edge_idx2[i]) in2 += x_2_vars[idx];    // Sum of incoming edges to node i for path 2
            for(int idx : out_edge_idx2[i]) out2 += x_2_vars[idx];  // Sum of outgoing edges from node i for path 2
            // Flow conservation constraints for each node in the graph
            if(i == 0) { model.addConstr(out1 - in1 == 1); model.addConstr(out2 - in2 == 1); }
            else if(i == sink) { model.addConstr(out1 - in1 == -1); model.addConstr(out2 - in2 == -1); }
            else { model.addConstr(out1 - in1 == 0); model.addConstr(out2 - in2 == 0); }
        }
        // Add recombination constraints if R is non-negative (i.e., if a recombination budget is specified), otherwise allow unlimited recombinations
        if (R >= 0) {
            GRBLinExpr nrec1 = 0, nrec2 = 0;
            for (size_t idx = 0; idx < edge_list.size(); ++idx) {
                if (edge_w[idx] == 0.0) continue;
                nrec1 += edge_w[idx] * x_1_vars[idx];
                nrec2 += edge_w[idx] * x_2_vars[idx];
            }
            model.addConstr(nrec1 <= R, "max_recomb_1");
            model.addConstr(nrec2 <= R, "max_recomb_2");
            std::cout << "Recombination cap: at most " << R << " per haplotype (W = 1/2 per junction edge)\n";
        }

        std::cout << "Optimizing IQP Model...\n";
        model.optimize();
        int status = model.get(GRB_IntAttr_Status);     // Get the optimization status

        if (status != GRB_OPTIMAL && status != GRB_SUBOPTIMAL) {
            std::cout << "No feasible solution found." << std::endl;
            exit_code = 1;
        } else {
            if (status == GRB_SUBOPTIMAL) std::cout << "Warning: Gurobi returned a suboptimal solution\n";
            double obj_val = model.get(GRB_DoubleAttr_ObjVal);

            std::cout << "\n========================================\n";
            std::cout << "Objective = " << obj_val << "\n";
            std::cout << "========================================\n";

            std::map<int,int> next_node_1;
            std::map<int,int> next_node_2;

            for (size_t e = 0; e < x_1_vars.size(); ++e) {
                if (x_1_vars[e].get(GRB_DoubleAttr_X) > 0.5) {
                    auto [u, v] = edge_list[e];
                    next_node_1[u] = v;
                }
            }

            for (size_t e = 0; e < x_2_vars.size(); ++e) {
                if (x_2_vars[e].get(GRB_DoubleAttr_X) > 0.5) {
                    auto [u, v] = edge_list[e];
                    next_node_2[u] = v;
                }
            }

            // Sequences of the segments on either path, read from the GFA in one pass.
            // Each path is walked with its own visited set: the two paths may share nodes.
            std::set<int> on_paths;
            for (const auto* next_node : {&next_node_1, &next_node_2}) {
                std::set<int> visited;
                for (int curr = 0; next_node->count(curr) && visited.insert(next_node->at(curr)).second; ) {
                    curr = next_node->at(curr);
                    if (curr != sink) on_paths.insert(curr);
                }
            }
            const std::map<int, std::string> path_seqs = loadSegmentSequences(GFA_FILE, gfa, on_paths);

            auto process_and_save_path = [&](const std::map<int,int>& next_node,
                                             const std::string& tag,
                                             const std::string& out_filename) -> int {
                std::cout << "\n--- Final Path Traversal (" << tag << ") ---\n";

                int src = 0;
                int t = sink;
                int curr = src;
                int n_recomb = 0;

                std::vector<int> path_nodes;
                std::set<int> visited;

                visited.insert(curr);
                path_nodes.push_back(curr);

                if (!next_node.count(src)) {
                    std::cout << "Error: No path from source\n";
                    return -1;
                }

                std::cout << curr;

                while (curr != t) {
                    auto it = next_node.find(curr);
                    if (it == next_node.end()) {
                        std::cout << " -> [BROKEN PATH]";
                        break;
                    }

                    int next = it->second;

                    if (visited.count(next)) {
                        std::cout << " -> " << next << " [CYCLE DETECTED]";
                        break;
                    }

                    std::cout << " -> " << next;

                    if (is_junction(next)) {
                        n_recomb++;
                    }

                    curr = next;
                    visited.insert(curr);
                    path_nodes.push_back(curr);
                }

                if (curr == t) {
                    std::cout << "\nPath verification: Reached Sink (" << t << ")";
                    std::cout << "\nNumber of recombinations: " << n_recomb << "\n";

                    // Debug check: every junction visit should be a real switch to a different row
                    int n_switches = 0, prev_row = -1;
                    for (int node_id : path_nodes) {
                        int row = gfa.hap[node_id];
                        if (row < 0) continue;
                        if (prev_row >= 0 && row != prev_row) n_switches++;
                        prev_row = row;
                    }
                    if (n_switches != n_recomb) {
                        std::cout << "Warning: " << n_recomb << " junction visits but only " << n_switches
                                  << " haplotype switches -- some recombinations return to the same row\n";
                    }

                    std::string sequence;
                    for (int node_id : path_nodes) {
                        if (node_id == src || node_id == t) continue; // source / sink have no sequence
                        auto it = path_seqs.find(node_id);
                        if (it == path_seqs.end())
                            throw std::runtime_error("No sequence loaded for segment " + gfa.name[node_id]);
                        sequence += it->second;
                    }

                    std::ofstream fout(out_filename);
                    if (fout.is_open()) {
                        fout << ">" << sample_name;
                        fout << "_" << tag << "_recomb" << n_recomb << "\n";

                        const size_t line_width = 60; // standard FASTA line wrap
                        for (size_t i = 0; i < sequence.size(); i += line_width) {
                            fout << sequence.substr(i, line_width) << "\n";
                        }
                        fout.close();
                        std::cout << "Saved path FASTA to: " << out_filename << "\n";
                    } else {
                        std::cerr << "Error writing to file: " << out_filename << "\n";
                    }
                    if (write_regions) {
                        std::string region_filename = out_filename;
                        size_t dot = region_filename.find_last_of('.');
                        if (dot != std::string::npos) region_filename = region_filename.substr(0, dot);
                        region_filename += "_haplotype_regions.tsv";

                        std::ofstream rout(region_filename);
                        if (rout.is_open()) {
                            rout << "start\tend\thaplotype\tgraph_row\tn_columns\n";

                            int seg_row = -2;      // -2 = no open segment yet
                            long long seg_start = -1, seg_end = -1;
                            int seg_cols = 0;

                            auto flush_segment = [&]() {
                                if (seg_row >= 0) {
                                    rout << seg_start << "\t" << seg_end << "\t" << gfa.hap_names[seg_row]
                                         << "\t" << seg_row << "\t" << seg_cols << "\n";
                                }
                            };

                            for (int node_id : path_nodes) {
                                int row = gfa.hap[node_id];
                                if (row < 0) continue; // source / sink / junction

                                long long node_start = gfa.start[node_id];
                                long long node_end   = gfa.end[node_id];
                                if (node_start < 0)
                                    throw std::runtime_error("Segment " + gfa.name[node_id] + " has no rr:Z: tag");

                                if (row != seg_row) {
                                    flush_segment();
                                    seg_row = row;
                                    seg_start = node_start;
                                    seg_cols = 0;
                                }
                                seg_end = node_end;
                                seg_cols++;
                            }
                            flush_segment();

                            rout.close();
                            std::cout << "Saved haplotype-region breakdown to: " << region_filename << "\n";
                        } else {
                            std::cerr << "Error writing to file: " << region_filename << "\n";
                        }
                    }

                    return n_recomb;
                } else {
                    std::cout << "\nWarning: Path did not reach Sink. File not written.\n";
                    return -1;
                }
            };

            std::filesystem::create_directories(out_dir);
            std::string file_path_1 = out_dir + "/" + sample_name + "_hap1.fa";
            std::string file_path_2 = out_dir + "/" + sample_name + "_hap2.fa";

            int n_recomb_1 = process_and_save_path(next_node_1, "1", file_path_1);
            int n_recomb_2 = process_and_save_path(next_node_2, "2", file_path_2);
            if (n_recomb_1 < 0 || n_recomb_2 < 0) exit_code = 1;

            double novel_cost = 0.0;
            for (size_t e = 0; e < x_1_vars.size(); ++e) {
                auto [u,v] = edge_list[e];
                if (x_1_vars[e].get(GRB_DoubleAttr_X) > 0.5) novel_cost += c2 * novel_kmer_count[v];
                if (x_2_vars[e].get(GRB_DoubleAttr_X) > 0.5) novel_cost += c2 * novel_kmer_count[v];
            }

            double missing_cost = 0.0;
            for (size_t i = 0; i < y_1_vars.size(); ++i) {
                bool y1 = y_1_vars[i].get(GRB_DoubleAttr_X) > 0.5;
                bool y2 = y_2_vars[i].get(GRB_DoubleAttr_X) > 0.5;
                if (!y1 && !y2) {
                    missing_cost += c1;
                }
            }

            std::cout << "\n-----------------------------------------\n";
            std::cout << "Recombinations (not in objective): " << n_recomb_1 << " + " << n_recomb_2 << "\n";
            std::cout << "Objective breakdown:\n";
            std::cout << "Novel kmer penalty: "      << novel_cost  << "\n";
            std::cout << "Missing kmer penalty: "    << missing_cost<< "\n";
            std::cout << "Total objective: "        << obj_val     << "\n";
        }

    } catch(GRBException& e){
        std::cerr << "Gurobi Error: " << e.getMessage() << std::endl;
        return 1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\nElapsed time: " << elapsed.count() << " seconds\n";

    return exit_code;
}

int main(int argc, char *argv[]){
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
