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

#include "gurobi_c++.h"
#include "../graph/graph_file.hpp"

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

// 2. Process Graph K-mer file in format: <node_id>\t<multiset of kmers>
void processGraphKmers(const std::string& graph_file,
                       const std::unordered_set<std::string>& Sk,
                       std::map<std::string, std::vector<int>>& targetKmerToNodes,
                       std::vector<uint32_t>& novel_kmer_count,
                       int num_nodes)
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

        int node_id = -1;
        try {
            node_id = std::stoi(line.substr(0, tab_pos));
        } catch (...) {
            continue;
        }

        if (node_id < 0 || node_id >= num_nodes) continue; 

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

// Maps a lattice node id back to (row, col). Returns (-1,-1) for the source,
// sink, or a junction node, none of which belong to a single haplotype.
std::pair<int,int> latticeRowCol(int node_id, int H, int r, int sink) {
    if (node_id == 0 || node_id == sink) return {-1, -1};
    int junction_start = H * (r + 1) + 1;
    int junction_end   = H * (r + 1) + r; // inclusive
    if (node_id >= junction_start && node_id <= junction_end) return {-1, -1};
    int idx = node_id - 1;
    int row = idx / (r + 1);
    int col = idx % (r + 1);
    if (row >= H) return {-1, -1};
    return {row, col};
}

int run(int argc, char *argv[]){
    auto start = std::chrono::high_resolution_clock::now();
    std::string GRAPH_KMERS_FILE = "";
    std::string SK_FILE          = "";
    // c1 = missing target k-mer penalty, c2 = novel k-mer penalty.
    // Recombinations carry no objective cost; they are only limited by the --R budget.
    double c1 = 1.0, c2 = 1.0;
    int R = -1;
    std::string GRAPH_FILE = "";
    std::string sample_name      = "Sample";
    bool write_regions = false; // off for now; --write_regions writes *_haplotype_regions.tsv

    std::string out_dir      = ".";
    std::string rel_weights  = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph_kmers" && i + 1 < argc) GRAPH_KMERS_FILE = argv[++i];
        else if (arg == "--sk" && i + 1 < argc)     SK_FILE = argv[++i];
        else if (arg == "--c1" && i + 1 < argc)     c1 = std::stod(argv[++i]);
        else if (arg == "--c2" && i + 1 < argc)     c2 = std::stod(argv[++i]);
        else if (arg == "--R" && i + 1 < argc)      R = std::stoi(argv[++i]);  
        else if (arg == "--graph" && i + 1 < argc)  GRAPH_FILE = argv[++i];
        else if (arg == "--sample" && i + 1 < argc) sample_name = argv[++i];
        else if (arg == "--write_regions")           write_regions = true;
        else if (arg == "--out_dir" && i + 1 < argc) out_dir = argv[++i];  
        else if (arg == "--rel_weights" && i + 1 < argc) rel_weights = argv[++i];
    }

    if (GRAPH_KMERS_FILE.empty() || SK_FILE.empty() || GRAPH_FILE.empty()) {
        std::cerr << "Usage: " << argv[0] << " --graph_kmers <file> --sk <file> --graph <graph_dir> [--R <int>] [options]\n";
        return 1;
    }

    const HaplotypeGraph graph = loadGraph(GRAPH_FILE);
    const int H = graph.H, r = graph.r;
    std::cout << "Graph: H=" << H << ", r=" << r << "\n";
    std::cout << "Graph read: (haplotypes: " << H << ", junctions: " << r << ")\n";

    const int n = H * (r + 1) + r + 2;
    int sink = H * (r + 1) + r + 1;
    Graph g(n);
    // Build Graph
    for(int row = 0; row < H; row++){
        g.addEdge(0, row*(r+1) + 1);    // source-to-first-column edges
        g.addEdge(row*(r+1) + r + 1, sink); // last-column-to-sink edges
    }
    for(int row = 0; row < H; row++){
        for(int col = 0; col < r; col++) g.addEdge(row*(r+1) + col + 1, row*(r+1) + col + 2); // zero-weight edges along each haplotype row
    }
    for(int recomb = 0; recomb < r; ++recomb){
        int w_node = H*(r+1) + recomb + 1;
        for(int row = 0; row < H; row++){
            g.addEdge(row*(r+1) + recomb + 1, w_node);  // haplotype-to-junction edges
            g.addEdge(w_node, row*(r+1) + recomb + 2);  // junction-to-haplotype edges
        }
    }

    // --- PIPELINE EXECUTION ---
    std::cout << "Loading Target Kmers...\n";
    auto Sk = loadSk(SK_FILE);

    std::cout << "Processing Graph Kmers in memory...\n";
    std::map<std::string, std::vector<int>> kmerToNodes;
    std::vector<uint32_t> novel_kmer_count(g.n, 0);

    processGraphKmers(GRAPH_KMERS_FILE, Sk, kmerToNodes, novel_kmer_count, g.n);

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
        auto is_junction = [&](int node) { return node > H*(r+1) && node != sink; };
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

                    if (next > H*(r+1) && next != t) {
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
                        int row = latticeRowCol(node_id, H, r, t).first;
                        if (row < 0) continue;
                        if (prev_row >= 0 && row != prev_row) n_switches++;
                        prev_row = row;
                    }
                    if (n_switches != n_recomb) {
                        std::cout << "Warning: " << n_recomb << " junction visits but only " << n_switches
                                  << " haplotype switches -- some recombinations return to the same row\n";
                    }

                    std::string sequence;
                    std::map<int, std::string> assemblies; // only the rows on this path
                    for (int node_id : path_nodes) {
                        auto [row, col] = latticeRowCol(node_id, H, r, t);
                        if (row < 0) continue; // source / sink / junction
                        const GraphHaplotype& hap = graph.haps[row];
                        auto it = assemblies.find(row);
                        if (it == assemblies.end())
                            it = assemblies.emplace(row, loadHaplotypeSequence(hap)).first;
                        sequence.append(it->second, hap.pos[col], hap.pos[col + 1] - hap.pos[col]);
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
                                    rout << seg_start << "\t" << seg_end << "\t" << graph.haps[seg_row].name
                                         << "\t" << seg_row << "\t" << seg_cols << "\n";
                                }
                            };

                            for (int node_id : path_nodes) {
                                auto [row, col] = latticeRowCol(node_id, H, r, t);
                                if (row < 0) continue; // source / sink / junction

                                long long node_start = graph.haps[row].pos[col];
                                long long node_end   = graph.haps[row].pos[col + 1];

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
            std::string file_path_1 = out_dir + "/" + sample_name + "_1_" + rel_weights + ".fa";
            std::string file_path_2 = out_dir + "/" + sample_name + "_2_" + rel_weights + ".fa";

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
