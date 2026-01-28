#include <io/netmp.h>
#include <graphdb/offline_evaluator.h>
#include <graphdb/online_evaluator.h>
#include <utils/circuit.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <omp.h>

#include "utils.h"
#include "graphutils.h"

using namespace graphdb;
using json = nlohmann::json;
namespace bpo = boost::program_options;

void printDaglistInfo(const DistributedDaglist& dist_daglist, const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << title << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int nC = dist_daglist.num_clients;
    int total_vertices = 0, total_edges = 0;
    int total_change_vertices = 0, total_change_edges = 0;
    
    for (int i = 0; i < nC; ++i) {
        total_vertices += dist_daglist.VSizes[i];
        total_edges += dist_daglist.ESizes[i];
        
        // Count changes
        for (size_t j = 0; j < dist_daglist.VSizes[i]; ++j) {
            if (dist_daglist.isChangeV[i][j] == Ring(1)) total_change_vertices++;
        }
        for (size_t j = 0; j < dist_daglist.ESizes[i]; ++j) {
            if (dist_daglist.isChangeE[i][j] == Ring(1)) total_change_edges++;
        }
    }
    
    std::cout << "Total Vertices: " << total_vertices << " (Marked for update: " << total_change_vertices << ")" << std::endl;
    std::cout << "Total Edges: " << total_edges << " (Marked for update: " << total_change_edges << ")" << std::endl;
    std::cout << "Total Entries: " << (total_vertices + total_edges) << std::endl;
    std::cout << std::endl;
    
    // Print distribution per client
    std::cout << "Distribution per Client:" << std::endl;
    for (int i = 0; i < nC; ++i) {
        std::cout << "  Client " << i << ": " << dist_daglist.VSizes[i] << " vertices, "
                  << dist_daglist.ESizes[i] << " edges" << std::endl;
    }
    std::cout << std::endl;
    
    // Print sample vertices
    std::cout << "Sample Vertices (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd | Chg | NewData" << std::endl;
    std::cout << "  " << std::string(80, '-') << std::endl;
    int count = 0;
    for (int c = 0; c < nC && count < 10; ++c) {
        for (size_t i = 0; i < dist_daglist.VSizes[c] && count < 10; ++i) {
            const auto& v = dist_daglist.VertexLists[c][i];
            std::cout << "  " << std::setw(2) << count << " | "
                      << std::setw(3) << v.src << " | "
                      << std::setw(3) << v.dst << " | "
                      << std::setw(4) << v.data << " | "
                      << std::setw(3) << v.isV << " | "
                      << std::setw(4) << v.sigv << " | "
                      << std::setw(4) << v.sigs << " | "
                      << std::setw(4) << v.sigd << " | "
                      << std::setw(3) << dist_daglist.isChangeV[c][i] << " | "
                      << std::setw(7) << dist_daglist.ChangeV[c][i] << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    
    // Print sample edges
    std::cout << "Sample Edges (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd | Chg | NewData" << std::endl;
    std::cout << "  " << std::string(80, '-') << std::endl;
    count = 0;
    for (int c = 0; c < nC && count < 10; ++c) {
        for (size_t i = 0; i < dist_daglist.ESizes[c] && count < 10; ++i) {
            const auto& e = dist_daglist.EdgeLists[c][i];
            std::cout << "  " << std::setw(2) << count << " | "
                      << std::setw(3) << e.src << " | "
                      << std::setw(3) << e.dst << " | "
                      << std::setw(4) << e.data << " | "
                      << std::setw(3) << e.isV << " | "
                      << std::setw(4) << e.sigv << " | "
                      << std::setw(4) << e.sigs << " | "
                      << std::setw(4) << e.sigd << " | "
                      << std::setw(3) << dist_daglist.isChangeE[c][i] << " | "
                      << std::setw(7) << dist_daglist.ChangeE[c][i] << std::endl;
            count++;
        }
    }
    std::cout << std::string(60, '=') << std::endl;
}

void printOutputs(const std::vector<Ring>& outputs, const DistributedDaglist& dist_daglist) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "OUTPUT: Graph After Data Updates" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int nC = dist_daglist.num_clients;
    size_t expected_outputs = 0;
    for (int c = 0; c < nC; ++c) {
        expected_outputs += dist_daglist.VSizes[c] * 7; // 7 fields per vertex
        expected_outputs += dist_daglist.ESizes[c] * 7; // 7 fields per edge
    }
    
    if (outputs.size() < expected_outputs) {
        std::cout << "Warning: Received " << outputs.size() << " outputs, expected " << expected_outputs << std::endl;
        return;
    }
    
    std::cout << "Updated Vertices (up to first 10 from each client):" << std::endl;
    std::cout << "  ID | Src | Dst | isV | Data | sigv | sigs | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    // Output structure: for each client, vertices come first, then edges
    // Client 0: all vertices, then all edges
    // Client 1: all vertices, then all edges, etc.
    size_t output_idx = 0;
    int display_count = 0;
    
    for (int c = 0; c < nC && display_count < 10; ++c) {
        size_t limit = std::min(static_cast<size_t>(dist_daglist.VSizes[c]), size_t(10 - display_count));
        for (size_t i = 0; i < limit && output_idx + 6 < outputs.size(); ++i) {
            std::cout << "  " << std::setw(2) << display_count << " | "
                      << std::setw(3) << outputs[output_idx + 0] << " | "
                      << std::setw(3) << outputs[output_idx + 1] << " | "
                      << std::setw(3) << outputs[output_idx + 2] << " | "
                      << std::setw(4) << outputs[output_idx + 3] << " | "
                      << std::setw(4) << outputs[output_idx + 4] << " | "
                      << std::setw(4) << outputs[output_idx + 5] << " | "
                      << std::setw(4) << outputs[output_idx + 6] << std::endl;
            output_idx += 7;
            display_count++;
        }
        // Skip remaining vertices for this client
        output_idx += (dist_daglist.VSizes[c] - limit) * 7;
        // Skip all edges for this client (since we're only showing vertices now)
        output_idx += dist_daglist.ESizes[c] * 7;
    }
    
    std::cout << "\nUpdated Edges (up to first 10 from each client):" << std::endl;
    std::cout << "  ID | Src | Dst | isV | Data | sigv | sigs | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    // Reset to start of outputs and iterate through clients for edges
    output_idx = 0;
    display_count = 0;
    for (int c = 0; c < nC && display_count < 10; ++c) {
        // Skip vertices for this client
        output_idx += dist_daglist.VSizes[c] * 7;
        
        // Now print edges
        size_t limit = std::min(static_cast<size_t>(dist_daglist.ESizes[c]), size_t(10 - display_count));
        for (size_t i = 0; i < limit && output_idx + 6 < outputs.size(); ++i) {
            std::cout << "  " << std::setw(2) << display_count << " | "
                      << std::setw(3) << outputs[output_idx + 0] << " | "
                      << std::setw(3) << outputs[output_idx + 1] << " | "
                      << std::setw(3) << outputs[output_idx + 2] << " | "
                      << std::setw(4) << outputs[output_idx + 3] << " | "
                      << std::setw(4) << outputs[output_idx + 4] << " | "
                      << std::setw(4) << outputs[output_idx + 5] << " | "
                      << std::setw(4) << outputs[output_idx + 6] << std::endl;
            output_idx += 7;
            display_count++;
        }
        // Skip remaining edges for this client
        output_idx += (dist_daglist.ESizes[c] - limit) * 7;
    }
    
    std::cout << std::string(60, '=') << std::endl;
}

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, DistributedDaglist dist_daglist) {

    int nC = dist_daglist.num_clients;
    auto VSizes = dist_daglist.VSizes;
    auto ESizes = dist_daglist.ESizes;

    std::cout << "Generating circuit" << std::endl;

    common::utils::Circuit<Ring> circ; 

    // Initialize all daglist field values
    std::vector<std::vector<wire_t>> vertex_src_values(nC);
    std::vector<std::vector<wire_t>> vertex_dst_values(nC);
    std::vector<std::vector<wire_t>> vertex_isV_values(nC);
    std::vector<std::vector<wire_t>> vertex_data_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigd_values(nC);
    std::vector<std::vector<wire_t>> vertex_changed(nC);
    std::vector<std::vector<wire_t>> vertex_changed_data(nC);

    std::vector<std::vector<wire_t>> edge_src_values(nC);
    std::vector<std::vector<wire_t>> edge_dst_values(nC);
    std::vector<std::vector<wire_t>> edge_isV_values(nC);
    std::vector<std::vector<wire_t>> edge_data_values(nC);
    std::vector<std::vector<wire_t>> edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> edge_sigv_values(nC);
    std::vector<std::vector<wire_t>> edge_sigd_values(nC);
    std::vector<std::vector<wire_t>> edge_changed(nC);
    std::vector<std::vector<wire_t>> edge_changed_data(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_vertex_src_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_dst_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_isV_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_data_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigs_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigv_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigd_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_changed(VSizes[i]);
        std::vector<wire_t> subg_vertex_changed_data(VSizes[i]);
        
        std::vector<wire_t> subg_edge_src_values(ESizes[i]);
        std::vector<wire_t> subg_edge_dst_values(ESizes[i]);
        std::vector<wire_t> subg_edge_isV_values(ESizes[i]);
        std::vector<wire_t> subg_edge_data_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigs_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigv_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigd_values(ESizes[i]);
        std::vector<wire_t> subg_edge_changed(ESizes[i]);
        std::vector<wire_t> subg_edge_changed_data(ESizes[i]);

        for (int j = 0; j < VSizes[i]; ++j){
            subg_vertex_src_values[j] = circ.newInputWire();
            subg_vertex_dst_values[j] = circ.newInputWire();
            subg_vertex_isV_values[j] = circ.newInputWire();
            subg_vertex_data_values[j] = circ.newInputWire();
            subg_vertex_sigs_values[j] = circ.newInputWire();
            subg_vertex_sigv_values[j] = circ.newInputWire();
            subg_vertex_sigd_values[j] = circ.newInputWire();
            subg_vertex_changed[j] = circ.newInputWire();
            subg_vertex_changed_data[j] = circ.newInputWire();
        }

        for (int j = 0; j < ESizes[i]; ++j){
            subg_edge_src_values[j] = circ.newInputWire();
            subg_edge_dst_values[j] = circ.newInputWire();
            subg_edge_isV_values[j] = circ.newInputWire();
            subg_edge_data_values[j] = circ.newInputWire();
            subg_edge_sigs_values[j] = circ.newInputWire();
            subg_edge_sigv_values[j] = circ.newInputWire();
            subg_edge_sigd_values[j] = circ.newInputWire();
            subg_edge_changed[j] = circ.newInputWire();
            subg_edge_changed_data[j] = circ.newInputWire();
        }

        vertex_src_values[i] = subg_vertex_src_values;
        vertex_dst_values[i] = subg_vertex_dst_values;
        vertex_isV_values[i] = subg_vertex_isV_values;
        vertex_data_values[i] = subg_vertex_data_values;
        vertex_sigs_values[i] = subg_vertex_sigs_values;
        vertex_sigv_values[i] = subg_vertex_sigv_values;
        vertex_sigd_values[i] = subg_vertex_sigd_values;
        vertex_changed[i] = subg_vertex_changed;
        vertex_changed_data[i] = subg_vertex_changed_data;

        edge_src_values[i] = subg_edge_src_values;
        edge_dst_values[i] = subg_edge_dst_values;
        edge_isV_values[i] = subg_edge_isV_values;
        edge_data_values[i] = subg_edge_data_values;
        edge_sigs_values[i] = subg_edge_sigs_values;
        edge_sigv_values[i] = subg_edge_sigv_values;
        edge_sigd_values[i] = subg_edge_sigd_values;
        edge_changed[i] = subg_edge_changed;
        edge_changed_data[i] = subg_edge_changed_data;
    
    }

    // Update vertex data
    std::vector<std::vector<wire_t>> updated_vertex_list(nC);
    std::vector<std::vector<wire_t>> updated_edge_list(nC);
    for (int i = 0; i < nC; ++i){
        updated_vertex_list[i].resize(VSizes[i]);
        updated_edge_list[i].resize(ESizes[i]);
    }

    for (int i = 0; i < nC; ++i){
        for (int j = 0; j < VSizes[i]; ++j){
            // Compute: updated_data = data + indicator * (new_data - data)
            // If indicator=0: updated_data = data
            // If indicator=1: updated_data = data + (new_data - data) = new_data
            auto diff = 
                circ.addGate(common::utils::GateType::kSub, vertex_changed_data[i][j], vertex_data_values[i][j]);
            auto change = 
                circ.addGate(common::utils::GateType::kMul, vertex_changed[i][j], diff);
            updated_vertex_list[i][j] =
                circ.addGate(common::utils::GateType::kAdd, vertex_data_values[i][j], change);
        }
    }

    for (int i = 0; i < nC; ++i){
        for (int j = 0; j < ESizes[i]; ++j){
            auto diff = 
                circ.addGate(common::utils::GateType::kSub, edge_changed_data[i][j], edge_data_values[i][j]);
            auto change = 
                circ.addGate(common::utils::GateType::kMul, edge_changed[i][j], diff);
            updated_edge_list[i][j] =
                circ.addGate(common::utils::GateType::kAdd, edge_data_values[i][j], change);
        }
    }

    for (int i = 0; i < nC; ++i){

        for (int j = 0; j < VSizes[i]; ++j){
            circ.setAsOutput(vertex_src_values[i][j]);
            circ.setAsOutput(vertex_dst_values[i][j]);
            circ.setAsOutput(vertex_isV_values[i][j]);
            circ.setAsOutput(updated_vertex_list[i][j]);
            circ.setAsOutput(vertex_sigv_values[i][j]);
            circ.setAsOutput(vertex_sigs_values[i][j]);
            circ.setAsOutput(vertex_sigd_values[i][j]);
        }

        for (int j = 0; j < ESizes[i]; ++j){
            circ.setAsOutput(edge_src_values[i][j]);
            circ.setAsOutput(edge_dst_values[i][j]);
            circ.setAsOutput(edge_isV_values[i][j]);
            circ.setAsOutput(updated_edge_list[i][j]);
            circ.setAsOutput(edge_sigv_values[i][j]);
            circ.setAsOutput(edge_sigs_values[i][j]);
            circ.setAsOutput(edge_sigd_values[i][j]);
        }
    }

    return circ;
}

void benchmark(const bpo::variables_map& opts) {

    bool save_output = false;
    std::string save_file;
    if (opts.count("output") != 0) {
        save_output = true;
        save_file = opts["output"].as<std::string>();
    }

    auto num_vert = opts["num-vert"].as<size_t>();
    auto num_edge = opts["num-edge"].as<size_t>();
    auto vec_size = num_vert + num_edge;
    auto nP = opts["num-parties"].as<int>();
    auto nC = opts["num-clients"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto use_pking = opts["use-pking"].as<bool>();
    auto random_inputs = opts["random-inputs"].as<bool>();

    omp_set_nested(1);
    if (nP < 10) { omp_set_num_threads(nP); }
    else { omp_set_num_threads(10); }

    std::cout << "Starting benchmarks" << std::endl;

    std::string net_config = opts.count("net-config") ? opts["net-config"].as<std::string>() : "";
    std::shared_ptr<io::NetIOMP> network = createNetwork(pid, nP, latency, port,
                                                          opts["localhost"].as<bool>(),
                                                          net_config);

    json output_data;
    output_data["details"] = {{"num_parties", nP},
                              {"num_clients", nC},
                              {"num_vert", num_vert},
                              {"num_vert", num_vert},
                              {"latency (ms)", latency},
                              {"pid", pid},
                              {"threads", threads},
                              {"seed", seed},
                              {"repeat", repeat}};
    output_data["benchmarks"] = json::array();

    std::cout << "--- Details ---" << std::endl;
    for (const auto& [key, value] : output_data["details"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    // Generate test graph
    Ring nV = static_cast<Ring>(num_vert);
    Ring nE = static_cast<Ring>(num_edge);
    
    DistributedDaglist dist_daglist;
    dist_daglist.num_clients = nC;
    dist_daglist.nV = nV;
    dist_daglist.nE = nE;
    
    if (!random_inputs) {
        std::cout << "============================\n" << std::endl;
        std::cout << "Generating random inputs " << std::endl;
        std::cout << "Generating scale-free graph with nV=" << nV << ", nE=" << nE << std::endl;
        auto edges = generate_scale_free(nV, nE);
        std::cout << "Generated " << edges.size() << " edges" << std::endl;
        
        std::cout << "Building daglist..." << std::endl;
        auto daglist = build_daglist(nV, edges);
        std::cout << "Built daglist with " << daglist.size() << " entries" << std::endl;
        
        // Distribute daglist across clients
        std::cout << "Distributing daglist across " << nC << " clients..." << std::endl;
        dist_daglist = distribute_daglist(daglist, nC);
        
        // Generate random data updates (change 20% of entries)
        Ring num_changes = static_cast<Ring>(vec_size * 0.05);
        std::cout << "Generating random data updates for " << num_changes << " entries..." << std::endl;
        dist_daglist = generate_random_data_updates(dist_daglist, num_changes, seed);
    } else {
        std::cout << "============================\n" << std::endl;
        std::cout << "Using random inputs" << std::endl;
        
        // Compute sizes for distribution
        dist_daglist.VSizes.resize(nC);
        dist_daglist.ESizes.resize(nC);
        
        int base_verts = nV / nC;
        int base_edges = nE / nC;
        int extra_verts = nV % nC;
        int extra_edges = nE % nC;
        
        for (int i = 0; i < nC; ++i) {
            dist_daglist.VSizes[i] = base_verts + (i < extra_verts ? 1 : 0);
            dist_daglist.ESizes[i] = base_edges + (i < extra_edges ? 1 : 0);
        }
        
        std::cout << "Computed distribution: " << nC << " clients, " << nV << " vertices, " << nE << " edges" << std::endl;
    }

    StatsPoint start(*network);
    network->sync();

    auto circ = generateCircuit(nP, pid, dist_daglist).orderGatesByLevel();
    network->sync();

    std::cout << "--- Circuit ---" << std::endl;
    std::cout << circ << std::endl;
    
    std::unordered_map<common::utils::wire_t, int> input_pid_map;
    for (const auto& g : circ.gates_by_level[0]) {
        if (g->type == common::utils::GateType::kInp) {
            input_pid_map[g->out] = 1;  // All inputs belong to party 1
        }
    }

    std::cout << "Starting preprocessing" << std::endl;
    StatsPoint preproc_start(*network);
    int latency_us = static_cast<int>(latency * 1000);  // Convert ms to microseconds
    OfflineEvaluator off_eval(nP, pid, network, circ, threads, seed, latency_us, use_pking);
    auto preproc = off_eval.run(input_pid_map);
    std::cout << "Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint preproc_end(*network);

    std::cout << "Setting inputs" << std::endl;
    // Why latency_us and use_pking?
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
    if (random_inputs) {
        // Use random inputs for benchmarking
        std::cout << "Using random inputs for party " << pid << std::endl;
        eval.setRandomInputs();
    } else {
        std::unordered_map<common::utils::wire_t, Ring> inputs;
        
        // Set inputs directly by wire ID, not through sorted array
        if (pid == 1) {

            // Print distribution info
            std::cout << "\n=== Daglist Distribution ===" << std::endl;
            for (int i = 0; i < nC; ++i) {
                std::cout << "Client " << i << ": " << dist_daglist.VSizes[i] << " vertices, "
                        << dist_daglist.ESizes[i] << " edges" << std::endl;
            }
            std::cout << "============================\n" << std::endl;
        
            // Set inputs in the SAME ORDER as circuit creates wires:
            // ALL vertices for ALL clients first, then ALL edges for ALL clients
            size_t wire_id = 0;
            
            // First pass: all vertices for all clients
            for (int c = 0; c < nC; ++c) {
                for (size_t i = 0; i < dist_daglist.VSizes[c]; ++i) {
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].src;
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].dst;
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].isV;
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].data;
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].sigs;
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].sigv;
                    inputs[wire_id++] = dist_daglist.VertexLists[c][i].sigd;
                    inputs[wire_id++] = dist_daglist.isChangeV[c][i];
                    inputs[wire_id++] = dist_daglist.ChangeV[c][i];
                }
            }

            // Second pass: all edges for all clients
            for (int c = 0; c < nC; ++c) {
                for (size_t i = 0; i < dist_daglist.ESizes[c]; ++i) {
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].src;
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].dst;
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].isV;
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].data;
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].sigs;
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].sigv;
                    inputs[wire_id++] = dist_daglist.EdgeLists[c][i].sigd;
                    inputs[wire_id++] = dist_daglist.isChangeE[c][i];
                    inputs[wire_id++] = dist_daglist.ChangeE[c][i];
                }
            }
        }

        std::cout << "Total inputs set by party " << pid << ": " << inputs.size() << std::endl;
        
        if (pid == 1) {
            std::cout << "Party 1 setting " << inputs.size() << " actual input values" << std::endl;
        } else {
            std::cout << "Party " << pid << " setting " << inputs.size() << " empty inputs (participant in MPC)" << std::endl;
        }
        
        eval.setInputs(inputs);
    }
    network->sync();
    
    std::cout << "Starting online evaluation" << std::endl;
    StatsPoint online_start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
        network->sync();
        network->flush();
    }
    network->sync();
    StatsPoint online_end(*network);
    StatsPoint end(*network);
    std::cout << "Online evaluation complete" << std::endl;

    std::cout << "Getting outputs..." << std::endl;
    network->flush();
    auto outputs = eval.getOutputs();
    network->sync();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print formatted inputs and outputs (skip if using random inputs)
    if (!random_inputs) {
        if (pid == 1) {
            printDaglistInfo(dist_daglist, "INPUT: Graph Before Data Updates");
        }
        
        // Print formatted outputs
        if (pid == 1 && outputs.size() > 0) {
            printOutputs(outputs, dist_daglist);
        }

            // Update the distributed daglist with output values
    // Outputs match setAsOutput order: for each client, all vertices then all edges
    // For each entry, 7 fields appear: src, dst, isV, data, sigv, sigs, sigd (note: sigv before sigs!)
    size_t output_idx = 0;
        
        for (int c = 0; c < nC; ++c) {
            // Parse vertex outputs for this client - fields are interleaved per vertex
            for (size_t i = 0; i < dist_daglist.VSizes[c]; ++i) {
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].src = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].dst = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].isV = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].data = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].sigv = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].sigs = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.VertexLists[c][i].sigd = outputs[output_idx++];
            }
            
            // Parse edge outputs for this client - fields are interleaved per edge
            for (size_t i = 0; i < dist_daglist.ESizes[c]; ++i) {
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].src = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].dst = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].isV = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].data = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].sigv = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].sigs = outputs[output_idx++];
                if (output_idx < outputs.size()) dist_daglist.EdgeLists[c][i].sigd = outputs[output_idx++];
            }
        }
    }



    

    auto preproc_rbench = preproc_end - preproc_start;
    auto online_rbench = online_end - online_start;
    auto total_rbench = end - start;
    output_data["benchmarks"].push_back(preproc_rbench);
    output_data["benchmarks"].push_back(online_rbench);
    output_data["benchmarks"].push_back(total_rbench);

    size_t pre_bytes_sent = 0;
    for (const auto& val : preproc_rbench["communication"]) {
        pre_bytes_sent += val.get<int64_t>();
    }
    size_t online_bytes_sent = 0;
    for (const auto& val : online_rbench["communication"]) {
        online_bytes_sent += val.get<int64_t>();
    }
    size_t total_bytes_sent = 0;
    for (const auto& val : total_rbench["communication"]) {
        total_bytes_sent += val.get<int64_t>();
    }

    std::cout << "preproc time: " << preproc_rbench["time"] << " ms" << std::endl;
    std::cout << "preproc sent: " << pre_bytes_sent << " bytes" << std::endl;
    std::cout << "online time: " << online_rbench["time"] << " ms" << std::endl;
    std::cout << "online sent: " << online_bytes_sent << " bytes" << std::endl;
    std::cout << "total time: " << total_rbench["time"] << " ms" << std::endl;
    std::cout << "total sent: " << total_bytes_sent << " bytes" << std::endl;
    std::cout << std::endl;

    output_data["stats"] = {{"peak_virtual_memory", peakVirtualMemory()},
                            {"peak_resident_set_size", peakResidentSetSize()}};

    std::cout << "--- Statistics ---" << std::endl;
    for (const auto& [key, value] : output_data["stats"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    if (save_output) {
        saveJson(output_data, save_file);
    }
}

// clang-format off
bpo::options_description programOptions() {
    bpo::options_description desc("Following options are supported by config file too.");
    desc.add_options()
        ("num-parties,n", bpo::value<int>()->required(), "Number of parties.")
        ("num-clients", bpo::value<int>()->default_value(2), "Number of parties.")
        ("num-vert", bpo::value<size_t>()->default_value(1000), "Number of vertices in the graph.")
        ("num-edge", bpo::value<size_t>()->default_value(4000), "Number of edges in the graph.")
        ("num-payloads", bpo::value<size_t>()->default_value(1), "Number of payload vectors.")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
        ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
        ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
        ("localhost", bpo::bool_switch(), "All parties are on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
        ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction (true) or direct reconstruction (false).")
        ("random-inputs", bpo::value<bool>()->default_value(false), "Use random inputs for benchmarking.");
  return desc;
}
// clang-format on

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark secure graph data update circuit.");
    cmdline.add(prog_opts);
    cmdline.add_options()(
      "config,c", bpo::value<std::string>(),
      "configuration file for easy specification of cmd line arguments")(
      "help,h", "produce help message");
    bpo::variables_map opts;
    bpo::store(bpo::command_line_parser(argc, argv).options(cmdline).run(), opts);
    if (opts.count("help") != 0) {
        std::cout << cmdline << std::endl;
        return 0;
    }
    if (opts.count("config") > 0) {
        std::string cpath(opts["config"].as<std::string>());
        std::ifstream fin(cpath.c_str());
        if (fin.fail()) {
            std::cerr << "Could not open configuration file at " << cpath << std::endl;
            return 1;
        }
        bpo::store(bpo::parse_config_file(fin, prog_opts), opts);
    }
    try {
        bpo::notify(opts);
        if (!opts["localhost"].as<bool>() && (opts.count("net-config") == 0)) {
            throw std::runtime_error("Expected one of 'localhost' or 'net-config'");
        }
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    try {
        benchmark(opts);
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << "\nFatal error" << std::endl;
        return 1;
    }
    return 0;
}

// usage: ./../run.sh data_change --num-parties 2 --num-clients 2 --num-vert 1000 --num-edge 4000