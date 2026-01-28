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
    int total_del_vertices = 0, total_del_edges = 0;
    
    for (int i = 0; i < nC; ++i) {
        total_vertices += dist_daglist.VSizes[i];
        total_edges += dist_daglist.ESizes[i];
        
        // Count deletions
        for (size_t j = 0; j < dist_daglist.VSizes[i]; ++j) {
            if (dist_daglist.isDelV[i][j] == Ring(1)) total_del_vertices++;
        }
        for (size_t j = 0; j < dist_daglist.ESizes[i]; ++j) {
            if (dist_daglist.isDelE[i][j] == Ring(1)) total_del_edges++;
        }
    }
    
    std::cout << "Total Vertices: " << total_vertices << " (Marked for deletion: " << total_del_vertices << ")" << std::endl;
    std::cout << "Total Edges: " << total_edges << " (Marked for deletion: " << total_del_edges << ")" << std::endl;
    std::cout << "Total Entries: " << (total_vertices + total_edges) << std::endl;
    std::cout << "Expected after deletion: " << (total_vertices + total_edges - total_del_vertices - total_del_edges) << std::endl;
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
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd | Del" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
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
                      << std::setw(3) << dist_daglist.isDelV[c][i] << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    
    // Print sample edges
    std::cout << "Sample Edges (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd | Del" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
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
                      << std::setw(3) << dist_daglist.isDelE[c][i] << std::endl;
            count++;
        }
    }
    std::cout << std::string(60, '=') << std::endl;
}

void printOutputs(const std::vector<Ring>& outputs, size_t vec_size) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "OUTPUT: Compacted Graph After Deletion" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    if (outputs.size() < 1) {
        std::cout << "No outputs received!" << std::endl;
        return;
    }
    
    Ring num_remaining = outputs[0];
    std::cout << "Number of entries remaining: " << num_remaining << std::endl;
    std::cout << "(out of " << vec_size << " original entries)" << std::endl;
    std::cout << "Entries deleted: " << (vec_size - num_remaining) << std::endl;
    std::cout << std::endl;
    
    size_t expected_outputs = 1 + vec_size * 7; // num_remaining + 7 fields per entry
    if (outputs.size() < expected_outputs) {
        std::cout << "Warning: Received " << outputs.size() << " outputs, expected " << expected_outputs << std::endl;
        return;
    }
    
    std::cout << "Compacted Entries (up to first 20):" << std::endl;
    std::cout << "  ID | Src | Dst | Data  | isV | sigv | sigs | sigd" << std::endl;
    std::cout << "  " << std::string(60, '-') << std::endl;
    
    size_t limit = std::min(static_cast<size_t>(num_remaining), size_t(20));
    for (size_t i = 0; i < limit; ++i) {
        size_t base = 1 + i * 7;
        if (base + 6 < outputs.size()) {
            std::cout << "  " << std::setw(2) << i << " | "
                      << std::setw(3) << outputs[base + 0] << " | "
                      << std::setw(3) << outputs[base + 1] << " | "
                      << std::setw(4) << outputs[base + 2] << " | "
                      << std::setw(4) << outputs[base + 3] << " | "
                      << std::setw(4) << outputs[base + 4] << " | "
                      << std::setw(4) << outputs[base + 5] << " | "
                      << std::setw(4) << outputs[base + 6] << std::endl;
        }
    }
    
    if (num_remaining > 20) {
        std::cout << "  ... (" << (num_remaining - 20) << " more entries)" << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
}

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, DistributedDaglist dist_daglist) {

    int nC = dist_daglist.num_clients;
    int nV = dist_daglist.nV;
    int nE = dist_daglist.nE;
    auto VSizes = dist_daglist.VSizes;
    auto ESizes = dist_daglist.ESizes;
    size_t vec_size = nV + nE;

    std::cout << "Generating circuit" << std::endl;

    common::utils::Circuit<Ring> circ; 

    // Reserve capacity for wire vectors upfront
    size_t total_vertex_wires = 0, total_edge_wires = 0;
    for (int i = 0; i < nC; ++i) {
        total_vertex_wires += VSizes[i];
        total_edge_wires += ESizes[i];
    }

    // Initialize all daglist field values
    std::vector<std::vector<wire_t>> vertex_src_values(nC);
    std::vector<std::vector<wire_t>> vertex_dst_values(nC);
    std::vector<std::vector<wire_t>> vertex_isV_values(nC);
    std::vector<std::vector<wire_t>> vertex_data_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigd_values(nC);
    std::vector<std::vector<wire_t>> vertex_deleted(nC);

    std::vector<std::vector<wire_t>> edge_src_values(nC);
    std::vector<std::vector<wire_t>> edge_dst_values(nC);
    std::vector<std::vector<wire_t>> edge_isV_values(nC);
    std::vector<std::vector<wire_t>> edge_data_values(nC);
    std::vector<std::vector<wire_t>> edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> edge_sigv_values(nC);
    std::vector<std::vector<wire_t>> edge_sigd_values(nC);
    std::vector<std::vector<wire_t>> edge_deleted(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_vertex_src_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_dst_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_isV_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_data_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigs_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigv_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigd_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_deleted(VSizes[i]);
        
        std::vector<wire_t> subg_edge_src_values(ESizes[i]);
        std::vector<wire_t> subg_edge_dst_values(ESizes[i]);
        std::vector<wire_t> subg_edge_isV_values(ESizes[i]);
        std::vector<wire_t> subg_edge_data_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigs_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigv_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigd_values(ESizes[i]);
        std::vector<wire_t> subg_edge_deleted(ESizes[i]);

        for (int j = 0; j < VSizes[i]; ++j){
            subg_vertex_src_values[j] = circ.newInputWire();
            subg_vertex_dst_values[j] = circ.newInputWire();
            subg_vertex_isV_values[j] = circ.newInputWire();
            subg_vertex_data_values[j] = circ.newInputWire();
            subg_vertex_sigs_values[j] = circ.newInputWire();
            subg_vertex_sigv_values[j] = circ.newInputWire();
            subg_vertex_sigd_values[j] = circ.newInputWire();
            subg_vertex_deleted[j] = circ.newInputWire();
        }

        for (int j = 0; j < ESizes[i]; ++j){
            subg_edge_src_values[j] = circ.newInputWire();
            subg_edge_dst_values[j] = circ.newInputWire();
            subg_edge_isV_values[j] = circ.newInputWire();
            subg_edge_data_values[j] = circ.newInputWire();
            subg_edge_sigs_values[j] = circ.newInputWire();
            subg_edge_sigv_values[j] = circ.newInputWire();
            subg_edge_sigd_values[j] = circ.newInputWire();
            subg_edge_deleted[j] = circ.newInputWire();
        }

        vertex_src_values[i] = std::move(subg_vertex_src_values);
        vertex_dst_values[i] = std::move(subg_vertex_dst_values);
        vertex_isV_values[i] = std::move(subg_vertex_isV_values);
        vertex_data_values[i] = std::move(subg_vertex_data_values);
        vertex_sigs_values[i] = std::move(subg_vertex_sigs_values);
        vertex_sigv_values[i] = std::move(subg_vertex_sigv_values);
        vertex_sigd_values[i] = std::move(subg_vertex_sigd_values);
        vertex_deleted[i] = std::move(subg_vertex_deleted);

        edge_src_values[i] = std::move(subg_edge_src_values);
        edge_dst_values[i] = std::move(subg_edge_dst_values);
        edge_isV_values[i] = std::move(subg_edge_isV_values);
        edge_data_values[i] = std::move(subg_edge_data_values);
        edge_sigs_values[i] = std::move(subg_edge_sigs_values);
        edge_sigv_values[i] = std::move(subg_edge_sigv_values);
        edge_sigd_values[i] = std::move(subg_edge_sigd_values);
        edge_deleted[i] = std::move(subg_edge_deleted);
    
    }

    // Generate permutation for shuffle
    // Here we just pass identity permutations
    std::vector<int> base_perm(nV + nE);
    for (size_t i = 0; i < nV + nE; ++i) {
        base_perm[i] = static_cast<int>(i);
    }
    std::vector<std::vector<int>> permutation;
    for (int i = 0; i < nP; ++i) {
        permutation.push_back(base_perm);
    }

    // Generate flat daglist
    auto zerowire = circ.addGate(kSub, vertex_src_values[0][0], vertex_src_values[0][0]);
    std::vector<wire_t> src;
    std::vector<wire_t> dst;
    std::vector<wire_t> isV;
    std::vector<wire_t> data;
    std::vector<wire_t> sigs;
    std::vector<wire_t> sigv;
    std::vector<wire_t> sigd;
    std::vector<wire_t> del_v;
    std::vector<wire_t> del_e;
    
    // Reserve capacity to avoid reallocations
    src.reserve(nV + nE);
    dst.reserve(nV + nE);
    isV.reserve(nV + nE);
    data.reserve(nV + nE);
    sigs.reserve(nV + nE);
    sigv.reserve(nV + nE);
    sigd.reserve(nV + nE);
    del_v.reserve(nV + nE);
    del_e.reserve(nV + nE);
    
    // First, push all vertices from all clients
    for (int i = 0; i < nC; ++i) {
        for (int j = 0; j < VSizes[i]; ++j) {
            src.push_back(vertex_src_values[i][j]);
            dst.push_back(vertex_dst_values[i][j]);
            isV.push_back(vertex_isV_values[i][j]);
            data.push_back(vertex_data_values[i][j]);
            sigs.push_back(vertex_sigs_values[i][j]);
            sigv.push_back(vertex_sigv_values[i][j]);
            sigd.push_back(vertex_sigd_values[i][j]);
            del_v.push_back(vertex_deleted[i][j]);
            del_e.push_back(zerowire);
        }
    }
    
    // Then, push all edges from all clients
    for (int i = 0; i < nC; ++i) {
        for (int j = 0; j < ESizes[i]; ++j) {
            src.push_back(edge_src_values[i][j]);
            dst.push_back(edge_dst_values[i][j]);
            isV.push_back(edge_isV_values[i][j]);
            data.push_back(edge_data_values[i][j]);
            sigs.push_back(edge_sigs_values[i][j]);
            sigv.push_back(edge_sigv_values[i][j]);
            sigd.push_back(edge_sigd_values[i][j]);
            del_v.push_back(zerowire);
            del_e.push_back(edge_deleted[i][j]);
        }
    }
    
    // Clear intermediate vectors to free memory
    vertex_src_values.clear();
    vertex_src_values.shrink_to_fit();
    vertex_dst_values.clear();
    vertex_dst_values.shrink_to_fit();
    vertex_isV_values.clear();
    vertex_isV_values.shrink_to_fit();
    vertex_data_values.clear();
    vertex_data_values.shrink_to_fit();
    vertex_sigs_values.clear();
    vertex_sigs_values.shrink_to_fit();
    vertex_sigv_values.clear();
    vertex_sigv_values.shrink_to_fit();
    vertex_sigd_values.clear();
    vertex_sigd_values.shrink_to_fit();
    vertex_deleted.clear();
    vertex_deleted.shrink_to_fit();
    
    edge_src_values.clear();
    edge_src_values.shrink_to_fit();
    edge_dst_values.clear();
    edge_dst_values.shrink_to_fit();
    edge_isV_values.clear();
    edge_isV_values.shrink_to_fit();
    edge_data_values.clear();
    edge_data_values.shrink_to_fit();
    edge_sigs_values.clear();
    edge_sigs_values.shrink_to_fit();
    edge_sigv_values.clear();
    edge_sigv_values.shrink_to_fit();
    edge_sigd_values.clear();
    edge_sigd_values.shrink_to_fit();
    edge_deleted.clear();
    edge_deleted.shrink_to_fit();

    // Propagate del tag to outgoing edges and reorder them back to vertex order
    auto del_S = circ.addSubCircPropagate(sigs, del_v, nV, permutation);
    // reorder del_S to vertex order
    auto sigs_to_sigv = circ.addSubCircPermList(sigs, {sigv}, permutation)[0];
    del_S = circ.addSubCircPermList(sigs_to_sigv, {del_S}, permutation)[0];

    // Propagate del tag to incoming edges and reorder them back to vertex order
    auto del_D = circ.addSubCircPropagate(sigd, del_v, nV, permutation, true);
    auto sigd_to_sigv = circ.addSubCircPermList(sigd, {sigv}, permutation)[0];
    del_D = circ.addSubCircPermList(sigd_to_sigv, {del_D}, permutation)[0];

    // Combine del tags
    std::vector<wire_t> del_final(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        auto temp = circ.addGate(common::utils::GateType::kAdd, del_S[i], del_D[i]);
        
        temp = circ.addGate(common::utils::GateType::kAdd, del_e[i], temp); 
        
        temp = circ.addGate(common::utils::GateType::kEqz, temp);
        del_final[i] = temp;
        temp = circ.addConstOpGate(common::utils::GateType::kConstMul, temp, Ring(-1));
        del_final[i] = circ.addConstOpGate(common::utils::GateType::kConstAdd, temp, Ring(1));
    }

    // // Update sigv
    std::vector<wire_t> updated_sigv;
    updated_sigv.reserve(vec_size);
    updated_sigv.push_back(sigv[0]);
    wire_t prefix_sum = del_final[0];
    for (size_t i = 1; i < vec_size; ++i) {
        updated_sigv.push_back(circ.addGate(common::utils::GateType::kSub, sigv[i], prefix_sum));
        prefix_sum = circ.addGate(common::utils::GateType::kAdd, prefix_sum, del_final[i]);
    }

    // // Combine graph vectors into single payload
    std::vector<std::vector<wire_t>> payload1;
    payload1.reserve(2);
    payload1.push_back(sigs);              
    payload1.push_back(del_final);

    std::vector<std::vector<wire_t>> payload2;
    payload2.reserve(2);
    payload2.push_back(sigd);              
    payload2.push_back(del_final);

    // // Reorder to source order and destination order
    auto payload_s = circ.addSubCircPermList(payload1[0], payload1, permutation);
    auto payload_d = circ.addSubCircPermList(payload2[0], payload2, permutation);

    // // Update sigs
    std::vector<wire_t> updated_sigs;
    updated_sigs.reserve(vec_size);
    auto& sigs_old = payload_s[0];
    auto& del_s = payload_s[1];
    wire_t prefix_sum_s;
    updated_sigs.push_back(sigs_old[0]);
    prefix_sum_s = del_s[0];
    
    for (size_t i = 1; i < vec_size; ++i) {
        updated_sigs.push_back(circ.addGate(common::utils::GateType::kSub, sigs_old[i], prefix_sum_s));
        prefix_sum_s = circ.addGate(common::utils::GateType::kAdd, prefix_sum_s, del_s[i]);
    }
    payload_s[0] = std::move(updated_sigs); 
    payload_s = circ.addSubCircPermList(sigs_to_sigv, payload_s, permutation);
    updated_sigs = std::move(payload_s[0]);

    // // Update sigd
    std::vector<wire_t> updated_sigd;
    updated_sigd.reserve(vec_size);
    auto& sigd_old = payload_d[0];
    auto& del_d = payload_d[1];
    wire_t prefix_sum_d;
    updated_sigd.push_back(sigd_old[0]);
    prefix_sum_d = del_d[0];
    for (size_t i = 1; i < vec_size; ++i) {
        updated_sigd.push_back(circ.addGate(common::utils::GateType::kSub, sigd_old[i], prefix_sum_d));
        prefix_sum_d = circ.addGate(common::utils::GateType::kAdd, prefix_sum_d, del_d[i]);
    }
    payload_d[0] = std::move(updated_sigd); 
    payload_d = circ.addSubCircPermList(sigd_to_sigv, payload_d, permutation);
    updated_sigd = std::move(payload_d[0]);
    
    // Clear del vectors to free memory
    del_v.clear();
    del_v.shrink_to_fit();
    del_e.clear();
    del_e.shrink_to_fit();

    payload1.resize(7);
    payload1[0] = std::move(src);
    payload1[1] = std::move(dst);
    payload1[2] = std::move(data);
    payload1[3] = std::move(isV);
    payload1[4] = std::move(updated_sigv);
    payload1[5] = std::move(updated_sigs);
    payload1[6] = std::move(updated_sigd);

    auto [num_remaining, payload1_deleted] = circ.addDeleteWiresGate(del_final, payload1, permutation);

    // Get compacted vectors (these have dynamic size = num_remaining)
    auto src_compacted = payload1_deleted[0];
    auto dst_compacted = payload1_deleted[1];
    auto data_compacted = payload1_deleted[2];
    auto isV_compacted = payload1_deleted[3];
    auto sigv_compacted = payload1_deleted[4];
    auto sigs_compacted = payload1_deleted[5];
    auto sigd_compacted = payload1_deleted[6];
    
    
    // Set outputs
    circ.setAsOutput(num_remaining);
    // Output all compacted entries (vertices + edges after deletion)
    // The compacted vectors maintain vertex-order: vertices first, then edges
    for (size_t i = 0; i < vec_size; ++i) {
        circ.setAsOutput(src_compacted[i]);
        circ.setAsOutput(dst_compacted[i]);
        circ.setAsOutput(data_compacted[i]);
        circ.setAsOutput(isV_compacted[i]);
        circ.setAsOutput(sigv_compacted[i]);
        circ.setAsOutput(sigs_compacted[i]);
        circ.setAsOutput(sigd_compacted[i]);
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
        std::cout << "Generating scale-free graph with nV=" << nV << ", nE=" << nE << " (seed=" << seed << ")" << std::endl;
        auto edges = generate_scale_free(nV, nE, seed);
        std::cout << "Generated " << edges.size() << " edges" << std::endl;
        
        std::cout << "Building daglist..." << std::endl;
        auto daglist = build_daglist(nV, edges);
        std::cout << "Built daglist with " << daglist.size() << " entries" << std::endl;
        
        // Distribute daglist across clients
        std::cout << "Distributing daglist across " << nC << " clients..." << std::endl;
        dist_daglist = distribute_daglist(daglist, nC);
        
        // Generate random deletion tags (delete 5% of entries)
        Ring num_deletes = static_cast<Ring>(vec_size * 0.05);
        std::cout << "Generating random deletion tags for " << num_deletes << " entries..." << std::endl;
        dist_daglist = generate_random_entry_deletes(dist_daglist, num_deletes, seed);

        // Print input daglist information
        if (pid == 1) {
            printDaglistInfo(dist_daglist, "INPUT: Graph Before Deletion");
        }
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
            input_pid_map[g->out] = 1;
        }
    }

    std::cout << "Starting preprocessing" << std::endl;
    StatsPoint preproc_start(*network);
    int latency_us = static_cast<int>(latency * 1000);
    OfflineEvaluator off_eval(nP, pid, network, circ, threads, seed, latency_us, use_pking);
    auto preproc = off_eval.run(input_pid_map);
    std::cout << "Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint preproc_end(*network);

    std::cout << "Setting inputs" << std::endl;
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
    if (random_inputs) {
        // Use random inputs for benchmarking
        std::cout << "Using random inputs for party " << pid << std::endl;
        eval.setRandomInputs();
    } else {
        std::unordered_map<common::utils::wire_t, Ring> inputs;
        
        // Collect all input wires owned by this party
        std::vector<common::utils::wire_t> input_wires;
        for (const auto& [wire, owner] : input_pid_map) {
            if (owner == static_cast<int>(pid)) {
                input_wires.push_back(wire);
            }
        }
        
        // Sort to ensure consistent ordering
        std::sort(input_wires.begin(), input_wires.end());
        
        std::cout << "Setting inputs for party " << pid << std::endl;
        
        // Only party 1 sets inputs
        std::vector<Ring> graph_input_values;

        if (pid == 1) {
            std::vector<Ring> all_input_values;
            
        // Collect all vertex and edge fields for all clients
        for (int c = 0; c < nC; ++c) {
            for (size_t i = 0; i < dist_daglist.VSizes[c]; ++i) {
                all_input_values.push_back(dist_daglist.VertexLists[c][i].src);
                all_input_values.push_back(dist_daglist.VertexLists[c][i].dst);
                all_input_values.push_back(dist_daglist.VertexLists[c][i].isV);
                all_input_values.push_back(dist_daglist.VertexLists[c][i].data);
                all_input_values.push_back(dist_daglist.VertexLists[c][i].sigs);
                all_input_values.push_back(dist_daglist.VertexLists[c][i].sigv);
                all_input_values.push_back(dist_daglist.VertexLists[c][i].sigd);
                all_input_values.push_back(dist_daglist.isDelV[c][i]);
            }

            for (size_t i = 0; i < dist_daglist.ESizes[c]; ++i) {
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].src);
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].dst);
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].isV);
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].data);
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigs);
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigv);
                all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigd);
                all_input_values.push_back(dist_daglist.isDelE[c][i]);
            }
        }
            
                    // Map collected values into circuit input wires (in order)
            size_t wire_idx = 0;
            for (size_t i = 0; i < all_input_values.size() && wire_idx < input_wires.size(); ++i) {
                inputs[input_wires[wire_idx++]] = all_input_values[i];
            }
            
            // Store for verification
            graph_input_values = all_input_values;
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
        // network->sync();
        // network->flush();
    }
    network->sync();
    StatsPoint end(*network);
    StatsPoint online_end(*network);
    std::cout << "Online evaluation complete" << std::endl;

    std::cout << "Getting outputs..." << std::endl;
    network->flush();
    auto outputs = eval.getOutputs();
    network->sync();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print formatted outputs (skip if using random inputs)
    if (!random_inputs && pid == 1 && outputs.size() > 0) {
        printOutputs(outputs, vec_size);
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
    bpo::options_description cmdline("Benchmark secure vertex/edge deletion circuit.");
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

// usage: ./../run.sh delete --num-parties 2 --num-clients 2 --num-vert 1000 --num-edge 4000