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
    int total_new_vertices = 0;
    
    for (int i = 0; i < nC; ++i) {
        total_vertices += dist_daglist.VSizes[i];
        total_edges += dist_daglist.ESizes[i];
        total_new_vertices += dist_daglist.InsertV[i].size();
    }
    
    std::cout << "Total Vertices: " << total_vertices << std::endl;
    std::cout << "Total Edges: " << total_edges << std::endl;
    std::cout << "New Vertices to Insert: " << total_new_vertices << std::endl;
    std::cout << "Total Entries: " << (total_vertices + total_edges) << std::endl;
    std::cout << std::endl;
    
    // Print distribution per client
    std::cout << "Distribution per Client:" << std::endl;
    for (int i = 0; i < nC; ++i) {
        std::cout << "  Client " << i << ": " << dist_daglist.VSizes[i] << " vertices, "
                  << dist_daglist.ESizes[i] << " edges, "
                  << dist_daglist.InsertV[i].size() << " new vertices" << std::endl;
    }
    std::cout << std::endl;
    
    // Print sample vertices
    std::cout << "Sample Existing Vertices (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd" << std::endl;
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
                      << std::setw(4) << v.sigd << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    
    // Print new vertices to insert
    std::cout << "New Vertices to Insert (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV" << std::endl;
    std::cout << "  " << std::string(40, '-') << std::endl;
    count = 0;
    for (int c = 0; c < nC && count < 10; ++c) {
        for (size_t i = 0; i < dist_daglist.InsertV[c].size() && count < 10; ++i) {
            const auto& v = dist_daglist.InsertV[c][i];
            std::cout << "  " << std::setw(2) << count << " | "
                      << std::setw(3) << v.src << " | "
                      << std::setw(3) << v.dst << " | "
                      << std::setw(4) << v.data << " | "
                      << std::setw(3) << v.isV << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    
    // Print sample edges
    std::cout << "Sample Edges (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd" << std::endl;
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
                      << std::setw(4) << e.sigd << std::endl;
            count++;
        }
    }
    std::cout << std::string(60, '=') << std::endl;
}

void printOutputs(const std::vector<Ring>& outputs, const DistributedDaglist& dist_daglist) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "OUTPUT: Graph After Vertex Insertion" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int nC = dist_daglist.num_clients;
    size_t expected_outputs = 0;
    for (int c = 0; c < nC; ++c) {
        expected_outputs += dist_daglist.VSizes[c] * 7; // existing vertices
        expected_outputs += dist_daglist.InsertV[c].size() * 7; // new vertices
        expected_outputs += dist_daglist.ESizes[c] * 7; // edges
    }
    
    if (outputs.size() < expected_outputs) {
        std::cout << "Warning: Received " << outputs.size() << " outputs, expected " << expected_outputs << std::endl;
        return;
    }
    
    std::cout << "Updated Vertices (up to first 10, including new vertices):" << std::endl;
    std::cout << "  ID | Src | Dst | isV | Data | sigs | sigv | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    int display_count = 0;
    size_t output_idx = 0;
    
    // Print per-client: existing vertices followed by new vertices
    for (int c = 0; c < nC && display_count < 10; ++c) {
        // Existing vertices for client c
        size_t limit = std::min(static_cast<size_t>(dist_daglist.VSizes[c]), size_t(10 - display_count));
        for (size_t i = 0; i < limit && output_idx + 6 < outputs.size(); ++i) {
            size_t idx = output_idx + i * 7;
            std::cout << "  " << std::setw(2) << display_count << " | "
                      << std::setw(3) << outputs[idx + 0] << " | "
                      << std::setw(3) << outputs[idx + 1] << " | "
                      << std::setw(3) << outputs[idx + 2] << " | "
                      << std::setw(4) << outputs[idx + 3] << " | "
                      << std::setw(4) << outputs[idx + 4] << " | "
                      << std::setw(4) << outputs[idx + 5] << " | "
                      << std::setw(4) << outputs[idx + 6] << std::endl;
            display_count++;
        }
        output_idx += dist_daglist.VSizes[c] * 7;

        // New vertices for client c
        if (display_count < 10) {
            limit = std::min(dist_daglist.InsertV[c].size(), size_t(10 - display_count));
            for (size_t i = 0; i < limit && output_idx + 6 < outputs.size(); ++i) {
                size_t idx = output_idx + i * 7;
                std::cout << "  " << std::setw(2) << display_count << " | "
                          << std::setw(3) << outputs[idx + 0] << " | "
                          << std::setw(3) << outputs[idx + 1] << " | "
                          << std::setw(3) << outputs[idx + 2] << " | "
                          << std::setw(4) << outputs[idx + 3] << " | "
                          << std::setw(4) << outputs[idx + 4] << " | "
                          << std::setw(4) << outputs[idx + 5] << " | "
                          << std::setw(4) << outputs[idx + 6] << " *NEW*" << std::endl;
                display_count++;
            }
        }
        output_idx += dist_daglist.InsertV[c].size() * 7;
    }
    
    std::cout << "\nUpdated Edges (up to first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | isV | Data | sigs | sigv | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    // Edges start after all vertex outputs
    size_t edges_start = output_idx;
    std::vector<size_t> edge_prefix(nC + 1, 0);
    for (int c = 0; c < nC; ++c) {
        edge_prefix[c + 1] = edge_prefix[c] + dist_daglist.ESizes[c] * 7;
    }
    display_count = 0;
    for (int c = 0; c < nC && display_count < 10; ++c) {
        size_t base = edges_start + edge_prefix[c];
        size_t limit = std::min(static_cast<size_t>(dist_daglist.ESizes[c]), size_t(10 - display_count));
        for (size_t i = 0; i < limit && base + i * 7 + 6 < outputs.size(); ++i) {
            size_t idx = base + i * 7;
            std::cout << "  " << std::setw(2) << display_count << " | "
                      << std::setw(3) << outputs[idx + 0] << " | "
                      << std::setw(3) << outputs[idx + 1] << " | "
                      << std::setw(3) << outputs[idx + 2] << " | "
                      << std::setw(4) << outputs[idx + 3] << " | "
                      << std::setw(4) << outputs[idx + 4] << " | "
                      << std::setw(4) << outputs[idx + 5] << " | "
                      << std::setw(4) << outputs[idx + 6] << std::endl;
            display_count++;
        }
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

    // Initialize all daglist field values
    std::vector<std::vector<wire_t>> vertex_src_values(nC);
    std::vector<std::vector<wire_t>> vertex_dst_values(nC);
    std::vector<std::vector<wire_t>> vertex_isV_values(nC);
    std::vector<std::vector<wire_t>> vertex_data_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigd_values(nC);

    std::vector<std::vector<wire_t>> edge_src_values(nC);
    std::vector<std::vector<wire_t>> edge_dst_values(nC);
    std::vector<std::vector<wire_t>> edge_isV_values(nC);
    std::vector<std::vector<wire_t>> edge_data_values(nC);
    std::vector<std::vector<wire_t>> edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> edge_sigv_values(nC);
    std::vector<std::vector<wire_t>> edge_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_vertex_src_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_dst_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_isV_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_data_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigs_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigv_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigd_values(VSizes[i]);
        
        std::vector<wire_t> subg_edge_src_values(ESizes[i]);
        std::vector<wire_t> subg_edge_dst_values(ESizes[i]);
        std::vector<wire_t> subg_edge_isV_values(ESizes[i]);
        std::vector<wire_t> subg_edge_data_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigs_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigv_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigd_values(ESizes[i]);

        for (int j = 0; j < VSizes[i]; ++j){
            subg_vertex_src_values[j] = circ.newInputWire();
            subg_vertex_dst_values[j] = circ.newInputWire();
            subg_vertex_isV_values[j] = circ.newInputWire();
            subg_vertex_data_values[j] = circ.newInputWire();
            subg_vertex_sigs_values[j] = circ.newInputWire();
            subg_vertex_sigv_values[j] = circ.newInputWire();
            subg_vertex_sigd_values[j] = circ.newInputWire();
        }

        for (int j = 0; j < ESizes[i]; ++j){
            subg_edge_src_values[j] = circ.newInputWire();
            subg_edge_dst_values[j] = circ.newInputWire();
            subg_edge_isV_values[j] = circ.newInputWire();
            subg_edge_data_values[j] = circ.newInputWire();
            subg_edge_sigs_values[j] = circ.newInputWire();
            subg_edge_sigv_values[j] = circ.newInputWire();
            subg_edge_sigd_values[j] = circ.newInputWire();
        }

        vertex_src_values[i] = subg_vertex_src_values;
        vertex_dst_values[i] = subg_vertex_dst_values;
        vertex_isV_values[i] = subg_vertex_isV_values;
        vertex_data_values[i] = subg_vertex_data_values;
        vertex_sigs_values[i] = subg_vertex_sigs_values;
        vertex_sigv_values[i] = subg_vertex_sigv_values;
        vertex_sigd_values[i] = subg_vertex_sigd_values;

        edge_src_values[i] = subg_edge_src_values;
        edge_dst_values[i] = subg_edge_dst_values;
        edge_isV_values[i] = subg_edge_isV_values;
        edge_data_values[i] = subg_edge_data_values;
        edge_sigs_values[i] = subg_edge_sigs_values;
        edge_sigv_values[i] = subg_edge_sigv_values;
        edge_sigd_values[i] = subg_edge_sigd_values;
    
    }

    int add_nV = 0;
    std::vector<int> add_VSizes(nC, 0);
    for (int i = 0; i < nC; ++i) {
        add_VSizes[i] = dist_daglist.InsertV[i].size();
        add_nV += add_VSizes[i];
    }
  
    // Initialize vertices to be added.
    std::vector<std::vector<wire_t>> new_vertex_src_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_dst_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_isV_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_data_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_new_vertex_src_values(add_VSizes[i]);
        std::vector<wire_t> subg_new_vertex_dst_values(add_VSizes[i]);
        std::vector<wire_t> subg_new_vertex_isV_values(add_VSizes[i]);
        std::vector<wire_t> subg_new_vertex_data_values(add_VSizes[i]);
        
        for (int j = 0; j < add_VSizes[i]; ++j){
            subg_new_vertex_src_values[j] = circ.newInputWire();
            subg_new_vertex_dst_values[j] = circ.newInputWire();
            subg_new_vertex_isV_values[j] = circ.newInputWire();
            subg_new_vertex_data_values[j] = circ.newInputWire();
        }

        new_vertex_src_values[i] = subg_new_vertex_src_values;
        new_vertex_dst_values[i] = subg_new_vertex_dst_values;
        new_vertex_isV_values[i] = subg_new_vertex_isV_values;
        new_vertex_data_values[i] = subg_new_vertex_data_values;

    }

    // zero wire for assigning constant value to wire
    wire_t zero_wire;
    if (nV > 0) {
        zero_wire = circ.addGate(kSub, vertex_src_values[0][0], vertex_src_values[0][0]);
    } else if (add_nV > 0) {
        zero_wire = circ.addGate(kSub, new_vertex_src_values[0][0], new_vertex_src_values[0][0]);
    } else {
        // No vertices at all, create a dummy input wire and zero it
        auto dummy = circ.newInputWire();
        zero_wire = circ.addGate(kSub, dummy, dummy);
    }

    // Compute per-client delta values (prefix sums of new vertices)
    std::vector<size_t> delta(nC);
    delta[0] = 0;
    for (int i = 1; i < nC; ++i) {
        delta[i] = delta[i-1] + add_VSizes[i-1];
    }

    // Assign positions for existing vertices
    
    std::vector<std::vector<wire_t>> updated_vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> updated_vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> updated_vertex_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        updated_vertex_sigs_values[i].resize(VSizes[i]);
        updated_vertex_sigv_values[i].resize(VSizes[i]);
        updated_vertex_sigd_values[i].resize(VSizes[i]);

        for (int k = 0; k < VSizes[i]; ++k) {
            updated_vertex_sigs_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    vertex_sigs_values[i][k], delta[i]);
            updated_vertex_sigv_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    vertex_sigv_values[i][k], delta[i]);
            updated_vertex_sigd_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    vertex_sigd_values[i][k], delta[i]);
        }
    }
    
    // Assign positions to new vertices
    wire_t prev_sigv_wire = zero_wire;
    wire_t prev_sigd_wire = zero_wire;
    for (int i = 0; i < nC; ++i) {
        new_vertex_sigs_values[i].resize(add_VSizes[i]);
        new_vertex_sigv_values[i].resize(add_VSizes[i]);
        new_vertex_sigd_values[i].resize(add_VSizes[i]);

        if (VSizes[i] > 0) {
            prev_sigv_wire = updated_vertex_sigv_values[i][VSizes[i]-1];
            prev_sigd_wire = updated_vertex_sigd_values[i][VSizes[i]-1];
        }

        for (int k = 0; k < add_VSizes[i]; ++k) {
            
            new_vertex_sigv_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    prev_sigv_wire, Ring(k+1));
            new_vertex_sigd_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    prev_sigd_wire, Ring(k+1));
        
        }        
    }


    wire_t next_sigs_wire = zero_wire;
    next_sigs_wire = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, Ring(vec_size + add_nV));
    for (int i = nC-1; i >= 0; --i) {
        for (int k = add_VSizes[i] - 1; k >= 0; --k) {
                new_vertex_sigs_values[i][k] = 
                    circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                        next_sigs_wire, Ring(k*-1 -1));
        }

        if (VSizes[i] > 0) {
            next_sigs_wire = updated_vertex_sigs_values[i][0];
        }
    }

    // Assign positions of edges
    // Assign delta values to wires
    std::vector<wire_t> delta_wires(vec_size);
    int index = 0;
    for (size_t i = 0 ; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            delta_wires[index] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, delta[i]);
            index++;
        }
    }
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < ESizes[i]; ++j) {
            // Initialize wires to zero
            delta_wires[index] = zero_wire;
            index++;
        }
    }

    // Flatten position maps
    std::vector<wire_t> sigv(vec_size);
    std::vector<wire_t> sigd(vec_size);
    index = 0;
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            sigv[index] = vertex_sigv_values[i][j];
            sigd[index] = vertex_sigd_values[i][j];
            index++;
        }
    }
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < ESizes[i]; ++j) {
            sigv[index] = edge_sigv_values[i][j];
            sigd[index] = edge_sigd_values[i][j];
            index++;
        }
    }

    // Generate permutation for shuffle
    // Here we just pass identity permutations
    std::vector<int> base_perm(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        base_perm[i] = static_cast<int>(i);
    }
    std::vector<std::vector<int>> permutation;
    permutation.push_back(base_perm);
    if (pid == 0) {
        for (int p = 1; p < nP; ++p) {
            permutation.push_back(base_perm);
        }
    }

    // Propagate (only if there are existing vertices/edges)
    std::vector<wire_t> prop_delta, sigv_d, delta_v;
    if (vec_size > 0) {
        if (nE > 0) {
            prop_delta = circ.addSubCircPropagate(sigd, delta_wires, nV, permutation, true);

            // Reorder sigv to destination order
            sigv_d = circ.addSubCircPermList(sigd, {sigv}, permutation)[0];

            // Reorder back to vertex order
            delta_v = circ.addSubCircPermList(sigv_d, {prop_delta}, permutation)[0];
        } else {
            // No edges, so propagation is identity
            prop_delta = delta_wires;
            sigv_d = sigv;
            delta_v = delta_wires;
        }
    } else {
        // No existing graph elements, create empty vectors
        prop_delta.resize(0);
        sigv_d.resize(0);
        delta_v.resize(0);
    }

    // Update sigd; vector includes dummy values for first nV entries (corresponding to vertices)
    std::vector<wire_t> updated_edge_sigd_values(vec_size);
    for (size_t i = nV; i < vec_size; ++i) {
        updated_edge_sigd_values[i] = circ.addGate(common::utils::GateType::kAdd, sigd[i], delta_v[i]);
    }

    // Update sigv and sigs
    std::vector<std::vector<wire_t>> updated_edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> updated_edge_sigv_values(nC);

    for (int i = 0; i < nC; ++i) {
        updated_edge_sigs_values[i].resize(ESizes[i]);
        updated_edge_sigv_values[i].resize(ESizes[i]);

        for (int k = 0; k < ESizes[i]; ++k) {
            updated_edge_sigs_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    edge_sigs_values[i][k], delta[i]);
            updated_edge_sigv_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    edge_sigv_values[i][k], add_nV);
        }
    }

    // Set outputs
    index = 0;
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            circ.setAsOutput(vertex_src_values[i][j]);
            circ.setAsOutput(vertex_dst_values[i][j]);
            circ.setAsOutput(vertex_isV_values[i][j]);
            circ.setAsOutput(vertex_data_values[i][j]);
            circ.setAsOutput(updated_vertex_sigs_values[i][j]);
            circ.setAsOutput(updated_vertex_sigv_values[i][j]);
            circ.setAsOutput(updated_vertex_sigd_values[i][j]);
            index++;
        }
        for (size_t j = 0; j < add_VSizes[i]; ++j) {
            circ.setAsOutput(new_vertex_src_values[i][j]);
            circ.setAsOutput(new_vertex_dst_values[i][j]);
            circ.setAsOutput(new_vertex_isV_values[i][j]);
            circ.setAsOutput(new_vertex_data_values[i][j]);
            circ.setAsOutput(new_vertex_sigs_values[i][j]);
            circ.setAsOutput(new_vertex_sigv_values[i][j]);
            circ.setAsOutput(new_vertex_sigd_values[i][j]);
        }
    }
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < ESizes[i]; ++j) {
            circ.setAsOutput(edge_src_values[i][j]);
            circ.setAsOutput(edge_dst_values[i][j]);
            circ.setAsOutput(edge_isV_values[i][j]);
            circ.setAsOutput(edge_data_values[i][j]);
            circ.setAsOutput(updated_edge_sigs_values[i][j]);
            circ.setAsOutput(updated_edge_sigv_values[i][j]);
            circ.setAsOutput(updated_edge_sigd_values[index]);
            index++;
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
    auto num_inserts_ = opts["num-inserts"].as<size_t>();
    auto vec_size = num_vert + num_edge + num_inserts_;
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
                              {"num_edge", num_edge},
                              {"num_inserts", num_inserts_},
                              {"operation", "insert_vertices"},
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
        
        // Generate random vertices to insert (insert 20% of current vertices)
        Ring num_inserts = static_cast<Ring>(num_inserts_);
        std::cout << "Generating random vertices to insert: " << num_inserts << " vertices..." << std::endl;
        dist_daglist = generate_random_vertices_to_insert(dist_daglist, num_inserts, seed);
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
        
        // Allocate InsertV sizes
        dist_daglist.InsertV.resize(nC);
        Ring num_inserts = static_cast<Ring>(num_inserts_);
        int base_inserts = num_inserts / nC;
        int extra_inserts = num_inserts % nC;
        
        for (int i = 0; i < nC; ++i) {
            int inserts = base_inserts + (i < extra_inserts ? 1 : 0);
            dist_daglist.InsertV[i].resize(inserts);
        }
        
        std::cout << "Computed distribution: " << nC << " clients, " << nV << " vertices, " << nE << " edges, " << num_inserts << " inserts" << std::endl;
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
    
    if (random_inputs) {
        // Use random inputs for benchmarking
        std::cout << "Using random inputs for party " << pid << std::endl;
        eval.setRandomInputs();
    } else {
        // Only party 1 sets inputs
        std::vector<Ring> graph_input_values;

        if (pid == 1) {

            // Print distribution info
            std::cout << "\n=== Daglist Distribution ===" << std::endl;
            for (int i = 0; i < nC; ++i) {
                std::cout << "Client " << i << ": " << dist_daglist.VSizes[i] << " vertices, "
                        << dist_daglist.ESizes[i] << " edges, " << dist_daglist.InsertV[i].size() << " to insert" << std::endl;
            }
            std::cout << "============================\n" << std::endl;
        
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
                }

                for (size_t i = 0; i < dist_daglist.ESizes[c]; ++i) {
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].src);
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].dst);
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].isV);
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].data);
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigs);
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigv);
                    all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigd);
                }
            }
            
            // Collect new vertex fields for all clients
            for (int c = 0; c < nC; ++c) {
                for (size_t i = 0; i < dist_daglist.InsertV[c].size(); ++i) {
                    all_input_values.push_back(dist_daglist.InsertV[c][i].src);
                    all_input_values.push_back(dist_daglist.InsertV[c][i].dst);
                    all_input_values.push_back(dist_daglist.InsertV[c][i].isV);
                    all_input_values.push_back(dist_daglist.InsertV[c][i].data);
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
    StatsPoint online_end(*network);
    std::cout << "Online evaluation complete" << std::endl;
    StatsPoint end(*network);

    std::cout << "Getting outputs..." << std::endl;
    network->flush();
    auto outputs = eval.getOutputs();
    network->sync();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print formatted inputs and outputs (skip if using random inputs)
    if (!random_inputs) {
        if (pid == 1) {
            printDaglistInfo(dist_daglist, "INPUT: Graph Before Vertex Insertion");
        }
        
        // Print formatted outputs
        if (pid == 1 && outputs.size() > 0) {
            printOutputs(outputs, dist_daglist);
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

bpo::options_description programOptions() {
    bpo::options_description desc("Following options are supported by config file too.");
    desc.add_options()
        ("num-parties,n", bpo::value<int>()->required(), "Number of parties.")
        ("num-clients", bpo::value<int>()->default_value(2), "Number of parties.")
        ("num-vert", bpo::value<size_t>()->default_value(1000), "Number of vertices in the graph.")
        ("num-edge", bpo::value<size_t>()->default_value(4000), "Number of edges in the graph.")
        ("num-inserts", bpo::value<size_t>(), "Number of vertex insertions (default: num-vert * 0.05).")
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

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark secure vertex insertion circuit.");
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
        // Set default value for num-inserts if not provided
        if (opts.count("num-inserts") == 0) {
            size_t num_vert = opts["num-vert"].as<size_t>();
            size_t default_inserts = static_cast<size_t>(num_vert * 0.05);
            opts.insert(std::make_pair("num-inserts", bpo::variable_value(default_inserts, false)));
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

// usage: ./../run.sh add_vertices --num-parties 2 --num-clients 2 --num-vert 1000 --num-edge 4000