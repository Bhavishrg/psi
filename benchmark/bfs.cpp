#include <io/netmp.h>
#include <graphdb/offline_evaluator.h>
#include <graphdb/online_evaluator.h>
#include <utils/circuit.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <chrono>
#include <iostream>
#include <memory>
#include <omp.h>
#include <fstream>
#include <random>
#include <vector>
#include <unordered_set>

#include "utils.h"
#include "graphutils.h"

using namespace graphdb;
using json = nlohmann::json;
namespace bpo = boost::program_options;

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, size_t vec_size, 
                                              size_t num_verts, int iterations,
                                              const Daglist& daglist) {

    std::cout << "Generating BFS circuit with vec_size=" << vec_size 
              << ", num_verts=" << num_verts 
              << ", iterations=" << iterations << std::endl;

    common::utils::Circuit<Ring> circ;

    // Generate identity permutation
    std::vector<int> base_perm(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        base_perm[i] = static_cast<int>(i);
    }
    std::vector<std::vector<int>> permutation;
    permutation.push_back(base_perm);
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            permutation.push_back(base_perm);
        }
    }

    // Create input wires for all daglist fields
    std::vector<common::utils::wire_t> src_wires(vec_size);
    std::vector<common::utils::wire_t> dst_wires(vec_size);
    std::vector<common::utils::wire_t> isV_wires(vec_size);
    std::vector<common::utils::wire_t> data_wires(vec_size);
    std::vector<common::utils::wire_t> sigs_wires(vec_size);
    std::vector<common::utils::wire_t> sigv_wires(vec_size);
    std::vector<common::utils::wire_t> sigd_wires(vec_size);
    
    std::generate(src_wires.begin(), src_wires.end(), [&]() { return circ.newInputWire(); });
    std::generate(dst_wires.begin(), dst_wires.end(), [&]() { return circ.newInputWire(); });
    std::generate(isV_wires.begin(), isV_wires.end(), [&]() { return circ.newInputWire(); });
    std::generate(data_wires.begin(), data_wires.end(), [&]() { return circ.newInputWire(); });
    std::generate(sigs_wires.begin(), sigs_wires.end(), [&]() { return circ.newInputWire(); });
    std::generate(sigv_wires.begin(), sigv_wires.end(), [&]() { return circ.newInputWire(); });
    std::generate(sigd_wires.begin(), sigd_wires.end(), [&]() { return circ.newInputWire(); });

    // BFS circuit: (propagate + permlist + gather) * iterations
    std::vector<common::utils::wire_t> current_data = data_wires;

    auto pos_map_v_to_s = sigs_wires;

    auto permuted_maps = circ.addSubCircPermList(pos_map_v_to_s, {sigd_wires, sigv_wires}, permutation);
    auto pos_map_s_to_d = permuted_maps[0];
    auto pos_map_d_to_v = circ.addSubCircPermList(pos_map_s_to_d, {permuted_maps[1]}, permutation)[0];

    for (int iter = 0; iter < iterations; ++iter) {
        std::cout << "  Adding BFS iteration " << (iter + 1) << std::endl;
        
        // Propagate step: forward message passing using sigs as position map
        auto propagated = circ.addSubCircPropagate(pos_map_v_to_s, current_data, num_verts, permutation, false);
        
        // Permutation gate: reorder propagated data using sigd as position map
        std::vector<std::vector<common::utils::wire_t>> payloads = {{propagated}};
        auto permuted = circ.addSubCircPermList(pos_map_s_to_d, payloads, permutation)[0];
        
        // Gather step: backward message passing using sigv as position map
        auto gathered = circ.addSubCircGather(pos_map_d_to_v, permuted, num_verts, permutation);
        
        // ApplyV step: update vertex data
        for (size_t i = 0; i < num_verts; ++i) {
            current_data[i] = circ.addGate(common::utils::GateType::kAdd, current_data[i], gathered[i]);
        }   
        for (size_t i = num_verts; i < vec_size; ++i) {
            current_data[i] = circ.addGate(common::utils::GateType::kSub, gathered[i], gathered[i]);
        }   
    }

    std::vector<common::utils::wire_t> final_data(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        final_data[i] = circ.addGate(common::utils::GateType::kEqz, current_data[i]);
    }
    // Create a constant-1 wire and set final outputs to (1 - final_data[i])
    auto zero_wire = circ.addGate(common::utils::GateType::kSub, final_data[0], final_data[0]);
    auto one_wire = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, static_cast<Ring>(1));
    for (size_t i = 0; i < vec_size; ++i) {
        auto out_wire = circ.addGate(common::utils::GateType::kSub, one_wire, final_data[i]);
        circ.setAsOutput(out_wire);
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

    auto nP = opts["num-parties"].as<int>();
    auto num_verts = opts["num-vert"].as<size_t>();
    auto num_edges = opts["num-edge"].as<size_t>();
    auto vec_size = num_verts + num_edges;
    auto iterations = opts["iterations"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto use_pking = opts["use-pking"].as<bool>();
    auto random_inputs = opts["random-inputs"].as<bool>();

    // Generate graph directly using num_verts and num_edges
    Ring nV = static_cast<Ring>(num_verts);
    Ring nE = static_cast<Ring>(num_edges);
    
    Daglist daglist;
    
    if (!random_inputs) {
        std::cout << "============================\n" << std::endl;
        std::cout << "Generating random inputs" << std::endl;
        std::cout << "Generating scale-free graph with nV=" << nV << ", nE=" << nE << " (seed=" << seed << ")" << std::endl;
        auto edges = generate_scale_free(nV, nE, seed);
        std::cout << "Generated " << edges.size() << " edges" << std::endl;
        
        std::cout << "Building daglist..." << std::endl;
        daglist = build_daglist(nV, edges);
        std::cout << "Built daglist with " << daglist.size() << " entries" << std::endl;
        
        // Update vec_size to match actual graph size
        vec_size = daglist.size();
    } else {
        std::cout << "============================\n" << std::endl;
        std::cout << "Using random inputs" << std::endl;
        // daglist will remain empty, inputs will be set randomly
    }

    // Print first 10 daglist entries
    if (!random_inputs) {
        std::cout << "\n=== First 10 Daglist Entries ===" << std::endl;
        size_t print_count = std::min(static_cast<size_t>(10), daglist.size());
        for (size_t i = 0; i < print_count; ++i) {
            const auto& entry = daglist.entries[i];
            std::cout << "Entry " << i << ": "
                      << "src=" << entry.src << ", "
                      << "dst=" << entry.dst << ", "
                      << "isV=" << entry.isV << ", "
                      << "data=" << entry.data << ", "
                      << "sigs=" << entry.sigs << ", "
                      << "sigv=" << entry.sigv << ", "
                      << "sigd=" << entry.sigd << std::endl;
        }
        std::cout << "================================\n" << std::endl;
    }

    std::cout << "Number of vertices: " << num_verts << std::endl;
    std::cout << "Vector size: " << vec_size << std::endl;

    omp_set_nested(1);
    if (nP < 10) { omp_set_num_threads(nP); }
    else { omp_set_num_threads(10); }
    std::cout << "Starting benchmarks" << std::endl;

    std::string net_config = opts.count("net-config") ? opts["net-config"].as<std::string>() : "";
    std::shared_ptr<io::NetIOMP> network = createNetwork(pid, nP, latency, port,
                                                          opts["localhost"].as<bool>(),
                                                          net_config);

    // Increase socket buffer sizes to prevent deadlocks with large messages
    increaseSocketBuffers(network.get(), 128 * 1024 * 1024);

    json output_data;
    output_data["details"] = {{"num_parties", nP},
                              {"vec_size", vec_size},
                              {"num_vertices", num_verts},
                              {"num_edges", num_edges},
                              {"iterations", iterations},
                              {"latency (ms)", latency},
                              {"pid", pid},
                              {"threads", threads},
                              {"seed", seed},
                              {"repeat", repeat},
                              {"random_inputs", random_inputs}};
    output_data["benchmarks"] = json::array();

    std::cout << "--- Details ---" << std::endl;
    for (const auto& [key, value] : output_data["details"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    StatsPoint start(*network);

    network->sync();

    auto circ = generateCircuit(nP, pid, vec_size, num_verts, iterations, daglist).orderGatesByLevel();
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
    int latency_us = static_cast<int>(latency * 1000);  // Convert ms to microseconds
    OfflineEvaluator off_eval(nP, pid, network, circ, threads, seed, latency_us);
    auto preproc = off_eval.run(input_pid_map);
    std::cout << "Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint preproc_end(*network);

    std::cout << "Setting inputs from daglist" << std::endl;
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
    std::unordered_map<common::utils::wire_t, Ring> inputs;
    std::vector<common::utils::wire_t> input_wires;
    input_wires.reserve(input_pid_map.size());
    for (const auto& [wire, owner] : input_pid_map) {
        if (owner == static_cast<int>(pid)) {
            input_wires.push_back(wire);
        }
    }
    std::sort(input_wires.begin(), input_wires.end());

    if (random_inputs) {
        // Use random inputs for benchmarking
        std::cout << "Using random inputs for party " << pid << std::endl;
        eval.setRandomInputs();
    } else {
        if (!input_wires.empty()) {
            std::cout << "\n=== SETTING BFS INPUTS FROM GENERATED GRAPH ===" << std::endl;
            std::cout << "Party " << pid << " setting inputs:" << std::endl;
            
            // Optionally modify daglist data before setting inputs (e.g., set starting vertex)
            // daglist.entries[0].data = 1;  // Set starting vertex data to 1
            
            // Use the helper function to set all daglist inputs
            set_daglist_inputs(daglist, input_wires, inputs);
            
            std::cout << "  Set " << (7 * vec_size) << " input values (7 fields × " << vec_size << " entries)" << std::endl;
            std::cout << "  BFS iterations: " << iterations << std::endl;
            std::cout << "=======================================\n" << std::endl;
        }
        eval.setInputs(inputs);
    }
    
    std::cout << "Starting online evaluation (BFS)" << std::endl;
    StatsPoint online_start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
    }

    auto outputs = eval.getOutputs();

    std::cout << "\n=== BFS RESULT ===" << std::endl;
    std::cout << "Party " << pid << " reconstructed outputs:" << std::endl;
    std::cout << "  Total number of outputs: " << outputs.size() << std::endl;
    
    // Display BFS outputs
    std::cout << "  BFS output (first 20): [";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), outputs.size()); ++i) {
        std::cout << outputs[i] << (i + 1 == std::min(static_cast<size_t>(20), outputs.size()) ? "" : ", ");
    }
    if (outputs.size() > 20) std::cout << ", ...";
    std::cout << "]" << std::endl;
    
    std::cout << "  ✓ BFS COMPLETE - Performed " << iterations 
              << " iterations of (propagate + gather)" << std::endl;
    std::cout << "==================\n" << std::endl;

    network->sync();
    StatsPoint online_end(*network);
    std::cout << "Online evaluation complete" << std::endl;

    StatsPoint end(*network);

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
        ("num-vert", bpo::value<size_t>()->default_value(1000), "Number of vertices in the graph.")
        ("num-edge", bpo::value<size_t>()->default_value(4000), "Number of edges in the graph.")
        ("iterations,i", bpo::value<int>()->required(), "Number of BFS iterations (propagate + gather).")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
        ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
        ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
        ("localhost", bpo::bool_switch(), "All parties are on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
        ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.")
        ("random-inputs", bpo::value<bool>()->default_value(false), "Use random inputs for benchmarking.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction (true) or direct reconstruction (false).");
  return desc;
}
// clang-format on

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark BFS using (propagate + gather) * iterations with generated scale-free graph.");
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

// usage: ./../run.sh bfs --num-parties 2 --num-vert 1000 --num-edge 4000 --iterations 3 --random-inputs 1