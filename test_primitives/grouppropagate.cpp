#include <io/netmp.h>
#include <graphdb/offline_evaluator.h>
#include <graphdb/online_evaluator.h>
#include <utils/circuit.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <omp.h>

#include "utils.h"

using namespace graphdb;
using json = nlohmann::json;
namespace bpo = boost::program_options;

common::utils::Circuit<Ring> generateGroupwisePropagateCircuit(int nP, int pid, 
                                                                size_t t1_vec_size, 
                                                                size_t t2_vec_size) {

    std::cout << "Generating Group-wise Propagate circuit for T1_size=" << t1_vec_size 
              << ", T2_size=" << t2_vec_size << std::endl;
    
    common::utils::Circuit<Ring> circ;

    // Input T1: key1 vector (binary) and v1 vector (data values)
    std::vector<common::utils::wire_t> key1_vector(t1_vec_size);
    std::vector<common::utils::wire_t> v1_vector(t1_vec_size);
    
    std::generate(key1_vector.begin(), key1_vector.end(), [&]() { return circ.newInputWire(); });
    std::generate(v1_vector.begin(), v1_vector.end(), [&]() { return circ.newInputWire(); });

    // Input T2: key2 vector (binary, group markers)
    std::vector<common::utils::wire_t> key2_vector(t2_vec_size);
    std::generate(key2_vector.begin(), key2_vector.end(), [&]() { return circ.newInputWire(); });

    // Prepare permutations for T1 compaction
    std::vector<std::vector<int>> permutation_t1;
    std::vector<int> t1_perm(t1_vec_size);
    for (size_t i = 0; i < t1_vec_size; ++i) {
        t1_perm[i] = i;
    }
    for (int p = 0; p < nP; ++p) {
        permutation_t1.push_back(t1_perm);
    }
    
    // Prepare permutations for T2 compaction
    std::vector<std::vector<int>> permutation_t2;
    std::vector<int> t2_perm(t2_vec_size);
    for (size_t i = 0; i < t2_vec_size; ++i) {
        t2_perm[i] = i;
    }
    for (int p = 0; p < nP; ++p) {
        permutation_t2.push_back(t2_perm);
    }
    
    // Prepare reverse permutation (same as t2 for identity permutation)
    std::vector<std::vector<int>> permutation_rev;
    for (int p = 0; p < nP; ++p) {
        permutation_rev.push_back(t2_perm);
    }

    // Use the groupwise propagate subcircuit
    auto [output_key2, output_v] = circ.addGroupwisePropagateSubcircuit(key1_vector, v1_vector, 
                                                                          key2_vector, permutation_t1, 
                                                                          permutation_t2, permutation_rev, pid);

    // Set outputs: key2 (restored) and v (propagated values)
    for (size_t i = 0; i < t2_vec_size; ++i) {
        circ.setAsOutput(output_key2[i]);  // T_O[i].key
    }

    for (size_t i = 0; i < t2_vec_size; ++i) {
        circ.setAsOutput(output_v[i]);     // T_O[i].v
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

    auto t1_size = opts["t1-size"].as<size_t>();
    auto t2_size = opts["t2-size"].as<size_t>();
    auto nP = opts["num-parties"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();

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
                              {"t1_size", t1_size},
                              {"t2_size", t2_size},
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

    StatsPoint start(*network);
    network->sync();

    auto circ = generateGroupwisePropagateCircuit(nP, pid, t1_size, t2_size).orderGatesByLevel();
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

    std::cout << "Starting online evaluation" << std::endl;
    StatsPoint online_start(*network);
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us);
    
    // Set inputs
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
    
    // Structure: first t1_size wires = key1, next t1_size wires = v1, next t2_size wires = key2
    std::vector<Ring> key1_input_values(t1_size);
    std::vector<Ring> v1_input_values(t1_size);
    std::vector<Ring> key2_input_values(t2_size);
    
    for (size_t i = 0; i < input_wires.size(); ++i) {
        auto wire = input_wires[i];
        if (i < t1_size) {
            // key1 vector: binary values (0 or 1), indicating data entries
            // For demonstration, let's have a few marked entries (key=1)
            Ring val = 0;
            if (i == 0 || i == 2 || i == 5) {  // Mark specific entries
                val = 1;
            }
            inputs[wire] = val;
            key1_input_values[i] = val;
        } else if (i < 2 * t1_size) {
            // v1 vector: data values to propagate
            Ring val = static_cast<Ring>(100 + (i - t1_size) * 10);  // Values: 100, 110, 120, ...
            inputs[wire] = val;
            v1_input_values[i - t1_size] = val;
        } else {
            // key2 vector: binary values (0 or 1), group markers
            size_t idx = i - 2 * t1_size;
            Ring val = 0;
            if (idx == 0 || idx % 3 == 0) {  // Start new group every 3 elements
                val = 1;
            }
            inputs[wire] = val;
            key2_input_values[idx] = val;
        }
    }

    // Print inputs in structured format (first 20 entries)
    std::cout << "\n=== INPUT VECTORS (Party " << pid << ") - First 20 entries ===" << std::endl;
    std::cout << "T1 - Key1 vector (1=marked entry): ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t1_size); ++i) {
        std::cout << key1_input_values[i] << " ";
    }
    if (t1_size > 20) std::cout << "...";
    std::cout << std::endl;
    
    std::cout << "T1 - Value vector (to propagate):  ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t1_size); ++i) {
        std::cout << v1_input_values[i] << " ";
    }
    if (t1_size > 20) std::cout << "...";
    std::cout << std::endl;
    
    std::cout << "T2 - Key2 vector (1=group start):  ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2_size); ++i) {
        std::cout << key2_input_values[i] << " ";
    }
    if (t2_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "Total inputs set: " << inputs.size() << std::endl;
    
    eval.setInputs(inputs);
    
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
        // Sync after each level to ensure all parties stay synchronized
        network->sync();
    }
    network->sync();
    
    auto outputs = eval.getOutputs();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print outputs in structured format (first 20 entries)
    // Output structure: t2_size for key2, t2_size for v (propagated)
    std::cout << "\n=== OUTPUT VECTORS (Party " << pid << ") - First 20 entries ===" << std::endl;
    std::cout << "Key2 (restored):             ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2_size); ++i) {
        std::cout << outputs[i] << " ";
    }
    if (t2_size > 20) std::cout << "...";
    std::cout << std::endl;
    
    std::cout << "Value (propagated to groups): ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2_size); ++i) {
        std::cout << outputs[t2_size + i] << " ";
    }
    if (t2_size > 20) std::cout << "...";
    std::cout << std::endl << std::endl;
    
    network->sync();
    StatsPoint online_end(*network);

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
        ("t1-size", bpo::value<size_t>()->required(), "Size of T1 list (key-value pairs).")
        ("t2-size", bpo::value<size_t>()->required(), "Size of T2 list (grouped list).")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
        ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
        ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
        ("localhost", bpo::bool_switch(), "All parties are on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
        ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.");
  return desc;
}
// clang-format on

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark Group-wise Propagate protocol circuit.");
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