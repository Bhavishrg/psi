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
                                                                size_t t1a_vec_size, 
                                                                size_t t2a_vec_size,
                                                                size_t t1b_vec_size,
                                                                size_t t2b_vec_size) {

    std::cout << "Generating Group-wise Propagate circuit for Gate1(T1_size=" << t1a_vec_size 
              << ", T2_size=" << t2a_vec_size << "), Gate2(T1_size=" << t1b_vec_size 
              << ", T2_size=" << t2b_vec_size << ")" << std::endl;
    
    common::utils::Circuit<Ring> circ;

    // Input for Gate 1: key1 vector (binary) and v1 vector (data values)
    std::vector<common::utils::wire_t> key11_vector(t1a_vec_size);
    std::vector<common::utils::wire_t> v1_vector(t1a_vec_size);
    std::vector<common::utils::wire_t> key12_vector(t2a_vec_size);

    // Input for Gate 2: key1 vector (binary) and v2 vector (data values)
    std::vector<common::utils::wire_t> key21_vector(t1b_vec_size);
    std::vector<common::utils::wire_t> v2_vector(t1b_vec_size);
    std::vector<common::utils::wire_t> key22_vector(t2b_vec_size);


    std::generate(key11_vector.begin(), key11_vector.end(), [&]() { return circ.newInputWire(); });
    std::generate(v1_vector.begin(), v1_vector.end(), [&]() { return circ.newInputWire(); });
    std::generate(key12_vector.begin(), key12_vector.end(), [&]() { return circ.newInputWire(); });

    std::generate(key21_vector.begin(), key21_vector.end(), [&]() { return circ.newInputWire(); });
    std::generate(v2_vector.begin(), v2_vector.end(), [&]() { return circ.newInputWire(); });
    std::generate(key22_vector.begin(), key22_vector.end(), [&]() { return circ.newInputWire(); });

    // Prepare permutations for Gate 1 - T1 compaction
    std::vector<std::vector<int>> permutation1_t1;
    std::vector<int> t1a_perm(t1a_vec_size);
    for (size_t i = 0; i < t1a_vec_size; ++i) {
        t1a_perm[i] = i;
    }
    for (int p = 0; p < nP; ++p) {
        permutation1_t1.push_back(t1a_perm);
    }
    
    // Prepare permutations for Gate 1 - T2 compaction
    std::vector<std::vector<int>> permutation1_t2;
    std::vector<int> t2a_perm(t2a_vec_size);
    for (size_t i = 0; i < t2a_vec_size; ++i) {
        t2a_perm[i] = i;
    }
    for (int p = 0; p < nP; ++p) {
        permutation1_t2.push_back(t2a_perm);
    }
    
    // Prepare reverse permutation for Gate 1
    std::vector<std::vector<int>> permutation1_rev;
    for (int p = 0; p < nP; ++p) {
        permutation1_rev.push_back(t2a_perm);
    }
    
    // Prepare permutations for Gate 2 - T1 compaction
    std::vector<std::vector<int>> permutation2_t1;
    std::vector<int> t1b_perm(t1b_vec_size);
    for (size_t i = 0; i < t1b_vec_size; ++i) {
        t1b_perm[i] = i;
    }
    for (int p = 0; p < nP; ++p) {
        permutation2_t1.push_back(t1b_perm);
    }
    
    // Prepare permutations for Gate 2 - T2 compaction
    std::vector<std::vector<int>> permutation2_t2;
    std::vector<int> t2b_perm(t2b_vec_size);
    for (size_t i = 0; i < t2b_vec_size; ++i) {
        t2b_perm[i] = i;
    }
    for (int p = 0; p < nP; ++p) {
        permutation2_t2.push_back(t2b_perm);
    }
    
    // Prepare reverse permutation for Gate 2
    std::vector<std::vector<int>> permutation2_rev;
    for (int p = 0; p < nP; ++p) {
        permutation2_rev.push_back(t2b_perm);
    }

    // Use the groupwise propagate subcircuit for both gates
    auto [output_key21, output_v1] = circ.addGroupwisePropagateSubcircuit(key11_vector, v1_vector, 
                                                                            key12_vector, permutation1_t1,
                                                                            permutation1_t2, permutation1_rev, pid);

    auto [output_key22, output_v2] = circ.addGroupwisePropagateSubcircuit(key21_vector, v2_vector, 
                                                                            key22_vector, permutation2_t1,
                                                                            permutation2_t2, permutation2_rev, pid);

    // Set outputs: key2 (restored) and v (propagated values)
    for (size_t i = 0; i < t2a_vec_size; ++i) {
        circ.setAsOutput(output_key21[i]);  // T_O[i].key
    }

    for (size_t i = 0; i < t2a_vec_size; ++i) {
        circ.setAsOutput(output_v1[i]);     // T_O[i].v
    }

    for (size_t i = 0; i < t2b_vec_size; ++i) {
        circ.setAsOutput(output_key22[i]);  // T_O[i].key
    }

    for (size_t i = 0; i < t2b_vec_size; ++i) {
        circ.setAsOutput(output_v2[i]);     // T_O[i].v
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

    auto t1a_size = opts["t1a-size"].as<size_t>();
    auto t2a_size = opts["t2a-size"].as<size_t>();
    auto t1b_size = opts["t1b-size"].as<size_t>();
    auto t2b_size = opts["t2b-size"].as<size_t>();
    auto nP = opts["num-parties"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto use_pking = opts["use-pking"].as<bool>();

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
                              {"t1a_size", t1a_size},
                              {"t2a_size", t2a_size},
                              {"t1b_size", t1b_size},
                              {"t2b_size", t2b_size},
                              {"latency (ms)", latency},
                              {"pid", pid},
                              {"threads", threads},
                              {"seed", seed},
                              {"repeat", repeat},
                              {"use_pking", use_pking}};
    output_data["benchmarks"] = json::array();

    std::cout << "--- Details ---" << std::endl;
    for (const auto& [key, value] : output_data["details"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    StatsPoint start(*network);
    network->sync();

    auto circ = generateGroupwisePropagateCircuit(nP, pid, t1a_size, t2a_size, t1b_size, t2b_size).orderGatesByLevel();
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
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
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
    
    // Structure: We have 2 gates with different sizes
    // Gate 1: key11 (t1a_size) + v1 (t1a_size) + key12 (t2a_size)
    // Gate 2: key21 (t1b_size) + v2 (t1b_size) + key22 (t2b_size)
    std::vector<Ring> key1a_input_values(t1a_size);
    std::vector<Ring> v1_input_values(t1a_size);
    std::vector<Ring> key2a_input_values(t2a_size);
    std::vector<Ring> key1b_input_values(t1b_size);
    std::vector<Ring> v2_input_values(t1b_size);
    std::vector<Ring> key2b_input_values(t2b_size);

    size_t gate1_inputs = 2 * t1a_size + t2a_size;
    size_t gate2_inputs = 2 * t1b_size + t2b_size;
    
    for (size_t i = 0; i < input_wires.size(); ++i) {
        auto wire = input_wires[i];
        
        if (i < gate1_inputs) {
            // Gate 1 inputs
            if (i < t1a_size) {
                // key1a vector
                Ring val = 0;
                if (i == 0 || i == 2 || (i == 5 && i < t1a_size)) {
                    val = 1;
                }
                inputs[wire] = val;
                key1a_input_values[i] = val;
            } else if (i < 2 * t1a_size) {
                // v1 vector
                size_t v_idx = i - t1a_size;
                Ring val = static_cast<Ring>(100 + v_idx * 10);
                inputs[wire] = val;
                v1_input_values[v_idx] = val;
            } else {
                // key2a vector
                size_t key2_idx = i - 2 * t1a_size;
                Ring val = 0;
                if (key2_idx == 0 || key2_idx % 3 == 0) {
                    val = 1;
                }
                inputs[wire] = val;
                key2a_input_values[key2_idx] = val;
            }
        } else {
            // Gate 2 inputs
            size_t gate2_offset = i - gate1_inputs;
            if (gate2_offset < t1b_size) {
                // key1b vector
                Ring val = 0;
                if (gate2_offset == 0 || gate2_offset == 2 || (gate2_offset == 5 && gate2_offset < t1b_size)) {
                    val = 1;
                }
                inputs[wire] = val;
                key1b_input_values[gate2_offset] = val;
            } else if (gate2_offset < 2 * t1b_size) {
                // v2 vector
                size_t v_idx = gate2_offset - t1b_size;
                Ring val = static_cast<Ring>(100 + v_idx * 10);
                inputs[wire] = val;
                v2_input_values[v_idx] = val;
            } else {
                // key2b vector
                size_t key2_idx = gate2_offset - 2 * t1b_size;
                Ring val = 0;
                if (key2_idx == 0 || key2_idx % 3 == 0) {
                    val = 1;
                }
                inputs[wire] = val;
                key2b_input_values[key2_idx] = val;
            }
        }
    }

    // Print inputs in structured format (first 20 entries)
    std::cout << "\n=== INPUT VECTORS (Party " << pid << ") - First 20 entries ===" << std::endl;
    std::cout << "Gate 1:" << std::endl;
    std::cout << "  T1a - Key1 vector (1=marked entry): ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t1a_size); ++i) {
        std::cout << key1a_input_values[i] << " ";
    }
    if (t1a_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "  T1a - Value vector (to propagate):  ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t1a_size); ++i) {
        std::cout << v1_input_values[i] << " ";
    }
    if (t1a_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "  T2a - Key2 vector (1=group start):  ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2a_size); ++i) {
        std::cout << key2a_input_values[i] << " ";
    }
    if (t2a_size > 20) std::cout << "...";
    std::cout << std::endl;
    
    std::cout << "Gate 2:" << std::endl;
    std::cout << "  T1b - Key1 vector (1=marked entry): ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t1b_size); ++i) {
        std::cout << key1b_input_values[i] << " ";
    }
    if (t1b_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "  T1b - Value vector (to propagate):  ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t1b_size); ++i) {
        std::cout << v2_input_values[i] << " ";
    }
    if (t1b_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "  T2b - Key2 vector (1=group start):  ";
    for (size_t i = 0; i < t2b_size; ++i) {
        std::cout << key2b_input_values[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Total inputs set: " << inputs.size() << std::endl;
    
    eval.setInputs(inputs);
    std::cout << "Total inputs setting done" << std::endl;
    
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
        // Sync after each level to ensure all parties stay synchronized
        network->sync();
    }
    network->sync();
    
    auto outputs = eval.getOutputs();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print outputs in structured format (first 20 entries)
    std::cout << "\n=== OUTPUT VECTORS (Party " << pid << ") - First 20 entries ===" << std::endl;
    std::cout << "Gate 1:" << std::endl;
    std::cout << "  Key2a (restored):             ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2a_size); ++i) {
        std::cout << outputs[i] << " ";
    }
    if (t2a_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "  Value1 (propagated to groups): ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2a_size); ++i) {
        std::cout << outputs[t2a_size + i] << " ";
    }
    if (t2a_size > 20) std::cout << "...";
    std::cout << std::endl;

    std::cout << "Gate 2:" << std::endl;
    std::cout << "  Key2b (restored):             ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2b_size); ++i) {
        std::cout << outputs[2*t2a_size + i] << " ";
    }
    if (t2b_size > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "  Value2 (propagated to groups): ";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), t2b_size); ++i) {
        std::cout << outputs[2*t2a_size + t2b_size + i] << " ";
    }
    if (t2b_size > 20) std::cout << "...";
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
        ("t1a-size", bpo::value<size_t>()->required(), "Size of T1 list for Gate 1 (key-value pairs).")
        ("t2a-size", bpo::value<size_t>()->required(), "Size of T2 list for Gate 1 (grouped list).")
        ("t1b-size", bpo::value<size_t>()->required(), "Size of T1 list for Gate 2 (key-value pairs).")
        ("t2b-size", bpo::value<size_t>()->required(), "Size of T2 list for Gate 2 (grouped list).")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
        ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction (true) or direct reconstruction (false).")
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
// ./../run.sh grouppropagate_parallel --num-parties 3 --t1a-size 10 --t2a-size 10 --t1b-size 10 --t2b-size 10