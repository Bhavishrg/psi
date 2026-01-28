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

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, size_t vec_size, size_t num_groups) {

    std::cout << "Generating propagate circuit with vec_size=" << vec_size 
              << ", num_groups=" << num_groups << std::endl;

    common::utils::Circuit<Ring> circ;

    // Generate permutation for shuffle
    // Here we just pass identity permutations
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

    // Create position map input wires (secret shares)
    std::vector<common::utils::wire_t> position_map_shares(vec_size);
    std::generate(position_map_shares.begin(), position_map_shares.end(), [&]() { return circ.newInputWire(); });

    // Create data values input wires (secret shares)
    std::vector<common::utils::wire_t> data_values(vec_size);
    std::generate(data_values.begin(), data_values.end(), [&]() { return circ.newInputWire(); });

    // Use utility function to add propagate sub-circuit
    auto prefix_sum = circ.addSubCircPropagate(position_map_shares, data_values, num_groups, permutation);

    // Set outputs
    for (size_t i = 0; i < vec_size; ++i) {
        circ.setAsOutput(prefix_sum[i]);
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
    auto vec_size = opts["vec-size"].as<size_t>();
    auto num_groups = opts["num-groups"].as<size_t>();
    auto iter = opts["iter"].as<int>();
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
                              {"vec_size", vec_size},
                              {"num_groups", num_groups},
                              {"iterations", iter},
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

    auto circ = generateCircuit(nP, pid, vec_size, num_groups).orderGatesByLevel();
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

    std::cout << "Setting inputs" << std::endl;
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

    if (!input_wires.empty()) {
        std::cout << "\n=== SETTING TEST INPUTS ===" << std::endl;
        std::cout << "Party " << pid << " setting inputs:" << std::endl;
        
        // Create a test position map (identity permutation for simplicity)
        std::vector<Ring> position_map_values(vec_size);

        position_map_values[0] = static_cast<Ring>(0);
        position_map_values[1] = static_cast<Ring>(3);
        position_map_values[2] = static_cast<Ring>(5);
        position_map_values[3] = static_cast<Ring>(9);
        position_map_values[4] = static_cast<Ring>(1);
        position_map_values[5] = static_cast<Ring>(8);
        position_map_values[6] = static_cast<Ring>(6);
        position_map_values[7] = static_cast<Ring>(7);
        position_map_values[8] = static_cast<Ring>(4);
        position_map_values[9] = static_cast<Ring>(2);

        for (size_t i = 10; i < vec_size; ++i) {
            position_map_values[i] = static_cast<Ring>(i);
        }
        
        // Create test data values (incremental values)
        std::vector<Ring> data_values(vec_size);
        for (size_t i = 0; i < vec_size; ++i) {
            data_values[0] = static_cast<Ring>(0);
            data_values[1] = static_cast<Ring>(1);
            data_values[2] = static_cast<Ring>(0);
            data_values[3] = static_cast<Ring>(0);
            data_values[4] = static_cast<Ring>(0);
            data_values[5] = static_cast<Ring>(1);
            data_values[6] = static_cast<Ring>(0);
            data_values[7] = static_cast<Ring>(0);
            data_values[8] = static_cast<Ring>(0);
            data_values[9] = static_cast<Ring>(0);
        }
        
        // Set position map inputs
        for (size_t idx = 0; idx < vec_size && idx < input_wires.size(); ++idx) {
            inputs[input_wires[idx]] = position_map_values[idx];
        }
        
        // Set data values inputs
        for (size_t idx = 0; idx < vec_size && (vec_size + idx) < input_wires.size(); ++idx) {
            inputs[input_wires[vec_size + idx]] = data_values[idx];
        }
        
        std::cout << "  Position map (first 10): [";
        for (size_t i = 0; i < std::min(static_cast<size_t>(10), vec_size); ++i) {
            std::cout << position_map_values[i] << (i + 1 == std::min(static_cast<size_t>(10), vec_size) ? "" : ", ");
        }
        if (vec_size > 10) std::cout << ", ...";
        std::cout << "]" << std::endl;
        
        std::cout << "  Data values (first 10): [";
        for (size_t i = 0; i < std::min(static_cast<size_t>(10), vec_size); ++i) {
            std::cout << data_values[i] << (i + 1 == std::min(static_cast<size_t>(10), vec_size) ? "" : ", ");
        }
        if (vec_size > 10) std::cout << ", ...";
        std::cout << "]" << std::endl;
        
        std::cout << "  Set " << input_wires.size() << " input values" << std::endl;
        std::cout << "  Num groups: " << num_groups << std::endl;
        std::cout << "========================\n" << std::endl;
    }

    eval.setInputs(inputs);
    
    std::cout << "Starting online evaluation" << std::endl;
    StatsPoint online_start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
    }

    auto outputs = eval.getOutputs();

    std::cout << "\n=== PROPAGATE RESULT ===" << std::endl;
    std::cout << "Party " << pid << " reconstructed outputs:" << std::endl;
    std::cout << "  Total number of outputs: " << outputs.size() << std::endl;
    
    // Display propagate outputs (prefix sums)
    std::cout << "  Prefix sum output (first 20): [";
    for (size_t i = 0; i < std::min(static_cast<size_t>(20), outputs.size()); ++i) {
        std::cout << outputs[i] << (i + 1 == std::min(static_cast<size_t>(20), outputs.size()) ? "" : ", ");
    }
    if (outputs.size() > 20) std::cout << ", ...";
    std::cout << "]" << std::endl;
    
    std::cout << "  ✓ PROPAGATE COMPLETE - Computed differences, shuffled, rewired, and prefix summed" << std::endl;
    std::cout << "============================\n" << std::endl;

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
        ("vec-size,v", bpo::value<size_t>()->required(), "Size of position map and data vectors.")
        ("num-groups,g", bpo::value<size_t>()->required(), "Number of groups (must be less than vec-size).")
        ("iter,i", bpo::value<int>()->default_value(1), "Number of iterations for message passing.")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
        ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
        ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
        ("localhost", bpo::bool_switch(), "All parties are on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
        ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction (true) or direct reconstruction (false).");
  return desc;
}
// clang-format on

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark online phase for propagate circuit.");
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
// usage: ./../run.sh propagate --num-parties 2 --vec-size 10 --num-groups 4