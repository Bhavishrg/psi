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

using namespace graphdb;
using json = nlohmann::json;
namespace bpo = boost::program_options;

common::utils::Circuit<Ring> generateSortCircuit(int nP, int pid, size_t vec_size, size_t input_size) {

    std::cout << "Generating sort circuit with vec_size=" << vec_size 
              << ", input_size=" << input_size << std::endl;
    std::cout << "Total wires: " << (vec_size * input_size) << std::endl;
    
    common::utils::Circuit<Ring> circ;

    // Input: vec_size * input_size wires (binary values 0 or 1)
    std::vector<common::utils::wire_t> input_wires(vec_size * input_size);
    
    for (size_t i = 0; i < vec_size * input_size; ++i) {
        input_wires[i] = circ.newInputWire();
    }

    // Generate identity permutation for preprocessing
    std::vector<std::vector<int>> permutations;
    std::vector<int> identity_perm(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        identity_perm[i] = i;
    }
    permutations.push_back(identity_perm);
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            permutations.push_back(identity_perm);
        }
    }

    // Use the sort subcircuit
    // Input wires are organized as: wire[i + j * input_size] for element j, component i
    // This matches the grouping expected by addSortSubcircuit
    auto sorted_indices = circ.addSortSubcircuit(input_wires, permutations, pid);
    
    // Set outputs: sorted indices
    for (size_t i = 0; i < sorted_indices.size(); ++i) {
        circ.setAsOutput(sorted_indices[i]);
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

    auto vec_size = opts["vec-size"].as<size_t>();
    auto input_size = opts["input-size"].as<size_t>();
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

    std::cout << "Starting sort benchmarks" << std::endl;

    std::string net_config = opts.count("net-config") ? opts["net-config"].as<std::string>() : "";
    std::shared_ptr<io::NetIOMP> network = createNetwork(pid, nP, latency, port,
                                                          opts["localhost"].as<bool>(),
                                                          net_config);

    // Increase socket buffer sizes to prevent deadlocks with large messages
    size_t expected_output_bytes = vec_size * sizeof(Ring);
    int buffer_size = std::max(128 * 1024 * 1024, 
                                static_cast<int>(expected_output_bytes * (nP - 1) * 3));
    increaseSocketBuffers(network.get(), buffer_size);

    json output_data;
    output_data["details"] = {{"num_parties", nP},
                              {"vec_size", vec_size},
                              {"input_size", input_size},
                              {"total_wires", vec_size * input_size},
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

    auto circ = generateSortCircuit(nP, pid, vec_size, input_size).orderGatesByLevel();
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
    OfflineEvaluator off_eval(nP, pid, network, circ, threads, seed, latency_us);
    auto preproc = off_eval.run(input_pid_map);
    std::cout << "Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint preproc_end(*network);

    std::cout << "Setting inputs" << std::endl;
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
    // Set inputs: binary values (0 or 1)
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
    
    // Generate random binary input values
    // Element j, component i: stored at wire[i + j * input_size]
    std::vector<std::vector<Ring>> input_values(vec_size, std::vector<Ring>(input_size));
    
    srand(seed);
    for (size_t elem = 0; elem < vec_size; ++elem) {
        for (size_t comp = 0; comp < input_size; ++comp) {
            // Binary values: 0 or 1
            input_values[elem][comp] = static_cast<Ring>(rand() % 2);
        }
    }
    
    // Set wire inputs: wire ordering is [comp0_elem0, comp1_elem0, ..., comp0_elem1, ...]
    size_t wire_idx = 0;
    for (size_t elem = 0; elem < vec_size && wire_idx < input_wires.size(); ++elem) {
        for (size_t comp = 0; comp < input_size && wire_idx < input_wires.size(); ++comp) {
            inputs[input_wires[wire_idx++]] = input_values[elem][comp];
        }
    }
    
    std::cout << "Total inputs set: " << inputs.size() << " wires" << std::endl;
    
    eval.setInputs(inputs);
    network->sync();
    
    std::cout << "Starting online evaluation" << std::endl;
    StatsPoint online_start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
    }
    network->flush();
    network->sync();
    StatsPoint online_end(*network);
    std::cout << "Online evaluation complete" << std::endl;

    std::cout << "Getting outputs..." << std::endl;
    network->flush();
    auto outputs = eval.getOutputs();
    network->sync();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Outputs: sorted indices [index_0, index_1, ..., index_{vec_size-1}]
    std::vector<Ring> sorted_indices(vec_size);
    
    if (outputs.size() >= vec_size) {
        for (size_t i = 0; i < vec_size; ++i) {
            sorted_indices[i] = outputs[i];
        }
    }
    
    // Print sorted output indices
    // Print first 10 input elements
    std::cout << "\n=== INPUT VALUES (Party " << pid << ") ===" << std::endl;
    std::cout << "Element | Components (binary) | Sort map" << std::endl;
    std::cout << "--------|---------------------|---------" << std::endl;
    // for (size_t c = 0; c < input_size; ++c) std::cout << "---";

    for (size_t elem = 0; elem < std::min(static_cast<size_t>(10), vec_size); ++elem) {
        std::cout << std::setw(7) << elem << " | ";
        for (size_t comp = 0; comp < input_size; ++comp) {
            std::cout << input_values[elem][comp];
            if (comp < input_size - 1) std::cout << " ";
        }
        std::cout << std::setw(7) << sorted_indices[elem] << std::endl;
    }
    std::cout << std::endl;
    
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

bpo::options_description programOptions() {
    bpo::options_description desc("Following options are supported by config file too.");
    desc.add_options()
        ("num-parties,n", bpo::value<int>()->required(), "Number of parties.")
        ("vec-size,v", bpo::value<size_t>()->required(), "Number of elements to sort.")
        ("input-size,i", bpo::value<size_t>()->required(), "Number of binary components per element.")
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

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark secure sort subcircuit.");
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
// usage: ./../run.sh sort --num-parties 2 --vec-size 100 --input-size 8
