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

common::utils::Circuit<Ring> generateCircuit1(size_t vec_size, int nP, int pid) {
    common::utils::Circuit<Ring> circ;
    
    // Create input wires
    std::vector<common::utils::wire_t> input_wires(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        input_wires[i] = circ.newInputWire();
    }
    
    // Create identity permutations
    std::vector<std::vector<int>> permutations;
    std::vector<int> identity(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        identity[i] = static_cast<int>(i);
    }
    permutations.push_back(identity);
    
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            permutations.push_back(identity);
        }
    }
    
    // Two parallel sort gates
    auto sorted1 = circ.addSortSubcircuit(input_wires, permutations, pid);
    auto sorted2 = circ.addSortSubcircuit(input_wires, permutations, pid);
    
    // Set outputs
    for (size_t i = 0; i < sorted1.size(); ++i) {
        circ.setAsOutput(sorted1[i]);
    }
    for (size_t i = 0; i < sorted2.size(); ++i) {
        circ.setAsOutput(sorted2[i]);
    }
    
    return circ;
}

common::utils::Circuit<Ring> generateCircuit2(size_t vec_size, int nP, int pid) {
    common::utils::Circuit<Ring> circ;
    
    // Create input wires
    std::vector<common::utils::wire_t> input_wires(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        input_wires[i] = circ.newInputWire();
    }
    
    // Create identity permutations
    std::vector<std::vector<int>> permutations;
    std::vector<int> identity(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        identity[i] = static_cast<int>(i);
    }
    permutations.push_back(identity);
    
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            permutations.push_back(identity);
        }
    }
    
    // Two sequential shuffles
    auto shuffled1 = circ.addMGate(common::utils::GateType::kShuffle, input_wires, permutations);
    auto shuffled2 = circ.addMGate(common::utils::GateType::kShuffle, shuffled1, permutations);
    
    // Two parallel reconstruction gates
    std::vector<common::utils::wire_t> reconstructed1(vec_size);
    std::vector<common::utils::wire_t> reconstructed2(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        reconstructed1[i] = circ.addGate(common::utils::GateType::kRec, shuffled1[i]);
        reconstructed2[i] = circ.addGate(common::utils::GateType::kRec, shuffled2[i]);
    }
    
    // Set outputs
    for (size_t i = 0; i < vec_size; ++i) {
        circ.setAsOutput(reconstructed1[i]);
        circ.setAsOutput(reconstructed2[i]);
    }
    
    return circ;
}

bpo::options_description programOptions() {
    bpo::options_description desc("Options");
    desc.add_options()
        ("num-parties,n", bpo::value<int>()->required(), "Number of parties.")
        ("num-vertices", bpo::value<size_t>()->required(), "Number of vertices.")
        ("num-edges", bpo::value<size_t>()->required(), "Number of edges.")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads.")
        ("seed", bpo::value<size_t>()->default_value(200), "Random seed.")
        ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
        ("localhost", bpo::bool_switch(), "All parties are on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
        ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction (true) or direct reconstruction (false).");
    return desc;
}

void benchmark(const bpo::variables_map& opts) {
    bool save_output = false;
    std::string save_file;
    if (opts.count("output") != 0) {
        save_output = true;
        save_file = opts["output"].as<std::string>();
    }

    auto num_parties = opts["num-parties"].as<int>();
    auto num_vertices = opts["num-vertices"].as<size_t>();
    auto num_edges = opts["num-edges"].as<size_t>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto use_pking = opts["use-pking"].as<bool>();

    size_t vec_size = num_vertices + num_edges;
    size_t log_factor = static_cast<size_t>(std::ceil(std::log2(std::max<size_t>(1, num_vertices))));

    omp_set_nested(1);
    if (num_parties < 10) { omp_set_num_threads(num_parties); }
    else { omp_set_num_threads(10); }

    std::cout << "Starting Graphiti init micro benchmarks" << std::endl;

    std::string net_config = opts.count("net-config") ? opts["net-config"].as<std::string>() : "";
    std::shared_ptr<io::NetIOMP> network = createNetwork(pid, num_parties, latency, port,
                                                          opts["localhost"].as<bool>(),
                                                          net_config);

    json output_data;
    output_data["details"] = {{"num_parties", num_parties},
                               {"num_vertices", num_vertices},
                               {"num_edges", num_edges},
                               {"vec_size", vec_size},
                               {"log_factor", log_factor},
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

    StatsPoint global_start(*network);
    network->sync();

    // ===== Circuit 1: Two parallel sort gates =====
    std::cout << "\n=== Evaluating Circuit 1 (Two parallel sorts) ===" << std::endl;
    auto circ1 = generateCircuit1(vec_size, num_parties, pid).orderGatesByLevel();
    network->sync();

    std::cout << "--- Circuit 1 ---" << std::endl;
    std::cout << circ1 << std::endl;

    std::unordered_map<common::utils::wire_t, int> input_pid_map1;
    for (const auto& g : circ1.gates_by_level[0]) {
        if (g->type == common::utils::GateType::kInp) {
            input_pid_map1[g->out] = 1;
        }
    }

    std::cout << "Circuit 1: Starting preprocessing" << std::endl;
    StatsPoint circ1_preproc_start(*network);
    int latency_us = static_cast<int>(latency * 1000);
    OfflineEvaluator off_eval1(num_parties, pid, network, circ1, threads, seed, latency_us, use_pking);
    auto preproc1 = off_eval1.run(input_pid_map1);
    std::cout << "Circuit 1: Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint circ1_preproc_end(*network);

    std::cout << "Circuit 1: Setting inputs" << std::endl;
    OnlineEvaluator eval1(num_parties, pid, network, std::move(preproc1), circ1, threads, seed, latency_us, use_pking);
    
    std::unordered_map<common::utils::wire_t, Ring> inputs1;
    std::vector<common::utils::wire_t> input_wires1;
    for (const auto& [wire, owner] : input_pid_map1) {
        if (owner == static_cast<int>(pid)) {
            input_wires1.push_back(wire);
        }
    }
    std::sort(input_wires1.begin(), input_wires1.end());
    
    srand(seed);
    for (size_t i = 0; i < input_wires1.size(); ++i) {
        inputs1[input_wires1[i]] = Ring(rand() % 100);
    }
    
    eval1.setInputs(inputs1);
    network->sync();
    
    std::cout << "Circuit 1: Starting online evaluation" << std::endl;
    StatsPoint circ1_online_start(*network);
    for (size_t i = 0; i < circ1.gates_by_level.size(); ++i) {
        eval1.evaluateGatesAtDepth(i);
    }
    network->flush();
    network->sync();
    StatsPoint circ1_online_end(*network);
    std::cout << "Circuit 1: Online evaluation complete" << std::endl;

    auto circ1_preproc_bench = circ1_preproc_end - circ1_preproc_start;
    auto circ1_online_bench = circ1_online_end - circ1_online_start;

    // ===== Circuit 2: Two shuffles + Two reconstructions =====
    std::cout << "\n=== Evaluating Circuit 2 (Two shuffles + Two reconstructions) ===" << std::endl;
    auto circ2 = generateCircuit2(vec_size, num_parties, pid).orderGatesByLevel();
    network->sync();

    std::cout << "--- Circuit 2 ---" << std::endl;
    std::cout << circ2 << std::endl;

    std::unordered_map<common::utils::wire_t, int> input_pid_map2;
    for (const auto& g : circ2.gates_by_level[0]) {
        if (g->type == common::utils::GateType::kInp) {
            input_pid_map2[g->out] = 1;
        }
    }

    std::cout << "Circuit 2: Starting preprocessing" << std::endl;
    StatsPoint circ2_preproc_start(*network);
    OfflineEvaluator off_eval2(num_parties, pid, network, circ2, threads, seed, latency_us, use_pking);
    auto preproc2 = off_eval2.run(input_pid_map2);
    std::cout << "Circuit 2: Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint circ2_preproc_end(*network);

    std::cout << "Circuit 2: Setting inputs" << std::endl;
    OnlineEvaluator eval2(num_parties, pid, network, std::move(preproc2), circ2, threads, seed, latency_us, use_pking);
    
    std::unordered_map<common::utils::wire_t, Ring> inputs2;
    std::vector<common::utils::wire_t> input_wires2;
    for (const auto& [wire, owner] : input_pid_map2) {
        if (owner == static_cast<int>(pid)) {
            input_wires2.push_back(wire);
        }
    }
    std::sort(input_wires2.begin(), input_wires2.end());
    
    for (size_t i = 0; i < input_wires2.size(); ++i) {
        inputs2[input_wires2[i]] = Ring(rand() % 100);
    }
    
    eval2.setInputs(inputs2);
    network->sync();
    
    std::cout << "Circuit 2: Starting online evaluation" << std::endl;
    StatsPoint circ2_online_start(*network);
    for (size_t i = 0; i < circ2.gates_by_level.size(); ++i) {
        eval2.evaluateGatesAtDepth(i);
    }
    network->flush();
    network->sync();
    StatsPoint circ2_online_end(*network);
    std::cout << "Circuit 2: Online evaluation complete" << std::endl;

    auto circ2_preproc_bench = circ2_preproc_end - circ2_preproc_start;
    auto circ2_online_bench = circ2_online_end - circ2_online_start;

    StatsPoint global_end(*network);
    auto total_bench = global_end - global_start;

    // Extract timing and communication data
    double circ1_preproc_time = circ1_preproc_bench["time"].get<double>();
    double circ1_online_time = circ1_online_bench["time"].get<double>();
    double circ2_preproc_time = circ2_preproc_bench["time"].get<double>();
    double circ2_online_time = circ2_online_bench["time"].get<double>();

    size_t circ1_preproc_bytes = 0;
    for (const auto& val : circ1_preproc_bench["communication"]) {
        circ1_preproc_bytes += val.get<int64_t>();
    }
    size_t circ1_online_bytes = 0;
    for (const auto& val : circ1_online_bench["communication"]) {
        circ1_online_bytes += val.get<int64_t>();
    }
    size_t circ2_preproc_bytes = 0;
    for (const auto& val : circ2_preproc_bench["communication"]) {
        circ2_preproc_bytes += val.get<int64_t>();
    }
    size_t circ2_online_bytes = 0;
    for (const auto& val : circ2_online_bench["communication"]) {
        circ2_online_bytes += val.get<int64_t>();
    }

    // Compute estimated cost with log factor
    double estimated_preproc_time = circ1_preproc_time * log_factor + circ2_preproc_time;
    double estimated_online_time = circ1_online_time * log_factor + circ2_online_time;
    double estimated_total_time = estimated_preproc_time + estimated_online_time;

    size_t estimated_preproc_bytes = circ1_preproc_bytes * log_factor + circ2_preproc_bytes;
    size_t estimated_online_bytes = circ1_online_bytes * log_factor + circ2_online_bytes;
    size_t estimated_total_bytes = estimated_preproc_bytes + estimated_online_bytes;

    std::cout << "\n=== RESULTS ===" << std::endl;
    std::cout << "Circuit 1 (two parallel sorts):" << std::endl;
    std::cout << "  preproc time: " << circ1_preproc_time << " ms" << std::endl;
    std::cout << "  preproc sent: " << circ1_preproc_bytes << " bytes" << std::endl;
    std::cout << "  online time: " << circ1_online_time << " ms" << std::endl;
    std::cout << "  online sent: " << circ1_online_bytes << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "Circuit 2 (two shuffles + two reconstructions):" << std::endl;
    std::cout << "  preproc time: " << circ2_preproc_time << " ms" << std::endl;
    std::cout << "  preproc sent: " << circ2_preproc_bytes << " bytes" << std::endl;
    std::cout << "  online time: " << circ2_online_time << " ms" << std::endl;
    std::cout << "  online sent: " << circ2_online_bytes << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "Estimated cost = Circuit1 * ceil(log_2(" << num_vertices << ")) + Circuit2" << std::endl;
    std::cout << "               = Circuit1 * " << log_factor << " + Circuit2" << std::endl;
    std::cout << "preproc time: " << estimated_preproc_time << " ms" << std::endl;
    std::cout << "preproc sent: " << estimated_preproc_bytes << " bytes" << std::endl;
    std::cout << "online time: " << estimated_online_time << " ms" << std::endl;
    std::cout << "online sent: " << estimated_online_bytes << " bytes" << std::endl;
    std::cout << "total time: " << estimated_total_time << " ms" << std::endl;
    std::cout << "total sent: " << estimated_total_bytes << " bytes" << std::endl;

    output_data["circuit1"] = {{"preproc_time_ms", circ1_preproc_time},
                                {"preproc_bytes", circ1_preproc_bytes},
                                {"online_time_ms", circ1_online_time},
                                {"online_bytes", circ1_online_bytes}};
    output_data["circuit2"] = {{"preproc_time_ms", circ2_preproc_time},
                                {"preproc_bytes", circ2_preproc_bytes},
                                {"online_time_ms", circ2_online_time},
                                {"online_bytes", circ2_online_bytes}};
    output_data["estimated_cost"] = {{"preproc_time_ms", estimated_preproc_time},
                                      {"preproc_bytes", estimated_preproc_bytes},
                                      {"online_time_ms", estimated_online_time},
                                      {"online_bytes", estimated_online_bytes},
                                      {"total_time_ms", estimated_total_time},
                                      {"total_bytes", estimated_total_bytes}};

    output_data["stats"] = {{"peak_virtual_memory", peakVirtualMemory()},
                            {"peak_resident_set_size", peakResidentSetSize()}};

    std::cout << "\n--- Statistics ---" << std::endl;
    for (const auto& [key, value] : output_data["stats"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    if (save_output) {
        saveJson(output_data, save_file);
    }
}

int main(int argc, char* argv[]) {
    auto prog_opts = programOptions();
    bpo::options_description cmdline("Graphiti init micro benchmark");
    cmdline.add(prog_opts);
    cmdline.add_options()
        ("config,c", bpo::value<std::string>(), "Configuration file")
        ("help,h", "Produce help message");

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
