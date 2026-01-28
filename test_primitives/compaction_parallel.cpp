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

common::utils::Circuit<Ring> generateParallelCompactionCircuit(int nP, int pid, size_t vec_size1, size_t vec_size2, size_t num_payloads) {

    std::cout << "Generating circuit with TWO parallel compaction gates:" << std::endl;
    std::cout << "  Gate 1: vec_size=" << vec_size1 << ", num_payloads=" << num_payloads << std::endl;
    std::cout << "  Gate 2: vec_size=" << vec_size2 << ", num_payloads=" << num_payloads << std::endl;
    
    common::utils::Circuit<Ring> circ;

    // Generate permutation for shuffle
    auto generatePermutation = [&](size_t size) {
        std::vector<std::vector<int>> permutation;
        std::vector<int> tmp_perm(size);
        for (size_t i = 0; i < size; ++i) {
            tmp_perm[i] = i;
        }
        permutation.push_back(tmp_perm);
        if (pid == 0) {
            for (int i = 1; i < nP; ++i) {
                permutation.push_back(tmp_perm);
            }
        }
        return permutation;
    };

    // First compaction gate with vec_size1
    std::vector<common::utils::wire_t> t_vector1(vec_size1);
    std::vector<std::vector<common::utils::wire_t>> p_vectors1(num_payloads, std::vector<common::utils::wire_t>(vec_size1));
    
    std::generate(t_vector1.begin(), t_vector1.end(), [&]() { return circ.newInputWire(); });
    for (size_t p = 0; p < num_payloads; ++p) {
        std::generate(p_vectors1[p].begin(), p_vectors1[p].end(), [&]() { return circ.newInputWire(); });
    }

    auto [t_compacted1, p_compacted1, label_reconstructed1] = circ.addCompactionSubcircuit(t_vector1, p_vectors1, generatePermutation(vec_size1), pid);

    // Second compaction gate with vec_size2
    std::vector<common::utils::wire_t> t_vector2(vec_size2);
    std::vector<std::vector<common::utils::wire_t>> p_vectors2(num_payloads, std::vector<common::utils::wire_t>(vec_size2));
    
    std::generate(t_vector2.begin(), t_vector2.end(), [&]() { return circ.newInputWire(); });
    for (size_t p = 0; p < num_payloads; ++p) {
        std::generate(p_vectors2[p].begin(), p_vectors2[p].end(), [&]() { return circ.newInputWire(); });
    }

    auto [t_compacted2, p_compacted2, label_reconstructed2] = circ.addCompactionSubcircuit(t_vector2, p_vectors2, generatePermutation(vec_size2), pid);
    
    // Set outputs from both compaction gates
    for (size_t i = 0; i < vec_size1; ++i) {
        circ.setAsOutput(t_compacted1[i]);
    }
    for (size_t p = 0; p < num_payloads; ++p) {
        for (size_t i = 0; i < vec_size1; ++i) {
            circ.setAsOutput(p_compacted1[p][i]);
        }
    }

    for (size_t i = 0; i < vec_size2; ++i) {
        circ.setAsOutput(t_compacted2[i]);
    }
    for (size_t p = 0; p < num_payloads; ++p) {
        for (size_t i = 0; i < vec_size2; ++i) {
            circ.setAsOutput(p_compacted2[p][i]);
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

    auto vec_size1 = opts["vec-size1"].as<size_t>();
    auto vec_size2 = opts["vec-size2"].as<size_t>();
    auto nP = opts["num-parties"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto num_payloads = opts["num-payloads"].as<size_t>();
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
                              {"vec_size1", vec_size1},
                              {"vec_size2", vec_size2},
                              {"num_payloads", num_payloads},
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

    auto circ = generateParallelCompactionCircuit(nP, pid, vec_size1, vec_size2, num_payloads).orderGatesByLevel();
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
    
    std::unordered_map<common::utils::wire_t, Ring> inputs;
    std::vector<common::utils::wire_t> input_wires;
    for (const auto& [wire, owner] : input_pid_map) {
        if (owner == static_cast<int>(pid)) {
            input_wires.push_back(wire);
        }
    }
    std::sort(input_wires.begin(), input_wires.end());
    
    // Set inputs for both gates
    size_t total_inputs_gate1 = vec_size1 * (1 + num_payloads);
    for (size_t i = 0; i < input_wires.size(); ++i) {
        auto wire = input_wires[i];
        if (i < total_inputs_gate1) {
            // Gate 1 inputs
            if (i < vec_size1) {
                inputs[wire] = static_cast<Ring>(rand() % 2);  // binary for t
            } else {
                inputs[wire] = static_cast<Ring>(rand() % 100);  // payload
            }
        } else {
            // Gate 2 inputs
            size_t gate2_idx = i - total_inputs_gate1;
            if (gate2_idx < vec_size2) {
                inputs[wire] = static_cast<Ring>(rand() % 2);  // binary for t
            } else {
                inputs[wire] = static_cast<Ring>(rand() % 100);  // payload
            }
        }
    }

    std::cout << "Total inputs set: " << inputs.size() << std::endl;
    
    // Print sample inputs for both gates (first 20 entries)
    if (input_wires.size() > 0) {
        std::cout << "\n=== GATE 1 Sample Inputs (first 20 elements) ===" << std::endl;
        std::cout << "t_vector1: ";
        for (size_t i = 0; i < std::min(std::min(size_t(20), vec_size1), input_wires.size()); ++i) {
            auto wire = input_wires[i];
            if (inputs.count(wire)) std::cout << inputs[wire] << " ";
        }
        if (vec_size1 > 20) std::cout << "...";
        std::cout << std::endl;
        std::cout << "p_vectors1[0]: ";
        for (size_t i = 0; i < std::min(size_t(20), vec_size1); ++i) {
            if (vec_size1 + i < input_wires.size()) {
                auto wire = input_wires[vec_size1 + i];
                if (inputs.count(wire)) std::cout << inputs[wire] << " ";
            }
        }
        if (vec_size1 > 20) std::cout << "...";
        std::cout << std::endl;
        
        size_t gate2_start = vec_size1 * (1 + num_payloads);
        if (gate2_start < input_wires.size()) {
            std::cout << "\n=== GATE 2 Sample Inputs (first 20 elements) ===" << std::endl;
            std::cout << "t_vector2: ";
            for (size_t i = 0; i < std::min(size_t(20), vec_size2); ++i) {
                if (gate2_start + i < input_wires.size()) {
                    auto wire = input_wires[gate2_start + i];
                    if (inputs.count(wire)) std::cout << inputs[wire] << " ";
                }
            }
            if (vec_size2 > 20) std::cout << "...";
            std::cout << std::endl;
            std::cout << "p_vectors2[0]: ";
            for (size_t i = 0; i < std::min(size_t(20), vec_size2); ++i) {
                if (gate2_start + vec_size2 + i < input_wires.size()) {
                    auto wire = input_wires[gate2_start + vec_size2 + i];
                    if (inputs.count(wire)) std::cout << inputs[wire] << " ";
                }
            }
            if (vec_size2 > 20) std::cout << "...";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    eval.setInputs(inputs);
    
    std::cout << "Starting online evaluation (testing parallel compaction gates)" << std::endl;
    StatsPoint online_start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
        network->sync();
        network->flush();
    }
    network->sync();
    StatsPoint online_end(*network);
    std::cout << "Online evaluation complete" << std::endl;

    auto outputs = eval.getOutputs();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print sample outputs for both gates (first 20 entries)
    std::vector<Ring> output_values = outputs;
    
    size_t gate1_outputs = vec_size1 * (1 + num_payloads);
    std::cout << "\n=== GATE 1 Sample Outputs (first 20 elements) ===" << std::endl;
    std::cout << "t_compacted1: ";
    for (size_t i = 0; i < std::min(size_t(20), vec_size1); ++i) {
        if (i < output_values.size()) std::cout << output_values[i] << " ";
    }
    if (vec_size1 > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "p_compacted1[0]: ";
    for (size_t i = 0; i < std::min(size_t(20), vec_size1); ++i) {
        size_t idx = vec_size1 + i;
        if (idx < output_values.size()) std::cout << output_values[idx] << " ";
    }
    if (vec_size1 > 20) std::cout << "...";
    std::cout << std::endl;
    
    std::cout << "\n=== GATE 2 Sample Outputs (first 20 elements) ===" << std::endl;
    std::cout << "t_compacted2: ";
    for (size_t i = 0; i < std::min(size_t(20), vec_size2); ++i) {
        size_t idx = gate1_outputs + i;
        if (idx < output_values.size()) std::cout << output_values[idx] << " ";
    }
    if (vec_size2 > 20) std::cout << "...";
    std::cout << std::endl;
    std::cout << "p_compacted2[0]: ";
    for (size_t i = 0; i < std::min(size_t(20), vec_size2); ++i) {
        size_t idx = gate1_outputs + vec_size2 + i;
        if (idx < output_values.size()) std::cout << output_values[idx] << " ";
    }
    if (vec_size2 > 20) std::cout << "...";
    std::cout << std::endl << std::endl;

    StatsPoint end(*network);

    auto preproc_rbench = preproc_end - preproc_start;
    auto online_rbench = online_end - online_start;
    auto total_rbench = end - start;
    output_data["benchmarks"].push_back(preproc_rbench);
    output_data["benchmarks"].push_back(online_rbench);
    output_data["benchmarks"].push_back(total_rbench);

    size_t online_bytes_sent = 0;
    for (const auto& val : online_rbench["communication"]) {
        online_bytes_sent += val.get<int64_t>();
    }

    std::cout << "preproc time: " << preproc_rbench["time"] << " ms" << std::endl;
    std::cout << "online time: " << online_rbench["time"] << " ms" << std::endl;
    std::cout << "online sent: " << online_bytes_sent << " bytes" << std::endl;
    std::cout << "total time: " << total_rbench["time"] << " ms" << std::endl;

    output_data["stats"] = {{"peak_virtual_memory", peakVirtualMemory()},
                            {"peak_resident_set_size", peakResidentSetSize()}};

    if (save_output) {
        saveJson(output_data, save_file);
    }
}

bpo::options_description programOptions() {
    bpo::options_description desc("Options for parallel compaction benchmark.");
    desc.add_options()
        ("num-parties,n", bpo::value<int>()->required(), "Number of parties.")
        ("vec-size1", bpo::value<size_t>()->required(), "Size of first compaction gate.")
        ("vec-size2", bpo::value<size_t>()->required(), "Size of second compaction gate.")
        ("num-payloads", bpo::value<size_t>()->default_value(1), "Number of payload vectors.")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads.")
        ("seed", bpo::value<size_t>()->default_value(200), "Random seed.")
        ("net-config", bpo::value<std::string>(), "Network config JSON file.")
        ("localhost", bpo::bool_switch(), "All parties on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port.")
        ("output,o", bpo::value<std::string>(), "Output file for benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of repetitions.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction.");
  return desc;
}

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark parallel compaction circuits.");
    cmdline.add(prog_opts);
    cmdline.add_options()("config,c", bpo::value<std::string>(), "Config file")("help,h", "Help message");
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
// usage: ./../run.sh compaction_parallel --num-parties 2 --vec-size1 10 --vec-size2 15