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

common::utils::Circuit<Ring> generateDeleteCircuit(int nP, int pid, size_t num_verts, size_t num_edges, size_t num_deletes) {
    
    size_t vec_size = num_verts + num_edges;
    
    std::cout << "Generating delete circuit with:" << std::endl;
    std::cout << "  num_verts=" << num_verts << std::endl;
    std::cout << "  num_edges=" << num_edges << std::endl;
    std::cout << "  vec_size=" << vec_size << " (num_verts + num_edges)" << std::endl;
    std::cout << "  num_deletes=" << num_deletes << std::endl;
    std::cout << "  Total daglist wires: " << (vec_size * 7) << " (7 fields per entry)" << std::endl;
    std::cout << "  Total delete wires: " << (num_deletes * 2) << " (src, dst per delete)" << std::endl;
    
    common::utils::Circuit<Ring> circ;

    // Input: vec_size * 7 wires for daglist entries (src, dst, isV, data, sigs, sigv, sigd)
    std::vector<common::utils::wire_t> daglist_src(vec_size);
    std::vector<common::utils::wire_t> daglist_dst(vec_size);
    std::vector<common::utils::wire_t> daglist_isV(vec_size);
    std::vector<common::utils::wire_t> daglist_data(vec_size);
    std::vector<common::utils::wire_t> daglist_sigs(vec_size);
    std::vector<common::utils::wire_t> daglist_sigv(vec_size);
    std::vector<common::utils::wire_t> daglist_sigd(vec_size);
    
    for (size_t i = 0; i < vec_size; ++i) {
        daglist_src[i] = circ.newInputWire();
        daglist_dst[i] = circ.newInputWire();
        daglist_isV[i] = circ.newInputWire();
        daglist_data[i] = circ.newInputWire();
        daglist_sigs[i] = circ.newInputWire();
        daglist_sigv[i] = circ.newInputWire();
        daglist_sigd[i] = circ.newInputWire();
    }

    // Input: num_deletes * 2 wires for delete entries (src, dst to match)
    std::vector<common::utils::wire_t> delete_src(num_deletes);
    std::vector<common::utils::wire_t> delete_dst(num_deletes);
    
    for (size_t i = 0; i < num_deletes; ++i) {
        delete_src[i] = circ.newInputWire();
        delete_dst[i] = circ.newInputWire();
    }

    // Create deletion tags for each daglist entry
    // tag[i] = 1 if entry i matches any delete entry, 0 otherwise
    std::vector<common::utils::wire_t> deletion_tags(vec_size);
    
    for (size_t i = 0; i < vec_size; ++i) {
        // Initialize tag as 0
        auto zero = circ.addGate(common::utils::GateType::kSub, daglist_src[i], daglist_src[i]);
        deletion_tags[i] = zero;
        
        // For each delete entry, check if it matches this daglist entry
        for (size_t d = 0; d < num_deletes; ++d) {
            // Equality check 1: delete.src == daglist.src
            auto diff_src = circ.addGate(common::utils::GateType::kSub, delete_src[d], daglist_src[i]);
            auto eq_src = circ.addGate(common::utils::GateType::kEqz, diff_src);
            
            // Equality check 2: delete.dst == daglist.dst
            auto diff_dst = circ.addGate(common::utils::GateType::kSub, delete_dst[d], daglist_dst[i]);
            auto eq_dst = circ.addGate(common::utils::GateType::kEqz, diff_dst);
            
            // Match = eq_src * eq_dst
            auto match = circ.addGate(common::utils::GateType::kMul, eq_src, eq_dst);
            
            // Update tag: tag = tag OR match (approximated as tag = tag + match - tag*match)
            // Or more simply: tag = max(tag, match), but we'll use tag = tag + match for simplicity
            deletion_tags[i] = circ.addGate(common::utils::GateType::kAdd, deletion_tags[i], match);
        }
    }

    // Generate identity permutation for compaction
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

    // Prepare payload vectors for compaction (7 fields per entry)
    std::vector<std::vector<common::utils::wire_t>> payload_vectors = {
        daglist_src,
        daglist_dst,
        daglist_isV,
        daglist_data,
        daglist_sigs,
        daglist_sigv,
        daglist_sigd
    };

    // Perform compaction: remove entries marked for deletion
    // Note: deletion_tags should be 1 for entries to KEEP, not delete (inverse of match)
    // So we need to invert deletion_tags first
    std::vector<common::utils::wire_t> keep_tags(vec_size);
    auto one = circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                     circ.addGate(common::utils::GateType::kSub, deletion_tags[0], deletion_tags[0]),
                                     static_cast<Ring>(1));
    for (size_t i = 0; i < vec_size; ++i) {
        keep_tags[i] = circ.addGate(common::utils::GateType::kSub, one, deletion_tags[i]);
    }

    auto [tags_compacted, payloads_compacted, labels] = 
        circ.addCompactionSubcircuit(keep_tags, payload_vectors, permutations, pid);

    // After compaction, we have fewer valid entries (vec_size - num_deletes approximately)
    // Extract the compacted daglist for further processing
    size_t remaining_entries = vec_size; // In practice this would be vec_size - num_deletes
    
    // Extract compacted src values for sorting (from first payload)
    std::vector<common::utils::wire_t> compacted_src_values = payloads_compacted[0];

    // For the remaining entries after deletion, perform two parallel sort subcircuits
    // Generate identity permutation for sorting
    std::vector<std::vector<int>> sort_permutations;
    std::vector<int> sort_identity_perm(remaining_entries);
    for (size_t i = 0; i < remaining_entries; ++i) {
        sort_identity_perm[i] = i;
    }
    sort_permutations.push_back(sort_identity_perm);
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            sort_permutations.push_back(sort_identity_perm);
        }
    }

    // First sort subcircuit on compacted daglist
    if (compacted_src_values.size() > 0) {
        auto sorted_indices_1 = circ.addSortSubcircuit(compacted_src_values, sort_permutations, pid);
        
        // Second parallel sort subcircuit (using same input)
        auto sorted_indices_2 = circ.addSortSubcircuit(compacted_src_values, sort_permutations, pid);
        
        // Set outputs: deletion tags + compacted fields + sorted indices from both sort gates
        for (size_t i = 0; i < vec_size; ++i) {
            circ.setAsOutput(deletion_tags[i]);
        }
        
        // Output all compacted payload fields
        for (const auto& payload : payloads_compacted) {
            for (const auto& wire : payload) {
                circ.setAsOutput(wire);
            }
        }
        
        for (const auto& idx : sorted_indices_1) {
            circ.setAsOutput(idx);
        }
        
        for (const auto& idx : sorted_indices_2) {
            circ.setAsOutput(idx);
        }
    } else {
        // Fallback: just output deletion tags and compacted fields
        for (size_t i = 0; i < vec_size; ++i) {
            circ.setAsOutput(deletion_tags[i]);
        }
        
        // Output all compacted payload fields
        for (const auto& payload : payloads_compacted) {
            for (const auto& wire : payload) {
                circ.setAsOutput(wire);
            }
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

    auto num_verts = opts["num-verts"].as<size_t>();
    auto num_edges = opts["num-edges"].as<size_t>();
    auto num_deletes = opts.count("num-deletes") ? opts["num-deletes"].as<size_t>() : static_cast<size_t>((num_verts + num_edges) * 0.05);
    auto nP = opts["num-parties"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto use_pking = opts["use-pking"].as<bool>();

    size_t vec_size = num_verts + num_edges;

    omp_set_nested(1);
    if (nP < 10) { omp_set_num_threads(nP); }
    else { omp_set_num_threads(10); }

    std::cout << "Starting delete benchmarks" << std::endl;

    std::string net_config = opts.count("net-config") ? opts["net-config"].as<std::string>() : "";
    std::shared_ptr<io::NetIOMP> network = createNetwork(pid, nP, latency, port,
                                                          opts["localhost"].as<bool>(),
                                                          net_config);

    // Increase socket buffer sizes to prevent deadlocks with large messages
    size_t expected_output_bytes = (vec_size * 7 + vec_size) * sizeof(Ring) * 2; // Compacted fields + tags + two sort gates
    int buffer_size = std::max(128 * 1024 * 1024, 
                                static_cast<int>(expected_output_bytes * (nP - 1) * 3));
    increaseSocketBuffers(network.get(), buffer_size);

    json output_data;
    output_data["details"] = {{"num_parties", nP},
                              {"num_verts", num_verts},
                              {"num_edges", num_edges},
                              {"num_deletes", num_deletes},
                              {"vec_size", vec_size},
                              {"total_daglist_wires", vec_size * 7},
                              {"total_delete_wires", num_deletes * 2},
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

    auto circ = generateDeleteCircuit(nP, pid, num_verts, num_edges, num_deletes).orderGatesByLevel();
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

    std::cout << "Setting random inputs" << std::endl;
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
    // Use random inputs for evaluation
    eval.setRandomInputs();
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
    
    // Print sample deletion tags and outputs
    std::cout << "\n=== OUTPUT: Deletion Tags (first 10 entries) ===" << std::endl;
    std::cout << "Entry | Tag" << std::endl;
    std::cout << "------|-----" << std::endl;

    size_t num_outputs = std::min(static_cast<size_t>(10), vec_size);
    for (size_t i = 0; i < num_outputs && i < outputs.size(); ++i) {
        std::cout << std::setw(5) << i << " | "
                  << std::setw(3) << outputs[i] << std::endl;
    }
    std::cout << std::endl;
    
    // Print sample compacted entries
    std::cout << "=== OUTPUT: Compacted Daglist (first 5 entries) ===" << std::endl;
    std::cout << "Entry | Src | Dst | isV | Data | sigs | sigv | sigd" << std::endl;
    std::cout << "------|-----|-----|-----|------|------|------|------" << std::endl;

    size_t compacted_start = vec_size;
    for (size_t i = 0; i < 5; ++i) {
        size_t idx = compacted_start + i * 7;
        if (idx + 6 < outputs.size()) {
            std::cout << std::setw(5) << i << " | "
                      << std::setw(3) << outputs[idx + 0] << " | "
                      << std::setw(3) << outputs[idx + 1] << " | "
                      << std::setw(3) << outputs[idx + 2] << " | "
                      << std::setw(4) << outputs[idx + 3] << " | "
                      << std::setw(4) << outputs[idx + 4] << " | "
                      << std::setw(4) << outputs[idx + 5] << " | "
                      << std::setw(4) << outputs[idx + 6] << std::endl;
        }
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
        ("num-verts", bpo::value<size_t>()->required(), "Number of vertices.")
        ("num-edges", bpo::value<size_t>()->required(), "Number of edges.")
        ("num-deletes", bpo::value<size_t>(), "Number of entries to delete (default: (num-verts + num-edges) * 0.05).")
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
    bpo::options_description cmdline("Benchmark delete circuit with compaction and sorting.");
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
// usage: ./../run.sh naive_delete --num-parties 2 --num-verts 100 --num-edges 200 --num-deletes 10
