#include <io/netmp.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <omp.h>

#include "utils.h"

using namespace common::utils;
using json = nlohmann::json;
namespace bpo = boost::program_options;

void benchmark(const bpo::variables_map& opts) {

    bool save_output = false;
    std::string save_file;
    if (opts.count("output") != 0) {
        save_output = true;
        save_file = opts["output"].as<std::string>();
    }

    auto vec_size = opts["vec-size"].as<int>();
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
    
    // Increase socket buffer sizes for large vectors
    // Default Linux buffer is ~200KB, increase to handle large messages
    int buffer_size = std::max(16 * 1024 * 1024, vec_size * static_cast<int>(sizeof(Ring)) * 2); // At least 16MB or 2x vector size
    increaseSocketBuffers(network.get(), buffer_size);

    json output_data;
    output_data["details"] = {{"num_parties", nP},
                              {"vec_size", vec_size},
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

    // Generate random vector with unique seed per party
    emp::PRG prg(&emp::zero_block, seed + pid);
    std::vector<Ring> random_vector(vec_size);
    for (int i = 0; i < vec_size; ++i) {
        prg.random_data(&random_vector[i], sizeof(Ring));
    }
    
    std::cout << "\n=== Party " << pid << " Random Vector (first 20 elements) ===" << std::endl;
    for (int i = 0; i < std::min(20, vec_size); ++i) {
        std::cout << "  Vector[" << i << "] = " << random_vector[i] << std::endl;
    }
    std::cout << "=============================================\n" << std::endl;
    
    network->sync();
    
    std::cout << "Starting reconstruction phase (sending to all parties)" << std::endl;
    StatsPoint recon_start(*network);
    
    // Each party sends their vector to all other parties
    // With increased buffer sizes, this should not deadlock
    for (int dest = 0; dest < nP; ++dest) {
        if (dest != static_cast<int>(pid)) {
            network->send(dest, random_vector.data(), vec_size * sizeof(Ring));
        }
    }
    network->flush();
    
    // Receive vectors from all other parties
    std::vector<std::vector<Ring>> received_vectors(nP);
    received_vectors[pid] = random_vector;  // Own vector
    
    for (int src = 0; src < nP; ++src) {
        if (src != static_cast<int>(pid)) {
            received_vectors[src].resize(vec_size);
            network->recv(src, received_vectors[src].data(), vec_size * sizeof(Ring));
        }
    }
    
    // Reconstruct by summing all shares
    std::vector<Ring> reconstructed_vector(vec_size);
    for (int i = 0; i < vec_size; ++i) {
        reconstructed_vector[i] = 0;
        for (int p = 0; p < nP; ++p) {
            reconstructed_vector[i] += received_vectors[p][i];
        }
    }
    
    network->sync();
    StatsPoint recon_end(*network);
    std::cout << "Reconstruction complete" << std::endl;
    
    std::cout << "\n=== Reconstruction Results ===" << std::endl;
    std::cout << "Party " << pid << " reconstructed " << reconstructed_vector.size() << " values" << std::endl;
    std::cout << "First 20 reconstructed values:" << std::endl;
    for (int i = 0; i < std::min(20, vec_size); ++i) {
        std::cout << "  Reconstructed[" << i << "]: " << reconstructed_vector[i] << std::endl;
    }
    std::cout << "================================\n" << std::endl;

    StatsPoint end(*network);

    auto recon_rbench = recon_end - recon_start;
    auto total_rbench = end - start;
    output_data["benchmarks"].push_back(recon_rbench);
    output_data["benchmarks"].push_back(total_rbench);

    size_t recon_bytes_sent = 0;
    for (const auto& val : recon_rbench["communication"]) {
        recon_bytes_sent += val.get<int64_t>();
    }
    size_t total_bytes_sent = 0;
    for (const auto& val : total_rbench["communication"]) {
        total_bytes_sent += val.get<int64_t>();
    }

    std::cout << "--- Benchmark Results ---" << std::endl;
    std::cout << "reconstruction time: " << recon_rbench["time"] << " ms" << std::endl;
    std::cout << "reconstruction sent: " << recon_bytes_sent << " bytes" << std::endl;
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
        ("vec-size,v", bpo::value<int>()->default_value(100), "Size of the vector to create and reconstruct.")
        ("latency,l", bpo::value<double>()->required(), "Network latency in ms.")
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
    bpo::options_description cmdline("Benchmark vector reconstruction using direct network send/recv: create random vector and reconstruct it.");
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
