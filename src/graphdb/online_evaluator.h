#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "../io/netmp.h"
#include "../utils/circuit.h"
#include "../utils/thread_pool.h"
#include "preproc.h"
#include "rand_gen_pool.h"
#include "sharing.h"
#include "../utils/types.h"

using namespace common::utils;

namespace graphdb {
  class OnlineEvaluator {
    int nP_;
    int id_;
    int latency_;  // Network latency in microseconds
    bool use_pking_;  // Use king party for reconstruction
    RandGenPool rgen_;
    std::shared_ptr<io::NetIOMP> network_;
    PreprocCircuit<Ring> preproc_;
    common::utils::LevelOrderedCircuit circ_;
    std::vector<Ring> wires_;
    std::shared_ptr<common::utils::ThreadPool> tpool_;

    // Helper function to reconstruct shares via king party or direct all-to-all
    static void reconstruct(int nP, int pid, std::shared_ptr<io::NetIOMP> network,
                           const std::vector<Ring>& shares_list,
                           std::vector<Ring>& reconstructed_list,
                           bool via_pking, int latency);

  public:
    OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                    PreprocCircuit<Ring> preproc,
                    common::utils::LevelOrderedCircuit circ,
                    int threads, int seed = 200, int latency = 100, bool use_pking = true);

    OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                    PreprocCircuit<Ring> preproc,
                    common::utils::LevelOrderedCircuit circ,
                    std::shared_ptr<common::utils::ThreadPool> tpool, int seed = 200, int latency = 100, bool use_pking = true);

    void setInputs(const std::unordered_map<common::utils::wire_t, Ring> &inputs);

    void setRandomInputs();

    // void evaluateGatesAtDepthPartySend(size_t depth, std::vector<Ring> &mult_vals);

    // void evaluateGatesAtDepthPartyRecv(size_t depth, std::vector<Ring> &mult_vals);
    
    void RecEvaluate(const std::vector<common::utils::FIn2Gate> &mult_gates);

    void multEvaluate(const std::vector<common::utils::FIn2Gate> &mult_gates);

    void eqzEvaluate(const std::vector<common::utils::FIn1Gate> &eqz_gates);

    void recEvaluate(const std::vector<common::utils::FIn1Gate> &rec_gates);

    void evaluateGatesAtDepth(size_t depth);

    void shuffleEvaluate(const std::vector<common::utils::SIMDOGate> &shuffle_gates);

    void permAndShEvaluate(const std::vector<common::utils::SIMDOGate> &permAndSh_gates);

    void compactEvaluate(const common::utils::SIMDOGate &compact_gate);

    void compactEvaluateParallel(const std::vector<common::utils::SIMDOGate> &compact_gates);

    void groupwiseIndexEvaluate(const common::utils::SIMDOGate &gi_gate);

    void groupwiseIndexEvaluateParallel(const std::vector<common::utils::SIMDOGate> &gi_gates);

    void groupwisePropagateEvaluate(const common::utils::SIMDOGate &gp_gate, int latency);

    void groupwisePropagateEvaluateParallel(const std::vector<common::utils::SIMDOGate> &gp_gates);

    void sortEvaluate(const std::vector<common::utils::SIMDOGate> &sort_gates);

    void publicPermEvaluate(const std::vector<common::utils::SIMDOGate> &public_perm_gates);

    void rewireEvaluate(const std::vector<common::utils::SIMDOGate> &rewire_gates);

    void deleteWiresEvaluate(const std::vector<common::utils::SIMDOGate> &delete_gates);

    std::vector<Ring> getOutputs();

    Ring reconstruct(AddShare<Ring> &shares);

    // Evaluate online phase for circuit
    std::vector<Ring> evaluateCircuit(const std::unordered_map<common::utils::wire_t, Ring> &inputs);
  };

}; // namespace graphdb
