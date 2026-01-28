#pragma once

#include "../utils/circuit.h"
#include "sharing.h"
#include "../utils/types.h"
#include <unordered_map>

using namespace common::utils;

namespace graphdb {
// Preprocessed data for a gate.
template <class R>
struct PreprocGate {
  PreprocGate() = default;
  virtual ~PreprocGate() = default;
};

template <class R>
using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

template <class R>
struct PreprocInput : public PreprocGate<R> {
  // ID of party providing input on wire.
  int pid{};
  PreprocInput() = default;
  PreprocInput(int pid) 
      : PreprocGate<R>(), pid(pid) {}
  PreprocInput(const PreprocInput<R>& pregate) 
      : PreprocGate<R>(), pid(pregate.pid) {}
};

template <class R>
struct PreprocRecGate : public PreprocGate<R> {
  bool Pking = false;
  PreprocRecGate() = default;
  PreprocRecGate(bool Pking)
    : PreprocGate<R>(), Pking(Pking) {}
};

template <class R>
struct PreprocMultGate : public PreprocGate<R> {
  // Secret shared product of inputs masks.
  AddShare<R> triple_a; // Holds one beaver triple share of a random value a
  AddShare<R> triple_b; // Holds one beaver triple share of a random value b
  AddShare<R> triple_c; // Holds one beaver triple share of c=a*b
  PreprocMultGate() = default;
  PreprocMultGate(const AddShare<R>& triple_a,
                  const AddShare<R>& triple_b,
                  const AddShare<R>& triple_c)
      : PreprocGate<R>(), triple_a(triple_a), triple_b(triple_b), triple_c(triple_c) {}
};

template <class R>
struct PreprocEqzGate : public PreprocGate<R> {
  AddShare<R> share_r1;
  AddShare<R> share_r2;
  std::vector<AddShare<R>> share_r1_bits;
  std::vector<AddShare<R>> share_r2_bits;
  PreprocEqzGate() = default;
  PreprocEqzGate(const AddShare<R> &share_r1, const AddShare<R> &share_r2,
                 const std::vector<AddShare<R>> &share_r1_bits,
                 const std::vector<AddShare<R>> &share_r2_bits)
    : PreprocGate<R>(), share_r1(share_r1), share_r2(share_r2), share_r1_bits(share_r1_bits), share_r2_bits(share_r2_bits) {}
};

template <class R>
struct PreprocShuffleGate : public PreprocGate<R> {
  std::vector<AddShare<R>> r_1; // Randomly sampled vector
  std::vector<AddShare<R>> r_2; // Randomly sampled vector
  std::vector<R> delta_r; // Randomly sampled vector
  std::vector<int> pi; // Randomly sampled permutation using HP
  std::vector<std::vector<int>> tp_pi_all; // Randomly sampled permutations of all parties using HP

  PreprocShuffleGate() = default;
  PreprocShuffleGate(const std::vector<AddShare<R>>& r_1,
                     const std::vector<AddShare<R>>& r_2,
                     const std::vector<R>& delta_r,
                     const std::vector<int>& pi, const std::vector<std::vector<int>>& tp_pi_all)
      : PreprocGate<R>(), r_1(r_1), r_2(r_2), delta_r(delta_r),
        pi(pi), tp_pi_all(tp_pi_all) {}
};

// Preprocessing for Sort gate
// Contains 32 compaction operations (one for each bit from MSB to LSB)
// Plus final shuffle for hiding the permutation before reconstruction
template <class R>
struct PreprocSortGate : public PreprocGate<R> {
  // Array of 32 compaction preprocessing structures (one per bit)
  // Index 0 corresponds to MSB (bit 31), index 31 corresponds to LSB (bit 0)
  std::vector<std::vector<AddShare<R>>> shuffle_a;  // [32][vec_size]
  std::vector<std::vector<TPShare<R>>> shuffle_tp_a;
  std::vector<std::vector<AddShare<R>>> shuffle_b;
  std::vector<std::vector<TPShare<R>>> shuffle_tp_b;
  std::vector<std::vector<AddShare<R>>> shuffle_c;
  std::vector<std::vector<TPShare<R>>> shuffle_tp_c;
  std::vector<std::vector<Ring>> shuffle_delta;  // [32][vec_size]
  std::vector<std::vector<int>> shuffle_pi;  // [32][vec_size]
  std::vector<std::vector<std::vector<int>>> shuffle_tp_pi_all;  // [32][nP][vec_size]
  
  // Multiplication triples for each of the 32 compact operations
  std::vector<std::vector<AddShare<R>>> mult_triple_a;  // [32][vec_size]
  std::vector<std::vector<TPShare<R>>> mult_tp_triple_a;
  std::vector<std::vector<AddShare<R>>> mult_triple_b;
  std::vector<std::vector<TPShare<R>>> mult_tp_triple_b;
  std::vector<std::vector<AddShare<R>>> mult_triple_c;
  std::vector<std::vector<TPShare<R>>> mult_tp_triple_c;
  
  // Final shuffle preprocessing for hiding the permutation before reconstruction
  std::vector<AddShare<R>> final_shuffle_a;
  std::vector<TPShare<R>> final_shuffle_tp_a;
  std::vector<AddShare<R>> final_shuffle_b;
  std::vector<TPShare<R>> final_shuffle_tp_b;
  std::vector<AddShare<R>> final_shuffle_c;
  std::vector<TPShare<R>> final_shuffle_tp_c;
  std::vector<Ring> final_shuffle_delta;
  std::vector<int> final_shuffle_pi;
  std::vector<std::vector<int>> final_shuffle_tp_pi_all;
  
  PreprocSortGate() = default;
};


template <class R>
struct PreprocPublicPermGate : public PreprocGate<R> {
  size_t vec_size;       // Size of position map and each payload vector
  size_t num_payloads;   // Number of payload vectors
  
  PreprocPublicPermGate() = default;
  PreprocPublicPermGate(size_t vec_size, size_t num_payloads)
      : PreprocGate<R>(), vec_size(vec_size), num_payloads(num_payloads) {}
};

template <class R>
struct PreprocRewireGate : public PreprocGate<R> {
  size_t vec_size;       // Size of position map and each payload vector
  size_t num_payloads;   // Number of payload vectors
  
  PreprocRewireGate() = default;
  PreprocRewireGate(size_t vec_size, size_t num_payloads)
      : PreprocGate<R>(), vec_size(vec_size), num_payloads(num_payloads) {}
};

// Preprocessing for Delete Wires gate
// Takes a delete mask `del` and multiple payload vectors
// Indices where `del[i] == 1` will be removed from the payloads
// 
// Gate logic:
//   1. Shuffle del and all payload vectors together
//   2. Reconstruct del to reveal which positions should be deleted
//   3. Compact/remove indices where del == 1 from all payloads
//
// This gate requires preprocessing for:
//   - Shuffle operation (for del + payloads)
//   - Reconstruction of del
template <class R>
struct PreprocDeleteWiresGate : public PreprocGate<R> {
  // Preprocessing for shuffle operation (del + all payloads)
  std::vector<AddShare<R>> shuffle_a;
  std::vector<TPShare<R>> shuffle_tp_a;
  std::vector<AddShare<R>> shuffle_b;
  std::vector<TPShare<R>> shuffle_tp_b;
  std::vector<AddShare<R>> shuffle_c;
  std::vector<TPShare<R>> shuffle_tp_c;
  std::vector<Ring> shuffle_delta;
  std::vector<int> shuffle_pi;
  std::vector<std::vector<int>> shuffle_tp_pi_all;
  
  PreprocDeleteWiresGate() = default;
  PreprocDeleteWiresGate(const std::vector<AddShare<R>>& shuffle_a, const std::vector<TPShare<R>>& shuffle_tp_a,
                         const std::vector<AddShare<R>>& shuffle_b, const std::vector<TPShare<R>>& shuffle_tp_b,
                         const std::vector<AddShare<R>>& shuffle_c, const std::vector<TPShare<R>>& shuffle_tp_c,
                         const std::vector<R>& shuffle_delta, const std::vector<int>& shuffle_pi,
                         const std::vector<std::vector<int>>& shuffle_tp_pi_all)
      : PreprocGate<R>(), shuffle_a(shuffle_a), shuffle_tp_a(shuffle_tp_a), shuffle_b(shuffle_b), shuffle_tp_b(shuffle_tp_b),
        shuffle_c(shuffle_c), shuffle_tp_c(shuffle_tp_c), shuffle_delta(shuffle_delta), shuffle_pi(shuffle_pi),
        shuffle_tp_pi_all(shuffle_tp_pi_all) {}
};

// Preprocessed data for the circuit.
template <class R>
struct PreprocCircuit {
  std::unordered_map<wire_t, preprocg_ptr_t<R>> gates;
  PreprocCircuit() = default;
};
};  // namespace graphdb
