#pragma once

#include <algorithm>
#include <array>
#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "helpers.h"
#include "types.h"

namespace common::utils {

using wire_t = size_t;

enum GateType {
  kInp,
  kRec,
  kAdd,
  kMul,
  kSub,
  kEqz,
  kConstAdd,
  kConstMul,
  kShuffle,
  kPermAndSh,
  kPublicPerm,
  kDeleteWires,
  kSort,
  kRewire,
  kInvalid,
  NumGates
};

std::ostream& operator<<(std::ostream& os, GateType type);

// Gates represent primitive operations.
// All gates have one output.
struct Gate {
  GateType type{GateType::kInvalid};
  int owner;
  wire_t out;
  std::vector<wire_t> outs;
  std::vector<std::vector<wire_t>> multi_outs;

  Gate() = default;
  Gate(GateType type, wire_t out);
  Gate(GateType type, int owner, wire_t out, std::vector<wire_t> outs);
  Gate(GateType type, int owner, wire_t out, std::vector<std::vector<wire_t>> multi_outs);

  virtual ~Gate() = default;
};

// Represents a gate with fan-in 2.
struct FIn2Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};

  FIn2Gate() = default;
  FIn2Gate(GateType type, wire_t in1, wire_t in2, wire_t out);
};

struct FIn3Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};
  wire_t in3{0};

  FIn3Gate() = default;
  FIn3Gate(GateType type, wire_t in1, wire_t in2, wire_t in3, wire_t out);
};

struct FIn4Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};
  wire_t in3{0};
  wire_t in4{0};

  FIn4Gate() = default;
  FIn4Gate(GateType type, wire_t in1, wire_t in2, wire_t in3, wire_t in4, wire_t out);
};

// Represents a gate with fan-in 1.
struct FIn1Gate : public Gate {
  wire_t in{0};

  FIn1Gate() = default;
  FIn1Gate(GateType type, wire_t in, wire_t out);
};

// Represents a gate used to denote SIMD operations.
// These type is used to represent operations that take vectors of inputs but
// might not necessarily be SIMD e.g., dot product.
struct SIMDGate : public Gate {
  std::vector<wire_t> in1{0};
  std::vector<wire_t> in2{0};

  SIMDGate() = default;
  SIMDGate(GateType type, std::vector<wire_t> in1, std::vector<wire_t> in2, wire_t out);
};

// Represents a gate used to denote SIMD operations.
// These type is used to represent operations that take vectors of inputs and give vector of output but
// might not necessarily be SIMD e.g., shuffle, permute+share.
struct SIMDOGate : public Gate {
  std::vector<wire_t> in{0};
  std::vector<std::vector<int>> permutation{0};
  size_t vec_size{0};  // Metadata: size of vector for SIMD operations
  int inv{0};  // For rewire gate: if 1, perform inverse rewiring

  SIMDOGate() = default;
  SIMDOGate(GateType type, int owner, std::vector<wire_t> in, std::vector<wire_t> out, std::vector<std::vector<int>> permutation, size_t vec_size = 0, int inv = 0);
};

// Represents a gate used to denote SIMD operations.
// These type is used to represent operations that take vectors of inputs and give 2D vector of output but
// might not necessarily be SIMD e.g., amortized permute+share.
struct SIMDMOGate : public Gate {
  std::vector<wire_t> in{0};
  std::vector<std::vector<int>> permutation{0};

  SIMDMOGate() = default;
  SIMDMOGate(GateType type, int owner, std::vector<wire_t> in, std::vector<std::vector<wire_t>> multi_outs,
             std::vector<std::vector<int>> permutation);
};

// Represents gates where one input is a constant.
template <class R>
struct ConstOpGate : public Gate {
  wire_t in{0};
  R cval;

  ConstOpGate() = default;
  ConstOpGate(GateType type, wire_t in, R cval, wire_t out)
      : Gate(type, out), in(in), cval(std::move(cval)) {}
};

using gate_ptr_t = std::shared_ptr<Gate>;

// Gates ordered by multiplicative depth.
//
// Addition gates are not considered to increase the depth.
// Moreover, if gates_by_level[l][i]'s output is input to gates_by_level[l][j]
// then i < j.
struct LevelOrderedCircuit {
  size_t num_gates;
  size_t num_wires;
  std::array<uint64_t, GateType::NumGates> count;
  std::vector<wire_t> outputs;
  std::vector<std::vector<gate_ptr_t>> gates_by_level;

  friend std::ostream& operator<<(std::ostream& os, const LevelOrderedCircuit& circ);
};

// Represents an arithmetic circuit.
template <class R>
class Circuit {
  std::vector<wire_t> outputs_;
  std::vector<gate_ptr_t> gates_;
  size_t num_wires;

  bool isWireValid(wire_t wid) { return wid < num_wires; }

 public:
  Circuit() : num_wires(0) {}

  // Methods to manually build a circuit.
  wire_t newInputWire() {
    wire_t wid = num_wires;
    gates_.push_back(std::make_shared<Gate>(GateType::kInp, wid));
    num_wires += 1;
    return wid;
  }

  void setAsOutput(wire_t wid) {
    if (!isWireValid(wid)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    outputs_.push_back(wid);
  }

  // Function to add a gate with fan-in 2.
  wire_t addGate(GateType type, wire_t input1, wire_t input2) {
    if (type != GateType::kAdd && type != GateType::kMul &&
        type != GateType::kSub) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(input1) || !isWireValid(input2)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = num_wires;
    gates_.push_back(std::make_shared<FIn2Gate>(type, input1, input2, output));
    num_wires += 1;

    return output;
  }

  // Function to add a gate with one input from a wire and a second constant
  // input.
  wire_t addConstOpGate(GateType type, wire_t wid, R cval) {
    if (type != kConstAdd && type != kConstMul) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(wid)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = num_wires;
    gates_.push_back(std::make_shared<ConstOpGate<R>>(type, wid, cval, output));
    num_wires += 1;

    return output;
  }

  // Function to add a single input gate.
  wire_t addGate(GateType type, wire_t input) {
    if (type != GateType::kEqz && type != GateType::kRec) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(input)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = num_wires;
    gates_.push_back(std::make_shared<FIn1Gate>(type, input, output));
    num_wires += 1;

    return output;
  }

  // Function to add a multiple in + out gate.
  std::vector<wire_t> addMGate(GateType type, const std::vector<wire_t>& input, const std::vector<std::vector<int>> &permutation,
                               int owner = 0) {
    if (type != GateType::kShuffle && type != GateType::kPermAndSh) {
      throw std::invalid_argument("Invalid gate type.");
    }

    for (size_t i = 0; i < input.size(); i++) {
      if (!isWireValid(input[i])) {
        throw std::invalid_argument("Invalid wire ID.");
      }
    }

    if (permutation.size() == 0) {
      throw std::invalid_argument("No permutation passed.");
    }

    for (size_t i = 0; i < permutation.size(); ++i) {
      if (input.size() != permutation[i].size()) {
        throw std::invalid_argument("Permutation size mismatch.");
      }
    }

    std::vector<wire_t> output(input.size());
    for (int i = 0; i < input.size(); i++) {
      output[i] = i + num_wires;
    }
    gates_.push_back(std::make_shared<SIMDOGate>(type, owner, input, output, permutation));
    num_wires += input.size();
    return output;
  }

  std::vector<wire_t> addConstOpMGate(GateType type, const std::vector<wire_t>& input, const std::vector<int> &permutation) {
    if (type != GateType::kPublicPerm) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (input.size() != permutation.size()) {
      throw std::invalid_argument("Permutation size mismatch.");
    }

    for (size_t i = 0; i < input.size(); i++) {
      if (!isWireValid(input[i])) {
        throw std::invalid_argument("Invalid wire ID.");
      }
    }

    std::vector<std::vector<int>> permutation_wrapper(1);
    permutation_wrapper[0] = std::move(permutation);

    std::vector<wire_t> output(input.size());
    for (int i = 0; i < input.size(); i++) {
      output[i] = i + num_wires;
    }
    gates_.push_back(std::make_shared<SIMDOGate>(type, 0, input, output, permutation_wrapper));
    num_wires += input.size();
    return output;
  }

    // Rewire gate that applies a public permutation based on a position map
  // Takes a position map vector and any number of payload vectors as input
  // Outputs permuted payload wires based on the position map
  // 
  // For each position i: if position_map[i] = idx_perm, then output[idx_perm] = payload[i]
  //
  // Input format: [pos_map_0, ..., pos_map_n, p1_0, ..., p1_n, p2_0, ..., p2_n, ...]
  // Output format: [p1_out_0, ..., p1_out_n, p2_out_0, ..., p2_out_n, ...]
  std::vector<std::vector<wire_t>> addPublicPerm(
      const std::vector<wire_t>& public_Perm,
      const std::vector<std::vector<wire_t>>& payload_vectors,
      int inv = 0) {
    
    size_t vec_size = public_Perm.size();
    size_t num_payloads = payload_vectors.size();
    
    if (num_payloads == 0) {
      throw std::invalid_argument("At least one payload vector is required.");
    }
    
    // Validate all payload vectors have same size as public_perm
    for (size_t p = 0; p < num_payloads; ++p) {
      if (payload_vectors[p].size() != vec_size) {
        throw std::invalid_argument("All payload vectors must have the same size as public_perm.");
      }
    }
    
    // Validate input wires
    for (size_t i = 0; i < vec_size; i++) {
      if (!isWireValid(public_Perm[i])) {
        throw std::invalid_argument("Invalid wire ID in position_map.");
      }
      for (size_t p = 0; p < num_payloads; ++p) {
        if (!isWireValid(payload_vectors[p][i])) {
          throw std::invalid_argument("Invalid wire ID in payload vector.");
        }
      }
    }
    
    // Create output wires: vec_size for each payload
    std::vector<std::vector<wire_t>> payload_outputs(num_payloads, std::vector<wire_t>(vec_size));
    
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        payload_outputs[p][i] = num_wires + vec_size * p + i;
      }
    }
    
    // Create input vector [pos_map_0, ..., pos_map_n, p1_0, ..., p1_n, p2_0, ...]
    std::vector<wire_t> input((1 + num_payloads) * vec_size);
    for (size_t i = 0; i < vec_size; i++) {
      input[i] = public_Perm[i];
    }
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        input[vec_size * (p + 1) + i] = payload_vectors[p][i];
      }
    }
    
    // Create output vector for the gate
    std::vector<wire_t> output(num_payloads * vec_size);
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        output[vec_size * p + i] = payload_outputs[p][i];
      }
    }
    
    // Empty permutation vector since the permutation is determined at runtime from position_map
    std::vector<std::vector<int>> empty_permutation;
    // Pass vec_size and inv as metadata to the gate constructor
    gates_.push_back(std::make_shared<SIMDOGate>(GateType::kPublicPerm, 0, input, output, empty_permutation, vec_size, inv));
    num_wires += num_payloads * vec_size;
    
    return payload_outputs;
  }

  // Rewire gate that applies a public permutation based on a position map
  // Takes a position map vector and any number of payload vectors as input
  // Outputs permuted payload wires based on the position map
  // 
  // For each position i: if position_map[i] = idx_perm, then output[idx_perm] = payload[i]
  //
  // Input format: [pos_map_0, ..., pos_map_n, p1_0, ..., p1_n, p2_0, ..., p2_n, ...]
  // Output format: [p1_out_0, ..., p1_out_n, p2_out_0, ..., p2_out_n, ...]
  std::vector<std::vector<wire_t>> addRewireGate(
      const std::vector<wire_t>& position_map,
      const std::vector<std::vector<wire_t>>& payload_vectors,
      int inv = 0) {
    
    size_t vec_size = position_map.size();
    size_t num_payloads = payload_vectors.size();
    
    if (num_payloads == 0) {
      throw std::invalid_argument("At least one payload vector is required.");
    }
    
    // Validate all payload vectors have same size as position_map
    for (size_t p = 0; p < num_payloads; ++p) {
      if (payload_vectors[p].size() != vec_size) {
        throw std::invalid_argument("All payload vectors must have the same size as position_map.");
      }
    }
    
    // Validate input wires
    for (size_t i = 0; i < vec_size; i++) {
      if (!isWireValid(position_map[i])) {
        throw std::invalid_argument("Invalid wire ID in position_map.");
      }
      for (size_t p = 0; p < num_payloads; ++p) {
        if (!isWireValid(payload_vectors[p][i])) {
          throw std::invalid_argument("Invalid wire ID in payload vector.");
        }
      }
    }
    
    // Create output wires: vec_size for each payload
    std::vector<std::vector<wire_t>> payload_outputs(num_payloads, std::vector<wire_t>(vec_size));
    
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        payload_outputs[p][i] = num_wires + vec_size * p + i;
      }
    }
    
    // Create input vector [pos_map_0, ..., pos_map_n, p1_0, ..., p1_n, p2_0, ...]
    std::vector<wire_t> input((1 + num_payloads) * vec_size);
    for (size_t i = 0; i < vec_size; i++) {
      input[i] = position_map[i];
    }
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        input[vec_size * (p + 1) + i] = payload_vectors[p][i];
      }
    }
    
    // Create output vector for the gate
    std::vector<wire_t> output(num_payloads * vec_size);
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        output[vec_size * p + i] = payload_outputs[p][i];
      }
    }
    
    // Empty permutation vector since the permutation is determined at runtime from position_map
    std::vector<std::vector<int>> empty_permutation;
    // Pass vec_size and inv as metadata to the gate constructor
    gates_.push_back(std::make_shared<SIMDOGate>(GateType::kRewire, 0, input, output, empty_permutation, vec_size, inv));
    num_wires += num_payloads * vec_size;
    
    return payload_outputs;
  }

  // Delete Wires gate that deletes specified wires and keeps others
  std::pair<wire_t, std::vector<std::vector<wire_t>>> addDeleteWiresGate(
      const std::vector<wire_t>& del_vector,
      const std::vector<std::vector<wire_t>>& payload_vectors,
      const std::vector<std::vector<int>>& permutation) {
    
    if (permutation.size() == 0) {
      throw std::invalid_argument("No permutation passed.");
    }

    size_t vec_size = permutation[0].size();
    size_t num_payloads = payload_vectors.size();

    if (num_payloads == 0) {
      throw std::invalid_argument("At least one payload vector is required.");
    }

    if (del_vector.size() != vec_size) {
      throw std::invalid_argument("del_vector size must match permutation size.");
    }

    for (size_t i = 0; i < permutation.size(); ++i) {
      if (permutation[i].size() != vec_size) {
        throw std::invalid_argument("Permutation size mismatch.");
      }
    }

    for (size_t p = 0; p < num_payloads; ++p) {
      if (payload_vectors[p].size() != vec_size) {
        throw std::invalid_argument("All payload vectors must have the same size as permutation.");
      }
    }

    for (size_t i = 0; i < vec_size; i++) {
      if (!isWireValid(del_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in del_vector.");
      }
      for (size_t p = 0; p < num_payloads; ++p) {
        if (!isWireValid(payload_vectors[p][i])) {
          throw std::invalid_argument("Invalid wire ID in payload vector.");
        }
      }
    }

    // Allocate output wires: 1 wire for keep_indices, vec_size for each payload
    wire_t keep_indices_wire = num_wires;
    
    std::vector<std::vector<wire_t>> payload_outputs(num_payloads, std::vector<wire_t>(vec_size));
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        payload_outputs[p][i] = num_wires + 1 + vec_size * p + i;
      }
    }

    // Create input vector [del_0,...,del_n, p1_0,...,p1_n, p2_0,...]
    std::vector<wire_t> input((1 + num_payloads) * vec_size);
    for (size_t i = 0; i < vec_size; i++) {
      input[i] = del_vector[i];
    }
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        input[vec_size * (p + 1) + i] = payload_vectors[p][i];
      }
    }

    // Flatten outputs into a single vector for the gate constructor
    // Output format: [keep_indices_wire, p1_0,...,p1_n, p2_0,...,p2_n, ...]
    std::vector<wire_t> output(1 + num_payloads * vec_size);
    output[0] = keep_indices_wire;
    for (size_t p = 0; p < num_payloads; ++p) {
      for (size_t i = 0; i < vec_size; i++) {
        output[1 + vec_size * p + i] = payload_outputs[p][i];
      }
    }

    gates_.push_back(std::make_shared<SIMDOGate>(GateType::kDeleteWires, 0, input, output, permutation, vec_size));
    num_wires += 1 + num_payloads * vec_size;

    return {keep_indices_wire, payload_outputs};
  }


  // ============================================================================
  // SUBCIRCUIT: Permlist
  // ============================================================================
  std::vector<std::vector<wire_t>> addSubCircPermList(
      const std::vector<wire_t>& position_map_shares,
      const std::vector<std::vector<wire_t>>& payloads,
      std::vector<std::vector<int>> permutation) {

    size_t vec_size = position_map_shares.size();

    // Validate payload sizes (each payload must match position map length)
    for (size_t p = 0; p < payloads.size(); ++p) {
      if (payloads[p].size() != vec_size) {
        throw std::invalid_argument("All payload vectors must have the same size as position_map_shares");
      }
    }

    // Shuffle position map
    auto shuffled_position_map = addMGate(GateType::kShuffle, position_map_shares, permutation);

    // Shuffle each payload vector
    std::vector<std::vector<wire_t>> shuffled_payloads;
    shuffled_payloads.reserve(payloads.size());
    for (const auto& pl : payloads) {
      shuffled_payloads.push_back(addMGate(GateType::kShuffle, pl, permutation));
    }

    // Reconstruct position map
    std::vector<wire_t> position_map_reconstructed(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
      position_map_reconstructed[i] = addGate(GateType::kRec, shuffled_position_map[i]);
    }

    // Rewire shuffled payloads using reconstructed position map
    auto rewired_outputs = addRewireGate(position_map_reconstructed, shuffled_payloads);

    return rewired_outputs;
  }

  // ============================================================================
  // SUBCIRCUIT: Propagate
  // ============================================================================
  std::vector<wire_t> addSubCircPropagate(
      const std::vector<wire_t>& position_map_shares,
      const std::vector<wire_t>& data_values,
      size_t num_groups,
      std::vector<std::vector<int>> permutation,
      bool in = false) {
    
    size_t vec_size = position_map_shares.size();
    
    // Validate num_groups < vec_size
    if (num_groups >= vec_size) {
      throw std::invalid_argument("num_groups must be less than vec_size");
    }

    // Step 1: Compute differences for group boundaries
    // data_values'[0] = data_values[0]
    // data_values'[i] = data_values[i] - data_values[i-1] for i = 1 to num_groups-1
    // data_values'[i] = 0 for i >= num_groups
    std::vector<wire_t> data_values_diff(vec_size);
    if (in) {
      for (int i = static_cast<int>(num_groups) - 1; i >= 0; --i) {
        if (i == static_cast<int>(num_groups) - 1) {
          // data_values'[num_groups-1] = data_values[num_groups-1]
          data_values_diff[i] = data_values[i];
        } else {
          // For i from 0 to num_groups-2: compute differences
          // data_values'[i] = data_values[i] - data_values[i+1]
          data_values_diff[i] = addGate(GateType::kSub, data_values[i], data_values[i + 1]);
        }
      }
      // Set remaining elements to 0
      for (size_t i = num_groups; i < vec_size; ++i) {
        data_values_diff[i] = addGate(GateType::kSub, data_values[i], data_values[i]);
      }
    } else {
      for (size_t i = 0; i < vec_size; ++i) {
        if (i == 0) {
          // data_values'[0] = data_values[0]
          data_values_diff[i] = data_values[i];
        } else if (i < num_groups) {
          // For i from 1 to num_groups-1: compute differences
          // data_values'[i] = data_values[i] - data_values[i-1]
          data_values_diff[i] = addGate(GateType::kSub, data_values[i], data_values[i - 1]);
        } else {
          // For i >= num_groups: set to 0
          // data_values'[i] = 0
          data_values_diff[i] = addGate(GateType::kSub, data_values[i], data_values[i]);
        }
      }
    }

    // Step 2-4: Shuffle, reconstruct, and rewire using addSubCircPermList
    std::vector<std::vector<wire_t>> payloads = {data_values_diff};
    auto rewired_outputs = addSubCircPermList(position_map_shares, payloads, permutation);
    auto reordered_data = rewired_outputs[0];

    // Step 5: Compute prefix sum of reordered data
    std::vector<wire_t> prefix_sum(vec_size);
    if (in) {
      prefix_sum[vec_size - 1] = reordered_data[vec_size - 1];
      for (int i = static_cast<int>(vec_size) - 2; i >= 0; --i) {
        prefix_sum[i] = addGate(GateType::kAdd, prefix_sum[i + 1], reordered_data[i]);
      }
    } else {
      prefix_sum[0] = reordered_data[0];
      for (size_t i = 1; i < vec_size; ++i) {
        prefix_sum[i] = addGate(GateType::kAdd, prefix_sum[i - 1], reordered_data[i]);
      }
    }    

    return prefix_sum;
  }


  // ============================================================================
  // SUBCIRCUIT: Gather
  // ============================================================================
  std::vector<wire_t> addSubCircGather(
      const std::vector<wire_t>& position_map_shares,
      const std::vector<wire_t>& data_values,
      size_t num_groups,
      std::vector<std::vector<int>> permutation) {
    
    size_t vec_size = position_map_shares.size();
    
    // Validate num_groups < vec_size
    if (num_groups >= vec_size) {
      throw std::invalid_argument("num_groups must be less than vec_size");
    }

    // Step 1: Compute prefix sum of data values
    std::vector<wire_t> prefix_sum(vec_size);
    prefix_sum[0] = data_values[0];
    for (size_t i = 1; i < vec_size; ++i) {
      prefix_sum[i] = addGate(GateType::kAdd, prefix_sum[i - 1], data_values[i]);
    }

    // Step 2-4: Shuffle, reconstruct, and rewire using addSubCircPermList
    std::vector<std::vector<wire_t>> payloads = {prefix_sum, data_values};
    auto rewired_outputs = addSubCircPermList(position_map_shares, payloads, permutation);
    auto reordered_data = rewired_outputs[0];
    auto rewired_data_values = rewired_outputs[1];

    // Step 5: Compute differences
    // data_values'[0] = reordered_data[0]
    // data_values'[i] = reordered_data[i] - reordered_data[i-1] - data_values[i] for i = 1 to num_groups-1
    // data_values'[i] = 0 for i >= num_groups
    std::vector<wire_t> data_values_diff(vec_size);
    
    for (size_t i = 0; i < vec_size; ++i) {
      if (i == 0) {
        // data_values'[0] = reordered_data[0] - data_values[0]
        data_values_diff[i] = addGate(GateType::kSub, reordered_data[i], rewired_data_values[i]);
      } else if (i < num_groups) {
        // For i from 1 to num_groups-1: compute differences
        // data_values'[i] = reordered_data[i] - reordered_data[i-1] - data_values[i]
        data_values_diff[i] = addGate(GateType::kSub, reordered_data[i], reordered_data[i - 1]);
        data_values_diff[i] = addGate(GateType::kSub, data_values_diff[i], rewired_data_values[i]);
      } else {
        // For i >= num_groups: set to 0
        // data_values'[i] = 0
        data_values_diff[i] = 0;
      }
    }

    return data_values_diff;
  }

  // ============================================================================
  // SUBCIRCUIT: Compaction
  // ============================================================================
  // Add a compaction subcircuit that uses basic gates to perform compaction
  // Takes t (tags) and multiple payload vectors
  // Returns: (t_compacted, p_compacted, reconstructed_labels) - compacted versions of inputs
  std::tuple<std::vector<wire_t>, std::vector<std::vector<wire_t>>, std::vector<wire_t>> addCompactionSubcircuit(
      const std::vector<wire_t>& t_vector,
      const std::vector<std::vector<wire_t>>& p_vectors,
      const std::vector<std::vector<int>>& permutation,
      int pid) {
    
    size_t vec_size = t_vector.size();
    size_t num_payloads = p_vectors.size();
    
    if (num_payloads == 0) {
      throw std::invalid_argument("At least one payload vector is required.");
    }
    
    // Validate all payload vectors have same size as t_vector
    for (size_t p = 0; p < num_payloads; ++p) {
      if (p_vectors[p].size() != vec_size) {
        throw std::invalid_argument("All payload vectors must have the same size as t_vector.");
      }
    }
    
    // Validate input wires
    for (size_t i = 0; i < vec_size; i++) {
      if (!isWireValid(t_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in t_vector.");
      }
      for (size_t p = 0; p < num_payloads; ++p) {
        if (!isWireValid(p_vectors[p][i])) {
          throw std::invalid_argument("Invalid wire ID in payload vector.");
        }
      }
    }
    
    // Step 1: Compute prefix sums c1 and c0 using addition gates
    // c1[j] = sum of t[0..j] (cumulative count of 1's)
    // c0[j] = j - c1[j] (cumulative count of 0's, with adjustment for party 1)
    std::vector<wire_t> c1(vec_size);
    std::vector<wire_t> c0(vec_size);

    // Create constant wires that all parties can use
    wire_t zero_wire = addGate(GateType::kSub, t_vector[0], t_vector[0]);
    wire_t one_wire = addConstOpGate(GateType::kConstAdd, zero_wire, static_cast<R>(1));
    
    // Pre-create constant wires for all needed values to avoid creating them in the loop
    std::vector<wire_t> const_wires(vec_size + 1);
    const_wires[0] = zero_wire;
    const_wires[1] = one_wire;
    for (size_t j = 2; j <= vec_size; ++j) {
      const_wires[j] = addConstOpGate(GateType::kConstAdd, zero_wire, static_cast<R>(j));
    }
    
    // c1[0] = t[0]
    c1[0] = t_vector[0];
    
    // c0[0] = 1 - c1[0]
    c0[0] = addGate(GateType::kSub, one_wire, c1[0]);
    
    // Compute remaining prefix sums
    for (size_t j = 1; j < vec_size; ++j) {
      // c1[j] = c1[j-1] + t[j]
      c1[j] = addGate(GateType::kAdd, c1[j-1], t_vector[j]);
      
      // c0[j] = (j+1) - c1[j]
      c0[j] = addGate(GateType::kSub, const_wires[j + 1], c1[j]);
    }
    
    // Step 2: Compute label shares using multiplication gates
    // label[j] = (c0[j] + c1[N-1] - c1[j]) * (1 - t[j]) + c1[j] - 1
    std::vector<wire_t> label(vec_size);
    wire_t c1_last = c1[vec_size - 1];
    
    for (size_t j = 0; j < vec_size; ++j) {
      // diff_term = c0[j] + c1_last - c1[j]
      wire_t temp1 = addGate(GateType::kAdd, c0[j], c1_last);
      wire_t diff_term = addGate(GateType::kSub, temp1, c1[j]);
      
      // one_minus_t = 1 - t[j]
      wire_t one_minus_t = addGate(GateType::kSub, one_wire, t_vector[j]);
      
      // mult_result = diff_term * one_minus_t
      wire_t mult_result = addGate(GateType::kMul, diff_term, one_minus_t);
      
      // label[j] = mult_result + c1[j] - 1
      wire_t temp2 = addGate(GateType::kAdd, mult_result, c1[j]);
      label[j] = addGate(GateType::kSub, temp2, one_wire);
    }
    
    // Step 3: Shuffle all vectors (t, payloads, label) using the same permutation
    std::vector<wire_t> t_shuffled = addMGate(GateType::kShuffle, t_vector, permutation);
    
    std::vector<std::vector<wire_t>> p_shuffled(num_payloads);
    for (size_t p = 0; p < num_payloads; ++p) {
      p_shuffled[p] = addMGate(GateType::kShuffle, p_vectors[p], permutation);
    }
    
    std::vector<wire_t> label_shuffled = addMGate(GateType::kShuffle, label, permutation);
    
    // Step 4: Reconstruct label_shuffled to get the position map
    // Each element needs to be reconstructed individually
    std::vector<wire_t> label_reconstructed(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
      label_reconstructed[i] = addGate(GateType::kRec, label_shuffled[i]);
    }
    
    // Step 5: Apply rewire gate using reconstructed labels as position map
    // Rewire both t and all payload vectors
    std::vector<std::vector<wire_t>> rewire_inputs;
    rewire_inputs.push_back(t_shuffled);
    for (size_t p = 0; p < num_payloads; ++p) {
      rewire_inputs.push_back(p_shuffled[p]);
    }
    
    std::vector<std::vector<wire_t>> rewired_outputs = addRewireGate(label_reconstructed, rewire_inputs);
    
    // Extract t_compacted and p_compacted from rewired outputs
    std::vector<wire_t> t_compacted = rewired_outputs[0];
    std::vector<std::vector<wire_t>> p_compacted(num_payloads);
    for (size_t p = 0; p < num_payloads; ++p) {
      p_compacted[p] = rewired_outputs[p + 1];
    }
    
    return {t_compacted, p_compacted, label_reconstructed};
  }

  // ============================================================================
  // SUBCIRCUIT: Groupwise Index
  // ============================================================================
  std::tuple<std::vector<wire_t>, std::vector<wire_t>, std::vector<wire_t>> addGroupwiseIndexSubcircuit(
      const std::vector<wire_t>& key_vector,
      const std::vector<wire_t>& v_vector,
      const std::vector<std::vector<int>>& permutation,
      int pid) {
    
    size_t vec_size = key_vector.size();
    
    if (v_vector.size() != vec_size) {
      throw std::invalid_argument("Value vector must have the same size as key vector.");
    }
    
    // Validate input wires
    for (size_t i = 0; i < vec_size; i++) {
      if (!isWireValid(key_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in key_vector.");
      }
      if (!isWireValid(v_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in v_vector.");
      }
    }
    
    // Step 1: Initialize ind vector with sequential indices [0, 1, 2, ..., N-1]
    // Create these as constant wires (party 1 holds the values, others hold 0)
    wire_t zero_wire = addGate(GateType::kSub, key_vector[0], key_vector[0]);
    std::vector<wire_t> ind_vector(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
      ind_vector[i] = addConstOpGate(GateType::kConstAdd, zero_wire, static_cast<R>(i));
    }
    
    // Step 2: Compact (key, ind) using compaction subcircuit
    // This gives us key_c (compacted keys) and ind_c (compacted indices)
    std::vector<std::vector<wire_t>> ind_payload = {ind_vector};
    auto [key_c, ind_c_vec, label_c] = addCompactionSubcircuit(key_vector, ind_payload, permutation, pid);
    std::vector<wire_t> ind_c = ind_c_vec[0];

    // key_c[i] = ind_compacted[i-1] - ind_compacted[i]
    // key_c[0] = ind_compacted[0] (no predecessor)
    std::vector<wire_t> key_c_prime(vec_size);
    key_c_prime[0] = zero_wire;
    for (size_t i = 1; i < vec_size; ++i){
      key_c_prime[i] = addGate(GateType::kSub, ind_c[i-1], ind_c[i]);
    }

    // Step 3: Compute ind-diff[i] = key_c[i] * key_c_prime[i]
    // This multiplication zeros out indices for non-group-start positions
    std::vector<wire_t> ind_diff(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
      ind_diff[i] = addGate(GateType::kMul, key_c[i], key_c_prime[i]);
    }
    
    // Step 4: Reverse compact to redistribute ind-diff back to original positions
    // We use the label_c (reconstructed labels from compaction) to reverse the rewire
    // First, we need to reverse shuffle ind_diff, then apply inverse rewire
    
    auto ind_diff_inv_rewired = addRewireGate(label_c, {ind_diff}, 1);
    std::vector<wire_t> ind_diff_rev = ind_diff_inv_rewired[0];

    // Compute final indices via prefix sum + index i
    wire_t cumulative_sum = zero_wire;
    std::vector<wire_t> ind_out(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
      cumulative_sum = addGate(GateType::kAdd, cumulative_sum, ind_diff_rev[i]);
      ind_out[i] = addConstOpGate(GateType::kConstAdd, cumulative_sum, static_cast<R>(i));
    }
    
    return {ind_out, key_vector, v_vector};
  }

  // ============================================================================
  // SUBCIRCUIT: Groupwise Propagate
  // ============================================================================
  std::pair<std::vector<wire_t>, std::vector<wire_t>> addGroupwisePropagateSubcircuit(
      const std::vector<wire_t>& key1_vector,
      const std::vector<wire_t>& v1_vector,
      const std::vector<wire_t>& key2_vector,
      const std::vector<std::vector<int>>& permutation_t1,
      const std::vector<std::vector<int>>& permutation_t2,
      const std::vector<std::vector<int>>& permutation_rev,
      int pid) {
    
    size_t t1_vec_size = key1_vector.size();
    size_t t2_vec_size = key2_vector.size();
    
    if (v1_vector.size() != t1_vec_size) {
      throw std::invalid_argument("Value vector must have the same size as key1 vector.");
    }
    
    // Validate input wires
    for (size_t i = 0; i < t1_vec_size; i++) {
      if (!isWireValid(key1_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in key1_vector.");
      }
      if (!isWireValid(v1_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in v1_vector.");
      }
    }
    for (size_t i = 0; i < t2_vec_size; i++) {
      if (!isWireValid(key2_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in key2_vector.");
      }
    }
    
    // Step 1: Compact T1 (key1, v1) and T2 (key2) in parallel using compaction subcircuit
    std::vector<std::vector<wire_t>> t1_payloads = {v1_vector};
    auto [key1_compacted, v1_compacted_vec, label_reconstructed_t1] = addCompactionSubcircuit(key1_vector, t1_payloads, permutation_t1, pid);
    std::vector<wire_t> v1_compacted = v1_compacted_vec[0];
    
    // For T2, we need to compact key2. We use key2 itself as both the selector and payload
    std::vector<std::vector<wire_t>> t2_payloads = {key2_vector};
    auto [key2_compacted, key2_payload_compacted, label_reconstructed_t2] = addCompactionSubcircuit(key2_vector, t2_payloads, permutation_t2, pid);
    
    // Step 2: Compute differences: diff[i] = v1_compacted[i] - v1_compacted[i-1]
    std::vector<wire_t> diff(t1_vec_size);
    
    // diff[0] = v1_compacted[0] (no subtraction for first element)
    diff[0] = v1_compacted[0];
    
    for (size_t i = 1; i < t1_vec_size; ++i) {
      // diff[i] = v1_compacted[i] - v1_compacted[i-1]
      diff[i] = addGate(GateType::kSub, v1_compacted[i], v1_compacted[i-1]);
    }
    
    // Step 3: Multiply: key_times_diff[i] = diff[i] * key1_compacted[i]
    std::vector<wire_t> key_times_diff(t1_vec_size);
    
    for (size_t i = 0; i < t1_vec_size; ++i) {
      key_times_diff[i] = addGate(GateType::kMul, diff[i], key1_compacted[i]);
    }
    
    // Step 4a: Extend key_times_diff to match T2 size by padding with zeros
    std::vector<wire_t> extended_diff(t2_vec_size);
    wire_t zero_wire = addGate(GateType::kSub, diff[0], diff[0]);
    
    for (size_t i = 0; i < t2_vec_size; ++i) {
      if (i < t1_vec_size) {
        extended_diff[i] = key_times_diff[i];
      } else {
        extended_diff[i] = zero_wire;
      }
    }
    
    // Step 4b: To reverse the compaction permutation, we need to:
    // 1. First reverse the rewire (public permutation) that was applied in compaction
    // 2. Then reverse the shuffle
    
    auto diff_rewired = addRewireGate(label_reconstructed_t2, {extended_diff}, 1);
    std::vector<wire_t> diff_reverse_shuffled = diff_rewired[0];
    
    // Step 5: Propagate values using prefix sum
    // propagated[i] = sum(diff_reverse_shuffled[0..i])
    std::vector<wire_t> propagated_values(t2_vec_size);
    
    propagated_values[0] = diff_reverse_shuffled[0];
    
    for (size_t i = 1; i < t2_vec_size; ++i) {
      propagated_values[i] = addGate(GateType::kAdd, propagated_values[i-1], diff_reverse_shuffled[i]);
    }
    
    return {key2_vector, propagated_values};
  }

  // ============================================================================
  // SUBCIRCUIT: Delete Wires
  // ============================================================================
  // Add a delete wires subcircuit using compaction
  // Takes del_vector as tags and payload vectors
  // Returns: (keep_indices, p_outputs) - keep indices and compacted payloads
  std::pair<wire_t, std::vector<std::vector<wire_t>>> addDeleteWiresSubcircuit(
      const std::vector<wire_t>& del_vector,
      const std::vector<std::vector<wire_t>>& payload_vectors,
      const std::vector<std::vector<int>>& permutation,
      int pid) {
    
    size_t vec_size = del_vector.size();
    size_t num_payloads = payload_vectors.size();
    
    if (num_payloads == 0) {
      throw std::invalid_argument("At least one payload vector is required.");
    }
    
    // Validate all payload vectors have same size as del_vector
    for (size_t p = 0; p < num_payloads; ++p) {
      if (payload_vectors[p].size() != vec_size) {
        throw std::invalid_argument("All payload vectors must have the same size as del_vector.");
      }
    }
    
    // Validate input wires
    for (size_t i = 0; i < vec_size; i++) {
      if (!isWireValid(del_vector[i])) {
        throw std::invalid_argument("Invalid wire ID in del_vector.");
      }
      for (size_t p = 0; p < num_payloads; ++p) {
        if (!isWireValid(payload_vectors[p][i])) {
          throw std::invalid_argument("Invalid wire ID in payload vector.");
        }
      }
    }
    
    // Use compaction subcircuit with del_vector as tags
    auto [t_compacted, p_compacted, label_reconstructed] = addCompactionSubcircuit(del_vector, payload_vectors, permutation, pid);
    
    // Return first reconstructed label as keep_indices (single wire) and compacted payloads
    return {p_compacted};
  }

  // ============================================================================
  // SUBCIRCUIT: Sorting
  // ============================================================================
  // Add a sorting subcircuit using iterative compaction.
  //
  // Input: 1. A vector of wires with dimension {vec_size * input_size}
  //        2. A permutation matrix of size [num_parties][vec_size]
  //
  // Circuit Logic:
  // 1. Group wires into input_size many groups as follows:
  //    - Group 0: {wire[0], wire[input_size], wire[2*input_size], ...}
  //    - Group 1: {wire[1], wire[input_size+1], wire[2*input_size+1], ...}
  //    - Group i: {wire[i], wire[input_size+i], wire[2*input_size+i], ...}
  //
  // 2. Initialize index wires [1, 2, 3, ..., vec_size]
  //
  // 3. For input_size many iterations (outer loop):
  //    - Take the corresponding group from the last iteration
  //    - Apply compaction with all wires in groups 0 to current + index wires as payload
  //    - Update all groups with compacted outputs
  //
  // 4. Output: index wires from the final iteration
  //
  std::vector<wire_t> addSortSubcircuit(
      const std::vector<wire_t>& wires,
      const std::vector<std::vector<int>>& permutation,
      int pid) {
    
    // Validate inputs
    if (permutation.size() == 0) {
      throw std::invalid_argument("Permutation matrix is empty.");
    }
    
    size_t vec_size = permutation[0].size();
    
    if (wires.size() % vec_size != 0) {
      throw std::invalid_argument("Total number of wires must be divisible by vec_size.");
    }
    
    size_t input_size = wires.size() / vec_size;
    
    // Validate all wires
    for (size_t i = 0; i < wires.size(); ++i) {
      if (!isWireValid(wires[i])) {
        throw std::invalid_argument("Invalid wire ID in wires vector.");
      }
    }
    
    // Create initial zero wire
    wire_t zero_wire = addGate(GateType::kSub, wires[0], wires[0]);
    
    // Step 1: Group wires into input_size groups
    std::vector<std::vector<wire_t>> groups(input_size);
    for (size_t i = 0; i < input_size; ++i) {
      groups[i].resize(vec_size);
      for (size_t j = 0; j < vec_size; ++j) {
        groups[i][j] = wires[i + j * input_size];
      }
    }
    
    // Step 2: Initialize index wires [0, 1, 2, ..., vec_size-1]
    std::vector<wire_t> sort_permutation(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
      sort_permutation[i] = addConstOpGate(GateType::kConstAdd, zero_wire, static_cast<R>(i));
    }


    // Step 3: Iterative compaction for input_size iterations
    // Use int instead of size_t to avoid underflow when decrementing from 0
    for (int iter = static_cast<int>(input_size) - 1; iter >= 0; --iter) {
      // Take the current group (sorting key) from the last iteration
      std::vector<wire_t>& current_key = groups[iter];
      
      // Prepare payload vectors: all groups from 0 to iter (inclusive) + index_wires
      std::vector<std::vector<wire_t>> payloads;
      
      for (int g = 0; g <= iter; ++g) {
        payloads.push_back(groups[g]);
      }
      payloads.push_back(sort_permutation);
      
      // Apply compaction with current_key as the sorting key
      auto [key_compacted, payloads_compacted, labels] = addCompactionSubcircuit(current_key, payloads, permutation, pid);

      // Update groups 0 to iter with compacted outputs
      for (int g = 0; g <= iter; ++g) {
        groups[g] = payloads_compacted[g];
      }
      
      // Update index_wires with compacted indices
      sort_permutation = payloads_compacted[iter + 1];
    }   
    
    return sort_permutation;
  }


  
  // Level ordered gates are helpful for evaluation.
  [[nodiscard]] LevelOrderedCircuit orderGatesByLevel() const {
    LevelOrderedCircuit res;
    res.outputs = outputs_;
    res.num_gates = gates_.size();
    res.num_wires = num_wires;

    // Map from output wire id to multiplicative depth/level.
    // Input gates have a depth of 0.
    std::vector<size_t> gate_level(num_wires, 0);
    size_t depth = 0;

    // This assumes that if gates_[i]'s output is input to gates_[j] then
    // i < j.
    for (const auto& gate : gates_) {
      switch (gate->type) {
        case GateType::kRec: {
          const auto* g = static_cast<FIn1Gate*>(gate.get());
          gate_level[g->out] = gate_level[g->in] + 1;
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }
        case GateType::kAdd:
        case GateType::kSub: {
          const auto* g = static_cast<FIn2Gate*>(gate.get());
          gate_level[g->out] = std::max(gate_level[g->in1], gate_level[g->in2]);
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kMul: {
          const auto* g = static_cast<FIn2Gate*>(gate.get());
          gate_level[g->out] = std::max(gate_level[g->in1], gate_level[g->in2]) + 1;
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kConstAdd:
        case GateType::kConstMul: {
          const auto* g = static_cast<ConstOpGate<R>*>(gate.get());
          gate_level[g->out] = gate_level[g->in];
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kEqz: {
          const auto* g = static_cast<FIn1Gate*>(gate.get());
          gate_level[g->out] = gate_level[g->in] + 1;
          depth = std::max(depth, gate_level[gate->out]);
          break;
        }

        case GateType::kShuffle:
        case GateType::kPermAndSh: {
          const auto* g = static_cast<SIMDOGate*>(gate.get());
          size_t gate_depth = 0;
          for (size_t i = 0; i < g->in.size(); i++) {
            gate_depth = std::max({gate_level[g->in[i]], gate_depth});
          }
          for (int i = 0; i < g->outs.size(); i++) {
            gate_level[g->outs[i]] = gate_depth + 1;
          }
          depth = std::max(depth, gate_level[gate->outs[0]]);
          break;
        }

        case GateType::kPublicPerm: {
          const auto* g = static_cast<SIMDOGate*>(gate.get());
          size_t gate_depth = 0;
          for (size_t i = 0; i < g->in.size(); i++) {
            gate_depth = std::max({gate_level[g->in[i]], gate_depth});
          }
          for (int i = 0; i < g->outs.size(); i++) {
            gate_level[g->outs[i]] = gate_depth;
          }
          depth = std::max(depth, gate_level[gate->outs[0]]);
          break;
        }

        case GateType::kDeleteWires: {
          const auto* g = static_cast<SIMDOGate*>(gate.get());
          size_t gate_depth = 0;
          for (size_t i = 0; i < g->in.size(); i++) {
            gate_depth = std::max({gate_level[g->in[i]], gate_depth});
          }
          for (int i = 0; i < g->outs.size(); i++) {
            gate_level[g->outs[i]] = gate_depth + 1;
          }
          depth = std::max(depth, gate_level[gate->outs[0]]);
          break;
        }


        case GateType::kSort: {
          const auto* g = static_cast<SIMDOGate*>(gate.get());
          size_t gate_depth = 0;
          for (size_t i = 0; i < g->in.size(); i++) {
            gate_depth = std::max({gate_level[g->in[i]], gate_depth});
          }
          for (int i = 0; i < g->outs.size(); i++) {
            gate_level[g->outs[i]] = gate_depth + 1;
          }
          depth = std::max(depth, gate_level[gate->outs[0]]);
          break;
        }

        case GateType::kRewire: {
          const auto* g = static_cast<SIMDOGate*>(gate.get());
          size_t gate_depth = 0;
          for (size_t i = 0; i < g->in.size(); i++) {
            gate_depth = std::max({gate_level[g->in[i]], gate_depth});
          }
          for (int i = 0; i < g->outs.size(); i++) {
            gate_level[g->outs[i]] = gate_depth;
          }
          depth = std::max(depth, gate_level[gate->outs[0]]);
          break;
        }

        default:
          break;
      }
    }

    std::fill(res.count.begin(), res.count.end(), 0);

    std::vector<std::vector<gate_ptr_t>> gates_by_level(depth + 1);
    for (const auto& gate : gates_) {
      res.count[gate->type]++;
      if (gate->type == GateType::kShuffle || gate->type == GateType::kPermAndSh || gate->type == GateType::kPublicPerm || gate->type == GateType::kDeleteWires || gate->type == GateType::kSort || gate->type == GateType::kRewire) {
        gates_by_level[gate_level[gate->outs[0]]].push_back(gate);
      } else {
        gates_by_level[gate_level[gate->out]].push_back(gate);
      } 
    }

    res.gates_by_level = std::move(gates_by_level);

    return res;
  }

  // Evaluate circuit on plaintext inputs.
  [[nodiscard]] std::vector<R> evaluate(const std::unordered_map<wire_t, R>& inputs) const {
    auto level_circ = orderGatesByLevel();
    std::vector<R> wires(level_circ.num_gates);

    auto num_inp_gates = level_circ.count[GateType::kInp];
    if (inputs.size() != num_inp_gates) {
      throw std::invalid_argument(boost::str(
          boost::format("Expected %1% inputs but received %2% inputs.") %
          num_inp_gates % inputs.size()));
    }

    for (const auto& level : level_circ.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            wires[gate->out] = inputs.at(gate->out);
            break;
          }

          case GateType::kMul: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] * wires[g->in2];
            break;
          }

          case GateType::kAdd: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] + wires[g->in2];
            break;
          }

          case GateType::kSub: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] - wires[g->in2];
            break;
          }

          case GateType::kConstAdd: {
            auto* g = static_cast<ConstOpGate<R>*>(gate.get());
            wires[g->out] = wires[g->in] + g->cval;
            break;
          }

          case GateType::kConstMul: {
            auto* g = static_cast<ConstOpGate<R>*>(gate.get());
            wires[g->out] = wires[g->in] * g->cval;
            break;
          }

          case GateType::kEqz: {
            auto* g = static_cast<FIn1Gate*>(gate.get());
            if (wires[g->in] == 0) {
              wires[g->out] = 1;
            }
            else {
              wires[g->out] = 0;
            }
            break;
          }

          default: {
            throw std::runtime_error("Invalid gate type.");
          }
        }
      }
    }

    std::vector<R> outputs;
    for (auto i : level_circ.outputs) {
      outputs.push_back(wires[i]);
    }

    return outputs;
  }

};
};  // namespace common::utils
