#include "circuit.h"

#include <stdexcept>

namespace common::utils {

Gate::Gate(GateType type, wire_t out) : type(type), out(out) {}

// Gate::Gate(GateType type, wire_t out, std::vector<wire_t> outs) : type(type), out(out), outs(outs) {}

Gate::Gate(GateType type, int owner, wire_t out, std::vector<wire_t> outs) : type(type), owner(owner), out(out), outs(outs) {}

Gate::Gate(GateType type, int owner, wire_t out, std::vector<std::vector<wire_t>> multi_outs)
    : type(type), owner(owner), out(out), multi_outs(multi_outs) {}

FIn2Gate::FIn2Gate(GateType type, wire_t in1, wire_t in2, wire_t out)
    : Gate(type, out), in1{in1}, in2{in2} {}

FIn3Gate::FIn3Gate(GateType type, wire_t in1, wire_t in2, wire_t in3, wire_t out)
    : Gate(type, out), in1{in1}, in2{in2}, in3{in3} {}

FIn4Gate::FIn4Gate(GateType type, wire_t in1, wire_t in2, wire_t in3, wire_t in4, wire_t out)
    : Gate(type, out), in1{in1}, in2{in2}, in3{in3}, in4{in4} {}

FIn1Gate::FIn1Gate(GateType type, wire_t in, wire_t out)
    : Gate(type, out), in{in} {}

SIMDGate::SIMDGate(GateType type, std::vector<wire_t> in1, std::vector<wire_t> in2, wire_t out)
    : Gate(type, out), in1(std::move(in1)), in2(std::move(in2)) {}

SIMDOGate::SIMDOGate(GateType type, int owner, std::vector<wire_t> in, std::vector<wire_t> outs, std::vector<std::vector<int>> permutation, size_t vec_size, int inv)
    : Gate(type, owner, outs[0], outs), in(std::move(in)), permutation(std::move(permutation)), vec_size(vec_size), inv(inv) {}

SIMDMOGate::SIMDMOGate(GateType type, int owner, std::vector<wire_t> in, std::vector<std::vector<wire_t>> multi_outs,
                       std::vector<std::vector<int>> permutation)
    : Gate(type, owner, multi_outs[0][0], multi_outs), in(std::move(in)), permutation(std::move(permutation)) {}

std::ostream& operator<<(std::ostream& os, GateType type) {
  switch (type) {
    case kInp:
      os << "Input";
      break;
    case kRec:
      os << "Reconstruction";
      break;

    case kAdd:
      os << "Addition";
      break;

    case kMul:
      os << "Multiplication";
      break;

    case kSub:
      os << "Subtraction";
      break;

    case kConstAdd:
      os << "Addition with constant";
      break;

    case kConstMul:
      os << "Multiplication with constant";
      break;
    
    case kEqz:
      os << "Equals to zero";
      break;

    case kShuffle:
      os << "Shuffle";
      break;

    case kPermAndSh:
      os << "Permute and Share";
      break;

    case kPublicPerm:
      os << "Public Permutation";
      break;

    case kDeleteWires:
      os << "Delete Wires";
      break;

    case kSort:
      os << "Sort";
      break;

    case kRewire:
      os << "Rewire";
      break;

    default:
      os << "Invalid";
      break;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const LevelOrderedCircuit& circ) {
  for (size_t i = 0; i < GateType::NumGates; ++i) {
    os << GateType(i) << ": " << circ.count[i] << std::endl;
  }
  os << "Total Gates: " << circ.num_gates << std::endl;
  os << "Total Wires: " << circ.num_wires << std::endl;
  os << "Depth: " << circ.gates_by_level.size() << std::endl;
  return os;
}
};  // namespace common::utils
