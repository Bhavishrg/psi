#include "sharing.h"

namespace graphdb {
//check the correctness of the following functions: 
template <>
void AddShare<BoolRing>::randomize(emp::PRG& prg) {
 bool data[1];
 prg.random_bool(static_cast<bool*>(data), 1);
 value_ = data[0];
}

//the following functions have dependencies on the number of parties
//How to handle the following functions
/*
template <>
TPShare<BoolRing>::TPShare(BoolRing secret, emp::PRG& prg) {
  bool values[5];
  prg.random_bool(static_cast<bool*>(values), 5);

  BoolRing sum;
  for (size_t i = 0; i < 5; ++i) {
    share_elements[i] = values[i];
    sum += share_elements[i];
  }
  share_elements[5] = secret - sum;
}

template <>
void TPShare<BoolRing>::randomize(emp::PRG& prg) {
  bool values[6];
  prg.random_bool(static_cast<bool*>(values), 6);

  for (size_t i = 0; i < 6; ++i) {
    share_elements[i] = values[i];
  }
}*/
};  // namespace graphdb
