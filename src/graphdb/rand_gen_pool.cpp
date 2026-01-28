#include "rand_gen_pool.h"

#include <algorithm>

#include "../utils/helpers.h"

namespace graphdb {

  RandGenPool::RandGenPool(int my_id, int num_parties, int num_clients, uint64_t seed) 
    : id_{my_id}, k_pi(num_parties + 1) { 
    
    
    auto seed_block = emp::makeBlock(seed, 0); 
    k_self.reseed(&seed_block, 0);
    k_all.reseed(&seed_block, 0);
    k_all_minus_0.reseed(&seed_block, 0);
    k_p0.reseed(&seed_block, 0);

    for (int i = 1; i <= num_parties + num_clients; i++) {
      k_pi[i].reseed(&seed_block, 0); 
    }
    
  }
  // all keys will be the same. for different keys look at emp toolkit

emp::PRG& RandGenPool::self() { return k_self; }

emp::PRG& RandGenPool::all() { return k_all; }

emp::PRG& RandGenPool::all_minus_0() { return k_all_minus_0; }

emp::PRG& RandGenPool::p0() { return k_p0; }

emp::PRG& RandGenPool::pi(int i) { return k_pi[i]; }
}  // namespace graphdb
