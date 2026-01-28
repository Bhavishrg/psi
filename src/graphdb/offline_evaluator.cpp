#include "offline_evaluator.h"

#include <NTL/BasicThreadPool.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>

// #include "../utils/helpers.h"

namespace graphdb {
OfflineEvaluator::OfflineEvaluator(int nP, int my_id,
                                   std::shared_ptr<io::NetIOMP> network,
                                   common::utils::LevelOrderedCircuit circ,
                                   int threads, int seed, int latency, bool use_pking)
    : nP_(nP),
      id_(my_id),
      latency_(latency),
      use_pking_(use_pking),
      rgen_(my_id, seed), 
      network_(std::move(network)),
      circ_(std::move(circ))
      // preproc_(circ.num_gates)

      { } // tpool_ = std::make_shared<ThreadPool>(threads); }



void OfflineEvaluator::randomShare(int nP, int pid, RandGenPool& rgen, AddShare<Ring>& share, TPShare<Ring>& tpShare) {
  Ring val = Ring(0);
  if (pid == 0) {
    share.pushValue(Ring(0));
    tpShare.pushValues(Ring(0));
    for (int i = 1; i <= nP; i++) {
      rgen.pi(i).random_data(&val, sizeof(Ring));
      tpShare.pushValues(val);
    }
  } else {
    rgen.p0().random_data(&val, sizeof(Ring));
    share.pushValue(val);
  }
}

void OfflineEvaluator::randomShare(int nP, int pid, RandGenPool& rgen, AddShare<Ring>& share) {
  Ring val = Ring(0);
  Ring sec = Ring(0);
  if (pid == 0) {
    share.pushValue(Ring(0));
    // tpShare.pushValues(Ring(0));
    for (int i = 1; i <= nP; i++) {
      rgen.pi(i).random_data(&val, sizeof(Ring));
      sec = sec + val; 
    }
    share.pushValue(sec);
  } else {
    rgen.p0().random_data(&val, sizeof(Ring));
    share.pushValue(val);
  }
}

void OfflineEvaluator::randomShareSecret(int nP, int pid, RandGenPool& rgen,
                                         AddShare<Ring>& share, TPShare<Ring>& tpShare, Ring secret,
                                         std::vector<Ring>& rand_sh_sec, size_t& idx_rand_sh_sec) {
  if (pid == 0) {
    Ring val = Ring(0);
    Ring valn = Ring(0);
    share.pushValue(Ring(0));
    tpShare.pushValues(Ring(0));
    for (int i = 1; i < nP; i++) {
      rgen.pi(i).random_data(&val, sizeof(Ring));
      tpShare.pushValues(val);
      valn += val;
    }
    valn = secret - valn;
    tpShare.pushValues(valn);
    rand_sh_sec.push_back(valn);
  } else {
    if (pid != nP) {
      Ring val;
      rgen.p0().random_data(&val, sizeof(Ring));
      share.pushValue(val);
    } else {
      share.pushValue(rand_sh_sec[idx_rand_sh_sec]);
      idx_rand_sh_sec++;
    }
  }
}

void OfflineEvaluator::randomShareSecret(int nP, int pid, RandGenPool& rgen,
                                         AddShare<Ring>& share, Ring secret,
                                         std::vector<Ring>& rand_sh_sec, size_t& idx_rand_sh_sec) {
  if (pid == 0) {
    Ring val = Ring(0);
    Ring valn = Ring(0);
    share.pushValue(secret);
    for (int i = 1; i < nP; i++) {
      rgen.pi(i).random_data(&val, sizeof(Ring));
      valn += val;
    }
    valn = secret - valn;
    rand_sh_sec.push_back(valn);
  } else {
    if (pid != nP) {
      Ring val;
      rgen.p0().random_data(&val, sizeof(Ring));
      share.pushValue(val);
    } else {
      share.pushValue(rand_sh_sec[idx_rand_sh_sec]);
      idx_rand_sh_sec++;
    }
  }
}

void OfflineEvaluator::randomPermutation(int nP, int pid, RandGenPool& rgen, std::vector<int>& pi, size_t& vec_size) {
  // Generate common random permutation between party 0 and party pid
  if (pid != 0) {
    pi.resize(vec_size);
    for (int i = 0; i < vec_size; ++i) {
      pi[i] = i;
    }
    // Fisher-Yates shuffle using common randomness with party 0
    for (int i = vec_size - 1; i > 0; --i) {
      uint32_t rand_val;
      rgen.p0().random_data(&rand_val, sizeof(uint32_t));
      int j = rand_val % (i + 1);
      std::swap(pi[i], pi[j]);
    }
  }
}

void OfflineEvaluator::generateShuffleDeltaVector(int nP, int pid, RandGenPool& rgen, std::vector<Ring>& delta,
                                                  std::vector<TPShare<Ring>>& tp_a, std::vector<TPShare<Ring>>& tp_b,
                                                  std::vector<TPShare<Ring>>& tp_c, std::vector<std::vector<int>>& tp_pi_all,
                                                  size_t& vec_size, std::vector<Ring>& rand_sh_sec, size_t& idx_rand_sh_sec) {
  if (pid == 0) {
    // Party 0 generates Δ for all shuffle gates and stores it in array
    // Δ = Π(a) + Σᵢ₌1 to ⁿ Πn·Πᵢ₋₁···Π₁(bn-ᵢ+1) + c
    // where Πᵢ is the permutation of party i
    // Π = Πₙ·Πₙ₋₁···Π₁ (composition of all party permutations)
    std::vector<Ring> deltan(vec_size);
    
    for (int i = 0; i < vec_size; ++i) {
      // Start with index i and apply all permutations from party 1 to nP: Π = Πₙ·Πₙ₋₁···Π₁
      int idx_perm = i;
      for (int j = 0; j < nP; ++j) {
        idx_perm = tp_pi_all[j][idx_perm];
      }
      
      // Compute Π(a[i])
      Ring pi_a = tp_a[idx_perm].secret();
      
      // Compute Σᵢ₌₁ⁿ Πn·Πᵢ₋₁···Π₁(bn-ᵢ+1)
      // For i from 1 to n:
      //   - Apply Π₁, Π₂, ..., Πᵢ₋₁ to position i to get intermediate index
      //   - Get b[intermediate_idx] for party (n-i+1)
      //   - Apply Πᵢ, Πᵢ₊₁, ..., Πn to intermediate_idx to get final position
      Ring sum_term = Ring(0);
      for (int i_term = 1; i_term <= nP; ++i_term) {
        // Party ID for this term (reversed: n-i+1)
        int party_id = nP - i_term + 1;
        
        // Apply Π₁, Π₂, ..., Πᵢ₋₁ to starting index i
        int intermediate_idx = i;
        for (int j = 0; j < i_term - 1; ++j) {
          intermediate_idx = tp_pi_all[j][intermediate_idx];
        }
        
        // Get b value for party_id at this intermediate position
        Ring b_value = tp_b[intermediate_idx][party_id];
        
        // Apply Πᵢ, Πᵢ₊₁, ..., Πn to find where this contributes in final output
        int final_idx = intermediate_idx;
        for (int j = i_term - 1; j < nP; ++j) {
          final_idx = tp_pi_all[j][final_idx];
        }
        
        // Accumulate at the final position
        sum_term += b_value;
      }
      
      // Compute c at final permuted position
      Ring c_final = tp_c[idx_perm].secret();
      
      // Final delta: Δ = Π(a) + Σᵢ₌₁ⁿ Πn·Πᵢ₋₁···Π₁(bn-ᵢ+1) - c
      deltan[idx_perm] = pi_a + sum_term + c_final;
    }    
    
    for (int i = 0; i < vec_size; ++i) {
      rand_sh_sec.push_back(deltan[i]);
    }
  } else if (pid == nP) {
    // Last party receives delta values from party 0
    delta.resize(vec_size);
    for (int i = 0; i < vec_size; ++i) {
      delta[i] = rand_sh_sec[idx_rand_sh_sec];
      idx_rand_sh_sec++;
    }
  
  }
  
}

void OfflineEvaluator::generateShuffleCorrelation(int nP, int pid, RandGenPool& rgen, std::vector<AddShare<Ring>>& r_1,
                                                 std::vector<AddShare<Ring>>& r_2, std::vector<Ring>& delta_r, std::vector<std::vector<int>>& tp_pi_all,
                                                  size_t& vec_size, std::vector<std::vector<Ring>>& shuffle_delta_sh, std::vector<size_t>& idx_shuffle_delta_sh) {
  if (pid == 0) {
    // Helper party (pid=0) generates intermediate random vectors and computes delta values
    // R[0] = r_1 (secret from shares)
    // R[i] = random vectors for i=1 to nP-1
    // R[nP] = r_2 (secret from shares)
    // delta_r for party i = tp_pi_all[i-1](R[i-1]) - R[i]
    
    // Initialize R with nP+1 vectors
    std::vector<std::vector<Ring>> R(nP + 1, std::vector<Ring>(vec_size));
    
    // R[0] = r_1 secret values
    for (int i = 0; i < vec_size; ++i) {
      R[0][i] = r_1[i].valueAt();
    }
    
    // Generate random intermediate vectors R[1] to R[nP-1]
    for (int j = 1; j < nP; ++j) {
      for (int i = 0; i < vec_size; ++i) {
        rgen.self().random_data(&R[j][i], sizeof(Ring));
      }
    }
    
    // R[nP] = r_2 secret values
    for (int i = 0; i < vec_size; ++i) {
      R[nP][i] = r_2[i].valueAt();
    }
    
    // Compute delta_r for each computing party (pid 1 to nP)
    // delta for party i (1-indexed) = tp_pi_all[i-1](R[i-1]) - R[i]
    // Online protocol: party p receives z[j], outputs z[pi[i]] + delta_r[i]
    // We want: (x[pi[i]] + R[p-1][pi[i]]) + delta = x[pi[i]] + R[p][i]
    // So delta = R[p][i] - R[p-1][pi[i]]
    for (int party_idx = 0; party_idx < nP; ++party_idx) {
      for (int i = 0; i < vec_size; ++i) {
        int idx_perm = tp_pi_all[party_idx][i];
        // For all parties: delta = R[party_idx+1][i] - R[party_idx][pi[i]]
        Ring delta = R[party_idx + 1][i] - R[party_idx][idx_perm];
        shuffle_delta_sh[party_idx].push_back(delta);
      }
    }

  } else {
    // Computing parties (pid 1 to nP) receive delta_r from shuffle_delta_sh
    for (int i = 0; i < vec_size; ++i) {
      delta_r[i] = shuffle_delta_sh[pid-1][idx_shuffle_delta_sh[pid-1]];
      idx_shuffle_delta_sh[pid-1]++;
    }
  }
}



void OfflineEvaluator::setWireMasksParty(const std::unordered_map<common::utils::wire_t, int>& input_pid_map, 
                                         std::vector<Ring>& rand_sh_sec,
                                         std::vector<std::vector<Ring>>& shuffle_delta_sh) {
  size_t idx_rand_sh_sec = 0;
  std::vector<size_t> idx_shuffle_delta_sh(nP_, 0);

  for (const auto& level : circ_.gates_by_level) {
    for (const auto& gate : level) {
      switch (gate->type) {
        case common::utils::GateType::kInp: {
          auto pregate = std::make_unique<PreprocInput<Ring>>();
          auto pid = input_pid_map.at(gate->out);
          pregate->pid = pid;
          preproc_.gates[gate->out] = std::move(pregate);
          break;
        }

        case common::utils::GateType::kRec: {
          auto pregate = std::make_unique<PreprocRecGate<Ring>>();
          // King party (party 1) receives the reconstructed value
          bool is_king = (id_ == 1);
          pregate->Pking = is_king;
          preproc_.gates[gate->out] = std::move(pregate);
          break;
        }

        case common::utils::GateType::kMul: {
          AddShare<Ring> triple_a; // Holds one beaver triple share of a random value a
          AddShare<Ring> triple_b; // Holds one beaver triple share of a random value b
          AddShare<Ring> triple_c; // Holds one beaver triple share of c=a*b
          // TPShare<Ring> tp_triple_c; // Holds all the beaver triple shares of c=a*b
          randomShare(nP_, id_, rgen_, triple_a);
          randomShare(nP_, id_, rgen_, triple_b);
          Ring tp_prod;
          if (id_ == 0) { tp_prod = triple_a.valueAt() * triple_b.valueAt(); }
          randomShareSecret(nP_, id_, rgen_, triple_c, tp_prod, rand_sh_sec, idx_rand_sh_sec);
          preproc_.gates[gate->out] =
              std::move(std::make_unique<PreprocMultGate<Ring>>(triple_a, triple_b, triple_c));
          break;
        }

         case common::utils::GateType::kEqz: {
          AddShare<Ring> share_r1;
          AddShare<Ring> share_r2;
          std::vector<AddShare<Ring>> share_r1_bits(RINGSIZEBITS);
          std::vector<AddShare<Ring>> share_r2_bits(RINGSIZEBITS);
          Ring r1 = Ring(0);
          Ring r2 = Ring(0);
          std::vector<Ring> r1_bits(RINGSIZEBITS);
          std::vector<Ring> r2_bits(RINGSIZEBITS);

          // sharing r1 and r1_bits
          randomShare(nP_, id_, rgen_, share_r1);
          
          if (id_ == 0) {
            r1 = share_r1.valueAt();
            r1_bits = bitDecomposeToInt(r1);
          }
          for (int i = 0; i < RINGSIZEBITS; ++i) {
              randomShareSecret(nP_, id_, rgen_, share_r1_bits[i], r1_bits[i],
                                                    rand_sh_sec, idx_rand_sh_sec);                                      
          }

          // sharing r2 and r2_bits
          if (id_ == 0) {
            rgen_.p0().random_data(&r2, sizeof(Ring));
            r2 = r2 % RINGSIZEBITS; // make sure r2 is in [0, RINGSIZEBITS-1]
          }
          randomShareSecret(nP_, id_, rgen_, share_r2, r2, rand_sh_sec, idx_rand_sh_sec);

          if (id_ == 0) {
            r2 = share_r2.valueAt();
            for (int i = 0; i < RINGSIZEBITS; ++i) {
              if (i == r2 % RINGSIZEBITS) {
                r2_bits[i] = 1;
              } else {
                r2_bits[i] = 0;
              }
            }
          }

          for (int i = 0; i < RINGSIZEBITS; ++i) {
            randomShareSecret(nP_, id_, rgen_, share_r2_bits[i], r2_bits[i],
                                                    rand_sh_sec, idx_rand_sh_sec);
          }
          preproc_.gates[gate->out] =
              std::make_unique<PreprocEqzGate<Ring>>(share_r1, share_r2, share_r1_bits, share_r2_bits);
          break;
        }

        case common::utils::GateType::kShuffle: {
          auto *shuffle_g = static_cast<common::utils::SIMDOGate *>(gate.get());
          auto vec_size = shuffle_g->in.size();

          std::vector<AddShare<Ring>> r_1(vec_size); // Randomly sampled vector for masking input
          std::vector<AddShare<Ring>> r_2(vec_size); // Randomly sampled vector for setting outputs
          std::vector<Ring> delta_r(vec_size); // delta_r[i] = pi_i(r[i]) - r[i+1]
          for (int i = 0; i < vec_size; i++) {
              randomShare(nP_, id_, rgen_, r_1[i]);
              randomShare(nP_, id_, rgen_, r_2[i]);
          }

          std::vector<int> pi; // Randomly sampled permutation of parties
          std::vector<std::vector<int>> tp_pi_all; // Randomly sampled permutations of all parties using HP
          if (id_ != 0) {
            pi = std::move(shuffle_g->permutation[0]);
          } else {
            tp_pi_all = std::move(shuffle_g->permutation);
          }

          generateShuffleCorrelation(nP_, id_, rgen_, r_1, r_2, delta_r, tp_pi_all, vec_size, shuffle_delta_sh, idx_shuffle_delta_sh);

          preproc_.gates[gate->out] =
              std::move(std::make_unique<PreprocShuffleGate<Ring>>(r_1, r_2, delta_r, pi, tp_pi_all));
          break;
        }

        case common::utils::GateType::kSort: {
          // Sort gate preprocessing: 32 compaction operations + final shuffle
          auto *sort_g = static_cast<common::utils::SIMDOGate *>(gate.get());
          // Input is bit-decomposed: [bit0_elem0, bit1_elem0, ..., bit31_elem0, bit0_elem1, ...]
          // Total size = 32 * vec_size
          auto total_size = sort_g->in.size();
          auto vec_size = total_size / 32;
          
          if (vec_size == 0 || total_size % 32 != 0) {
            throw std::runtime_error("Sort gate input size must be divisible by 32");
          }
          
          auto preproc_sort = std::make_unique<PreprocSortGate<Ring>>();
          
          // Initialize 2D vectors for 32 compaction operations
          preproc_sort->shuffle_a.resize(32);
          preproc_sort->shuffle_tp_a.resize(32);
          preproc_sort->shuffle_b.resize(32);
          preproc_sort->shuffle_tp_b.resize(32);
          preproc_sort->shuffle_c.resize(32);
          preproc_sort->shuffle_tp_c.resize(32);
          preproc_sort->shuffle_delta.resize(32);
          preproc_sort->shuffle_pi.resize(32);
          preproc_sort->shuffle_tp_pi_all.resize(32);
          
          preproc_sort->mult_triple_a.resize(32);
          preproc_sort->mult_tp_triple_a.resize(32);
          preproc_sort->mult_triple_b.resize(32);
          preproc_sort->mult_tp_triple_b.resize(32);
          preproc_sort->mult_triple_c.resize(32);
          preproc_sort->mult_tp_triple_c.resize(32);
          
          // Generate preprocessing for each of the 32 bit positions (MSB to LSB)
          for (int bit = 0; bit < 32; ++bit) {
            // Shuffle preprocessing for this bit's compaction
            preproc_sort->shuffle_a[bit].resize(vec_size);
            preproc_sort->shuffle_tp_a[bit].resize(vec_size);
            preproc_sort->shuffle_b[bit].resize(vec_size);
            preproc_sort->shuffle_tp_b[bit].resize(vec_size);
            preproc_sort->shuffle_c[bit].resize(vec_size);
            preproc_sort->shuffle_tp_c[bit].resize(vec_size);
            
            for (int i = 0; i < vec_size; i++) {
              randomShare(nP_, id_, rgen_, preproc_sort->shuffle_a[bit][i], 
                         preproc_sort->shuffle_tp_a[bit][i]);
              randomShare(nP_, id_, rgen_, preproc_sort->shuffle_b[bit][i], 
                         preproc_sort->shuffle_tp_b[bit][i]);
              randomShare(nP_, id_, rgen_, preproc_sort->shuffle_c[bit][i], 
                         preproc_sort->shuffle_tp_c[bit][i]);
            }
            
            // Generate random permutation for this bit's compaction
            if (id_ != 0) {
              preproc_sort->shuffle_pi[bit].resize(vec_size);
              for (int i = 0; i < vec_size; i++) {
                preproc_sort->shuffle_pi[bit][i] = i;
              }
              // Identity permutation
              for (size_t i = 0; i < vec_size; ++i) {
                preproc_sort->shuffle_pi[bit][i] = i;
              }
            } else {
              preproc_sort->shuffle_tp_pi_all[bit].resize(nP_);
              for (int p = 0; p < nP_; ++p) {
                preproc_sort->shuffle_tp_pi_all[bit][p].resize(vec_size);
                for (int i = 0; i < vec_size; i++) {
                  preproc_sort->shuffle_tp_pi_all[bit][p][i] = i;
                }
              }
            }
            
            preproc_sort->shuffle_delta[bit].resize(vec_size);
            generateShuffleDeltaVector(nP_, id_, rgen_, preproc_sort->shuffle_delta[bit], 
                                      preproc_sort->shuffle_tp_a[bit], 
                                      preproc_sort->shuffle_tp_b[bit], 
                                      preproc_sort->shuffle_tp_c[bit],
                                      preproc_sort->shuffle_tp_pi_all[bit], vec_size, 
                                      rand_sh_sec, idx_rand_sh_sec);
            
            // Multiplication triples for label computation in this compaction
            preproc_sort->mult_triple_a[bit].resize(vec_size);
            preproc_sort->mult_tp_triple_a[bit].resize(vec_size);
            preproc_sort->mult_triple_b[bit].resize(vec_size);
            preproc_sort->mult_tp_triple_b[bit].resize(vec_size);
            preproc_sort->mult_triple_c[bit].resize(vec_size);
            preproc_sort->mult_tp_triple_c[bit].resize(vec_size);
            
            for (int i = 0; i < vec_size; i++) {
              randomShare(nP_, id_, rgen_, preproc_sort->mult_triple_a[bit][i], 
                         preproc_sort->mult_tp_triple_a[bit][i]);
              randomShare(nP_, id_, rgen_, preproc_sort->mult_triple_b[bit][i], 
                         preproc_sort->mult_tp_triple_b[bit][i]);
              Ring tp_prod;
              if (id_ == 0) { 
                tp_prod = preproc_sort->mult_tp_triple_a[bit][i].secret() * 
                         preproc_sort->mult_tp_triple_b[bit][i].secret(); 
              }
              randomShareSecret(nP_, id_, rgen_, preproc_sort->mult_triple_c[bit][i], 
                               preproc_sort->mult_tp_triple_c[bit][i], tp_prod, 
                               rand_sh_sec, idx_rand_sh_sec);
            }
          }
          
          // Final shuffle preprocessing for hiding the sorting permutation before reconstruction
          preproc_sort->final_shuffle_a.resize(vec_size);
          preproc_sort->final_shuffle_tp_a.resize(vec_size);
          preproc_sort->final_shuffle_b.resize(vec_size);
          preproc_sort->final_shuffle_tp_b.resize(vec_size);
          preproc_sort->final_shuffle_c.resize(vec_size);
          preproc_sort->final_shuffle_tp_c.resize(vec_size);
          
          for (int i = 0; i < vec_size; i++) {
            randomShare(nP_, id_, rgen_, preproc_sort->final_shuffle_a[i], 
                       preproc_sort->final_shuffle_tp_a[i]);
            randomShare(nP_, id_, rgen_, preproc_sort->final_shuffle_b[i], 
                       preproc_sort->final_shuffle_tp_b[i]);
            randomShare(nP_, id_, rgen_, preproc_sort->final_shuffle_c[i], 
                       preproc_sort->final_shuffle_tp_c[i]);
          }
          
          // Generate random permutation for final shuffle
          if (id_ != 0) {
            preproc_sort->final_shuffle_pi.resize(vec_size);
            for (int i = 0; i < vec_size; i++) {
              preproc_sort->final_shuffle_pi[i] = i;
            }
            // Identity permutation
            for (size_t i = 0; i < vec_size; ++i) {
              preproc_sort->final_shuffle_pi[i] = i;
            }
          } else {
            preproc_sort->final_shuffle_tp_pi_all.resize(nP_);
            for (int p = 0; p < nP_; ++p) {
              preproc_sort->final_shuffle_tp_pi_all[p].resize(vec_size);
              for (int i = 0; i < vec_size; i++) {
                preproc_sort->final_shuffle_tp_pi_all[p][i] = i;
              }
            }
          }
          
          preproc_sort->final_shuffle_delta.resize(vec_size);
          generateShuffleDeltaVector(nP_, id_, rgen_, preproc_sort->final_shuffle_delta, 
                                    preproc_sort->final_shuffle_tp_a, 
                                    preproc_sort->final_shuffle_tp_b, 
                                    preproc_sort->final_shuffle_tp_c,
                                    preproc_sort->final_shuffle_tp_pi_all, vec_size, 
                                    rand_sh_sec, idx_rand_sh_sec);
          
          preproc_.gates[gate->out] = std::move(preproc_sort);
          break;
        }

        case common::utils::GateType::kRewire: {
          // Rewire gate requires no preprocessing since it applies a public permutation
          // The permutation is determined from the position_map wire values
          auto *rewire_g = static_cast<common::utils::SIMDOGate *>(gate.get());
          
          // Extract vec_size from gate metadata and calculate num_payloads
          size_t vec_size = rewire_g->vec_size;
          size_t total_size = rewire_g->in.size();
          size_t output_size = rewire_g->outs.size();
          
          // Calculate num_payloads: total_size = vec_size * (1 + num_payloads)
          size_t num_payloads = (total_size / vec_size) - 1;
          
          auto pregate = std::make_unique<PreprocRewireGate<Ring>>(vec_size, num_payloads);
          preproc_.gates[gate->out] = std::move(pregate);
          break;
        }

        case common::utils::GateType::kPublicPerm: {
          // The permutation is determined from the position_map wire values
          auto *publicperm_g = static_cast<common::utils::SIMDOGate *>(gate.get());
          
          // Extract vec_size from gate metadata and calculate num_payloads
          size_t vec_size = publicperm_g->vec_size;
          size_t total_size = publicperm_g->in.size();
          size_t output_size = publicperm_g->outs.size();
          
          // Calculate num_payloads: total_size = vec_size * (1 + num_payloads)
          size_t num_payloads = (total_size / vec_size) - 1;
          
          auto pregate = std::make_unique<PreprocPublicPermGate<Ring>>(vec_size, num_payloads);
          preproc_.gates[gate->out] = std::move(pregate);
          break;
        }

        case common::utils::GateType::kDeleteWires: {
          // Delete wires gate preprocessing: shuffle for del + payloads, then reconstruction
          auto *delete_g = static_cast<common::utils::SIMDOGate *>(gate.get());
          // Input format: [del_0,...,del_n, p1_0,...,p1_n, p2_0,...,p2_n, ...]
          // Output format: [p1_out_0,...,p1_out_n, p2_out_0,...,p2_out_n, ...] (compacted)
          
          size_t vec_size = delete_g->vec_size;
          size_t total_size = delete_g->in.size();
          size_t num_payloads = (total_size / vec_size) - 1;
          
          // Shuffle operates on vec_size elements
          // Preprocessing for shuffle operation (del + all payloads together)
          std::vector<AddShare<Ring>> shuffle_a(vec_size);
          std::vector<TPShare<Ring>> shuffle_tp_a(vec_size);
          std::vector<AddShare<Ring>> shuffle_b(vec_size);
          std::vector<TPShare<Ring>> shuffle_tp_b(vec_size);
          std::vector<AddShare<Ring>> shuffle_c(vec_size);
          std::vector<TPShare<Ring>> shuffle_tp_c(vec_size);
          
          for (size_t i = 0; i < vec_size; i++) {
            randomShare(nP_, id_, rgen_, shuffle_a[i], shuffle_tp_a[i]);
            randomShare(nP_, id_, rgen_, shuffle_b[i], shuffle_tp_b[i]);
            randomShare(nP_, id_, rgen_, shuffle_c[i], shuffle_tp_c[i]);
          }
          
          std::vector<int> shuffle_pi;
          std::vector<std::vector<int>> shuffle_tp_pi_all;
          if (id_ != 0) {
            shuffle_pi = std::move(delete_g->permutation[0]);
          } else {
            shuffle_tp_pi_all = std::move(delete_g->permutation);
          }
          
          std::vector<Ring> shuffle_delta(vec_size);
          generateShuffleDeltaVector(nP_, id_, rgen_, shuffle_delta, shuffle_tp_a, shuffle_tp_b, shuffle_tp_c,
                                    shuffle_tp_pi_all, vec_size, rand_sh_sec, idx_rand_sh_sec);
          
          preproc_.gates[gate->out] = std::move(std::make_unique<PreprocDeleteWiresGate<Ring>>(
              shuffle_a, shuffle_tp_a, shuffle_b, shuffle_tp_b, shuffle_c, shuffle_tp_c,
              shuffle_delta, shuffle_pi, shuffle_tp_pi_all));
          break;
        }

        default: {
          break;
        }
      }
    }
  }
}


void OfflineEvaluator::setWireMasks(const std::unordered_map<common::utils::wire_t, int>& input_pid_map) {
  std::vector<Ring> rand_sh_sec;
  std::vector<std::vector<Ring>> shuffle_delta_sh(nP_, std::vector<Ring>());

  if (id_ == 0) {
    setWireMasksParty(input_pid_map, rand_sh_sec, shuffle_delta_sh);
    
    size_t rand_sh_sec_num = rand_sh_sec.size();

    // Send shuffle_delta_sh to each computing party (1 to nP)
    // Party i receives shuffle_delta_sh[i-1] which contains its delta values
    for (int pid = 1; pid <= nP_; ++pid) {
      size_t shuffle_delta_sh_num = shuffle_delta_sh[pid - 1].size();
      network_->send(pid, &shuffle_delta_sh_num, sizeof(size_t));
      network_->send(pid, shuffle_delta_sh[pid - 1].data(), shuffle_delta_sh_num * sizeof(Ring));
    }

    // Send rand_sh_sec only to party nP
    network_->send(nP_, &rand_sh_sec_num, sizeof(size_t));
    network_->send(nP_, rand_sh_sec.data(), sizeof(Ring) * rand_sh_sec_num);

  } else if (id_ != nP_) {
    // Parties 1 to nP-1 receive their shuffle_delta_sh values
    size_t shuffle_delta_sh_num;
    usleep(latency_);
    network_->recv(0, &shuffle_delta_sh_num, sizeof(size_t));
    std::vector<Ring> shuffle_delta_sh_flat(shuffle_delta_sh_num);
    network_->recv(0, shuffle_delta_sh_flat.data(), shuffle_delta_sh_num * sizeof(Ring));
    
    std::vector<std::vector<Ring>> shuffle_delta_sh(nP_, std::vector<Ring>());
    shuffle_delta_sh[id_ - 1] = shuffle_delta_sh_flat;
    setWireMasksParty(input_pid_map, rand_sh_sec, shuffle_delta_sh);

  } else {
    // Party nP receives both shuffle_delta_sh and rand_sh_sec
    size_t shuffle_delta_sh_num;
    usleep(latency_);
    network_->recv(0, &shuffle_delta_sh_num, sizeof(size_t));
    std::vector<Ring> shuffle_delta_sh_flat(shuffle_delta_sh_num);
    network_->recv(0, shuffle_delta_sh_flat.data(), shuffle_delta_sh_num * sizeof(Ring));

    size_t rand_sh_sec_num;
    usleep(latency_);
    network_->recv(0, &rand_sh_sec_num, sizeof(size_t));
    rand_sh_sec.resize(rand_sh_sec_num);
    network_->recv(0, rand_sh_sec.data(), sizeof(Ring) * rand_sh_sec_num);

    std::vector<std::vector<Ring>> shuffle_delta_sh(nP_, std::vector<Ring>());
    shuffle_delta_sh[nP_ - 1] = shuffle_delta_sh_flat;
    setWireMasksParty(input_pid_map, rand_sh_sec, shuffle_delta_sh);
  }
}

PreprocCircuit<Ring> OfflineEvaluator::getPreproc() {
  return std::move(preproc_);
}

PreprocCircuit<Ring> OfflineEvaluator::run(const std::unordered_map<common::utils::wire_t, int>& input_pid_map) {
  setWireMasks(input_pid_map);
  return std::move(preproc_);
}

};  // namespace graphdb
