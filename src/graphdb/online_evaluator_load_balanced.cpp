#include "online_evaluator.h"

#include "../utils/helpers.h"

namespace graphdb
{
    // Helper function to reconstruct shares via king party or direct all-to-all
    void OnlineEvaluator::reconstruct(int nP, int pid, std::shared_ptr<io::NetIOMP> network,
                                     const std::vector<Ring>& shares_list,
                                     std::vector<Ring>& reconstructed_list,
                                     bool via_pking, int latency) {
        int pKing = 1;
        size_t num_shares = shares_list.size();
        reconstructed_list.resize(num_shares, 0);
        
        if (via_pking) {
            // Reconstruction via king party
            if (pid != pKing) {
                network->send(pKing, shares_list.data(), shares_list.size() * sizeof(Ring));
                usleep(latency);
                network->recv(pKing, reconstructed_list.data(), reconstructed_list.size() * sizeof(Ring));
            } else {
                std::vector<std::vector<Ring>> share_recv(nP);
                share_recv[pKing - 1] = shares_list;
                usleep(latency);
                
                // Receive from all parties (not parallelized as recv is blocking)
                for (int p = 1; p <= nP; ++p) {
                    if (p != pKing) {
                        share_recv[p - 1].resize(num_shares);
                        network->recv(p, share_recv[p - 1].data(), share_recv[p - 1].size() * sizeof(Ring));
                    }
                }
                
                // Aggregate shares
                for (int p = 0; p < nP; ++p) {
                    for (size_t i = 0; i < num_shares; ++i) {
                        reconstructed_list[i] += share_recv[p][i];
                    }
                }
                
                // Send result to all parties (sequential for now to avoid race conditions)
                for (int p = 1; p <= nP; ++p) {
                    if (p != pKing) {
                        network->send(p, reconstructed_list.data(), reconstructed_list.size() * sizeof(Ring));
                        network->flush(p);
                    }
                }
            }
        } else {
            // Direct reconstruction (all parties exchange shares)
            std::vector<std::vector<Ring>> share_recv(nP);
            share_recv[pid - 1] = shares_list;
            
            // Send to all parties (sequential for now to avoid race conditions)
            for (int p = 1; p <= nP; ++p) {
                if (p != pid) {
                    network->send(p, shares_list.data(), shares_list.size() * sizeof(Ring));
                }
            }
            
            usleep(latency);
            
            // Receive from all parties
            for (int p = 1; p <= nP; ++p) {
                if (p != pid) {
                    share_recv[p - 1].resize(num_shares);
                    network->recv(p, share_recv[p - 1].data(), share_recv[p - 1].size() * sizeof(Ring));
                }
            }
            
            // Aggregate shares
            for (int p = 0; p < nP; ++p) {
                for (size_t i = 0; i < num_shares; ++i) {
                    reconstructed_list[i] += share_recv[p][i];
                }
            }
        }
    }

    OnlineEvaluator::OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                                     PreprocCircuit<Ring> preproc,
                                     common::utils::LevelOrderedCircuit circ,
                                     int threads, int seed, int latency, bool use_pking)
        : nP_(nP),
          id_(id),
          latency_(latency),
          use_pking_(use_pking),
          rgen_(id, seed),
          network_(std::move(network)),
          preproc_(std::move(preproc)),
          circ_(std::move(circ)),
          wires_(circ.num_wires)
    {
        tpool_ = std::make_shared<common::utils::ThreadPool>(threads);
    }

    OnlineEvaluator::OnlineEvaluator(int nP, int id, std::shared_ptr<io::NetIOMP> network,
                                     PreprocCircuit<Ring> preproc,
                                     common::utils::LevelOrderedCircuit circ,
                                     std::shared_ptr<common::utils::ThreadPool> tpool, int seed, int latency, bool use_pking)
        : nP_(nP),
          id_(id),
          latency_(latency),
          use_pking_(use_pking),
          rgen_(id, seed),
          network_(std::move(network)),
          preproc_(std::move(preproc)),
          circ_(std::move(circ)),
          tpool_(std::move(tpool)),
          wires_(circ.num_wires) {}

    
          void OnlineEvaluator::setInputs(const std::unordered_map<common::utils::wire_t, Ring> &inputs) {
        // Input gates have depth 0
        for (auto &g : circ_.gates_by_level[0]) {
            if (g->type == common::utils::GateType::kInp) {
                auto *pre_input = static_cast<PreprocInput<Ring> *>(preproc_.gates[g->out].get());
                auto pid = pre_input->pid;
                if (id_ != 0) {
                    if (pid == id_) {
                        Ring accumulated_val = Ring(0);
                        for (size_t i = 1; i <= nP_; i++) {
                            if (i != pid) {
                                Ring rand_sh;
                                rgen_.pi(i).random_data(&rand_sh, sizeof(Ring));
                                accumulated_val += rand_sh;
                            }
                        }
                        wires_[g->out] = inputs.at(g->out) - accumulated_val;
                    } else {
                        rgen_.pi(id_).random_data(&wires_[g->out], sizeof(Ring));
                    }
                }
            }
        }
    }

    void OnlineEvaluator::setRandomInputs() {
        // Input gates have depth 0.
        for (auto &g : circ_.gates_by_level[0]) {
            if (g->type == common::utils::GateType::kInp) {
                rgen_.pi(id_).random_data(&wires_[g->out], sizeof(Ring));
            }
        }
    }

    void OnlineEvaluator::eqzEvaluate(const std::vector<common::utils::FIn1Gate> &eqz_gates) {
        if (id_ == 0) { return; }
        int pKing = 1; // Designated king party
        size_t num_eqz_gates = eqz_gates.size();
        std::vector<Ring> r1_send(num_eqz_gates);
        std::vector<Ring> r2_send(num_eqz_gates);

        // Compute share of m1 = input + random_value r1
        #pragma omp parallel for
        for (size_t i = 0; i < num_eqz_gates; ++i) {
            auto &eqz_gate = eqz_gates[i];
            auto *pre_eqz = static_cast<PreprocEqzGate<Ring> *>(preproc_.gates[eqz_gate.out].get());
            Ring share_m1 = wires_[eqz_gate.in] + pre_eqz->share_r1.valueAt();
            r1_send[i] = share_m1;
        }

        // Reconstruct the masked input m1
        std::vector<Ring> recon_m1(num_eqz_gates, 0);
        reconstruct(nP_, id_, network_, r1_send, recon_m1, use_pking_, latency_);

        // Compute hamming distance between bits of m1 and bits of r1
        std::vector<Ring> share_m2(num_eqz_gates, 0);
        #pragma omp parallel for
        for (int i = 0; i < num_eqz_gates; ++i) {
            auto *pre_eqz = static_cast<PreprocEqzGate<Ring> *>(preproc_.gates[eqz_gates[i].out].get());
            std::vector<Ring> m1_bits(RINGSIZEBITS);
            m1_bits = bitDecomposeToInt(recon_m1[i]);
            std::vector<Ring> r1_bits(RINGSIZEBITS);
            for (int j = 0; j < RINGSIZEBITS; ++j) {
                r1_bits[j] = pre_eqz->share_r1_bits[j].valueAt();
            }
            if (id_ == 1) {
                for (int j = 0; j < RINGSIZEBITS; ++j) {
                    share_m2[i] += m1_bits[j] + r1_bits[j] - 2 * m1_bits[j] * r1_bits[j];
                }
                share_m2[i] += pre_eqz->share_r2.valueAt();
            }
            else{
                for (int j = 0; j < RINGSIZEBITS; ++j) {
                    share_m2[i] += r1_bits[j] - 2 * m1_bits[j] * r1_bits[j];
                }
                share_m2[i] += pre_eqz->share_r2.valueAt();
            }
        }


        // Reconstruct the masked input m2
        std::vector<Ring> recon_m2(num_eqz_gates, 0);
        reconstruct(nP_, id_, network_, share_m2, recon_m2, use_pking_, latency_);
 

        // Compute final output share
        std::vector<Ring> recon_out(num_eqz_gates, 0);
        #pragma omp parallel for
        for (int i = 0; i < num_eqz_gates; ++i) {
            auto *pre_eqz = static_cast<PreprocEqzGate<Ring> *>(preproc_.gates[eqz_gates[i].out].get());
            std::vector<Ring> r2_bits(RINGSIZEBITS);
            for (int j = 0; j < RINGSIZEBITS; ++j) {
                r2_bits[j] = pre_eqz->share_r2_bits[j].valueAt();
            }
            recon_out[i] = r2_bits[recon_m2[i]%RINGSIZEBITS];
            wires_[eqz_gates[i].out] = recon_out[i]; // Reconstructed output
        }
    }

    void OnlineEvaluator::recEvaluate(const std::vector<common::utils::FIn1Gate> &rec_gates) {
        if (id_ == 0) { return; }
        size_t num_rec_gates = rec_gates.size();
        std::vector<Ring> shares_to_send(num_rec_gates);

        // Gather shares to reconstruct
        #pragma omp parallel for
        for (size_t i = 0; i < num_rec_gates; ++i) {
            auto &rec_gate = rec_gates[i];
            shares_to_send[i] = wires_[rec_gate.in];
        }

        // Reconstruct the values
        std::vector<Ring> reconstructed(num_rec_gates, 0);
        
        reconstruct(nP_, id_, network_, shares_to_send, reconstructed, use_pking_, latency_);

        // Store reconstructed values in output wires
        #pragma omp parallel for
        for (size_t i = 0; i < num_rec_gates; ++i) {
            wires_[rec_gates[i].out] = reconstructed[i];
        }
    }


    void OnlineEvaluator::shuffleEvaluate(const std::vector<common::utils::SIMDOGate> &shuffle_gates) {
        if (id_ == 0 || shuffle_gates.empty()) { return; }

        std::vector<Ring> z_all;
        std::vector<std::vector<Ring>> z_sum;
        size_t total_comm = 0;

        for (const auto &gate : shuffle_gates) {
            auto *pre_shuffle = static_cast<PreprocShuffleGate<Ring> *>(preproc_.gates[gate.out].get());
            size_t vec_size = gate.in.size();
            total_comm += vec_size;
            std::vector<Ring> z(vec_size, 0);
            
            for (int i = 0; i < vec_size; ++i) {
                z[i] = wires_[gate.in[i]] + pre_shuffle->r_1[i].valueAt();
            }
            z_sum.push_back(z);
        }
        
        if (id_ == 1) {
            // Party 1 collects z values from all parties
            std::vector<std::vector<Ring>> z_recv_all(nP_);
            
            // Sequential receive instead of parallel to avoid deadlock
            for (int pid = 2; pid <= nP_; ++pid) {
                z_recv_all[pid - 1] = std::vector<Ring>(total_comm);
                network_->recv(pid, z_recv_all[pid - 1].data(), z_recv_all[pid - 1].size() * sizeof(Ring));
            }
            usleep(latency_);

            // Aggregate received z values from other parties
            for (int pid = 2; pid <= nP_; ++pid) {
                size_t idx_vec = 0;
                for (int idx_gate = 0; idx_gate < shuffle_gates.size(); ++idx_gate) {
                    size_t vec_size = shuffle_gates[idx_gate].in.size();
                    std::vector<Ring> z(z_recv_all[pid - 1].begin() + idx_vec, z_recv_all[pid - 1].begin() + idx_vec + vec_size);
                    for (int i = 0; i < vec_size; ++i) {
                        z_sum[idx_gate][i] += z[i];
                    }
                    idx_vec += vec_size;
                }
            }


            z_all.reserve(total_comm);

            for (int idx_gate = 0; idx_gate < shuffle_gates.size(); ++idx_gate) {
                auto *pre_shuffle = static_cast<PreprocShuffleGate<Ring> *>(preproc_.gates[shuffle_gates[idx_gate].out].get());
                size_t vec_size = shuffle_gates[idx_gate].in.size();
                std::vector<Ring> z(vec_size);

                #pragma omp parallel for
                for (int i = 0; i < vec_size; ++i) {
                    // Party 1: apply permutation, add delta
                    z[i] = z_sum[idx_gate][pre_shuffle->pi[i]] + pre_shuffle->delta_r[i];
                    // Set output as share of r_2
                    wires_[shuffle_gates[idx_gate].outs[i]] = -pre_shuffle->r_2[i].valueAt();
                }
                z_all.insert(z_all.end(), z.begin(), z.end());
            }
            network_->send(2, z_all.data(), z_all.size() * sizeof(Ring));

        } else {

            // Flatten z_sum into z_all before sending
            z_all.clear();
            z_all.reserve(total_comm);
            for (const auto& z_vec : z_sum) {
                z_all.insert(z_all.end(), z_vec.begin(), z_vec.end());
            }
            
            network_->send(1, z_all.data(), z_all.size() * sizeof(Ring));
            network_->flush(1);


            z_all.clear();
            z_all.resize(total_comm);
            network_->recv(id_ - 1, z_all.data(), z_all.size() * sizeof(Ring));
            usleep(latency_);

            for (int idx_gate = 0, idx_vec = 0; idx_gate < shuffle_gates.size(); ++idx_gate) {
                auto *pre_shuffle = static_cast<PreprocShuffleGate<Ring> *>(preproc_.gates[shuffle_gates[idx_gate].out].get());
                size_t vec_size = shuffle_gates[idx_gate].in.size();
                std::vector<Ring> z(z_all.begin() + idx_vec, z_all.begin() + idx_vec + vec_size);
                
                for (int i = 0; i < vec_size; ++i) {
                    if (id_ != nP_) {
                        // Parties 2 to nP-1: apply permutation, add delta, forward
                        Ring z_send = z[pre_shuffle->pi[i]] + pre_shuffle->delta_r[i];
                        // Set output as share of r_2
                        wires_[shuffle_gates[idx_gate].outs[i]] = -pre_shuffle->r_2[i].valueAt();
                        z_all[idx_vec++] = z_send;
                    } else {
                        // Party nP: apply permutation, add delta, add r_2 share
                        Ring permuted_val = z[pre_shuffle->pi[i]];
                        wires_[shuffle_gates[idx_gate].outs[i]] = permuted_val + pre_shuffle->delta_r[i] - pre_shuffle->r_2[i].valueAt();
                        idx_vec++;
                    }
                }
            }
            if (id_ != nP_) {
                network_->send(id_ + 1, z_all.data(), z_all.size() * sizeof(Ring));
            }
        }
    }


    void OnlineEvaluator::multEvaluate(const std::vector<common::utils::FIn2Gate> &mult_gates) {
        if (id_ == 0) { return; }
        
        size_t num_mult_gates = mult_gates.size();
        if (num_mult_gates == 0) { return; }
        
        // Step 1: Compute shares of u and v for each multiplication gate
        // For multiplication z = x * y using Beaver triples (a, b, c) where c = a * b:
        // - u = x - a (each party computes their share)
        // - v = y - b (each party computes their share)
        // - Reconstruct u and v to get actual values
        // - Compute z = u*v + u*b + v*a + c
        
        std::vector<Ring> u_shares(num_mult_gates);
        std::vector<Ring> v_shares(num_mult_gates);
        
        #pragma omp parallel for
        for (size_t i = 0; i < num_mult_gates; ++i) {
            auto &mult_gate = mult_gates[i];
            auto *pre_out = static_cast<PreprocMultGate<Ring> *>(preproc_.gates[mult_gate.out].get());
            
            // Compute this party's share of u and v
            u_shares[i] = wires_[mult_gate.in1] - pre_out->triple_a.valueAt();
            v_shares[i] = wires_[mult_gate.in2] - pre_out->triple_b.valueAt();
        }
        
        // Step 2: Reconstruct u and v values
        std::vector<Ring> u_reconstructed(num_mult_gates, 0);
        std::vector<Ring> v_reconstructed(num_mult_gates, 0);
        
        // Prepare shares to send (interleaved: u0, v0, u1, v1, ...)
        std::vector<Ring> shares_to_send(2 * num_mult_gates);
        for (size_t i = 0; i < num_mult_gates; ++i) {
            shares_to_send[2*i] = u_shares[i];
            shares_to_send[2*i + 1] = v_shares[i];
        }
        
        std::vector<Ring> reconstructed(2 * num_mult_gates, 0);
        reconstruct(nP_, id_, network_, shares_to_send, reconstructed, use_pking_, latency_);
        
        // Unpack reconstructed values
        for (size_t i = 0; i < num_mult_gates; ++i) {
            u_reconstructed[i] = reconstructed[2*i];
            v_reconstructed[i] = reconstructed[2*i + 1];
        }
        
        // Step 3: Compute multiplication result using Beaver triple formula
        // z = u*v + u*b + v*a + c (where only c is a share, rest are public values)
        #pragma omp parallel for
        for (size_t i = 0; i < num_mult_gates; ++i) {
            auto &mult_gate = mult_gates[i];
            auto *pre_out = static_cast<PreprocMultGate<Ring> *>(preproc_.gates[mult_gate.out].get());
            
            Ring u = u_reconstructed[i];
            Ring v = v_reconstructed[i];
            Ring a = pre_out->triple_a.valueAt();
            Ring b = pre_out->triple_b.valueAt();
            Ring c = pre_out->triple_c.valueAt();
            
            // Beaver triple formula: z = u*v + u*b + v*a + c
            if (id_ == 1){
                wires_[mult_gate.out] = u * v + u * b + v * a + c;
            }
            else{
                wires_[mult_gate.out] = u * b + v * a + c;
            }
        }
    }

    void OnlineEvaluator::sortEvaluate(const std::vector<common::utils::SIMDOGate> &sort_gates) {
        if (id_ == 0 || sort_gates.empty()) { return; }
        
        for (const auto& sort_gate : sort_gates) {
            auto *pre_sort = static_cast<PreprocSortGate<Ring> *>(preproc_.gates[sort_gate.out].get());
            
            // Input: [bit0_elem0, bit1_elem0, ..., bit31_elem0, bit0_elem1, ...]
            size_t total_size = sort_gate.in.size();
            size_t vec_size = total_size / 32;
            
            if (vec_size == 0 || total_size % 32 != 0) {
                throw std::runtime_error("Sort gate input size must be divisible by 32");
            }
            
            // Initialize payload as identity permutation [0, 1, 2, ..., vec_size-1] (0-indexed)
            std::vector<Ring> payload_shares(vec_size);
            for (size_t i = 0; i < vec_size; ++i) {
                if (id_ == 1) {
                    payload_shares[i] = static_cast<Ring>(i);
                } else {
                    payload_shares[i] = Ring(0);
                }
            }
            
            // Storage for bits being processed
            // bits_current[i] = shares of bit i for all elements (after compaction)
            std::vector<std::vector<Ring>> bits_current(32, std::vector<Ring>(vec_size));
            
            // Initialize: extract bit-decomposed input
            // Input format: consecutive 32 bits per element
            for (size_t elem = 0; elem < vec_size; ++elem) {
                for (size_t bit = 0; bit < 32; ++bit) {
                    bits_current[bit][elem] = wires_[sort_gate.in[elem * 32 + bit]];
                }
            }
            
            // Main loop: process each bit from MSB (bit 31) to LSB (bit 0)
            // In our storage: bit 31 is MSB, bit 0 is LSB
            for (int bit_idx = 31; bit_idx >= 0; --bit_idx) {
                // For preprocessing, we use array index (0 = first iteration = MSB)
                int preproc_idx = 31 - bit_idx;
                
                // Extract current bit as key vector
                std::vector<Ring> key_shares = bits_current[bit_idx];
                
                // Prepare payloads: all higher bits + current payload
                std::vector<std::vector<Ring>> compact_payloads;
                for (int higher_bit = 31; higher_bit > bit_idx; --higher_bit) {
                    compact_payloads.push_back(bits_current[higher_bit]);
                }
                compact_payloads.push_back(payload_shares);
                
                size_t num_payloads = compact_payloads.size();
                
                // Perform compaction using the preprocessing for this bit
                // This is similar to compactEvaluate but using pre_sort data
                
                // Step 1: Compute prefix sums c0 and c1 locally
                std::vector<Ring> c1_shares(vec_size);
                std::vector<Ring> c0_shares(vec_size);
                
                c1_shares[0] = key_shares[0];
                if (id_ == 1) {
                    c0_shares[0] = Ring(1) - c1_shares[0];
                } else {
                    c0_shares[0] = -c1_shares[0];
                }
                
                for (size_t j = 1; j < vec_size; ++j) {
                    c1_shares[j] = c1_shares[j-1] + key_shares[j];
                    if (id_ == 1) {
                        c0_shares[j] = Ring(j+1) - c1_shares[j];
                    } else {
                        c0_shares[j] = -c1_shares[j];
                    }
                }
                
                // Step 2: Compute label shares using multiplications
                std::vector<Ring> label_shares(vec_size);
                Ring c1_last = c1_shares[vec_size - 1];
                
                std::vector<Ring> u_shares(vec_size);
                std::vector<Ring> v_shares(vec_size);
                
                for (size_t j = 0; j < vec_size; ++j) {
                    Ring diff_term = c0_shares[j] + c1_last - c1_shares[j];
                    Ring one_minus_t;
                    if (id_ == 1) {
                        one_minus_t = Ring(1) - key_shares[j];
                    } else {
                        one_minus_t = -key_shares[j];
                    }
                    
                    u_shares[j] = diff_term - pre_sort->mult_triple_a[preproc_idx][j].valueAt();
                    v_shares[j] = one_minus_t - pre_sort->mult_triple_b[preproc_idx][j].valueAt();
                }
                
                std::vector<Ring> u_reconstructed(vec_size, 0);
                std::vector<Ring> v_reconstructed(vec_size, 0);
                
                std::vector<Ring> shares_to_send(2 * vec_size);
                for (size_t i = 0; i < vec_size; ++i) {
                    shares_to_send[2*i] = u_shares[i];
                    shares_to_send[2*i + 1] = v_shares[i];
                }
                
                std::vector<Ring> reconstructed(2 * vec_size, 0);
                reconstruct(nP_, id_, network_, shares_to_send, reconstructed, use_pking_, latency_);
                
                for (size_t i = 0; i < vec_size; ++i) {
                    u_reconstructed[i] = reconstructed[2*i];
                    v_reconstructed[i] = reconstructed[2*i + 1];
                }
                
                for (size_t j = 0; j < vec_size; ++j) {
                    Ring u = u_reconstructed[j];
                    Ring v = v_reconstructed[j];
                    Ring a = pre_sort->mult_triple_a[preproc_idx][j].valueAt();
                    Ring b = pre_sort->mult_triple_b[preproc_idx][j].valueAt();
                    Ring c = pre_sort->mult_triple_c[preproc_idx][j].valueAt();
                    
                    Ring mult_result;
                    if (id_ == 1) {
                        mult_result = u * v + u * b + v * a + c;
                        label_shares[j] = mult_result + c1_shares[j] - Ring(1);
                    } else {
                        mult_result = u * b + v * a + c;
                        label_shares[j] = mult_result + c1_shares[j];
                    }
                }
                
                // Step 3: Shuffle key, all payloads, and label
                std::vector<Ring> key_shuffled(vec_size);
                std::vector<std::vector<Ring>> payloads_shuffled(num_payloads, std::vector<Ring>(vec_size));
                std::vector<Ring> label_shuffled(vec_size);
                
                // Compute z = input + a for shuffle
                std::vector<Ring> z_key(vec_size);
                std::vector<std::vector<Ring>> z_payloads(num_payloads, std::vector<Ring>(vec_size));
                std::vector<Ring> z_label(vec_size);
                
                for (size_t i = 0; i < vec_size; ++i) {
                    z_key[i] = key_shares[i] + pre_sort->shuffle_a[preproc_idx][i].valueAt();
                    z_label[i] = label_shares[i] + pre_sort->shuffle_a[preproc_idx][i].valueAt();
                }
                for (size_t p = 0; p < num_payloads; ++p) {
                    for (size_t i = 0; i < vec_size; ++i) {
                        z_payloads[p][i] = compact_payloads[p][i] + pre_sort->shuffle_a[preproc_idx][i].valueAt();
                    }
                }
                
                // Shuffle protocol (same as in compactEvaluate)
                if (id_ == 1) {
                    std::vector<std::vector<Ring>> z_key_recv(nP_);
                    std::vector<std::vector<std::vector<Ring>>> z_payloads_recv(nP_, std::vector<std::vector<Ring>>(num_payloads));
                    std::vector<std::vector<Ring>> z_label_recv(nP_);
                    
                    z_key_recv[0] = z_key;
                    z_payloads_recv[0] = z_payloads;
                    z_label_recv[0] = z_label;
                    
                    #pragma omp parallel for
                    for (int pid = 2; pid <= nP_; ++pid) {
                        z_key_recv[pid - 1].resize(vec_size);
                        z_label_recv[pid - 1].resize(vec_size);
                        network_->recv(pid, z_key_recv[pid - 1].data(), vec_size * sizeof(Ring));
                        for (size_t p = 0; p < num_payloads; ++p) {
                            z_payloads_recv[pid - 1][p].resize(vec_size);
                            network_->recv(pid, z_payloads_recv[pid - 1][p].data(), vec_size * sizeof(Ring));
                        }
                        network_->recv(pid, z_label_recv[pid - 1].data(), vec_size * sizeof(Ring));
                    }
                    usleep(latency_);
                    
                    std::vector<Ring> z_key_sum(vec_size, 0);
                    std::vector<std::vector<Ring>> z_payloads_sum(num_payloads, std::vector<Ring>(vec_size, 0));
                    std::vector<Ring> z_label_sum(vec_size, 0);
                    
                    for (int pid = 0; pid < nP_; ++pid) {
                        for (size_t i = 0; i < vec_size; ++i) {
                            z_key_sum[i] += z_key_recv[pid][i];
                            z_label_sum[i] += z_label_recv[pid][i];
                        }
                        for (size_t p = 0; p < num_payloads; ++p) {
                            for (size_t i = 0; i < vec_size; ++i) {
                                z_payloads_sum[p][i] += z_payloads_recv[pid][p][i];
                            }
                        }
                    }
                    
                    std::vector<Ring> z_key_perm(vec_size);
                    std::vector<std::vector<Ring>> z_payloads_perm(num_payloads, std::vector<Ring>(vec_size));
                    std::vector<Ring> z_label_perm(vec_size);
                    
                    for (size_t i = 0; i < vec_size; ++i) {
                        int pi_i = pre_sort->shuffle_pi[preproc_idx][i];
                        z_key_perm[i] = z_key_sum[pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt();
                        z_label_perm[i] = z_label_sum[pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt();
                        
                        key_shuffled[i] = pre_sort->shuffle_c[preproc_idx][i].valueAt();
                        label_shuffled[i] = pre_sort->shuffle_c[preproc_idx][i].valueAt();
                    }
                    for (size_t p = 0; p < num_payloads; ++p) {
                        for (size_t i = 0; i < vec_size; ++i) {
                            int pi_i = pre_sort->shuffle_pi[preproc_idx][i];
                            z_payloads_perm[p][i] = z_payloads_sum[p][pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt();
                            payloads_shuffled[p][i] = pre_sort->shuffle_c[preproc_idx][i].valueAt();
                        }
                    }
                    
                    network_->send(2, z_key_perm.data(), vec_size * sizeof(Ring));
                    for (size_t p = 0; p < num_payloads; ++p) {
                        network_->send(2, z_payloads_perm[p].data(), vec_size * sizeof(Ring));
                    }
                    network_->send(2, z_label_perm.data(), vec_size * sizeof(Ring));
                    
                } else {
                    network_->send(1, z_key.data(), vec_size * sizeof(Ring));
                    for (size_t p = 0; p < num_payloads; ++p) {
                        network_->send(1, z_payloads[p].data(), vec_size * sizeof(Ring));
                    }
                    network_->send(1, z_label.data(), vec_size * sizeof(Ring));
                    network_->flush(1);
                    
                    std::vector<Ring> z_key_recv(vec_size);
                    std::vector<std::vector<Ring>> z_payloads_recv(num_payloads, std::vector<Ring>(vec_size));
                    std::vector<Ring> z_label_recv(vec_size);
                    
                    network_->recv(id_ - 1, z_key_recv.data(), vec_size * sizeof(Ring));
                    for (size_t p = 0; p < num_payloads; ++p) {
                        network_->recv(id_ - 1, z_payloads_recv[p].data(), vec_size * sizeof(Ring));
                    }
                    network_->recv(id_ - 1, z_label_recv.data(), vec_size * sizeof(Ring));
                    usleep(latency_);
                    
                    std::vector<Ring> z_key_send(vec_size);
                    std::vector<std::vector<Ring>> z_payloads_send(num_payloads, std::vector<Ring>(vec_size));
                    std::vector<Ring> z_label_send(vec_size);
                    
                    for (size_t i = 0; i < vec_size; ++i) {
                        int pi_i = pre_sort->shuffle_pi[preproc_idx][i];
                        
                        if (id_ != nP_) {
                            z_key_send[i] = z_key_recv[pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt();
                            z_label_send[i] = z_label_recv[pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt();
                            
                            key_shuffled[i] = pre_sort->shuffle_c[preproc_idx][i].valueAt();
                            label_shuffled[i] = pre_sort->shuffle_c[preproc_idx][i].valueAt();
                        } else {
                            z_key_send[i] = z_key_recv[pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt() - 
                                           pre_sort->shuffle_delta[preproc_idx][i];
                            z_label_send[i] = z_label_recv[pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt() - 
                                             pre_sort->shuffle_delta[preproc_idx][i];
                            
                            key_shuffled[i] = z_key_send[i] + pre_sort->shuffle_c[preproc_idx][i].valueAt();
                            label_shuffled[i] = z_label_send[i] + pre_sort->shuffle_c[preproc_idx][i].valueAt();
                        }
                    }
                    
                    for (size_t p = 0; p < num_payloads; ++p) {
                        for (size_t i = 0; i < vec_size; ++i) {
                            int pi_i = pre_sort->shuffle_pi[preproc_idx][i];
                            
                            if (id_ != nP_) {
                                z_payloads_send[p][i] = z_payloads_recv[p][pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt();
                                payloads_shuffled[p][i] = pre_sort->shuffle_c[preproc_idx][i].valueAt();
                            } else {
                                z_payloads_send[p][i] = z_payloads_recv[p][pi_i] + pre_sort->shuffle_b[preproc_idx][pi_i].valueAt() - 
                                                       pre_sort->shuffle_delta[preproc_idx][i];
                                payloads_shuffled[p][i] = z_payloads_send[p][i] + pre_sort->shuffle_c[preproc_idx][i].valueAt();
                            }
                        }
                    }
                    
                    if (id_ != nP_) {
                        network_->send(id_ + 1, z_key_send.data(), vec_size * sizeof(Ring));
                        for (size_t p = 0; p < num_payloads; ++p) {
                            network_->send(id_ + 1, z_payloads_send[p].data(), vec_size * sizeof(Ring));
                        }
                        network_->send(id_ + 1, z_label_send.data(), vec_size * sizeof(Ring));
                    }
                }
                
                // Step 4: Reconstruct label (don't need to reconstruct key or payloads)
                std::vector<Ring> label_reconstructed(vec_size, 0);
                reconstruct(nP_, id_, network_, label_shuffled, label_reconstructed, use_pking_, latency_);
                
                // Step 5: Apply permutation based on reconstructed labels
                // Reorder key, payloads based on label indices
                std::vector<Ring> key_compacted(vec_size);
                std::vector<std::vector<Ring>> payloads_compacted(num_payloads, std::vector<Ring>(vec_size));
                
                for (size_t i = 0; i < vec_size; ++i) {
                    int idx_perm = static_cast<int>(label_reconstructed[i]);
                    if (idx_perm >= 0 && idx_perm < vec_size) {
                        key_compacted[idx_perm] = key_shuffled[i];
                        for (size_t p = 0; p < num_payloads; ++p) {
                            payloads_compacted[p][idx_perm] = payloads_shuffled[p][i];
                        }
                    }
                }
                
                // Update for next iteration
                bits_current[bit_idx] = key_compacted;
                for (int higher_bit = 31, p = 0; higher_bit > bit_idx; --higher_bit, ++p) {
                    bits_current[higher_bit] = payloads_compacted[p];
                }
                payload_shares = payloads_compacted[num_payloads - 1];
            }
            
            // After all 32 iterations, payload_shares contains the sorting permutation (0-indexed)
            // Now shuffle and reveal the permutation
            std::vector<Ring> z_perm(vec_size);
            for (size_t i = 0; i < vec_size; ++i) {
                z_perm[i] = payload_shares[i] + pre_sort->final_shuffle_a[i].valueAt();
            }
            
            std::vector<Ring> payload_shuffled(vec_size);
            
            // Final shuffle protocol
            if (id_ == 1) {
                std::vector<std::vector<Ring>> z_perm_recv(nP_);
                z_perm_recv[0] = z_perm;
                
                #pragma omp parallel for
                for (int pid = 2; pid <= nP_; ++pid) {
                    z_perm_recv[pid - 1].resize(vec_size);
                    network_->recv(pid, z_perm_recv[pid - 1].data(), vec_size * sizeof(Ring));
                }
                usleep(latency_);
                
                std::vector<Ring> z_perm_sum(vec_size, 0);
                for (int pid = 0; pid < nP_; ++pid) {
                    for (size_t i = 0; i < vec_size; ++i) {
                        z_perm_sum[i] += z_perm_recv[pid][i];
                    }
                }
                
                std::vector<Ring> z_perm_perm(vec_size);
                for (size_t i = 0; i < vec_size; ++i) {
                    int pi_i = pre_sort->final_shuffle_pi[i];
                    z_perm_perm[i] = z_perm_sum[pi_i] + pre_sort->final_shuffle_b[pi_i].valueAt();
                    payload_shuffled[i] = pre_sort->final_shuffle_c[i].valueAt();
                }
                
                network_->send(2, z_perm_perm.data(), vec_size * sizeof(Ring));
                
            } else {
                network_->send(1, z_perm.data(), vec_size * sizeof(Ring));
                network_->flush(1);
                
                std::vector<Ring> z_perm_recv(vec_size);
                network_->recv(id_ - 1, z_perm_recv.data(), vec_size * sizeof(Ring));
                usleep(latency_);
                
                std::vector<Ring> z_perm_send(vec_size);
                for (size_t i = 0; i < vec_size; ++i) {
                    int pi_i = pre_sort->final_shuffle_pi[i];
                    
                    if (id_ != nP_) {
                        z_perm_send[i] = z_perm_recv[pi_i] + pre_sort->final_shuffle_b[pi_i].valueAt();
                        payload_shuffled[i] = pre_sort->final_shuffle_c[i].valueAt();
                    } else {
                        z_perm_send[i] = z_perm_recv[pi_i] + pre_sort->final_shuffle_b[pi_i].valueAt() - 
                                        pre_sort->final_shuffle_delta[i];
                        payload_shuffled[i] = z_perm_send[i] + pre_sort->final_shuffle_c[i].valueAt();
                    }
                }
                
                if (id_ != nP_) {
                    network_->send(id_ + 1, z_perm_send.data(), vec_size * sizeof(Ring));
                }
            }
            
            // Reconstruct the permutation
            std::vector<Ring> perm_reconstructed(vec_size, 0);
            reconstruct(nP_, id_, network_, payload_shuffled, perm_reconstructed, use_pking_, latency_);
            
            // Final step: Rearrange all input wires according to revealed permutation
            // For each element, all 32 bits move together
            for (size_t elem = 0; elem < vec_size; ++elem) {
                int perm_idx = static_cast<int>(perm_reconstructed[elem]);
                if (perm_idx >= 0 && perm_idx < vec_size) {
                    // Copy all 32 bits from input element perm_idx to output element elem
                    for (size_t bit = 0; bit < 32; ++bit) {
                        wires_[sort_gate.outs[elem * 32 + bit]] = wires_[sort_gate.in[perm_idx * 32 + bit]];
                    }
                }
            }
        }
    }

    void OnlineEvaluator::deleteWiresEvaluate(const std::vector<common::utils::SIMDOGate> &delete_gates) {
        if (id_ == 0) { return; }
        if (delete_gates.empty()) { return; }

        for (const auto& delete_gate : delete_gates) {
            auto* pre_delete = static_cast<PreprocDeleteWiresGate<Ring>*>(preproc_.gates[delete_gate.out].get());
            
            // Input format: [del_0,...,del_n, p1_0,...,p1_n, p2_0,...,p2_n, ...]
            // Output format: [p1_out_0,...,p1_out_n, p2_out_0,...,p2_out_n, ...] (compacted)
            
            size_t vec_size = delete_gate.vec_size;  // From gate metadata
            size_t total_inputs = delete_gate.in.size();
            size_t num_payloads = (total_inputs / vec_size) - 1;
            
            // Extract del and payload shares
            std::vector<Ring> del_shares(vec_size);
            std::vector<std::vector<Ring>> payload_shares(num_payloads, std::vector<Ring>(vec_size));
            
            for (size_t i = 0; i < vec_size; ++i) {
                del_shares[i] = wires_[delete_gate.in[i]];
            }
            for (size_t p = 0; p < num_payloads; ++p) {
                for (size_t i = 0; i < vec_size; ++i) {
                    payload_shares[p][i] = wires_[delete_gate.in[vec_size * (p + 1) + i]];
                }
            }
            
            // Step 1: Shuffle del and all payloads using the same permutation
            std::vector<Ring> shuffled_del(vec_size);
            std::vector<std::vector<Ring>> shuffled_payloads(num_payloads, std::vector<Ring>(vec_size));
            
            // Shuffle logic - mask with 'a' values
            std::vector<Ring> z_del(vec_size);
            std::vector<std::vector<Ring>> z_payloads(num_payloads, std::vector<Ring>(vec_size));
            
            for (size_t i = 0; i < vec_size; ++i) {
                z_del[i] = del_shares[i] + pre_delete->shuffle_a[i].valueAt();
            }
            for (size_t p = 0; p < num_payloads; ++p) {
                for (size_t i = 0; i < vec_size; ++i) {
                    z_payloads[p][i] = payload_shares[p][i] + pre_delete->shuffle_a[i].valueAt();
                }
            }
            
            // Shuffle protocol for del and payloads using same permutation
            if (id_ == 1) {
                // Party 1 collects from all parties
                std::vector<Ring> z_del_sum = z_del;
                std::vector<std::vector<Ring>> z_payloads_sum = z_payloads;
                
                for (int pid = 2; pid <= nP_; ++pid) {
                    std::vector<Ring> z_del_recv(vec_size);
                    network_->recv(pid, z_del_recv.data(), vec_size * sizeof(Ring));
                    for (size_t i = 0; i < vec_size; ++i) {
                        z_del_sum[i] += z_del_recv[i];
                    }
                    
                    for (size_t p = 0; p < num_payloads; ++p) {
                        std::vector<Ring> z_payload_recv(vec_size);
                        network_->recv(pid, z_payload_recv.data(), vec_size * sizeof(Ring));
                        for (size_t i = 0; i < vec_size; ++i) {
                            z_payloads_sum[p][i] += z_payload_recv[i];
                        }
                    }
                }
                usleep(latency_);
                
                // Apply permutation and mask for del
                std::vector<Ring> z_del_permuted(vec_size);
                for (size_t i = 0; i < vec_size; ++i) {
                    z_del_permuted[i] = z_del_sum[pre_delete->shuffle_pi[i]] + pre_delete->shuffle_b[pre_delete->shuffle_pi[i]].valueAt();
                    shuffled_del[i] = pre_delete->shuffle_c[i].valueAt();
                }
                network_->send(2, z_del_permuted.data(), vec_size * sizeof(Ring));
                
                // Apply permutation and mask for payloads
                for (size_t p = 0; p < num_payloads; ++p) {
                    std::vector<Ring> z_payload_permuted(vec_size);
                    for (size_t i = 0; i < vec_size; ++i) {
                        z_payload_permuted[i] = z_payloads_sum[p][pre_delete->shuffle_pi[i]] + pre_delete->shuffle_b[pre_delete->shuffle_pi[i]].valueAt();
                        shuffled_payloads[p][i] = pre_delete->shuffle_c[i].valueAt();
                    }
                    network_->send(2, z_payload_permuted.data(), vec_size * sizeof(Ring));
                }
            } else {
                // Send to party 1
                network_->send(1, z_del.data(), vec_size * sizeof(Ring));
                for (size_t p = 0; p < num_payloads; ++p) {
                    network_->send(1, z_payloads[p].data(), vec_size * sizeof(Ring));
                }
                network_->flush(1);
                
                // Receive from previous party
                std::vector<Ring> z_del_permuted(vec_size);
                network_->recv(id_ - 1, z_del_permuted.data(), vec_size * sizeof(Ring));
                usleep(latency_);
                
                std::vector<std::vector<Ring>> z_payloads_permuted(num_payloads, std::vector<Ring>(vec_size));
                for (size_t p = 0; p < num_payloads; ++p) {
                    network_->recv(id_ - 1, z_payloads_permuted[p].data(), vec_size * sizeof(Ring));
                }
                
                // Apply own permutation for del
                std::vector<Ring> z_del_next(vec_size);
                for (size_t i = 0; i < vec_size; ++i) {
                    if (id_ != nP_) {
                        z_del_next[i] = z_del_permuted[pre_delete->shuffle_pi[i]] + pre_delete->shuffle_b[pre_delete->shuffle_pi[i]].valueAt();
                        shuffled_del[i] = pre_delete->shuffle_c[i].valueAt();
                    } else {
                        z_del_next[i] = z_del_permuted[pre_delete->shuffle_pi[i]] + pre_delete->shuffle_b[pre_delete->shuffle_pi[i]].valueAt() - pre_delete->shuffle_delta[i];
                        shuffled_del[i] = z_del_next[i] + pre_delete->shuffle_c[i].valueAt();
                    }
                }
                
                if (id_ != nP_) {
                    network_->send(id_ + 1, z_del_next.data(), vec_size * sizeof(Ring));
                }
                
                // Apply own permutation for payloads
                for (size_t p = 0; p < num_payloads; ++p) {
                    std::vector<Ring> z_payload_next(vec_size);
                    for (size_t i = 0; i < vec_size; ++i) {
                        if (id_ != nP_) {
                            z_payload_next[i] = z_payloads_permuted[p][pre_delete->shuffle_pi[i]] + pre_delete->shuffle_b[pre_delete->shuffle_pi[i]].valueAt();
                            shuffled_payloads[p][i] = pre_delete->shuffle_c[i].valueAt();
                        } else {
                            z_payload_next[i] = z_payloads_permuted[p][pre_delete->shuffle_pi[i]] + pre_delete->shuffle_b[pre_delete->shuffle_pi[i]].valueAt() - pre_delete->shuffle_delta[i];
                            shuffled_payloads[p][i] = z_payload_next[i] + pre_delete->shuffle_c[i].valueAt();
                        }
                    }
                    
                    if (id_ != nP_) {
                        network_->send(id_ + 1, z_payload_next.data(), vec_size * sizeof(Ring));
                    }
                }
            }
            
            // Step 2: Reconstruct del to reveal deletion positions
            std::vector<Ring> del_reconstructed(vec_size, 0);
            reconstruct(nP_, id_, network_, shuffled_del, del_reconstructed, use_pking_, latency_);
            
            // Step 3: Compact payloads by removing indices where del == 1
            // Count non-deleted elements
            std::vector<size_t> keep_indices;
            for (size_t i = 0; i < vec_size; ++i) {
                if (del_reconstructed[i] == 0) {
                    keep_indices.push_back(i);
                }
            }
            
            // Write keep_indices to the first output wire (only party 1 has the actual value)
            // We'll encode the count as the wire value for party 1
            if (id_ == 1) {
                wires_[delete_gate.outs[0]] = static_cast<Ring>(keep_indices.size());
            } else {
                wires_[delete_gate.outs[0]] = Ring(0);
            }
            
            // Write compacted outputs to wires (front-packed)
            // Output format: [keep_indices_wire, p1_0,...,p1_n, p2_0,...,p2_n, ...]
            for (size_t p = 0; p < num_payloads; ++p) {
                for (size_t i = 0; i < vec_size; ++i) {
                    if (i < keep_indices.size()) {
                        wires_[delete_gate.outs[1 + vec_size * p + i]] = shuffled_payloads[p][keep_indices[i]];
                    } else {
                        // Fill remaining slots with zero (or leave as-is)
                        wires_[delete_gate.outs[1 + vec_size * p + i]] = Ring(0);
                    }
                }
            }
        }
    }

    void OnlineEvaluator::rewireEvaluate(const std::vector<common::utils::SIMDOGate> &rewire_gates) {
        if (id_ == 0) { return; }
        if (rewire_gates.empty()) { return; }

        size_t num_gates = rewire_gates.size();
        
        // Extract metadata and prepare data structures for all gates in parallel
        std::vector<size_t> vec_sizes(num_gates);
        std::vector<size_t> num_payloads_vec(num_gates);
        std::vector<std::vector<Ring>> all_position_maps(num_gates);
        std::vector<std::vector<std::vector<Ring>>> all_payload_shares(num_gates);
        
        #pragma omp parallel for
        for (size_t g = 0; g < num_gates; ++g) {
            const auto& rewire_gate = rewire_gates[g];
            
            // Input format: [pos_map_0, ..., pos_map_n, p1_0, ..., p1_n, p2_0, ..., p2_n, ...]
            // Output format: [p1_out_0, ..., p1_out_n, p2_out_0, ..., p2_out_n, ...]
            // Note: position_map wires already hold reconstructed values (public)
            
            // Get vec_size and num_payloads from preprocessing
            auto* preproc_rewire = static_cast<PreprocRewireGate<Ring>*>(preproc_.gates[rewire_gate.out].get());
            size_t vec_size = preproc_rewire->vec_size;
            size_t num_payloads = preproc_rewire->num_payloads;
            
            vec_sizes[g] = vec_size;
            num_payloads_vec[g] = num_payloads;
            
            // Extract position_map values (already reconstructed/public)
            all_position_maps[g].resize(vec_size);
            for (size_t i = 0; i < vec_size; ++i) {
                all_position_maps[g][i] = wires_[rewire_gate.in[i]];
            }
            
            // Extract payload shares
            all_payload_shares[g].resize(num_payloads, std::vector<Ring>(vec_size));
            for (size_t p = 0; p < num_payloads; ++p) {
                for (size_t i = 0; i < vec_size; ++i) {
                    all_payload_shares[g][p][i] = wires_[rewire_gate.in[vec_size * (p + 1) + i]];
                }
            }
        }
        
        // Apply permutations for all gates in parallel
        std::vector<std::vector<std::unordered_map<common::utils::wire_t, Ring>>> all_temp_outputs(num_gates);
        
        #pragma omp parallel for
        for (size_t g = 0; g < num_gates; ++g) {
            const auto& rewire_gate = rewire_gates[g];
            size_t vec_size = vec_sizes[g];
            size_t num_payloads = num_payloads_vec[g];
            
            all_temp_outputs[g].resize(num_payloads);
            
            // Apply public permutation based on position_map
            // If inv=0: For each position i: if position_map[i] = idx_perm, then output[idx_perm] = payload[i]
            // If inv=1: For each position i: if position_map[i] = idx_perm, then output[i] = payload[idx_perm] (inverse)
            if (rewire_gate.inv == 0) {
                // Normal rewiring
                for (size_t i = 0; i < vec_size; ++i) {
                    int idx_perm = static_cast<int>(all_position_maps[g][i]);
                    if (idx_perm >= 0 && idx_perm < vec_size) {
                        for (size_t p = 0; p < num_payloads; ++p) {
                            all_temp_outputs[g][p][rewire_gate.outs[vec_size * p + idx_perm]] = all_payload_shares[g][p][i];
                        }
                    }
                }
            } else {
                // Inverse rewiring
                for (size_t i = 0; i < vec_size; ++i) {
                    int idx_perm = static_cast<int>(all_position_maps[g][i]);
                    if (idx_perm >= 0 && idx_perm < vec_size) {
                        for (size_t p = 0; p < num_payloads; ++p) {
                            all_temp_outputs[g][p][rewire_gate.outs[vec_size * p + i]] = all_payload_shares[g][p][idx_perm];
                        }
                    }
                }
            }
        }
        
        // Write outputs to wires (sequential to avoid race conditions)
        for (size_t g = 0; g < num_gates; ++g) {
            size_t num_payloads = num_payloads_vec[g];
            for (size_t p = 0; p < num_payloads; ++p) {
                for (const auto& [wire_id, value] : all_temp_outputs[g][p]) {
                    wires_[wire_id] = value;
                }
            }
        }
    }

    void OnlineEvaluator::publicPermEvaluate(const std::vector<common::utils::SIMDOGate> &public_perm_gates) {
        if (id_ == 0) { return; }
        if (public_perm_gates.empty()) { return; }

        for (const auto& gate : public_perm_gates) {
            // Check if this is a kPublicPerm gate with stored permutation or position map in wires
            if (!gate.permutation.empty() && !gate.permutation[0].empty()) {
                // Old-style: permutation stored in gate (addConstOpMGate)
                size_t vec_size = gate.in.size();
                const auto& perm = gate.permutation[0];
                
                
                // Apply permutation: output[perm[i]] = input[i]
                for (size_t i = 0; i < vec_size; ++i) {
                    size_t perm_idx = perm[i];
                    if (perm_idx >= vec_size) {
                        continue;
                    }
                    wires_[gate.outs[perm_idx]] = wires_[gate.in[i]];
                }
            } else {
                // New-style: position map in first vec_size input wires (addPublicPerm)
                size_t vec_size = gate.vec_size;
                size_t total_inputs = gate.in.size();
                                
                size_t num_payloads = (total_inputs - vec_size) / vec_size;

                // Read position map from first vec_size wires (reconstructed values)
                std::vector<size_t> position_map(vec_size);
                for (size_t i = 0; i < vec_size; ++i) {
                    Ring pos_value = wires_[gate.in[i]];
                    position_map[i] = static_cast<size_t>(pos_value);

                }
                std::cout << std::endl;
                
                
              
                // Apply permutation to each payload
                for (size_t p = 0; p < num_payloads; ++p) {
                    if (gate.inv == 0) {
                        // Forward permutation: output[position_map[i]] = input[i]
                        for (size_t i = 0; i < vec_size; ++i) {
                            size_t idx_perm = position_map[i];
                            
                            size_t out_wire = gate.outs[vec_size * p + idx_perm];
                            size_t in_wire = gate.in[vec_size * (p + 1) + i];
                            
                            wires_[out_wire] = wires_[in_wire];
                        }
                    } else {
                        // Inverse permutation: output[i] = input[position_map[i]]
                        for (size_t i = 0; i < vec_size; ++i) {
                            size_t idx_perm = position_map[i];

                            
                            size_t out_wire = gate.outs[vec_size * p + i];
                            size_t in_wire = gate.in[vec_size * (p + 1) + idx_perm];
                                                                                    
                            wires_[out_wire] = wires_[in_wire];
                        }
                    }
                }
            }
        }
    }

    void OnlineEvaluator::evaluateGatesAtDepth(size_t depth) {
        if (id_ == 0) { return; }

        std::vector<common::utils::FIn2Gate> mult_gates;
        std::vector<common::utils::FIn1Gate> eqz_gates;
        std::vector<common::utils::FIn1Gate> rec_gates;
        std::vector<common::utils::SIMDOGate> shuffle_gates;
        std::vector<common::utils::SIMDOGate> sort_gates;
        std::vector<common::utils::SIMDOGate> public_perm_gates;
        std::vector<common::utils::SIMDOGate> rewire_gates;
        std::vector<common::utils::SIMDOGate> delete_gates;

        // First pass: collect the multi-round gates so their batched handlers can run once.
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case common::utils::GateType::kMul: {
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    mult_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kEqz: {
                    auto *g = static_cast<common::utils::FIn1Gate *>(gate.get());
                    eqz_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kRec: {
                    auto *g = static_cast<common::utils::FIn1Gate *>(gate.get());
                    rec_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kShuffle: {
                    auto *g = static_cast<common::utils::SIMDOGate *>(gate.get());
                    shuffle_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kSort: {
                    auto *g = static_cast<common::utils::SIMDOGate *>(gate.get());
                    sort_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kPublicPerm: {
                    auto *g = static_cast<common::utils::SIMDOGate *>(gate.get());
                    public_perm_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kRewire: {
                    auto *g = static_cast<common::utils::SIMDOGate *>(gate.get());
                    rewire_gates.push_back(*g);
                    break;
                }
                case common::utils::GateType::kDeleteWires: {
                    auto *g = static_cast<common::utils::SIMDOGate *>(gate.get());
                    // Delete wires gates are processed immediately
                    delete_gates.push_back(*g);
                    break;
                }
                default:
                    break;
            }
        }

        if (!mult_gates.empty()) { multEvaluate(mult_gates); }
        if (!eqz_gates.empty()) { eqzEvaluate(eqz_gates); }
        if (!rec_gates.empty()) { recEvaluate(rec_gates); }
        if (!shuffle_gates.empty()) { shuffleEvaluate(shuffle_gates); }
        if (!sort_gates.empty()) { sortEvaluate(sort_gates); }
        if (!public_perm_gates.empty()) { publicPermEvaluate(public_perm_gates); }
        if (!rewire_gates.empty()) { rewireEvaluate(rewire_gates); }
        if (!delete_gates.empty()) { deleteWiresEvaluate(delete_gates); }
        
        // Second pass: handle locally evaluable gates.
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case common::utils::GateType::kAdd: {
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = wires_[g->in1] + wires_[g->in2];
                    break;
                }
                case common::utils::GateType::kSub: {
                    auto *g = static_cast<common::utils::FIn2Gate *>(gate.get());
                    wires_[g->out] = wires_[g->in1] - wires_[g->in2];
                    break;
                }
                case common::utils::GateType::kConstAdd: {
                    auto *g = static_cast<common::utils::ConstOpGate<Ring> *>(gate.get());
                    if (id_ == 1) {
                        wires_[g->out] = wires_[g->in] + g->cval;
                    } else {
                        wires_[g->out] = wires_[g->in];
                    }
                    break;
                }
                case common::utils::GateType::kConstMul: {
                    auto *g = static_cast<common::utils::ConstOpGate<Ring> *>(gate.get());
                    wires_[g->out] = wires_[g->in] * g->cval;
                    break;
                }

                default:
                    break;
            }
        }
    }

    std::vector<Ring> OnlineEvaluator::getOutputs() {
        std::vector<Ring> outvals(circ_.outputs.size());
        if (circ_.outputs.empty()) {
            return outvals;
        }
        if (id_ == 0) { return outvals; }

        // Use the shared reconstruct helper to aggregate output shares across parties.
        std::vector<Ring> shares_to_send(circ_.outputs.size());
        for (size_t i = 0; i < circ_.outputs.size(); ++i) {
            shares_to_send[i] = wires_[circ_.outputs[i]];
        }

        reconstruct(nP_, id_, network_, shares_to_send, outvals, 1, latency_);
        return outvals;
    }

    std::vector<Ring> OnlineEvaluator::evaluateCircuit(const std::unordered_map<common::utils::wire_t, Ring> &inputs) {
        setInputs(inputs);
        for (size_t i = 0; i < circ_.gates_by_level.size(); ++i) {
            evaluateGatesAtDepth(i);
        }
        return getOutputs();
    }

}; // namespace graphdb
