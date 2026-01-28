#ifndef GRAPH_GEN_H
#define GRAPH_GEN_H

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <mutex>
#include <random>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../src/utils/types.h"

using namespace std;
using common::utils::Ring;
using common::utils::wire_t;

// Struct to represent a single daglist entry (vertex or edge)
struct DagEntry {
    Ring src;   // Source vertex ID
    Ring dst;   // Destination vertex ID
    Ring isV;   // 1 if vertex entry, 0 if edge entry
    Ring data;  // Data payload
    Ring sigs;  // Position in source-ordered index
    Ring sigv;  // Position in vertex-ordered index
    Ring sigd;  // Position in destination-ordered index
    
    // Constructor
    DagEntry(Ring src_ = 0, Ring dst_ = 0, Ring isV_ = 0, Ring data_ = 0,
             Ring sigs_ = 0, Ring sigv_ = 0, Ring sigd_ = 0)
        : src(src_), dst(dst_), isV(isV_), data(data_), 
          sigs(sigs_), sigv(sigv_), sigd(sigd_) {}
    
    // Convert to vector<Ring> format for compatibility
    vector<Ring> toVector() const {
        return {src, dst, isV, data, sigs, sigv, sigd};
    }
    
    // Create from vector<Ring>
    static DagEntry fromVector(const vector<Ring>& v) {
        if (v.size() < 7) return DagEntry();
        return DagEntry(v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
    }
};

// Struct to represent a complete daglist graph
struct Daglist {
    vector<DagEntry> entries;
    Ring nV;  // Number of vertices
    Ring nE;  // Number of edges
    
    Daglist() : nV(0), nE(0) {}
    
    Daglist(const vector<DagEntry>& entries_) : entries(entries_) {
        nV = 0;
        nE = 0;
        for (const auto& e : entries) {
            if (e.isV == 1) nV++;
            else nE++;
        }
    }
    
    size_t size() const { return entries.size(); }
    bool empty() const { return entries.empty(); }
    
    // Convert to vector<vector<Ring>> format for compatibility
    vector<vector<Ring>> toVectorFormat() const {
        vector<vector<Ring>> result;
        result.reserve(entries.size());
        for (const auto& e : entries) {
            result.push_back(e.toVector());
        }
        return result;
    }
    
    // Create from vector<vector<Ring>>
    static Daglist fromVectorFormat(const vector<vector<Ring>>& vecs) {
        vector<DagEntry> entries;
        entries.reserve(vecs.size());
        for (const auto& v : vecs) {
            entries.push_back(DagEntry::fromVector(v));
        }
        return Daglist(entries);
    }
};

// Struct to represent distributed daglist across multiple parties
struct DistributedDaglist {
    int num_clients;
    int nV;  // Total number of vertices
    int nE;  // Total number of edges
    vector<vector<DagEntry>> VertexLists;  // VertexLists[i] = vertices owned by party i
    vector<vector<DagEntry>> EdgeLists;    // EdgeLists[i] = edges owned by party i
    vector<Ring> VSizes;  // Number of vertices in each VertexList
    vector<Ring> ESizes;  // Number of edges in each EdgeList
    
    // Operation lists for graph modifications
    vector<vector<DagEntry>> InsertV;  // InsertV[i] = vertices to insert for party i
    vector<vector<DagEntry>> InsertE;  // InsertE[i] = edges to insert for party i
    vector<vector<Ring>> VIn;       // VIn[i][j]  1/0 indicating if j-th vertex of party i has an in-edge added
    vector<vector<Ring>> VOut;      // VOut[i][j] 1/0 indicating if j-th vertex of party i has an out-edge added
    vector<vector<Ring>> ChangeV;   // ChangeV[i][j] = data change for j-th vertex of party i
    vector<vector<Ring>> isChangeV; // ChangeV[i][j] = 1/0 indicating if j-th vertex of party i is changed
    vector<vector<Ring>> ChangeE;   // ChangeE[i][j] = data change for j-th edge of party i
    vector<vector<Ring>> isChangeE; // ChangeE[i][j] = 1/0 indicating if j-th edge of party i is changed
    vector<vector<Ring>> isDelV; // ChangeE[i][j] = 1/0 indicating if j-th edge of party i is changed
    vector<vector<Ring>> isDelE; // ChangeE[i][j] = 1/0 indicating if j-th edge of party i is changed

    DistributedDaglist() : num_clients(0), nV(0), nE(0) {}
    
    DistributedDaglist(int np) : num_clients(np), nV(0), nE(0) {
        VertexLists.resize(np);
        EdgeLists.resize(np);
        VSizes.resize(np, 0);
        ESizes.resize(np, 0);
        InsertV.resize(np);
        InsertE.resize(np);
        VIn.resize(np);
        VOut.resize(np);
        ChangeV.resize(np);
        ChangeE.resize(np);
        isChangeV.resize(np);
        isChangeE.resize(np);
        isDelV.resize(np);
        isDelE.resize(np);
    }
    
    DistributedDaglist(int np, int total_nV, int total_nE) 
        : num_clients(np), nV(total_nV), nE(total_nE) {
        VertexLists.resize(np);
        EdgeLists.resize(np);
        VSizes.resize(np, 0);
        ESizes.resize(np, 0);
        InsertV.resize(np);
        InsertE.resize(np);
        VIn.resize(np);
        VOut.resize(np);
        ChangeV.resize(np);
        ChangeE.resize(np);
        isChangeV.resize(np);
        isChangeE.resize(np);
        isDelV.resize(np);
        isDelE.resize(np);
    }
    
    // Get all entries (vertices + edges) for a specific party
    vector<DagEntry> getPartyEntries(int party_id) const {
        vector<DagEntry> result;
        result.reserve(VSizes[party_id] + ESizes[party_id]);
        result.insert(result.end(), VertexLists[party_id].begin(), VertexLists[party_id].end());
        result.insert(result.end(), EdgeLists[party_id].begin(), EdgeLists[party_id].end());
        return result;
    }
    
    // Get total entries across all parties
    size_t totalSize() const {
        return nV + nE;
    }
};

// Struct to represent a single daglist entry with secret-shared wires (for circuit construction)
struct SSDagEntry {
    wire_t src;   // Wire for source vertex ID
    wire_t dst;   // Wire for destination vertex ID
    wire_t isV;   // Wire for isVertex flag (1 if vertex entry, 0 if edge entry)
    wire_t data;  // Wire for data payload
    wire_t sigs;  // Wire for position in source-ordered index
    wire_t sigv;  // Wire for position in vertex-ordered index
    wire_t sigd;  // Wire for position in destination-ordered index
    
    // Constructor
    SSDagEntry(wire_t src_ = 0, wire_t dst_ = 0, wire_t isV_ = 0, wire_t data_ = 0,
               wire_t sigs_ = 0, wire_t sigv_ = 0, wire_t sigd_ = 0)
        : src(src_), dst(dst_), isV(isV_), data(data_), 
          sigs(sigs_), sigv(sigv_), sigd(sigd_) {}
};

// Struct to represent a complete daglist graph with secret-shared wires
struct SSDaglist {
    vector<SSDagEntry> entries;
    Ring nV;  // Number of vertices (plaintext metadata)
    Ring nE;  // Number of edges (plaintext metadata)
    
    SSDaglist() : nV(0), nE(0) {}
    
    SSDaglist(const vector<SSDagEntry>& entries_, Ring nV_, Ring nE_) 
        : entries(entries_), nV(nV_), nE(nE_) {}
    
    size_t size() const { return entries.size(); }
    bool empty() const { return entries.empty(); }
};

// Struct to represent distributed daglist with secret-shared wires across multiple parties
struct SSDistributedDaglist {
    int num_clients;
    int nV;  // Total number of vertices (plaintext metadata)
    int nE;  // Total number of edges (plaintext metadata)
    vector<vector<SSDagEntry>> VertexLists;  // VertexLists[i] = secret-shared vertices owned by party i
    vector<vector<SSDagEntry>> EdgeLists;    // EdgeLists[i] = secret-shared edges owned by party i
    vector<Ring> VSizes;  // Number of vertices in each VertexList
    vector<Ring> ESizes;  // Number of edges in each EdgeList
    
    SSDistributedDaglist() : num_clients(0), nV(0), nE(0) {}
    
    SSDistributedDaglist(int np) : num_clients(np), nV(0), nE(0) {
        VertexLists.resize(np);
        EdgeLists.resize(np);
        VSizes.resize(np, 0);
        ESizes.resize(np, 0);
    }
    
    SSDistributedDaglist(int np, int total_nV, int total_nE) 
        : num_clients(np), nV(total_nV), nE(total_nE) {
        VertexLists.resize(np);
        EdgeLists.resize(np);
        VSizes.resize(np, 0);
        ESizes.resize(np, 0);
    }
    
    // Get all entries (vertices + edges) for a specific party
    vector<SSDagEntry> getPartyEntries(int party_id) const {
        vector<SSDagEntry> result;
        result.reserve(VSizes[party_id] + ESizes[party_id]);
        result.insert(result.end(), VertexLists[party_id].begin(), VertexLists[party_id].end());
        result.insert(result.end(), EdgeLists[party_id].begin(), EdgeLists[party_id].end());
        return result;
    }
    
    // Get total entries across all parties
    size_t totalSize() const {
        return nV + nE;
    }
};

static inline Ring pack_pair(Ring a, Ring b) {
  return (static_cast<uint64_t>(a) << 32) | b;
}

inline vector<pair<Ring, Ring>> generate_scale_free(Ring nV, Ring nE, Ring fixed_seed = 42) {
  // Deterministic (single-threaded) version: previously parallel with shared
  // mutable state could reorder edge generation and break reproducibility.
  vector<pair<Ring, Ring>> edges;
  if (nV == 0 || nE == 0) return edges;

  // degree list for preferential attachment (destination-preferential)
  vector<Ring> degreeList;
  degreeList.reserve(nV * 4 + 10);

  // initialize degreeList with each node once to avoid zero-probability
  for (Ring i = 0; i < nV; ++i) degreeList.push_back(i);

  unordered_set<Ring> seen;
  seen.reserve(nE * 2 + 10);

  std::mt19937_64 rng(fixed_seed);
  std::uniform_int_distribution<Ring> uniV(0, nV - 1);

  while (edges.size() < static_cast<size_t>(nE)) {
    Ring src = uniV(rng);
    if (degreeList.empty()) break;
    Ring dst = degreeList[rng() % degreeList.size()];

    if (src == dst) continue; // avoid self-loops

    Ring key = pack_pair(src, dst);
    if (seen.find(key) != seen.end()) continue; // avoid duplicate

    edges.emplace_back(src, dst);
    seen.insert(key);

    // update degreeList to increase attachment probability
    degreeList.push_back(dst);
    degreeList.push_back(src);
  }

  return edges;
}

inline Daglist build_daglist(Ring nV, const vector<pair<Ring, Ring>>& edges) {
  Ring nE = edges.size();
  Ring total = nV + nE;

  // Handle empty graph case
  if (total == 0) {
    return Daglist();
  }

  // Build in vertex order: vertices, then edges grouped by source
  vector<DagEntry> entries;
  entries.reserve(total);

  // First, add all vertex entries (parallelizable but overhead likely not worth it for small nV)
  for (Ring v = 0; v < nV; ++v) {
    entries.emplace_back(v, v, 1, 0, 0, 0, 0);
  }

  // Group edges by source vertex
  vector<vector<Ring>> outEdges(nV);
  
  // Step 1: Count edges per source (parallel)
  vector<Ring> outDegrees(nV, 0);
  #pragma omp parallel for
  for (Ring i = 0; i < nE; ++i) {
    #pragma omp atomic
    outDegrees[edges[i].first]++;
  }
  
  // Step 2: Reserve space for each vertex's edges
  for (Ring v = 0; v < nV; ++v) {
    outEdges[v].reserve(outDegrees[v]);
  }
  
  // Step 3: Fill edges
  for (Ring i = 0; i < nE; ++i) {
    outEdges[edges[i].first].push_back(i);
  }

  // Compute edge start positions for each vertex
  vector<Ring> edgeStartPos(nV);
  edgeStartPos[0] = nV;
  for (Ring v = 1; v < nV; ++v) {
    edgeStartPos[v] = edgeStartPos[v-1] + outEdges[v-1].size();
  }
  
  // Resize entries to final size so we can fill in parallel
  entries.resize(total);
  
  // Add edges grouped by source - parallel by vertex
  #pragma omp parallel for schedule(dynamic)
  for (Ring v = 0; v < nV; ++v) {
    Ring writePos = edgeStartPos[v];
    for (Ring edgeIdx : outEdges[v]) {
      entries[writePos++] = DagEntry(edges[edgeIdx].first, edges[edgeIdx].second, 0, 0, 0, 0, 0);
    }
  }

  // Build source-ordered index - parallel construction
  vector<Ring> srcOrder(total);
  #pragma omp parallel for schedule(dynamic)
  for (Ring v = 0; v < nV; ++v) {
    Ring writePos = v * (1 + static_cast<Ring>(outEdges[v].size()));
    // Adjust writePos to account for vertices with different out-degrees
    writePos = v; // Start with vertex position
    for (Ring u = 0; u < v; ++u) {
      writePos += outEdges[u].size();
    }
    
    srcOrder[writePos] = v;
    Ring edgeStart = edgeStartPos[v];
    for (Ring i = 0; i < outEdges[v].size(); ++i) {
      srcOrder[writePos + i + 1] = edgeStart + i;
    }
  }

  // Build dest-ordered index - first group edges by destination (parallel)
  vector<vector<Ring>> inEdges(nV);
  
  // Count in-degrees (parallel)
  vector<Ring> inDegrees(nV, 0);
  #pragma omp parallel for
  for (Ring i = nV; i < total; ++i) {
    Ring dst = entries[i].dst;
    #pragma omp atomic
    inDegrees[dst]++;
  }
  
  // Reserve space
  for (Ring v = 0; v < nV; ++v) {
    inEdges[v].reserve(inDegrees[v]);
  }
  
  // Fill in-edges (sequential to maintain order)
  for (Ring edgePos = nV; edgePos < total; ++edgePos) {
    Ring dst = entries[edgePos].dst;
    inEdges[dst].push_back(edgePos);
  }

  // Build destination-ordered array - parallel construction
  vector<Ring> dstOrder(total);
  #pragma omp parallel for schedule(dynamic)
  for (Ring v = 0; v < nV; ++v) {
    // Calculate write position for vertex v in dstOrder
    Ring writePos = v;
    for (Ring u = 0; u < v; ++u) {
      writePos += inEdges[u].size();
    }
    
    // Write in-edges first
    for (Ring i = 0; i < inEdges[v].size(); ++i) {
      dstOrder[writePos + i] = inEdges[v][i];
    }
    // Then write vertex entry
    dstOrder[writePos + inEdges[v].size()] = v;
  }

  // Build reverse maps: position of each index in each ordering (parallel)
  vector<Ring> pos_in_src(total), pos_in_dst(total), pos_in_vert(total);
  
  #pragma omp parallel for
  for (Ring i = 0; i < total; ++i) {
    pos_in_src[srcOrder[i]] = i;
    pos_in_dst[dstOrder[i]] = i;
    pos_in_vert[i] = i;
  }

  // Fill sigs, sigv, sigd fields (parallel)
  #pragma omp parallel for
  for (Ring idx = 0; idx < total; ++idx) {
    entries[idx].sigs = pos_in_src[idx];
    entries[idx].sigv = pos_in_vert[idx];
    entries[idx].sigd = pos_in_dst[idx];
    if (idx == 0) {
      entries[idx].data = 1; // Vertex data initialized to 1 for first vertex
    }
  }

  return Daglist(entries);
}

// Distribute daglist_graph across np parties
// Input: daglist_graph (nV vertices + nE edges), np (number of parties)
// Output: DistributedDaglist with np party graphs, where party i owns V/np vertices and their outgoing edges
inline DistributedDaglist distribute_daglist(const Daglist& daglist_graph, int np) {
  if (np <= 0) {
    throw std::invalid_argument("Number of parties must be positive");
  }
  
  Ring nV = daglist_graph.nV;
  Ring nE = daglist_graph.nE;
  
  // Handle empty graph case
  if (nV == 0 && nE == 0) {
    return DistributedDaglist(np, 0, 0);
  }

  if (nV == 0) {
    return DistributedDaglist(np, 0, nE);
  }

  // Calculate vertices per party
  Ring verts_per_party = (nV + np - 1) / np;  // Ceiling division
  
  // Initialize result
  DistributedDaglist result(np, nV, nE);
  
  // Build a map from vertex ID to party assignment
  // Party i owns vertices [i * verts_per_party, (i+1) * verts_per_party)
  auto get_party = [&](Ring vertex_id) -> int {
    Ring party = vertex_id / verts_per_party;
    return std::min(static_cast<Ring>(party), static_cast<Ring>(np - 1));
  };

  // Distribute entries
  for (const auto& entry : daglist_graph.entries) {
    if (entry.isV == 1) {
      // This is a vertex entry - assign to the party that owns this vertex
      int party = get_party(entry.src);
      result.VertexLists[party].push_back(entry);
      result.VSizes[party]++;
    } else {
      // This is an edge entry - assign to the party that owns the source vertex
      int party = get_party(entry.src);
      result.EdgeLists[party].push_back(entry);
      result.ESizes[party]++;
    }
  }

  return result;
}

// Build a complete Daglist from a DistributedDaglist
// Reconstructs the full graph by merging all party entries
inline Daglist build_daglist_from_distributed(const DistributedDaglist& dist_daglist) {
  vector<DagEntry> all_entries;
  all_entries.reserve(dist_daglist.nV + dist_daglist.nE);
  
  // Collect all vertices first, then all edges (maintaining vertex-order structure)
  for (int i = 0; i < dist_daglist.num_clients; ++i) {
    all_entries.insert(all_entries.end(), 
                      dist_daglist.VertexLists[i].begin(), 
                      dist_daglist.VertexLists[i].end());
  }
  
  for (int i = 0; i < dist_daglist.num_clients; ++i) {
    all_entries.insert(all_entries.end(), 
                      dist_daglist.EdgeLists[i].begin(), 
                      dist_daglist.EdgeLists[i].end());
  }
  
  return Daglist(all_entries);
}

// Build a complete SSDaglist from an SSDistributedDaglist
// Reconstructs the full graph with secret-shared wires by merging all party entries
inline SSDaglist build_ssdaglist_from_distributed(const SSDistributedDaglist& dist_daglist) {
  vector<SSDagEntry> all_entries;
  all_entries.reserve(dist_daglist.nV + dist_daglist.nE);
  
  // Collect all vertices first, then all edges (maintaining vertex-order structure)
  for (int i = 0; i < dist_daglist.num_clients; ++i) {
    all_entries.insert(all_entries.end(), 
                      dist_daglist.VertexLists[i].begin(), 
                      dist_daglist.VertexLists[i].end());
  }
  
  for (int i = 0; i < dist_daglist.num_clients; ++i) {
    all_entries.insert(all_entries.end(), 
                      dist_daglist.EdgeLists[i].begin(), 
                      dist_daglist.EdgeLists[i].end());
  }
  
  return SSDaglist(all_entries, dist_daglist.nV, dist_daglist.nE);
}

// Create circuit inputs from a Daglist
// Populates the inputs map with values from the daglist entries based on wire assignments
// wire_idx tracks the current position in input_wires
inline void set_daglist_inputs(
    const Daglist& daglist,
    const std::vector<wire_t>& input_wires,
    std::unordered_map<wire_t, Ring>& inputs,
    size_t& wire_idx) {
  
  size_t vec_size = daglist.size();
  
  // Initialize all daglist field values
  std::vector<Ring> src_values(vec_size);
  std::vector<Ring> dst_values(vec_size);
  std::vector<Ring> isV_values(vec_size);
  std::vector<Ring> data_values(vec_size);
  std::vector<Ring> sigs_values(vec_size);
  std::vector<Ring> sigv_values(vec_size);
  std::vector<Ring> sigd_values(vec_size);
  
  for (size_t i = 0; i < vec_size; ++i) {
    src_values[i] = daglist.entries[i].src;
    dst_values[i] = daglist.entries[i].dst;
    isV_values[i] = daglist.entries[i].isV;
    data_values[i] = daglist.entries[i].data;
    sigs_values[i] = daglist.entries[i].sigs;
    sigv_values[i] = daglist.entries[i].sigv;
    sigd_values[i] = daglist.entries[i].sigd;
  }
  
  // Set all daglist fields as inputs (7 * vec_size inputs total)
  // The order must match the circuit wire creation order
  
  // src inputs
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = src_values[i];
  }
  
  // dst inputs
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = dst_values[i];
  }
  
  // isV inputs
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = isV_values[i];
  }
  
  // data inputs
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = data_values[i];
  }
  
  // sigs inputs (position map for propagate)
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = sigs_values[i];
  }
  
  // sigv inputs (position map for gather)
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = sigv_values[i];
  }
  
  // sigd inputs (position map for permlist)
  for (size_t i = 0; i < vec_size && wire_idx < input_wires.size(); ++i, ++wire_idx) {
    inputs[input_wires[wire_idx]] = sigd_values[i];
  }
}

// Create circuit inputs from a Daglist (overload with initial wire_idx = 0)
inline void set_daglist_inputs(
    const Daglist& daglist,
    const std::vector<wire_t>& input_wires,
    std::unordered_map<wire_t, Ring>& inputs) {
  size_t wire_idx = 0;
  set_daglist_inputs(daglist, input_wires, inputs, wire_idx);
}


// Generate random inputs for changing data of vertices and edges in a DistributedDaglist
// Returns a DistributedDaglist with ChangeV and ChangeE populated
// ChangeV[i][j] = new_data for j-th vertex of party i (0 if no change)
// ChangeE[i][j] = new_data for j-th edge of party i (0 if no change)
// Input: dist_daglist - the distributed daglist
//        num_changes - total number of entries to modify across all parties
//        seed - random seed for reproducibility
inline DistributedDaglist generate_random_data_updates(const DistributedDaglist& dist_daglist, 
                                                       Ring num_changes, 
                                                       Ring seed = 42) {
  
  // Copy the input distributed daglist
  DistributedDaglist result = dist_daglist;
  
  // Initialize ChangeV and ChangeE with zeros
  for (int p = 0; p < result.num_clients; ++p) {
    result.ChangeV[p].resize(result.VSizes[p], 0);
    result.ChangeE[p].resize(result.ESizes[p], 0);
    // Initialize isChange vectors (0 = no change, 1 = changed)
    result.isChangeV[p].resize(result.VSizes[p], 0);
    result.isChangeE[p].resize(result.ESizes[p], 0);
  }
  
  if (num_changes == 0 || dist_daglist.totalSize() == 0) {
    return result;
  }
  
  std::mt19937_64 rng(seed);
  
  // Collect all possible entries to update (party_id, is_vertex, local_idx)
  struct EntryLocation {
    int party_id;
    bool is_vertex;
    Ring local_idx;
  };
  
  std::vector<EntryLocation> all_entries;
  for (int p = 0; p < dist_daglist.num_clients; ++p) {
    for (Ring i = 0; i < dist_daglist.VSizes[p]; ++i) {
      all_entries.push_back({p, true, i});
    }
    for (Ring i = 0; i < dist_daglist.ESizes[p]; ++i) {
      all_entries.push_back({p, false, i});
    }
  }
  
  if (all_entries.empty()) {
    return result;
  }
  
  // Sample num_changes entries randomly without replacement
  Ring actual_changes = std::min(num_changes, static_cast<Ring>(all_entries.size()));
  std::shuffle(all_entries.begin(), all_entries.end(), rng);
  
  // Generate random data values and populate ChangeV/ChangeE
  std::uniform_int_distribution<Ring> data_dist(1, 1000);
  
  for (Ring i = 0; i < actual_changes; ++i) {
    const auto& entry = all_entries[i];
    Ring new_data = data_dist(rng);
    
    if (entry.is_vertex) {
      result.ChangeV[entry.party_id][entry.local_idx] = new_data;
      result.isChangeV[entry.party_id][entry.local_idx] = 1;
    } else {
      result.ChangeE[entry.party_id][entry.local_idx] = new_data;
      result.isChangeE[entry.party_id][entry.local_idx] = 1;
    }
  }
  
  return result;
}

// Generate random deletion tags for vertices and edges in a DistributedDaglist
// Returns a DistributedDaglist with isDelV and isDelE populated
// isDelV[i][j] = 1 if j-th vertex of party i should be deleted, 0 otherwise
// isDelE[i][j] = 1 if j-th edge of party i should be deleted, 0 otherwise
// Input: dist_daglist - the distributed daglist
//        num_deletes - total number of entries to delete across all parties
//        seed - random seed for reproducibility
inline DistributedDaglist generate_random_entry_deletes(const DistributedDaglist& dist_daglist, 
                                                        Ring num_deletes, 
                                                        Ring seed = 42) {
  
  // Copy the input distributed daglist
  DistributedDaglist result = dist_daglist;
  
  // Initialize isDelV and isDelE with zeros
  for (int p = 0; p < result.num_clients; ++p) {
    result.isDelV[p].resize(result.VSizes[p], 0);
    result.isDelE[p].resize(result.ESizes[p], 0);
  }
  
  if (num_deletes == 0 || dist_daglist.totalSize() == 0) {
    return result;
  }
  
  std::mt19937_64 rng(seed);
  
  // Collect all possible entries to delete (party_id, is_vertex, local_idx)
  struct EntryLocation {
    int party_id;
    bool is_vertex;
    Ring local_idx;
  };
  
  std::vector<EntryLocation> all_entries;
  for (int p = 0; p < dist_daglist.num_clients; ++p) {
    for (Ring i = 0; i < dist_daglist.VSizes[p]; ++i) {
      all_entries.push_back({p, true, i});
    }
    for (Ring i = 0; i < dist_daglist.ESizes[p]; ++i) {
      all_entries.push_back({p, false, i});
    }
  }
  
  if (all_entries.empty()) {
    return result;
  }
  
  // Sample num_deletes entries randomly without replacement
  Ring actual_deletes = std::min(num_deletes, static_cast<Ring>(all_entries.size()));
  std::shuffle(all_entries.begin(), all_entries.end(), rng);
  
  // Mark entries for deletion
  for (Ring i = 0; i < actual_deletes; ++i) {
    const auto& entry = all_entries[i];
    
    if (entry.is_vertex) {
      result.isDelV[entry.party_id][entry.local_idx] = 1;
    } else {
      result.isDelE[entry.party_id][entry.local_idx] = 1;
    }
  }
  
  return result;
}

// Generate random vertices to be inserted into a DistributedDaglist
// Returns a new DistributedDaglist with InsertV populated with random vertex entries
// Vertices are distributed across clients proportionally to their current sizes
// Input: dist_daglist - the distributed daglist (used for distribution reference)
//        num_inserts - total number of vertices to insert across all clients
//        seed - random seed for reproducibility
inline DistributedDaglist generate_random_vertices_to_insert(const DistributedDaglist& dist_daglist, 
                                                             Ring num_inserts, 
                                                             Ring seed = 42) {
  
  // Create a new DistributedDaglist with same structure but populated InsertV
  DistributedDaglist result(dist_daglist.num_clients, dist_daglist.nV, dist_daglist.nE);
  
  // Initialize InsertV for each client
  for (int c = 0; c < result.num_clients; ++c) {
    result.InsertV[c].clear();
  }
  
  if (num_inserts == 0) {
    return result;
  }
  
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<Ring> data_dist(1, 1000);
  
  // Distribute num_inserts across clients proportionally to their current vertex sizes
  std::vector<Ring> inserts_per_client(result.num_clients, 0);
  Ring total_verts = 0;
  for (int c = 0; c < result.num_clients; ++c) {
    total_verts += dist_daglist.VSizes[c];
  }
  
  // Assign vertices to clients proportionally
  Ring remaining = num_inserts;
  for (int c = 0; c < result.num_clients; ++c) {
    if (total_verts > 0) {
      Ring share = (num_inserts * dist_daglist.VSizes[c]) / total_verts;
      inserts_per_client[c] = share;
      remaining -= share;
    }
  }
  
  // Distribute any remaining vertices to the first client
  if (remaining > 0 && result.num_clients > 0) {
    inserts_per_client[0] += remaining;
  } else if (remaining > 0 && result.num_clients == 0) {
    // No clients available - shouldn't happen in normal usage
    throw std::runtime_error("Cannot insert vertices - no clients available");
  }
  
  // Generate random vertex entries for each client
  Ring global_vertex_id = dist_daglist.nV;  // New vertices start after existing ones
  
  for (int c = 0; c < result.num_clients; ++c) {
    for (Ring i = 0; i < inserts_per_client[c]; ++i) {
      Ring src = global_vertex_id;
      Ring dst = global_vertex_id;
      Ring isV = 1;  // Mark as vertex entry
      Ring data = data_dist(rng);
      
      // sigs, sigv, sigd will be computed during circuit generation
      DagEntry new_vertex(src, dst, isV, data, 0, 0, 0);
      result.InsertV[c].push_back(new_vertex);
      
      global_vertex_id++;
    }
  }
  
  // Copy over existing graph info for reference
  for (int c = 0; c < result.num_clients; ++c) {
    result.VertexLists[c] = dist_daglist.VertexLists[c];
    result.EdgeLists[c] = dist_daglist.EdgeLists[c];
    result.VSizes[c] = dist_daglist.VSizes[c];
    result.ESizes[c] = dist_daglist.ESizes[c];
  }
  
  return result;
}

// Generate random edges to be inserted into a DistributedDaglist
// Returns a DistributedDaglist with InsertE populated with random edge entries
// Edges are distributed across clients proportionally to their current sizes
// Input: dist_daglist - the distributed daglist (used for distribution reference)
//        num_inserts - total number of edges to insert across all clients
//        seed - random seed for reproducibility
inline DistributedDaglist generate_random_edges_to_insert(const DistributedDaglist& dist_daglist, 
                                                          Ring num_inserts, 
                                                          Ring seed = 42) {
  
  // Create a new DistributedDaglist with same structure but populated InsertE
  DistributedDaglist result(dist_daglist.num_clients, dist_daglist.nV, dist_daglist.nE);
  
  // Initialize InsertE for each client
  for (int c = 0; c < result.num_clients; ++c) {
    result.InsertE[c].clear();
  }
  
  if (num_inserts == 0 || dist_daglist.nV == 0) {
    // Cannot insert edges if no vertices exist
    return result;
  }
  
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<Ring> vertex_dist(0, dist_daglist.nV - 1);
  std::uniform_int_distribution<Ring> data_dist(1, 1000);
  
  // Distribute num_inserts across clients proportionally to their current edge sizes
  std::vector<Ring> inserts_per_client(result.num_clients, 0);
  Ring total_edges = 0;
  for (int c = 0; c < result.num_clients; ++c) {
    total_edges += dist_daglist.ESizes[c];
  }
  
  // Assign edges to clients proportionally
  Ring remaining = num_inserts;
  for (int c = 0; c < result.num_clients; ++c) {
    if (total_edges > 0) {
      Ring share = (num_inserts * dist_daglist.ESizes[c]) / total_edges;
      inserts_per_client[c] = share;
      remaining -= share;
    } else if (result.num_clients > 0 && c == 0) {
      // If no edges exist, give all to first client
      inserts_per_client[c] = num_inserts;
      remaining = 0;
    }
  }
  
  // Distribute any remaining edges to the first client
  if (remaining > 0 && result.num_clients > 0) {
    inserts_per_client[0] += remaining;
  }
  
  // Initialize VIn and VOut with zeros (each entry has size nV for all vertices)
  for (int c = 0; c < result.num_clients; ++c) {
    result.VIn[c].resize(dist_daglist.nV, 0);
    result.VOut[c].resize(dist_daglist.nV, 0);
  }
  
  // Generate random edge entries for each client
  unordered_set<Ring> seen;
  
  for (int c = 0; c < result.num_clients; ++c) {
    for (Ring i = 0; i < inserts_per_client[c]; ++i) {
      Ring src, dst;
      // Generate random edge without duplicates
      do {
        src = vertex_dist(rng);
        dst = vertex_dist(rng);
      } while (src == dst || seen.find(pack_pair(src, dst)) != seen.end());
      
      seen.insert(pack_pair(src, dst));
      Ring isV = 0;  // Mark as edge entry
      Ring data = data_dist(rng);
      
      // sigs, sigv, sigd will be computed during circuit generation
      DagEntry new_edge(src, dst, isV, data, 0, 0, 0);
      result.InsertE[c].push_back(new_edge);
      
      // Mark source vertex as having an outgoing edge
      result.VOut[c][src] = 1;
      
      // Mark destination vertex as having an incoming edge
      result.VIn[c][dst] = 1;
    }
  }
  
  // Copy over existing graph info for reference
  for (int c = 0; c < result.num_clients; ++c) {
    result.VertexLists[c] = dist_daglist.VertexLists[c];
    result.EdgeLists[c] = dist_daglist.EdgeLists[c];
    result.VSizes[c] = dist_daglist.VSizes[c];
    result.ESizes[c] = dist_daglist.ESizes[c];
  }
  
  return result;
}





#endif // GRAPH_GEN_H
