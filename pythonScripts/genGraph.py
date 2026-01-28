import networkx as nx
import random
import itertools
import os
import sys

def powerlaw_directed_graph(n_nodes, m_edges, exponent):
    
    # Initialize an empty directed graph
    graph = nx.DiGraph()
    graph.add_nodes_from(range(n_nodes))


    # Add edges based on a power-law distribution
    edges_added = 0
    while edges_added < m_edges:
        # Choose two nodes randomly
        node_u = random.randint(0, n_nodes - 1)
        node_v = random.randint(0, n_nodes - 1)

        if node_u != node_v and not graph.has_edge(node_u, node_v):
          # Add edge with probability proportional to the degree of node_v raised to the power of -exponent
          probability = 0.75**(-exponent) # Use out-degree to fit powerlaw
          if random.random() < probability:
            graph.add_edge(node_u, node_v)
            edges_added += 1

    return graph


def get_neighbors(graph, vertices):

    neighbors = []
    for v in vertices:
        for u in graph.predecessors(v):
          if u not in vertices:
            neighbors.append(u)
    neighbors = list(set(neighbors))
    neighbors.sort()
    return neighbors



def generate_split_graph(num_parties, graph_size):
    # num_parties = 7
    # graph_size = 100000
    num_vert = graph_size//10
    num_edge = graph_size - num_vert

    subg_num_vert = []
    subg_num_edge = []

    G = powerlaw_directed_graph(num_vert, num_edge, 2.5)

    for i in range(num_parties):
        if i != num_parties-1:
                subg_num_vert.append(num_vert//num_parties)
                subg_num_edge.append(num_edge//num_parties)
        else:
                subg_num_vert.append(num_vert//num_parties + num_vert % num_parties)
                subg_num_edge.append(num_edge//num_parties + num_edge % num_parties)

    subg_idx_vert = [0] + list(itertools.accumulate(subg_num_vert))

    subg_vert_set = [list(G.nodes())[subg_idx_vert[i]:subg_idx_vert[i+1]] for i in range(len(subg_idx_vert)-1)]

    subg_vert_set_bar = []
    subg_edge_set = []

    for i in range(num_parties):
        neighbours = get_neighbors(G, subg_vert_set[i])
  
        subg_vert_set_bar.append(neighbours)
        
        subg_edge_set.append(list(G.in_edges(subg_vert_set[i])))
        
          
    return G, num_vert, num_edge, subg_num_vert, subg_num_edge, subg_idx_vert, subg_vert_set, subg_vert_set_bar, subg_edge_set



def inverse_permutation(permutation):
    n = len(permutation)
    inverse = [0] * n
    for i in range(n):
        if 0 <= permutation[i] < n:
            inverse[permutation[i]] = i
        else:
            return None  # Invalid permutation
    return inverse

  
  
def compose_permutations(perm1, perm2):
    if len(perm1) != len(perm2):
        return None  # Permutations must have the same length

    composed_perm = [0] * len(perm1)
    for i in range(len(perm1)):
        composed_perm[i] = perm2[perm1[i]]
    return composed_perm



def apply_permutation(data, permutation):
    if len(data) != len(permutation):
        return None  # Invalid permutation

    permuted_data = [None] * len(data)
    for i in range(len(data)):
        if 0 <= permutation[i] < len(data):
            permuted_data[i] = data[permutation[i]]
        else:
            return None  # Invalid permutation index

    return permuted_data


def src_val(tup):
  if tup[0] == tup[1]:
    return tup[0] - 0.1
  else:
    return tup[0]

def dst_val(tup):
  if tup[0] == tup[1]:
    return tup[1] + 0.1
  else:
    return tup[1]


def find_source_permutation(dag):
    if not dag:
        return []

    # Create a list of (destination, index) pairs to sort by destination
    source_index_pairs = [(source, dest, index) for index, (source, dest) in enumerate(dag)]

    # Sort by destination
    source_index_pairs.sort(key=src_val)

    # Extract the sorted indices
    sorted_indices = [index for source, dest, index in source_index_pairs]
    
    return sorted_indices


def find_dest_permutation(dag):
    if not dag:
        return []

    # Create a list of (destination, index) pairs to sort by destination
    dest_index_pairs = [(source, dest, index) for index, (source, dest) in enumerate(dag)]

    # Sort by destination
    dest_index_pairs.sort(key=dst_val)

    # Extract the sorted indices
    dest_indices = [index for source, dest, index in dest_index_pairs]
    
    return dest_indices


def genSubgPerm(G, num_vert, num_edge, subg_vert_set, subg_vert_set_bar, subg_edge_set):
    perm_subg =[]
    subg_perms = []
    subg_permd = []
    subg_permv = []
    for i in range(len(subg_vert_set_bar)):
        vert_set = subg_vert_set[i] + subg_vert_set_bar[i]
        perm = vert_set + [j for j in range(num_vert) if j not in vert_set]
        perm_subg.append(perm)


    


    for s_id in range(len(subg_vert_set)):

        vert = subg_vert_set[s_id]
        vert_bar = subg_vert_set_bar[s_id]
        edge_set = subg_edge_set[s_id]

        n = len(vert) + len(vert_bar)
        idx = [i for i in range(n)]
        
        map = dict(zip(vert+vert_bar, idx))

        vert_map = [(map[u], map[u]) for u in vert]
        vert_bar_map = [(map[u], map[u]) for u in vert_bar]
        edge_map = [(map[u], map[v]) for u,v in edge_set]


        dag = vert_map + vert_bar_map + edge_map
        perms = find_source_permutation(dag)
        subg_perms.append(perms)


        dag1 = apply_permutation(dag, perms)

        permd = find_dest_permutation(dag1)
        subg_permd.append(permd)

        dag2 = apply_permutation(dag1, permd)

        inv_permd = inverse_permutation(permd)
        inv_perms = inverse_permutation(perms)

        permv = compose_permutations(inv_perms, inv_permd)
        subg_permv.append(permv)
    
    return perm_subg, subg_perms, subg_permd, subg_permv




def genRandGraph(num_parties, graph_size):
    G, num_vert, num_edge, subg_num_vert, subg_num_edge, subg_idx_vert, subg_vert_set, subg_vert_set_bar, subg_edge_set = generate_split_graph(num_parties, graph_size)
    
    perm_subg, subg_perms, subg_permd, subg_permv = genSubgPerm(G, num_vert, num_edge, subg_vert_set, subg_vert_set_bar, subg_edge_set)
    
    
    
    f = open("graphData.txt", "w")
    
    
    f.write("num parties: "+ str(num_parties)+"\n")
    f.write("graph size: "+ str(graph_size))
    f.write("\n\n\n")
    
    # f.write("graph info:\n")
    # f.write(nx.info(G))
    # f.write("\n\n\n")
    
    f.write("number of vertices in each subgraph:")
    f.write(str(subg_num_vert))
    f.write("\n")
    f.write("number of edges in each subgraph:")
    f.write(str(subg_num_edge))
    f.write("\n\n\n")
    
    f.write("permutations for amortised PnS:")
    f.write(str(perm_subg))
    f.write("\n\n\n")
    
    f.write("permutations for source sort:")
    f.write(str(subg_perms))
    f.write("\n\n\n")
    
    f.write("permutations for dest sort:")
    f.write(str(subg_permd))
    f.write("\n\n\n")
    
    f.write("permutations for vert sort:")
    f.write(str(subg_permv))
    f.write("\n\n\n")
    
    return True



num_parties = int(sys.argv[1])
graph_size = int(sys.argv[2])
genRandGraph(num_parties, graph_size)