import networkx as nx
import sys
from queue import Queue
import ast
from time import perf_counter
import copy

import hierarchical_clusters as h


def LCA(u, v, agglom_nodes, agglom_tree, agglom_tree_level):
    pu = agglom_nodes[u]
    pv = agglom_nodes[v]
    while pu != pv:
        if agglom_tree_level[pu] > agglom_tree_level[pv]:
            pu = agglom_tree[pu]
        elif agglom_tree_level[pu] < agglom_tree_level[pv]:
            pv = agglom_tree[pv]
        else:
            pu = agglom_tree[pu]
            pv = agglom_tree[pv]
    return pv



def agglom_tree_size_func(agglom_tree, top_down_tree, leaf_inner_indicator):
    agglom_tree_size = [0 for _ in range(len(agglom_tree))]
    peak = 0
    for i in range(len(agglom_tree)):
        if agglom_tree[i] == -1:
            peak = i
            break

    stk = list()
    stk.append(peak)
    while stk:
        p = stk[-1]
        if agglom_tree_size[p] >= 2:
            stk.pop()
            continue
        else:
            child0 = top_down_tree[p][0]
            child1 = top_down_tree[p][1]
            if leaf_inner_indicator[p][0] == 0 and leaf_inner_indicator[p][1] == 0:
                agglom_tree_size[p] = 2
            elif leaf_inner_indicator[p][0] == 0 and leaf_inner_indicator[p][1] == 1:
                if agglom_tree_size[child1] >= 2:
                    agglom_tree_size[p] = agglom_tree_size[child1] + 1
                else:
                    stk.append(child1)
            elif leaf_inner_indicator[p][0] == 1 and leaf_inner_indicator[p][1] == 0:
                if agglom_tree_size[child0] >= 2:
                    agglom_tree_size[p] = agglom_tree_size[child0] + 1
                else:
                    stk.append(child0)
            else:
                if agglom_tree_size[child0] >= 2 and agglom_tree_size[child1] >= 2:
                    agglom_tree_size[p] = agglom_tree_size[child0] + agglom_tree_size[child1]
                else:
                    if agglom_tree_size[child0] < 2:
                        stk.append(child0)
                    if agglom_tree_size[child1] < 2:
                        stk.append(child1)
    return agglom_tree_size


def agglom_level_func(agglom_tree, top_down_tree, leaf_inner_indicator):
    agglom_node_level = [0 for _ in range(len(agglom_tree) + 1)]
    agglom_tree_level = [0 for _ in range(len(agglom_tree))]
    stk = list()
    stk.append(len(agglom_tree) - 1)
    while stk:
        p = stk.pop()
        plevel = agglom_tree_level[p]
        child0 = top_down_tree[p][0]
        child1 = top_down_tree[p][1]
        if leaf_inner_indicator[p][0] == 0:
            agglom_node_level[child0] = plevel + 1
        else:
            agglom_tree_level[child0] = plevel + 1
            stk.append(child0)
        if leaf_inner_indicator[p][1] == 0:
            agglom_node_level[child1] = plevel + 1
        else:
            agglom_tree_level[child1] = plevel + 1
            stk.append(child1)
    return agglom_node_level, agglom_tree_level


def topo_sim(n0, n1, F):
    w = F[n0][n1]['weight']
    size0 = F.nodes[n0]['size']
    size1 = F.nodes[n1]['size']
    return w / (size0 * size1)


def graph_weighting(linkage, g, feats):
    if linkage == 'label_edge':
        for e in g.edges:
            n0 = e[0]
            n1 = e[1]
            if feats[n0] and feats[n1]:
                g[n0][n1]['weight'] = h.LAMDA
            else:
                g[n0][n1]['weight'] = (1 - h.LAMDA)
    elif linkage == 'unlabel':
        for e in g.edges:
            n0 = e[0]
            n1 = e[1]
            g[n0][n1]['weight'] = 1.0
    else:
        print("wrong linkage\n")
        sys.exit()


def nn_find(chain_top, F, diff, candidates=None):
    next_nn = -1
    sim_max = -1.0
    if candidates is None:
        for c in F:
            if c == chain_top:
                continue
            sim = topo_sim(chain_top, c, F)
            if sim_max < (sim - diff):
                next_nn = c
                sim_max = sim
            elif (sim - diff) <= sim_max <= (sim + diff) and c < next_nn:
                next_nn = c
    else:
        for c in candidates:
            if c == chain_top:
                continue
            sim = topo_sim(chain_top, c, F)
            if sim_max < (sim - diff):
                next_nn = c
                sim_max = sim
            elif (sim - diff) <= sim_max <= (sim + diff) and c < next_nn:
                next_nn = c
    return next_nn, sim_max


def nnc_with_candidates(g, feats, linkage):
    seed = 0
    seed_merged = [False for _ in range(2 * g.number_of_nodes() - 1)]

    num = g.number_of_nodes()
    F = g.copy()
    graph_weighting(linkage, F, feats)

    max_size = 1
    max2_size = 1

    for n in F:
        F.nodes[n]['size'] = 1.0

    agglom_nodes = [-1 for _ in range(g.number_of_nodes())]
    agglom_tree = []
    top_down_tree = []
    leaf_inner_indicator = []

    S = []
    while F.number_of_nodes() > 1:
        if len(S) == 0:
            S.append(seed)
        chain_head = S[len(S) - 1]

        # generate the candidate nodes for nearest neighbor
        candidates = []
        for nbr in F[chain_head]:
            candidates.append(nbr)

        diff = 0.09 / ((max_size ** 2) * (max2_size ** 2))
        # diff is set to 9e-2 for LAMDA with one decimal place

        next_nn, _ = nn_find(chain_head, F, diff, candidates)

        if len(S) == 1 or next_nn != S[len(S) - 2]:
            S.append(next_nn)
        else:
            # merge next_nn and S[len(S)-2]
            nn0 = S[len(S) - 1]
            nn1 = S[len(S) - 2]

            agglom_tree.append(-1)
            if nn0 < g.number_of_nodes():
                agglom_nodes[nn0] = num - g.number_of_nodes()
            else:
                agglom_tree[nn0 - g.number_of_nodes()] = num - g.number_of_nodes()
            if nn1 < g.number_of_nodes():
                agglom_nodes[nn1] = num - g.number_of_nodes()
            else:
                agglom_tree[nn1 - g.number_of_nodes()] = num - g.number_of_nodes()

            if nn0 < g.number_of_nodes() and nn1 < g.number_of_nodes():
                top_down_tree.append([nn1, nn0])
                leaf_inner_indicator.append([0, 0])
            elif nn0 < g.number_of_nodes() <= nn1:
                top_down_tree.append([nn1 - g.number_of_nodes(), nn0])
                leaf_inner_indicator.append([1, 0])
            elif nn1 < g.number_of_nodes() <= nn0:
                top_down_tree.append([nn1, nn0 - g.number_of_nodes()])
                leaf_inner_indicator.append([0, 1])
            else:
                top_down_tree.append([nn1 - g.number_of_nodes(), nn0 - g.number_of_nodes()])
                leaf_inner_indicator.append([1, 1])

            nbr_dict = dict()
            for nbr in F[nn0]:
                if nbr != nn1:
                    if nbr in nbr_dict:
                        nbr_dict[nbr] += F[nn0][nbr]['weight']
                    else:
                        nbr_dict[nbr] = F[nn0][nbr]['weight']
            for nbr in F[nn1]:
                if nbr != nn0:
                    if nbr in nbr_dict:
                        nbr_dict[nbr] += F[nn1][nbr]['weight']
                    else:
                        nbr_dict[nbr] = F[nn1][nbr]['weight']

            merge_size = F.nodes[nn0]['size'] + F.nodes[nn1]['size']

            if merge_size > max_size:
                max2_size = max_size
                max_size = merge_size
            elif merge_size == max_size:
                max2_size = merge_size
            if max_size + max2_size > g.number_of_nodes():
                max2_size = g.number_of_nodes() - max_size

            F.add_node(num, size=merge_size)
            F.remove_node(nn0)
            F.remove_node(nn1)

            for nbr in nbr_dict:
                F.add_edge(num, nbr, weight=nbr_dict[nbr])

            S.pop()
            S.pop()
            num += 1

            seed_merged[nn0] = True
            seed_merged[nn1] = True
            while seed_merged[seed]:
                seed += 1
    return agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator


def LocalHAC_L(G, q, feats, agglom_nodes_ul, agglom_tree_ul, top_down_tree_ul, leaf_inner_indicator_ul, L):
    tn = agglom_nodes_ul[q]
    pre_tn = q
    pre_tn_leaf = 0
    com = [q]

    mid_tn = []

    agglom_nodes = copy.deepcopy(agglom_nodes_ul)
    agglom_tree = copy.deepcopy(agglom_tree_ul)
    top_down_tree = copy.deepcopy(top_down_tree_ul)
    leaf_inner_indicator = copy.deepcopy(leaf_inner_indicator_ul)

    for _ in range(L):
        if tn == -1:
            break
        tn_leaves = []
        child0 = top_down_tree_ul[tn][0]
        child1 = top_down_tree_ul[tn][1]
        if child0 != pre_tn or pre_tn_leaf != leaf_inner_indicator_ul[tn][0]:
            if leaf_inner_indicator_ul[tn][0] == 0:
                tn_leaves.append(child0)
            else:
                stk = [child0]
                while stk:
                    p = stk.pop()
                    mid_tn.append(p)
                    c0 = top_down_tree_ul[p][0]
                    c1 = top_down_tree_ul[p][1]
                    if leaf_inner_indicator_ul[p][0] == 0:
                        tn_leaves.append(c0)
                    else:
                        stk.append(c0)
                    if leaf_inner_indicator_ul[p][1] == 0:
                        tn_leaves.append(c1)
                    else:
                        stk.append(c1)
        else:
            if leaf_inner_indicator_ul[tn][1] == 0:
                tn_leaves.append(child1)
            else:
                stk = [child1]
                while stk:
                    p = stk.pop()
                    mid_tn.append(p)
                    c0 = top_down_tree_ul[p][0]
                    c1 = top_down_tree_ul[p][1]
                    if leaf_inner_indicator_ul[p][0] == 0:
                        tn_leaves.append(c0)
                    else:
                        stk.append(c0)
                    if leaf_inner_indicator_ul[p][1] == 0:
                        tn_leaves.append(c1)
                    else:
                        stk.append(c1)
        com += tn_leaves
        mid_tn.append(tn)
        pre_tn = tn
        tn = agglom_tree_ul[tn]
        pre_tn_leaf = 1

    tn = pre_tn

    # re-hierarchical-clustering the sub-tree under tn
    S = G.subgraph(com)

    com_map = dict()
    reverse_com_map = [0 for _ in range(len(com))]

    cnt = 0
    for n in com:
        com_map[n] = cnt
        reverse_com_map[cnt] = n
        cnt += 1
    Subg = nx.Graph()
    for e in S.edges:
        Subg.add_edge(com_map[e[0]], com_map[e[1]])

    subfeats = [False for _ in range(len(com))]
    for n in com:
        if feats[n]:
            subfeats[com_map[n]] = True
    sub_agglom_nodes, sub_agglom_tree, sub_top_down_tree, sub_leaf_inner_indicator = \
        nnc_with_candidates(Subg, subfeats, 'label_edge')

    for cnt in range(len(sub_agglom_nodes)):  # cnt denotes nodes in the subgraph
        node = reverse_com_map[cnt]  # node in the original graph
        if sub_agglom_tree[sub_agglom_nodes[cnt]] == -1:
            par = tn
        else:
            par = mid_tn[sub_agglom_nodes[cnt]]
        agglom_nodes[node] = par

    for cnt in range(len(sub_agglom_tree)):  # cnt denotes a mid node in the sub-tree
        p_mid = sub_agglom_tree[cnt]
        if p_mid == -1:
            mid = tn
        else:
            mid = mid_tn[cnt]
        # mid denote the mid-node in the original graph

        if sub_agglom_tree[cnt] >= 0:
            if sub_agglom_tree[sub_agglom_tree[cnt]] == -1:  # mid-node "cnt" under "tn"
                agglom_tree[mid] = tn
            else:
                agglom_tree[mid] = mid_tn[sub_agglom_tree[cnt]]

        leaf_inner_indicator[mid][0] = sub_leaf_inner_indicator[cnt][0]
        leaf_inner_indicator[mid][1] = sub_leaf_inner_indicator[cnt][1]

        if leaf_inner_indicator[mid][0] == 0:
            top_down_tree[mid][0] = reverse_com_map[sub_top_down_tree[cnt][0]]
        else:
            top_down_tree[mid][0] = mid_tn[sub_top_down_tree[cnt][0]]

        if leaf_inner_indicator[mid][1] == 0:
            top_down_tree[mid][1] = reverse_com_map[sub_top_down_tree[cnt][1]]
        else:
            top_down_tree[mid][1] = mid_tn[sub_top_down_tree[cnt][1]]
    return agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator


def LocalHAC(G, q, feats, agglom_nodes_ul, agglom_tree_ul, top_down_tree_ul, leaf_inner_indicator_ul):
    agglom_tree_size = agglom_tree_size_func(agglom_tree_ul, top_down_tree_ul, leaf_inner_indicator_ul)
    _, agglom_tree_level = agglom_level_func(agglom_tree_ul, top_down_tree_ul, leaf_inner_indicator_ul)

    agglom_tree_Lcut = [0.0 for _ in range(len(agglom_tree_ul))]
    for u in G:
        for v in G[u]:
            if v > u and feats[u] and feats[v]:
                cp = LCA(u, v, agglom_nodes_ul, agglom_tree_ul, agglom_tree_level)
                agglom_tree_Lcut[cp] += 1
    L = 1
    temp_L = 1
    temp_tn = agglom_nodes_ul[q]
    tn = agglom_nodes_ul[q]

    max_label_edge_cutr = 0
    sum_label_edge_level = 0
    while temp_tn != -1:
        sum_label_edge_level += agglom_tree_level[temp_tn] * agglom_tree_Lcut[temp_tn]
        temp_label_edge_cutr = sum_label_edge_level / agglom_tree_size[temp_tn]
        if temp_label_edge_cutr > max_label_edge_cutr:
            max_label_edge_cutr = temp_label_edge_cutr
            L = temp_L
            tn = temp_tn
        temp_tn = agglom_tree_ul[temp_tn]
        temp_L += 1
    
    if L > 1:
        agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator = \
            LocalHAC_L(G, q, feats, agglom_nodes_ul, agglom_tree_ul,
                       top_down_tree_ul, leaf_inner_indicator_ul, L)
    else:
        agglom_nodes = agglom_nodes_ul
        agglom_tree = agglom_tree_ul
        top_down_tree = top_down_tree_ul
        leaf_inner_indicator = leaf_inner_indicator_ul
    return agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator, tn
