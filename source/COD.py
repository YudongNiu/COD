import networkx as nx
import sys
import random
import ast
from time import perf_counter
import math
import numpy as np

import hierarchical_clusters as h
import lore
import himor


NUM_RRSET = 10


def ocs_ic_cpr_prop_n(g, qn, K, c_indicator):
    influence_matrix = [[0 for _ in range(c_indicator[n] + 1)] for n in range(g.number_of_nodes())]
    # influence_matrix[n][l] denote the estimated influence of node n at level l
    activated_levels = [-1 for _ in range(g.number_of_nodes())]
    for rrsource in g:
        if c_indicator[rrsource] < 0:
            pass
        for _ in range(NUM_RRSET):
            new_active = [rrsource]
            activated_levels[rrsource] = c_indicator[rrsource]
            activated_nodes = [rrsource]
            level_que = [c_indicator[rrsource]]
            while new_active:
                u = new_active.pop(0)
                lu = level_que.pop(0)
                p = (1.0 / g.degree[u])
                for w in g[u]:
                    if c_indicator[w] >= 0:
                        l = min(lu, c_indicator[w])
                        if l > activated_levels[w]:
                            rand = random.random()
                            if rand < p:
                                new_active.append(w)
                                level_que.append(l)
                                if activated_levels[w] >= 0:
                                    influence_matrix[w][activated_levels[w]] -= 1
                                influence_matrix[w][l] += 1
                                activated_levels[w] = l
                                activated_nodes.append(w)
            for anode in activated_nodes:
                activated_levels[anode] = -1
    
    for i in range(len(influence_matrix)):
        sum = 0
        for j in range(len(influence_matrix[i]) - 1, -1, -1):
            sum += influence_matrix[i][j]
            influence_matrix[i][j] = sum

    depth = max(c_indicator) + 1
    for level in range(depth):
        candidates = []
        for n in range(len(c_indicator)):
            if c_indicator[n] >= level:
                candidates.append(n)
        rank = 0
        qn_influence = influence_matrix[qn][level]
        for n in candidates:
            if influence_matrix[n][level] >= qn_influence:
                rank += 1
        if rank <= K:
            return candidates
    return []


def ocs_ic_prop_n(g, qn, K, c_indicator):
    depth = max(c_indicator) + 1
    for level in range(depth):
        influences = [0 for _ in range(g.number_of_nodes())]
        candidates = []
        for n in range(len(c_indicator)):
            if c_indicator[n] >= level:
                candidates.append(n)

        activated = [False for _ in range(g.number_of_nodes())]
        for rrsource in candidates:
            for _ in range(NUM_RRSET):
                new_active = [rrsource]
                activated[rrsource] = True
                activated_nodes = [rrsource]

                while new_active:
                    u = new_active.pop(0)
                    p = (1.0 / g.degree[u])
                    for w in g[u]:
                        if c_indicator[w] >= level and (not activated[w]):
                            rand = random.random()
                            if rand < p:
                                new_active.append(w)
                                influences[w] += 1
                                activated[w] = True
                                activated_nodes.append(w)
                for anode in activated_nodes:
                    activated[anode] = False

        rank = 0
        qn_influence = influences[qn]
        for n in candidates:
            if influences[n] >= qn_influence:
                rank += 1
        if rank <= K:
            return candidates
    return []


if __name__ == "__main__":
    dataset = sys.argv[1]
    labeled = int(sys.argv[2])
    K = int(sys.argv[3])
    cpr = int(sys.argv[4])  # compress
    loc = int(sys.argv[5])  # if local recluster over $C_\ell$
    indexed = int(sys.argv[6])  # if use the himor index

    if labeled == 0:
        linkage = 'unlabel'
    else:
        linkage = 'label_edge'

    # Load Graph
    path = 'labeled_graph/' + dataset + '/' + dataset
    G = nx.Graph()
    with open(path + '.txt', 'r') as txt:
        line = txt.readline()
        while line:
            line = line.split(' ')
            n0 = int(line[0])
            n1 = int(line[1])
            if n0 != n1:
                G.add_edge(n0, n1)
            line = txt.readline()

    # Load himor-index
    index_file = path + '.himor'
    index = himor.index_load(G, index_file)

    # Load the unlabel HAC
    hac_path = path + '_unlabel.hac'
    with open(hac_path, 'r') as f:
        line = f.readline()
        agglom_nodes_ul = ast.literal_eval(line.split(':')[1])
        line = f.readline()
        agglom_tree_ul = ast.literal_eval(line.split(':')[1])
        line = f.readline()
        top_down_tree_ul = ast.literal_eval(line.split(':')[1])
        line = f.readline()
        leaf_inner_indicator_ul = ast.literal_eval(line.split(':')[1])

    # Load Queries
    qnodes = []
    qfeats = []
    with open(path + '.query_nodes', 'r') as f:
        line = f.readline()
        while line:
            line_split = line.split(' ')
            qnodes.append(int(line_split[0]))
            qfeats.append(int(line_split[1]))
            line = f.readline()

    # Load feats
    feats = []
    feat_path = path + '.feat'
    with open(feat_path, 'r') as f:
        line = f.readline()
        while line:
            feat = ast.literal_eval(line)
            if type(feat) is int:
                feats.append([feat])
            else:
                if 0 in feat:
                    feat.remove(0)
                feats.append(feat)
            line = f.readline()

    # HUS
    sizes = []
    topo_dens = []
    attr_dens = []
    c_indicator_times = []
    ocs_times = []
    for q_count in range(len(qnodes)):
        q = qnodes[q_count]
        qfeat = qfeats[q_count]
        bin_feats = []
        for feat_n in feats:
            if qfeat in feat_n:
                bin_feats.append(True)
            else:
                bin_feats.append(False)

        # computing the c_indicator
        starttime = perf_counter()
        c_indicator = []
        if linkage == 'unlabel':
            c_indicator = h.getShell(q, agglom_nodes_ul, agglom_tree_ul, top_down_tree_ul, leaf_inner_indicator_ul)
        elif linkage == 'label_edge':
            if loc == 1:
                agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator, local_tn = \
                    lore.LocalHAC(G, q, bin_feats, agglom_nodes_ul, agglom_tree_ul, top_down_tree_ul,
                                   leaf_inner_indicator_ul)
                if indexed == 1:
                    c_indicator = himor.shell_from_index_rb(q, agglom_nodes, agglom_tree, top_down_tree,
                                                           leaf_inner_indicator, local_tn, index, K)
                else:
                    c_indicator = h.getShell(q, agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator)
            else:
                agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator = \
                    lore.nnc_with_candidates(G, bin_feats, linkage)
                c_indicator = h.getShell(q, agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator)
        endtime = perf_counter()
        c_indicator_time = endtime - starttime
        c_indicator_times.append(c_indicator_time)

        # computing the outstanding community
        starttime = perf_counter()
        com = []
        if cpr == 1:
            com = ocs_ic_cpr_prop_n(G, q, K, c_indicator)
        elif cpr == 0:
            com = ocs_ic_prop_n(G, q, K, c_indicator)
        endtime = perf_counter()
        ocs_time = endtime - starttime
        ocs_times.append(ocs_time)

        if len(com) > 0:
            sizes.append(len(com))
            F = G.subgraph(com).copy()
            topo_dens.append(nx.density(F))
            temp_attr_dens = 0
            for n in com:
                if qfeat in feats[n]:
                    temp_attr_dens += 1.0
            temp_attr_dens /= len(com)
            attr_dens.append(temp_attr_dens)
        else:
            sizes.append(0)
            topo_dens.append(0)
            attr_dens.append(0)

    if linkage == 'unlabel':
        out_path = path + '_K' + str(K) + '.ocs_ic' + '_unlabel_cpr' + str(cpr)
    elif linkage == 'label_edge':
        if loc == 1:
            if indexed == 1:
                out_path = path + '_K' + str(K) + '.ocs_ic' + '_' + linkage + str(h.LAMDA) + '_loc_himor_cpr' + str(cpr)
            else:
                out_path = path + '_K' + str(K) + '.ocs_ic' + '_' + linkage + str(h.LAMDA) + '_loc_cpr' + str(cpr)
        else:
            out_path = path + '_K' + str(K) + '.ocs_ic' + '_' + linkage + str(h.LAMDA) + '_cpr' + str(cpr)

    with open(out_path, 'w') as f:
        for q_count in range(len(qnodes)):
            f.write('size:'+str(sizes[q_count])+'\n')
            f.write('topo_dens:'+str(topo_dens[q_count])+'\n')
            f.write('attr_dens:'+str(attr_dens[q_count])+'\n')
            f.write('shell_time:'+str(c_indicator_times[q_count])+'\n')
            f.write('ocs_time:'+str(ocs_times[q_count])+'\n')
            f.write('\n')
        f.write('_avg_size:'+str(np.mean(sizes))+'\n')
        f.write('_avg_topo_dens:'+str(np.mean(topo_dens))+'\n')
        f.write('_avg_attr_dens:'+str(np.mean(attr_dens))+'\n')
        f.write('_avg_shell_time:'+str(np.mean(c_indicator_times))+'\n')
        f.write('_avg_ocs_time:'+str(np.mean(ocs_times))+'\n')
