import networkx as nx
import sys
import ast
from time import perf_counter
import math

import lore
import hierarchical_clusters as h


# ************************
# himor-index construction
# ************************
def DupRank(items):
    ranks = dict()
    pointer = 0
    pre_inf = items[pointer][1]
    while pointer < len(items):
        item_inf = items[pointer][1]
        if item_inf < pre_inf:
            back_pointer = pointer - 1
            while items[back_pointer][1] == pre_inf:
                ranks[items[back_pointer][0]] = pointer
                back_pointer -= 1
            pre_inf = item_inf
        pointer += 1
    pointer = len(items) - 1
    item_inf = items[pointer][1]
    while items[pointer][1] == item_inf and pointer >= 0:
        ranks[items[pointer][0]] = len(items)
        pointer -= 1
    return ranks


def influence_pairs_cpr_prop_n(g, agglom_nodes, agglom_tree, tn_level, infs):
    tn = len(agglom_tree) - 1
    activated_levels = [-1 for _ in range(g.number_of_nodes())]
    for rrsource in g:
        rrsource_ancs = set()

        temp_p = agglom_nodes[rrsource]
        while 0 <= temp_p <= tn:
            rrsource_ancs.add(temp_p)
            temp_p = agglom_tree[temp_p]

        for _ in range(hus.NUM_RRSET):
            new_active = [rrsource]
            activated_levels[rrsource] = tn_level[agglom_nodes[rrsource]] - tn_level[tn]
            activated_nodes = [rrsource]
            level_que = [activated_levels[rrsource]]

            while new_active:
                u = new_active.pop(0)
                lu = level_que.pop(0)
                p = (1.0 / g.degree[u])
                for w in g[u]:
                    pp = agglom_nodes[w]
                    while not (pp in rrsource_ancs):
                        pp = agglom_tree[pp]
                    c_indicator_w = tn_level[pp] - tn_level[tn]
                    l = min(lu, c_indicator_w)
                    if l > activated_levels[w]:
                        rand = random.random()
                        if rand < p:
                            new_active.append(w)
                            level_que.append(l)
                            for tt in range(activated_levels[w] + 1, l + 1):
                                infs[w][tt + tn_level[tn]] += 1
                            activated_levels[w] = l
                            activated_nodes.append(w)
            for anode in activated_nodes:
                activated_levels[anode] = -1


def rub_cpr_prop_n(g, agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator, rub_file):
    start = perf_counter()

    index = dict()
    for n in g:
        index[n] = []

    agglom_node_level, tn_level = lore.agglom_level_func(agglom_tree, top_down_tree, leaf_inner_indicator)
    infs = [[0 for _ in range(agglom_node_level[n])] for n in range(len(agglom_tree) + 1)]

    influence_pairs_cpr_prop_n(g, agglom_nodes, agglom_tree, tn_level, infs)

    que = [len(agglom_tree) - 1]
    while que:
        temp_tn = que.pop(0)
        cc0 = top_down_tree[temp_tn][0]
        cc1 = top_down_tree[temp_tn][1]
        if leaf_inner_indicator[temp_tn][0] == 1:
            que.append(cc0)
        if leaf_inner_indicator[temp_tn][1] == 1:
            que.append(cc1)

        com = []
        stk = [temp_tn]
        while stk:
            pp = stk.pop(0)
            cc0 = top_down_tree[pp][0]
            cc1 = top_down_tree[pp][1]
            if leaf_inner_indicator[pp][0] == 1:
                stk.append(cc0)
            else:
                com.append([cc0, infs[cc0][tn_level[temp_tn]]])
            if leaf_inner_indicator[pp][1] == 1:
                stk.append(cc1)
            else:
                com.append([cc1, infs[cc1][tn_level[temp_tn]]])
        sort_com = sorted(com, key=lambda kv: (kv[1], kv[0]), reverse=True)
        ranks = DupRank(sort_com)

        for n in ranks:
            n_inf = infs[n][tn_level[temp_tn]] / hus.NUM_RRSET
            if len(index[n]) == 0:
                index[n].append([temp_tn, ranks[n], n_inf])
            else:
                r_gap = index[n][-1][1] - ranks[n]
                i_gap = index[n][-1][2] - n_inf
                if i_gap < 0:
                    index[n][-1][2] = n_inf
                if r_gap > 0:
                    index[n].append([temp_tn, ranks[n], n_inf])

    for n in range(len(index)):
        max_inf = index[n][len(index[n]) - 1][2]
        for t in range(len(index[n]) - 1, -1, -1):
            if index[n][t][2] < max_inf:
                index[n][t][2] = max_inf
            else:
                max_inf = index[n][t][2]

    end = perf_counter()
    indexing_time = end - start

    with open(rub_file, 'w') as f:
        for n in range(g.number_of_nodes()):
            index_n = index[n]
            for index_item in index_n:
                f.write(str(index_item[0]) + ',' + str(index_item[1]) + ',' + str(index_item[2]) + ' ')
            f.write('\n')
        f.write("# index_time:" + str(indexing_time) + '\n')
        f.write("# RRSET_NUM:" + str(hus.NUM_RRSET) + '\n')


# *******************
# himor-index loading
# *******************
def index_load(g, index_file, T=-1):
    index = dict()
    for n in range(g.number_of_nodes()):
        index[n] = []
    with open(index_file, 'r') as f:
        line = f.readline()
        node = 0
        while line:
            if line.startswith('#'):
                line = f.readline()
                continue
            line = line[:-2]
            line_split = line.split(' ')

            interval = 0
            if T > 0:
                interval = math.ceil(len(line_split) * T)
            if interval < 1:
                interval = 1
            for item_count in range(len(line_split)):
                item = line_split[item_count]
                item_split = item.split(',')
                tn = int(item_split[0])
                r = int(item_split[1])
                if item_count % interval == 0:
                    index[node].append((tn, r))
            node += 1
            line = f.readline()
    return index


def shell_from_index_rb(q, agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator, local_tn, index, K):
    # c_indicator generation based on ranks in index
    q_index = index[q]
    q_pointer = len(q_index) - 1

    while q_pointer > 0 and (q_index[q_pointer][1] < K or q_index[q_pointer][0] < local_tn):
        q_pointer -= 1
    tn = q_index[q_pointer][0]
    c_indicator = h.getShell(q, agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator, tn)
    return c_indicator


if __name__ == "__main__":
    dataset = sys.argv[1]

    path = 'labeled_graph/' + dataset + '/' + dataset
    g = nx.Graph()
    with open(path + '.txt', 'r') as txt:
        line = txt.readline()
        while line:
            line = line.split(' ')
            n0 = int(line[0])
            n1 = int(line[1])
            if n0 != n1:
                g.add_edge(n0, n1)
            line = txt.readline()

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
    index_file = path + '.himor_index'
    rub_cpr_prop_n(g, agglom_nodes_ul, agglom_tree_ul, top_down_tree_ul, leaf_inner_indicator_ul, index_file)
