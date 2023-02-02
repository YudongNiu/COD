import networkx as nx
import sys
from time import perf_counter
import ast

import lore
import hierarchical_clusters as h


if __name__ == "__main__":
    dataset = sys.argv[1]
    linkage = sys.argv[2]  # [1]label_edge [2]unlabel
    qfeat = int(sys.argv[3])

    data_path = 'labeled_graph/' + dataset + '/' + dataset
    g = nx.Graph()
    with open('labeled_graph/' + dataset + '/' + dataset + '.txt', 'r') as txt:
        line = txt.readline()
        while line:
            line = line.split(' ')
            n0 = int(line[0])
            n1 = int(line[1])
            if n0 != n1:
                g.add_edge(n0, n1)
            line = txt.readline()

    feats = []
    feat_path = data_path + '.feat'
    with open(feat_path, 'r') as gf:
        line = gf.readline()
        while line:
            feat = ast.literal_eval(line)
            if type(feat) is int:
                if qfeat == feat:
                    feats.append(True)
                else:
                    feats.append(False)
            elif type(feat) is list:
                if qfeat in feat:
                    feats.append(True)
                else:
                    feats.append(False)
            line = gf.readline()

    if linkage == 'unlabel':
        out_file = data_path + '_unlabel.hac'
    else:
        out_file = data_path + '_' + linkage + str(h.LAMDA) + '_qlabel' + str(qfeat) + '.hac'


    start = perf_counter()
    agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator = lore.nnc_with_candidates(g, feats, linkage)
    end = perf_counter()
    time = end - start

    with open(out_file, 'w') as f:
        f.write("agglom_nodes:" + str(agglom_nodes) + "\n")
        f.write("agglom_tree:" + str(agglom_tree) + "\n")
        f.write("top_down_tree:" + str(top_down_tree) + "\n")
        f.write("leaf_inner_indicator:" + str(leaf_inner_indicator) + "\n")
        f.write("# time:" + str(time) + "\n")
