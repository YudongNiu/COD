import sys
import ast
from queue import Queue
from time import perf_counter

LAMDA = 0.9


def getShell(v, agglom_nodes, agglom_tree, top_down_tree, leaf_inner_indicator, tn=-1):
    v_path = [0 for _ in range(len(agglom_tree))]
    p = agglom_nodes[v]
    if tn >= 0:
        tn = agglom_tree[tn]
    while agglom_tree[p] != tn:
        v_path[p] = 1
        p = agglom_tree[p]

    c_indicator = [-1 for _ in range(len(agglom_nodes))]
    tnode_queue = Queue()
    core_value_queue = Queue()
    tnode_queue.put(p)
    core_value_queue.put(0)
    while not tnode_queue.empty():
        tnode = tnode_queue.get()
        core_value = core_value_queue.get()
        child0 = top_down_tree[tnode][0]
        child1 = top_down_tree[tnode][1]
        if leaf_inner_indicator[tnode][0] == 0:
            c_indicator[child0] = core_value
        else:
            tnode_queue.put(child0)
            if v_path[child0] == 1:
                core_value_queue.put(core_value + 1)
            else:
                core_value_queue.put(core_value)

        if leaf_inner_indicator[tnode][1] == 0:
            c_indicator[child1] = core_value
        else:
            tnode_queue.put(child1)
            if v_path[child1] == 1:
                core_value_queue.put(core_value + 1)
            else:
                core_value_queue.put(core_value)
    return c_indicator
