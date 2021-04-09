import csv
import time
import os.path

import numpy as np

import parameters
from utils import (
    generate_cycle_types,
    create_new_directory_in_cycles_info,
    create_new_directory_for_logging_experiments,
)


def read_parameters_from_console():
    print("Enter number of fragile regions")
    n = int(input())

    print("Enter number of steps")
    k = int(input())

    print("Enter probability of rearrangement between A-type fragile regions:")
    p_aa = float(input())

    print("Enter probability of rearrangement between B-type fragile regions:")
    p_bb = float(input())

    # Probability of rearrangement between A-type and D-type fragile regions.
    p_ab = 1 - p_aa - p_bb

    assert p_aa >= 0 and p_bb >= 0 and p_ab >= 0

    return n, k, p_aa, p_bb, p_ab


class EdgeInfo:
    to = 0
    edge_type = "A"

    def __init__(self, to, edge_type):
        self.to = to
        self.edge_type = edge_type


def split_fragile_edges(n, a_type_edges_proportion):
    a_type = []
    b_type = []
    edges = []

    will_be_a_type = np.random.choice(
        range(0, n), int(a_type_edges_proportion * n), replace=False
    )
    is_a_type = [False] * n
    for e in will_be_a_type:
        is_a_type[e] = True

    for i in range(n):
        edge = [2 * i, 2 * i + 1]
        if is_a_type[i]:
            a_type.append(edge)
            edges.append(EdgeInfo(2 * i + 1, "A"))
            edges.append(EdgeInfo(2 * i, "A"))
        else:
            b_type.append(edge)
            edges.append(EdgeInfo(2 * i + 1, "B"))
            edges.append(EdgeInfo(2 * i, "B"))

    return a_type, b_type, edges


def is_in_cycle(x, edges, y):
    cur = edges[x].to

    if y == cur:
        return True

    while cur != x:
        if cur % 2 == 0:
            cur = cur + 1
        else:
            cur = cur - 1

        if y == cur:
            return True

        if cur == x:
            break

        cur = edges[cur].to

        if y == cur:
            return True

    return False


def set_color(x, edges, colors, color):
    colors[x] = color
    cur = edges[x].to

    while cur != x:
        colors[cur] = color

        if cur % 2 == 0:
            cur = cur + 1
        else:
            cur = cur - 1

        if cur == x:
            break

        colors[cur] = color
        cur = edges[cur].to


def update_cycles(edge1, edge2, joining_type, edges, colors, last_color):
    x, y = edge1
    u, v = edge2

    if joining_type == "xu":
        new_edge1_type, new_edge2_type = np.random.choice(
            [edges[x].edge_type, edges[u].edge_type], 2, replace=False
        )

        edges[x], edges[u] = EdgeInfo(u, new_edge1_type), EdgeInfo(x, new_edge1_type)
        edges[y], edges[v] = EdgeInfo(v, new_edge2_type), EdgeInfo(y, new_edge2_type)

        new_edges = [[x, u], [y, v]]
        if new_edge1_type != new_edge2_type and new_edge1_type == "B":
            new_edges.reverse()
    else:
        new_edge1_type, new_edge2_type = np.random.choice(
            [edges[x].edge_type, edges[v].edge_type], 2, replace=False
        )

        edges[x], edges[v] = EdgeInfo(v, new_edge1_type), EdgeInfo(x, new_edge1_type)
        edges[y], edges[u] = EdgeInfo(u, new_edge2_type), EdgeInfo(y, new_edge2_type)

        new_edges = [[x, v], [y, u]]
        if new_edge1_type != new_edge2_type and new_edge1_type == "B":
            new_edges.reverse()

    if colors[x] == colors[u]:
        if is_in_cycle(x, edges, y):
            new_color = last_color
        else:
            set_color(x, edges, colors, last_color)
            new_color = last_color + 1
    else:
        set_color(x, edges, colors, colors[x])
        new_color = last_color

    return new_color, new_edges


class CyclesInfo:
    num_all_cycles = 0
    # cycles with length from 1 to 5 with different types
    num_n_cycles = {}
    max_len = 0

    def __init__(self, num_all_cycles, num_n_cycles, max_len):
        self.num_all_cycles = num_all_cycles
        self.num_n_cycles = num_n_cycles
        self.max_len = max_len


# Возвращает:
# количество циклов;
# [число циклов длины 1: (с хрупким ребром A; с хрупким ребром B),
#  длины 2: (с хрупкими ребрами A, A; с хрупкими ребрами A, B; с хрупкими ребрами B, B),
#  и так далее до длины 5.];
# длина максимального цикла.
# Длина цикла измеряется в количестве хрупких ребер в нем.
def compute_cycles_info(colors, edges):
    cycles = {}
    for i, c in enumerate(colors):
        if c in cycles.keys():
            cycles[c].append(edges[i].edge_type)
        else:
            cycles[c] = [edges[i].edge_type]

    cycles_info = {}
    max_len = 0
    for c in cycles.keys():
        cur_len = len(cycles[c]) // 2
        if cur_len < 6:
            cycle_type = to_cycle_type(cycles[c])
            if cycle_type in cycles_info:
                cycles_info[cycle_type] += 1
            else:
                cycles_info[cycle_type] = 1
        max_len = max(max_len, cur_len)

    return CyclesInfo(len(cycles), cycles_info, max_len)


def to_cycle_type(edges_type):
    edges_type.sort()
    cycle_type = ""
    for i in range(0, len(edges_type), 2):
        cycle_type += edges_type[i]
    return cycle_type


def markov_process(n, k, p_aa, p_bb, p_ab, a_type, b_type, edges):
    steps_cycles_info = []
    len_a = len(a_type)
    len_b = len(b_type)

    probabilities = [p_aa, p_bb, p_ab]

    colors = [0] * (2 * n)
    for i in range(n):
        colors[2 * i] = i
        colors[2 * i + 1] = i

    new_color = n

    steps_cycles_info.append(compute_cycles_info(colors, edges))

    for i in range(k):
        dcj_type = np.random.choice(["aa", "bb", "ab"], 1, p=probabilities)[0]
        joining_type = np.random.choice(["xu", "xv"], 1)[0]

        if dcj_type == "aa":
            edge_ind1, edge_ind2 = np.random.choice(len(a_type), 2, replace=False)
            edge1 = a_type[edge_ind1]
            edge2 = a_type[edge_ind2]

            new_color, new_edges = update_cycles(
                edge1, edge2, joining_type, edges, colors, new_color
            )
            a_type[edge_ind1], a_type[edge_ind2] = new_edges
        elif dcj_type == "bb":
            edge_ind1, edge_ind2 = np.random.choice(len(b_type), 2, replace=False)
            edge1 = b_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_color, new_edges = update_cycles(
                edge1, edge2, joining_type, edges, colors, new_color
            )
            b_type[edge_ind1], b_type[edge_ind2] = new_edges
        else:
            edge_ind1 = np.random.choice(len(a_type), 1)[0]
            edge_ind2 = np.random.choice(len(b_type), 1)[0]

            edge1 = a_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_color, new_edges = update_cycles(
                edge1, edge2, joining_type, edges, colors, new_color
            )
            a_type[edge_ind1], b_type[edge_ind2] = new_edges

        steps_cycles_info.append(compute_cycles_info(colors, edges))

        assert len_a == len(a_type) and len_b == len(b_type)

    return steps_cycles_info


def log_experiments(experiments, file, max_cycle_len=5):
    file_exists = os.path.isfile(file)

    with open(file, "a", newline="") as f_log:
        cycle_types = generate_cycle_types(1, max_cycle_len)
        fieldnames = ["k", "all"] + cycle_types + ["max_cycle_len"]

        log_cycles = csv.DictWriter(f_log, fieldnames=fieldnames)
        if not file_exists:
            log_cycles.writeheader()

        for experiment in experiments:
            for step, info in enumerate(experiment):
                cur_result = {
                    "k": step,
                    "all": info.num_all_cycles,
                    "max_cycle_len": info.max_len,
                }
                for cycle_type in cycle_types:
                    if cycle_type in info.num_n_cycles:
                        cur_result[cycle_type] = info.num_n_cycles[cycle_type]
                    else:
                        cur_result[cycle_type] = 0
                log_cycles.writerow(cur_result)
        f_log.close()


def main():
    # n, k, p_aa, p_bb, p_ab = read_parameters_from_console()

    start_time = time.time()

    n = parameters.NUMBER_OF_FRAGILE_EDGES
    k = n
    experiments_log_path = create_new_directory_for_logging_experiments()

    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        string_parameters, p_aa, p_bb, a_type_edges_proportion = parameter
        file = string_parameters + ".csv"
        p_ab = 1 - p_aa - p_bb

        start_i = 0
        if string_parameters == "paa0,2_pbb0,6_alpha0,9":
            start_i = parameters.DIFFERENT_FRAGILE_EDGES_SPLITS
        elif string_parameters == "paa0,4_pbb0,35_alpha0,7":
            start_i = 8699

        experiments = []
        print("parameters:", string_parameters)
        for i in range(start_i, parameters.DIFFERENT_FRAGILE_EDGES_SPLITS):
            a_type, b_type, edges = split_fragile_edges(n, a_type_edges_proportion)
            experiments.append(
                markov_process(
                    n,
                    k,
                    p_aa,
                    p_bb,
                    p_ab,
                    a_type.copy(),
                    b_type.copy(),
                    edges.copy(),
                )
            )

            if len(experiments) == 100:
                print("i:", i, ", time: ", time.time() - start_time, " s.")
                log_experiments(experiments, experiments_log_path + file)
                experiments = []

    print(time.time() - start_time)


if __name__ == "__main__":
    main()
